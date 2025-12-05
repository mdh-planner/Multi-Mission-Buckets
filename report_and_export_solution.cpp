// ------------------------------------------------------------
// Report + export solved plan
// ------------------------------------------------------------
static void reportAndExportSolution(
    const ProblemData& D,
    const SolveOptions& opt,
    IloCplex& cplex,
    const IloArray<IloBoolVarArray>& xDep,
    const IloArray< IloArray<IloBoolVarArray> >& xTask,
    const IloArray<IloBoolVarArray>& xEnd,
    const IloArray<IloNumVarArray>& t,
    const IloNumVarArray& T,
    // pass pointers; set to nullptr if not used in this baseline
    const IloArray< IloArray<IloBoolVarArray> >* yPtr,
    const IloArray< IloArray<IloBoolVarArray> >* zPtr,
    bool credit_spillover = true
) {
    const int R = D.R, N = D.N, Mmissions = D.Mmissions, H = D.H;

    const bool activeBuckets = hasActiveBucketAccounting(D);
    const bool useY = (yPtr != nullptr) &&
        opt.useBucketLayer && opt.useQuotas && opt.useReservations && activeBuckets;
    const bool useZ = (zPtr != nullptr) &&
        opt.useBucketLayer && opt.useQuotas && !opt.useReservations && activeBuckets;

    auto kpi = evaluateApplesToApples(
        D, cplex, xDep, xTask, xEnd, t, T, credit_spillover,
        useY ? yPtr : nullptr, useZ ? zPtr : nullptr
    );

    std::cout << "\n[Baseline] "
        << (opt.baseline == Baseline::OURS ? "OURS" :
            opt.baseline == Baseline::MONOLITHIC ? "MONOLITHIC" :
            opt.baseline == Baseline::MONO_SOFT_BUCKETS ? "MONO_SOFT_BUCKETS" :
            opt.baseline == Baseline::WEIGHTS_ONLY ? "WEIGHTS_ONLY" : "DECOUPLED")
        << "\n";

    std::cout << "\n--- Apples-to-Apples Evaluation ---\n";
    std::cout << "Raw CPLEX objective     : " << std::fixed << std::setprecision(3) << kpi.rawObj << "\n";
    std::cout << "A2A score (common)      : " << kpi.score << "\n";
    std::cout << "Sum w[m] T[m]           : "
        << [&]() { double s = 0; for (size_t m = 0; m < kpi.T.size(); ++m) s += D.w[m] * kpi.T[m]; return s; }() << "\n";
    std::cout << "  Travel (weighted)     : " << D.lambdaTravel * kpi.travel << "  (travel=" << kpi.travel << ")\n";
    std::cout << "  Switches (weighted)   : " << D.etaSwitch * D.switchPenalty * kpi.switches
        << "  (#switches=" << kpi.switches << ")\n";
    std::cout << "  Quota shortfall term  : " << D.rhoQ * kpi.shortQ << "  (shortQ=" << kpi.shortQ << ")\n";
    std::cout << "  Rolling total term    : " << D.rhoR * kpi.shortR << "  (shortR=" << kpi.shortR << ")\n";

    // Per-bucket table
    if (!D.q.empty()) {
        std::cout << "\nPer-bucket robots reserved (inferred) vs quota:\n";
        for (int h = 0; h < D.H; ++h) {
            std::cout << "  Bucket " << h << ": ";
            for (int m = 0; m < D.Mmissions; ++m) {
                int got = kpi.cap[m][h];
                int req = D.q[m][h];
                std::cout << "M" << m << "=" << got << "/" << req << "  ";
            }
            std::cout << "\n";
        }
    }

    // Print y (OURS)
    if (useY) {
        const auto& y = *yPtr;
        std::cout << "\nReservations y[i][h][m]=1:\n";
        for (int h = 0; h < H; ++h) {
            std::cout << "  Bucket " << h << ": ";
            bool any = false;
            for (int i = 0; i < R; ++i)
                for (int m = 0; m < Mmissions; ++m) {
                    const IloBoolVar& v = y[i][h][m];
                    if (cplex.isExtracted(v) && cplex.getValue(v) > 0.5) {
                        std::cout << "[i=" << i << ",m=" << m << "] ";
                        any = true;
                    }
                }
            if (!any) std::cout << "(none)";
            std::cout << "\n";
        }
    }

    // Print z (soft mono)
    if (useZ) {
        const auto& z = *zPtr;
        std::cout << "\nPer-bucket robots (z, soft mono) vs quota:\n";
        for (int h = 0; h < H; ++h) {
            std::cout << "  Bucket " << h << ": ";
            for (int m = 0; m < Mmissions; ++m) {
                int got = 0;
                for (int i = 0; i < R; ++i) {
                    const IloBoolVar& v = z[i][h][m];
                    if (cplex.isExtracted(v) && cplex.getValue(v) > 0.5) ++got;
                }
                std::cout << "M" << m << "=" << got << "/" << (D.q.empty() ? 0 : D.q[m][h]) << "  ";
            }
            std::cout << "\n";
        }
    }

    // ===== Per-robot travel legs (with earliest arrival & waiting) =====
    struct LegRow {
        int  i;
        std::string from;
        std::string to;
        int  fromTask;     // -1 if depot
        int  toTask;       // -1 if end/depot
        double travel;
        double departEarliest;
        double arriveEarliest;
        double startAtDest;
        double waitAtDest;
        bool switchMission;
        // geometry
        double fromX, fromY;
        double toX, toY;
        int    fromMission; // -1 if depot
        int    toMission;   // -1 if depot/end
    };

    std::vector<LegRow> legs;
    legs.reserve(R * (N + 1));

    auto fmtNode = [&](int i, int j, bool isDepot, bool isEnd)->std::string {
        if (isDepot) return "depot_" + std::to_string(i);
        if (isEnd)   return "end_" + std::to_string(i);
        return "task " + std::to_string(j);
    };

    for (int i = 0; i < R; ++i) {
        int cur = -1;
        for (int j = 0; j < N; ++j) if (cplex.getValue(xDep[i][j]) > 0.5) { cur = j; break; }
        if (cur == -1) continue;

        // depot -> first task
        {
            double travel = D.cDep[i][cur];
            double departEarliest = 0.0;
            double arriveEarliest = departEarliest + travel;
            double startAtDest = cplex.getValue(t[i][cur]);
            double waitAtDest = std::max(0.0, startAtDest - arriveEarliest);
            double fx = D.depotXY[i].x, fy = D.depotXY[i].y;
            double tx = D.taskXY[cur].x, ty = D.taskXY[cur].y;
            int fm = -1, tm = D.missionOf[cur];

            legs.push_back({ i, fmtNode(i,-1,true,false), fmtNode(i,cur,false,false),
                             -1, cur, travel, departEarliest, arriveEarliest,
                             startAtDest, waitAtDest, false,
                             fx, fy, tx, ty, fm, tm });
        }

        while (true) {
            int nxt = -1;
            for (int k = 0; k < N; ++k)
                if (k != cur && cplex.getValue(xTask[i][cur][k]) > 0.5) { nxt = k; break; }

            if (nxt != -1) {
                double travel = D.cTask[cur][nxt];
                double departEarliest = cplex.getValue(t[i][cur]) + D.s[cur];
                double arriveEarliest = departEarliest + travel;
                double startAtDest = cplex.getValue(t[i][nxt]);
                double waitAtDest = std::max(0.0, startAtDest - arriveEarliest);
                bool   switchM = (D.missionOf[cur] != D.missionOf[nxt]);
                double fx = D.taskXY[cur].x, fy = D.taskXY[cur].y;
                double tx = D.taskXY[nxt].x, ty = D.taskXY[nxt].y;
                int fm = D.missionOf[cur], tm = D.missionOf[nxt];

                legs.push_back({ i, fmtNode(i,cur,false,false), fmtNode(i,nxt,false,false),
                                 cur, nxt, travel, departEarliest, arriveEarliest,
                                 startAtDest, waitAtDest, switchM,
                                 fx, fy, tx, ty, fm, tm });

                cur = nxt;
                continue;
            }

            if (cplex.getValue(xEnd[i][cur]) > 0.5) {
                double travel = D.cEnd[i][cur];
                double departEarliest = cplex.getValue(t[i][cur]) + D.s[cur];
                double arriveEarliest = departEarliest + travel;
                double fx = D.taskXY[cur].x, fy = D.taskXY[cur].y;
                double tx = D.depotXY[i].x, ty = D.depotXY[i].y;
                int fm = D.missionOf[cur], tm = -1;

                legs.push_back({ i, fmtNode(i,cur,false,false), fmtNode(i,-1,false,true),
                                 cur, -1, travel, departEarliest, arriveEarliest,
                                 std::numeric_limits<double>::quiet_NaN(), 0.0, false,
                                 fx, fy, tx, ty, fm, tm });
            }
            break;
        }
    }

    // Print per robot
    std::cout << "\n=== Travel Legs (per robot) ===\n";
    std::cout << std::fixed << std::setprecision(2);
    for (int i = 0; i < R; ++i) {
        double totalTravel = 0.0;
        bool any = false;
        for (const auto& L : legs) if (L.i == i) { any = true; break; }
        if (!any) continue;

        std::cout << "\nRobot " << i << ":\n";
        std::cout << "  From        -> To          Travel  Depart*  Arrive*  Start@To  Wait  Switch\n";
        for (const auto& L : legs) if (L.i == i) {
            totalTravel += L.travel;
            std::cout << "  " << std::setw(10) << L.from
                << " -> " << std::setw(10) << L.to
                << "  " << std::setw(6) << L.travel
                << "  " << std::setw(7) << L.departEarliest
                << "  " << std::setw(7) << L.arriveEarliest
                << "  " << std::setw(8) << (L.toTask >= 0 ? std::to_string(L.startAtDest) : "-")
                << "  " << std::setw(5) << (L.toTask >= 0 ? std::to_string(L.waitAtDest) : "0")
                << "  " << (L.switchMission ? "yes" : "no") << "\n";
        }
        std::cout << "  Total travel time = " << totalTravel << "\n";
    }
    std::cout << "  (* earliest times assuming no waiting before departure)\n";

    // ===== EXPORT: JSON + CSV =====
    const std::string baselineTag =
        (opt.baseline == Baseline::OURS) ? "OURS" :
        (opt.baseline == Baseline::MONOLITHIC) ? "MONOLITHIC" :
        (opt.baseline == Baseline::MONO_SOFT_BUCKETS) ? "MONO_SOFT_BUCKETS" :
        (opt.baseline == Baseline::WEIGHTS_ONLY) ? "WEIGHTS_ONLY" : "DECOUPLED";

    fs::create_directories(baselineTag);

    auto path = [&](const std::string& fname) {
        return (fs::path(baselineTag) / fname).string();
    };

    // 1) meta
    {
        std::ofstream jf(path("sol_meta.json"));
        jf << std::fixed << std::setprecision(6);
        jf << "{\n";
        jf << "  \"R\": " << D.R << ",\n";
        jf << "  \"N\": " << D.N << ",\n";
        jf << "  \"M\": " << D.Mmissions << ",\n";
        jf << "  \"H\": " << D.H << ",\n";
        jf << "  \"Delta\": " << D.Delta << ",\n";
        jf << "  \"objective\": " << cplex.getObjValue() << ",\n";
        jf << "  \"baseline\": \"" << baselineTag << "\"\n";
        jf << "}\n";
    }

    // 2) schedule
    {
        std::ofstream sch(path("sol_schedule.csv"));
        sch << std::fixed << std::setprecision(6);
        sch << "robot,task,mission,start,end,duration,bucket_start,bucket_end\n";

        for (int i = 0; i < R; ++i) {
            int cur = -1;
            for (int j = 0; j < N; ++j)
                if (cplex.getValue(xDep[i][j]) > 0.5) { cur = j; break; }

            while (cur != -1) {
                const double start = cplex.getValue(t[i][cur]);
                const double end = start + D.s[cur];
                const int    m = D.missionOf[cur];
                const int    bh = (int)std::floor(start / D.Delta);
                const int    be = (int)std::floor((end - 1e-9) / D.Delta);

                sch << i << "," << cur << "," << m << ","
                    << start << "," << end << "," << (end - start) << ","
                    << bh << "," << be << "\n";

                int nxt = -1;
                for (int k = 0; k < N; ++k)
                    if (k != cur && cplex.getValue(xTask[i][cur][k]) > 0.5) { nxt = k; break; }
                cur = nxt;
            }
        }
    }

    // 3) bucket cap vs quota
    {
        std::ofstream bc(path("sol_buckets.csv"));
        bc << "bucket,mission,cap,quota\n";
        for (int h = 0; h < H; ++h)
            for (int m = 0; m < Mmissions; ++m)
                bc << h << "," << m << "," << kpi.cap[m][h] << "," << D.q[m][h] << "\n";
    }

    // 4) reservation/credit matrices
    if (useY) {
        const auto& y = *yPtr;
        std::ofstream yz(path("sol_y.csv"));
        yz << "robot,bucket,mission,value\n";
        for (int i = 0; i < R; ++i)
            for (int h = 0; h < H; ++h)
                for (int m = 0; m < Mmissions; ++m)
                    yz << i << "," << h << "," << m << ","
                    << ((cplex.getValue(y[i][h][m]) > 0.5) ? 1 : 0) << "\n";
    }
    if (useZ) {
        const auto& z = *zPtr;
        std::ofstream yz(path("sol_z.csv"));
        yz << "robot,bucket,mission,value\n";
        for (int i = 0; i < R; ++i)
            for (int h = 0; h < H; ++h)
                for (int m = 0; m < Mmissions; ++m)
                    yz << i << "," << h << "," << m << ","
                    << ((cplex.getValue(z[i][h][m]) > 0.5) ? 1 : 0) << "\n";
    }

    // 5) legs
    {
        std::ofstream lc(path("sol_legs.csv"));
        lc << std::fixed << std::setprecision(6);
        lc << "robot,from,to,from_task,to_task,travel,depart_earliest,arrive_earliest,start_at_to,wait_at_to,switch_mission\n";
        for (const auto& L : legs) {
            lc << L.i << "," << L.from << "," << L.to << ","
                << L.fromTask << "," << L.toTask << ","
                << L.travel << "," << L.departEarliest << "," << L.arriveEarliest << ",";
            if (L.toTask >= 0) lc << L.startAtDest << "," << L.waitAtDest;
            else               lc << "," << 0.0;
            lc << "," << (L.switchMission ? 1 : 0) << "\n";
        }
    }

    // 6) routes
    {
        std::ofstream rcsv(path("sol_routes.csv"));
        rcsv << "robot,seq,from_label,to_label,from_task,to_task,"
            "from_mission,to_mission,from_x,from_y,to_x,to_y,travel\n";
        int lastRobot = -1, seq = 0;
        for (const auto& L : legs) {
            if (L.i != lastRobot) { lastRobot = L.i; seq = 0; }
            rcsv << L.i << "," << seq++ << ","
                << L.from << "," << L.to << ","
                << L.fromTask << "," << L.toTask << ","
                << L.fromMission << "," << L.toMission << ","
                << std::setprecision(10) << L.fromX << "," << L.fromY << ","
                << L.toX << "," << L.toY << ","
                << L.travel << "\n";
        }
    }
}
