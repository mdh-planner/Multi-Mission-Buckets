// build.cpp (replacement)
// Enforces compatibility via setUB(0) at var creation;
// uses big-M timing links + tightening cuts for stronger LP bounds.

#include "build.h"
#include "a2a.h"     // hasActiveBucketAccounting if you need the flag
#include <iomanip>
#include <algorithm>
ILOSTLBEGIN

// ----------------- helper: routing + timing + buckets (no indicators) -----------------
static void addRoutingTimingAndBuckets(
    const ProblemData& D, const SolveOptions& opt,
    int R, int N, int M, int H, double Delta, double eps,
    IloModel& model,
    const IloArray<IloBoolVarArray>& xDep,
    const IloArray<IloBoolVarArray>& xEnd,
    const IloArray<IloArray<IloBoolVarArray>>& xTask,
    const IloArray<IloNumVarArray>& t,
    const IloNumVarArray& T,
    const IloNumVarArray& rEnd,
    // optional bucket vars (can be null if feature off)
    const IloArray<IloArray<IloBoolVarArray>>* y,
    const IloArray<IloArray<IloBoolVarArray>>* b,
    const IloArray<IloArray<IloBoolVarArray>>* z,
    IloArray<IloNumVarArray>* sQ,
    IloNumVarArray* sR
    
) {
    auto env = model.getEnv();

    // ---------- visit each task exactly once ----------
    for (int j = 0; j < N; ++j) {
        IloExpr pred(env), succ(env);
        for (int i = 0; i < R; ++i) {
            pred += xDep[i][j];
            succ += xEnd[i][j];
            for (int h2 = 0; h2 < N; ++h2) pred += xTask[i][h2][j];
            for (int k = 0; k < N; ++k)  succ += xTask[i][j][k];
        }
        model.add(pred == 1); pred.end();
        model.add(succ == 1); succ.end();
    }

    // ---------- flow conservation & single route per robot ----------
    for (int i = 0; i < R; ++i) {
        IloExpr dep(env), fin(env);
        for (int j = 0; j < N; ++j) { dep += xDep[i][j]; fin += xEnd[i][j]; }
        model.add(dep == fin);
        model.add(dep <= 1);
        dep.end(); fin.end();

        for (int j = 0; j < N; ++j) {
            IloExpr in(env), out(env);
            in += xDep[i][j];  for (int h2 = 0; h2 < N; ++h2) in += xTask[i][h2][j];
            out += xEnd[i][j];  for (int k = 0; k < N; ++k)  out += xTask[i][j][k];
            model.add(in == out); in.end(); out.end();
        }
    }

    // ---------- timing: big-M links (strong LP), + tightening cuts ----------
    const double horizonUB = H * Delta;

    auto M_dep = [&](int i, int j) {
        double v = horizonUB - D.cDep[i][j];
        return (v > 0.0 ? v : 0.0);
    };
    auto M_arc = [&](int j, int k) {
        double v = horizonUB - (D.s[j] + D.cTask[j][k]);
        return (v > 0.0 ? v : 0.0);
    };
    auto M_miss = [&]() {
        double maxS = 0.0; for (double sj : D.s) if (sj > maxS) maxS = sj;
        return horizonUB + maxS;
    };

    // t[i][j] >= cDep[i][j] - M*(1 - xDep[i][j])
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < N; ++j)
            model.add(t[i][j] >= D.cDep[i][j] - M_dep(i, j) * (1 - xDep[i][j]));

    // t[i][k] >= t[i][j] + s[j] + c[j][k] - M*(1 - xTask[i][j][k])
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) if (j != k)
                model.add(t[i][k] >= t[i][j] + D.s[j] + D.cTask[j][k]
                    - M_arc(j, k) * (1 - xTask[i][j][k]));

    // T[m] >= t[i][j] + s[j] - M*(1 - act_{i,j}), where act_{i,j} = xDep[i][j] + sum_u xTask[i][u][j]
    for (int m = 0; m < M; ++m) {
        for (int idx = 0; idx < (int)D.tasksOfMission[m].size(); ++idx) {
            int j = D.tasksOfMission[m][idx];
            for (int i = 0; i < R; ++i) {
                IloExpr act(env);
                act += xDep[i][j];
                for (int u = 0; u < N; ++u) if (u != j) act += xTask[i][u][j];
                model.add(T[m] >= t[i][j] + D.s[j] - M_miss() * (1 - act));
                act.end();
            }
        }
    }

    // Tightening: if task j not used by robot i, keep t[i][j] small
    // t[i][j] <= horizonUB * (xDep[i][j] + sum_u xTask[i][u][j])
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < N; ++j) {
            IloExpr used(env);
            used += xDep[i][j];
            for (int u = 0; u < N; ++u) if (u != j) used += xTask[i][u][j];
            model.add(t[i][j] <= horizonUB * used);
            used.end();
        }

    // ---------- robot end-time constraints (return-to-depot) ----------
    {
        auto M_end = [&](int i, int j) {
            double v = horizonUB - (D.s[j] + D.cEnd[i][j]);
            return (v > 0.0 ? v : 0.0);
        };

        for (int i = 0; i < R; ++i) {
            for (int j = 0; j < N; ++j) {
                // If xEnd[i][j] = 1, robot i finishes at t[i][j] + s[j] + cEnd[i][j].
                // Otherwise, constraint is relaxed by M_end.
                model.add(
                    rEnd[i] >=
                    t[i][j] + D.s[j] + D.cEnd[i][j]
                    - M_end(i, j) * (1 - xEnd[i][j])
                );
            }
        }
    }

    // ---------- buckets ----------
    if (!opt.useBucketLayer) return;

    // F1: one mission per (i,h) if reservations
    if (opt.useReservations && y && opt.enforceOneMissionPerBucket) {
        for (int i = 0; i < R; ++i) for (int h = 0; h < H; ++h) {
            IloExpr sumM(env); for (int m = 0; m < M; ++m) sumM += (*y)[i][h][m];
            model.add(sumM <= 1); sumM.end();
        }
    }

    // F2: ≤K starts per (i,h)
    if (opt.enforceUpperLimitTaskStartPerBucket && b) {
        for (int i = 0; i < R; ++i) for (int h = 0; h < H; ++h) {
            IloExpr sumB(env); for (int j = 0; j < N; ++j) sumB += (*b)[i][j][h];
            model.add(sumB <= opt.K); sumB.end();
        }
    }

    // F3: start-bucket equals “i does j”
    if (b) {
        for (int i = 0; i < R; ++i) for (int j = 0; j < N; ++j) {
            IloExpr lhs(env), rhs(env);
            for (int h = 0; h < H; ++h) lhs += (*b)[i][j][h];
            rhs += xDep[i][j]; for (int u = 0; u < N; ++u) rhs += xTask[i][u][j];
            model.add(lhs == rhs); lhs.end(); rhs.end();
        }
    }

    // F4: respect reservation mission
    if (opt.useReservations && y && b) {
        for (int i = 0; i < R; ++i) for (int j = 0; j < N; ++j) {
            int m = D.missionOf[j];
            for (int h = 0; h < H; ++h) model.add((*b)[i][j][h] <= (*y)[i][h][m]);
        }
    }

 // F4': if y[i,h,m] = 1, then robot i must start at least one task of mission m
// in bucket h:  sum_{j: missionOf[j] = m} b[i][j][h] >= y[i][h][m]
    //if (opt.useReservations && y && b) {
    //    for (int i = 0; i < R; ++i) {
    //        for (int h = 0; h < H; ++h) {
    //            for (int m = 0; m < M; ++m) {
    //                IloExpr sumB(env);
    //                for (int j = 0; j < N; ++j) {
    //                    if (D.missionOf[j] == m) {
    //                        sumB += (*b)[i][j][h];
    //                    }
    //                }
    //                // (*y)[i][h][m] is a single IloBoolVar
    //                model.add(sumB >= (*y)[i][h][m]);
    //                sumB.end();
    //            }
    //        }
    //    }
    //}


    //// F4b: reservation excludes other missions within (i,h)
    //if (opt.useReservations && y && b) {
    //    const int K_M = std::max(1, opt.K);
    //    for (int i = 0; i < R; ++i)
    //        for (int h = 0; h < H; ++h)
    //            for (int m = 0; m < M; ++m) {
    //                IloExpr other(env);
    //                for (int j = 0; j < N; ++j)
    //                    if (D.missionOf[j] != m) other += (*b)[i][j][h];
    //                model.add(other <= K_M * (1 - (*y)[i][h][m]));
    //                other.end();
    //            }
    //}

    // F5: pin start time within bucket window when selected
    if (b) {
        for (int i = 0; i < R; ++i) for (int j = 0; j < N; ++j) for (int h = 0; h < H; ++h) {
            model.add(IloIfThen(env, ((*b)[i][j][h] == 1), (t[i][j] >= (h * Delta))));
            model.add(IloIfThen(env, ((*b)[i][j][h] == 1), (t[i][j] + D.s[j] <= (h + 1) * Delta - eps)));
        }
    }

    // Z: soft-mono crediting (no reservations)
    if (!opt.useReservations && z && b) {
        for (int i = 0; i < R; ++i)
            for (int h = 0; h < H; ++h)
                for (int m = 0; m < M; ++m) {
                    for (int j = 0; j < N; ++j)
                        if (D.missionOf[j] == m) model.add((*z)[i][h][m] >= (*b)[i][j][h]);
                    IloExpr sumB(env);
                    for (int j = 0; j < N; ++j)
                        if (D.missionOf[j] == m) sumB += (*b)[i][j][h];
                    model.add((*z)[i][h][m] <= sumB);
                    sumB.end();
                }
        // at most one mission per (i,h)
        for (int i = 0; i < R; ++i) for (int h = 0; h < H; ++h) {
            IloExpr sumZ(env); for (int m = 0; m < M; ++m) sumZ += (*z)[i][h][m];
            model.add(sumZ <= 1); sumZ.end();
        }
    }

    // F6: per-bucket minima (quota) and F7: rolling totals
    if (opt.useQuotas) {
        for (int m = 0; m < M; ++m) {
            for (int h = 0; h < H; ++h) {
                int req = (D.q.empty() ? 0 : D.q[m][h]);
                /*if (req <= 0)  continue;

                IloExpr cap(env);
                
                if (opt.useReservations && y) {
                    for (int i = 0; i < R; ++i) {
                        cap += (*y)[i][h][m];
                    }
                }
                else if (z) {
                    for (int i = 0; i < R; ++i) {
                        cap += (*z)[i][h][m];
                    }
                }
                
                if (D.rhoQ > 0.0) model.add(cap + (*sQ)[m][h] >= req);
                else              model.add(cap >= req);
                cap.end();*/

                // Hard quotas:
                if (req > 0) {
                    IloExpr cap(env);
                    for (int i = 0; i < R; ++i) cap += (*y)[i][h][m];
                    model.add(cap >= req);   // no + sQ
                    cap.end();
                }

            }
        }
        if (!D.Rtot.empty()) {
            const int Ms = std::min(M, (int)D.Rtot.size());
            for (int m = 0; m < Ms; ++m) if (D.Rtot[m] > 0.0) {
                IloExpr tot(env);
                if (opt.useReservations && y)
                    for (int i = 0; i < R; ++i) for (int h = 0; h < H; ++h) tot += Delta * (*y)[i][h][m];
                else if (z)
                    for (int i = 0; i < R; ++i) for (int h = 0; h < H; ++h) tot += Delta * (*z)[i][h][m];
                if (D.rhoR > 0.0) model.add(tot + (*sR)[m] >= D.Rtot[m]);
                else              model.add(tot >= D.Rtot[m]);
                tot.end();
            }
        }
    }
}

// ----------------- buildModel (alloc vars with compat setUB(0); objective; MTZ opt) -----------------
BuildArtifacts buildModel(const ProblemData& D, const SolveOptions& opt, IloEnv& env)
{
    const int R = D.R, N = D.N, M = D.Mmissions, H = D.H;
    const double Delta = D.Delta, eps = D.eps;

    BuildArtifacts A(env, R, N, M, H);
    A.useY = (opt.useBucketLayer && opt.useQuotas && opt.useReservations);
    A.useZ = (opt.useBucketLayer && opt.useQuotas && !opt.useReservations);

    // ---------- allocate variables ----------
    // buckets
    if (opt.useBucketLayer) {
        if (opt.useReservations) {
            for (int i = 0; i < R; ++i) {
                A.y[i] = IloArray<IloBoolVarArray>(env, H);
                for (int h = 0; h < H; ++h) {
                    A.y[i][h] = IloBoolVarArray(env, M);
                    for (int m = 0; m < M; ++m) A.y[i][h][m] = IloBoolVar(env);
                }
            }
        }
        for (int i = 0; i < R; ++i) {
            A.b[i] = IloArray<IloBoolVarArray>(env, N);
            for (int j = 0; j < N; ++j) {
                A.b[i][j] = IloBoolVarArray(env, H);
                for (int h = 0; h < H; ++h) {
                    A.b[i][j][h] = IloBoolVar(env);
                    // Incompatible robot–task => cannot start in any bucket
                    if (!D.compat[i][j]) A.b[i][j][h].setUB(0);
                }
            }
        }
        if (!opt.useReservations) {
            for (int i = 0; i < R; ++i) {
                A.z[i] = IloArray<IloBoolVarArray>(env, H);
                for (int h = 0; h < H; ++h) {
                    A.z[i][h] = IloBoolVarArray(env, M);
                    for (int m = 0; m < M; ++m) A.z[i][h][m] = IloBoolVar(env);
                }
            }
        }
        if (opt.useQuotas) {
            for (int m = 0; m < M; ++m) A.sQ[m] = IloNumVarArray(env, H, 0.0, IloInfinity, ILOFLOAT);
            // sR constructed in BuildArtifacts ctor
        }
    }

    // routing & times
    const double horizonUB = H * Delta;
    double maxS = 0.0; for (double sj : D.s) if (sj > maxS) maxS = sj;

    for (int i = 0; i < R; ++i) {
        A.xDep[i] = IloBoolVarArray(env, N);
        A.xEnd[i] = IloBoolVarArray(env, N);
        A.t[i] = IloNumVarArray(env, N, 0.0, IloInfinity, ILOFLOAT);

        for (int j = 0; j < N; ++j) {
            A.xDep[i][j] = IloBoolVar(env);
            A.xEnd[i][j] = IloBoolVar(env);
            A.xDep[i][j].setName(("xDep_" + std::to_string(i) + "_" + std::to_string(j)).c_str());
            A.xEnd[i][j].setName(("xEnd_" + std::to_string(i) + "_" + std::to_string(j)).c_str());

            A.t[i][j].setName(("t_" + std::to_string(i) + "_" + std::to_string(j)).c_str());
            A.t[i][j].setUB(horizonUB);

            // Incompatible robot–task: forbid start/end by fixing to 0
            if (!D.compat[i][j]) {
                A.xDep[i][j].setUB(0);
                A.xEnd[i][j].setUB(0);
            }
        }

        A.xTask[i] = IloArray<IloBoolVarArray>(env, N);
        for (int j = 0; j < N; ++j) {
            A.xTask[i][j] = IloBoolVarArray(env, N);
            for (int k = 0; k < N; ++k) {
                A.xTask[i][j][k] = IloBoolVar(env);
                // No self-loops
                if (j == k) A.xTask[i][j][k].setUB(0);
                // Any arc touching incompatible endpoint is forbidden
                if (!D.compat[i][j] || !D.compat[i][k]) A.xTask[i][j][k].setUB(0);
            }
        }
    }

    for (int m = 0; m < M; ++m) {
        A.T[m] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
        A.T[m].setName(("T_" + std::to_string(m)).c_str());
        A.T[m].setUB(horizonUB + maxS);
    }

    // Allocate robot end-time variables
    A.rEnd = IloNumVarArray(env, R, 0.0, IloInfinity, ILOFLOAT);
    for (int i = 0; i < R; ++i) {
        A.rEnd[i].setName(("rEnd_" + std::to_string(i)).c_str());
        A.rEnd[i].setUB(horizonUB + maxS); // safe upper bound
    }
    // ---------- constraints ----------
    addRoutingTimingAndBuckets(
        D, opt, R, N, M, H, Delta, eps,
        A.model, A.xDep, A.xEnd, A.xTask, A.t, A.T, A.rEnd,
        (opt.useBucketLayer && opt.useReservations) ? &A.y : nullptr,
        (opt.useBucketLayer) ? &A.b : nullptr,
        (opt.useBucketLayer && !opt.useReservations) ? &A.z : nullptr,
        (opt.useBucketLayer && opt.useQuotas) ? &A.sQ : nullptr,
        (opt.useBucketLayer && opt.useQuotas) ? &A.sR : nullptr
    );

    // --- (optional) MTZ strengthening when buckets OFF ---
    if (!opt.useBucketLayer) {
        IloArray<IloNumVarArray> u(env, R);
        for (int i = 0; i < R; ++i) u[i] = IloNumVarArray(env, N, 0.0, N, ILOFLOAT);
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k) if (j != k) {
                    // u_k >= u_j + 1 - N*(1 - x_ijk)
                    A.model.add(u[i][k] >= u[i][j] + 1 - N * (1 - A.xTask[i][j][k]));
                }
    }
    // Create the auxiliary variable for the makespan
    IloNumVar Tmax(env, 0.0, IloInfinity, ILOFLOAT, "Tmax");

    // Tmax >= T[m] for all missions m
    for (int m = 0; m < M; ++m) {
        A.model.add(Tmax >= A.T[m]);
    }

    // ---------- objective ----------
    // Objective: minimize Tmax (i.e., minimize max_m T[m])
    IloObjective obj = IloMinimize(env, Tmax);

    /*IloExpr obj(env);
    for (int m = 0; m < M; ++m) {
        
        obj += A.T[m];
    }*/
    

    /*for (int i = 0; i < R; ++i)
        obj += A.rEnd[i];*/

    // small travel weight helps LB (optional; only if provided)
  /*  if (D.lambdaTravel > 0.0) {
        for (int i = 0; i < R; ++i) {
            for (int j = 0; j < N; ++j) obj += D.lambdaTravel * D.cDep[i][j] * A.xDep[i][j];
            for (int j = 0; j < N; ++j) for (int k = 0; k < N; ++k) if (j != k)
                obj += D.lambdaTravel * D.cTask[j][k] * A.xTask[i][j][k];
            for (int j = 0; j < N; ++j) obj += D.lambdaTravel * D.cEnd[i][j] * A.xEnd[i][j];
        }
    }*/

 /*   if (D.etaSwitch > 0.0 && D.switchPenalty > 0.0) {
        for (int i = 0; i < R; ++i) for (int j = 0; j < N; ++j) for (int k = 0; k < N; ++k) if (j != k) {
            double psi = (D.missionOf[j] == D.missionOf[k]) ? 0.0 : D.switchPenalty;
            if (psi > 0.0) obj += D.etaSwitch * psi * A.xTask[i][j][k];
        }
    }*/

    A.model.add(IloMinimize(env, obj));
    obj.end();

    return A;
}
