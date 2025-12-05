// a2a.cpp
#include "a2a.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

// ---------- helpers ----------
static inline double overlapLen(double a0, double a1, double b0, double b1) {
    double L = std::max(a0, b0), R = std::min(a1, b1);
    return (R > L) ? (R - L) : 0.0;
}

// Build intervals per robot of [start,end) with mission id
std::vector<std::vector<std::tuple<double, double, int>>> buildRobotIntervals(
    const ProblemData& D,
    const IloCplex& cplex,
    const IloArray<IloBoolVarArray>& xDep,
    const IloArray< IloArray<IloBoolVarArray> >& xTask,
    const IloArray<IloBoolVarArray>& /*xEnd*/,
    const IloArray<IloNumVarArray>& t
) {
    const int R = D.R, N = D.N;
    std::vector<std::vector<std::tuple<double, double, int>>> intervals(R);

    for (int i = 0; i < R; ++i) {
        // pick first task after depot (if any)
        int cur = -1;
        for (int j = 0; j < N; ++j) {
            if (cplex.getValue(xDep[i][j]) > 0.5) { cur = j; break; }
        }
        if (cur < 0) continue; // robot unused

        // walk route and stash [start,end) intervals w/ mission id
        while (true) {
            double tj = cplex.getValue(t[i][cur]);
            double ej = tj + D.s[cur];
            int    mj = D.missionOf[cur];
            intervals[i].push_back(std::make_tuple(tj, ej, mj));

            int nxt = -1;
            for (int k = 0; k < N; ++k) {
                if (k != cur && cplex.getValue(xTask[i][cur][k]) > 0.5) { nxt = k; break; }
            }
            if (nxt != -1) { cur = nxt; continue; }
            break;
        }
        std::sort(intervals[i].begin(), intervals[i].end(),
            [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
    }
    return intervals;
}

bool hasActiveBucketAccounting(const ProblemData& D) {
    bool anyQ = false;
    if (!D.q.empty()) {
        for (int m = 0; m < (int)D.q.size() && !anyQ; ++m)
            for (int h = 0; h < (int)D.q[m].size(); ++h)
                if (D.q[m][h] > 0) { anyQ = true; break; }
    }
    bool anyR = false;
    if (!D.Rtot.empty()) for (double v : D.Rtot) if (v > 0.0) { anyR = true; break; }
    return anyQ || anyR;
}

// ---------- main KPI evaluator ----------
EvalKPI evaluateApplesToApples(
    const ProblemData& D, const IloCplex& cplex,
    const IloArray<IloBoolVarArray>& xDep,
    const IloArray< IloArray<IloBoolVarArray> >& xTask,
    const IloArray<IloBoolVarArray>& xEnd,
    const IloArray<IloNumVarArray>& t,
    const IloNumVarArray& Tvars,
    bool credit_spillover,
    const IloArray< IloArray<IloBoolVarArray> >* y_ptr,
    const IloArray< IloArray<IloBoolVarArray> >* z_ptr
) {
    const int R = D.R, N = D.N, M = D.Mmissions, H = D.H;
    const double Delta = D.Delta;

    EvalKPI kpi;
    kpi.rawObj = cplex.getObjValue();
    kpi.T.assign(M, 0.0);
    for (int m = 0; m < M; ++m) kpi.T[m] = cplex.getValue(Tvars[m]);

    // Travel
    double travel = 0.0;
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < N; ++j) travel += D.cDep[i][j] * (cplex.getValue(xDep[i][j]) > 0.5 ? 1.0 : 0.0);
        for (int j = 0; j < N; ++j) for (int k = 0; k < N; ++k) if (j != k)
            travel += D.cTask[j][k] * (cplex.getValue(xTask[i][j][k]) > 0.5 ? 1.0 : 0.0);
        for (int j = 0; j < N; ++j) travel += D.cEnd[i][j] * (cplex.getValue(xEnd[i][j]) > 0.5 ? 1.0 : 0.0);
    }
    kpi.travel = travel;

    // Mission switches
    int switches = 0;
    for (int i = 0; i < R; ++i) for (int j = 0; j < N; ++j) for (int k = 0; k < N; ++k) if (j != k) {
        if (cplex.getValue(xTask[i][j][k]) > 0.5 && D.missionOf[j] != D.missionOf[k]) ++switches;
    }
    kpi.switches = switches;

    // Per-robot intervals
    auto intervals = buildRobotIntervals(D, cplex, xDep, xTask, xEnd, t);

    // Infer per-bucket capacity cap[m][h]
    kpi.cap.assign(D.Mmissions, std::vector<int>(D.H, 0));

    if (y_ptr && y_ptr->getSize() == D.R) {
        // reservations in OURS
        const auto& yy = *y_ptr;
        for (int i = 0; i < D.R; ++i)
            for (int h = 0; h < D.H; ++h)
                for (int m = 0; m < D.Mmissions; ++m) {
                    const IloBoolVar& v = yy[i][h][m];
                    if (cplex.isExtracted(v) && cplex.getValue(v) > 0.5) kpi.cap[m][h] += 1;
                }
    }
    else if (z_ptr && z_ptr->getSize() == D.R) {
        // starts-per-mission in soft-mono baseline
        const auto& zz = *z_ptr;
        for (int i = 0; i < D.R; ++i)
            for (int h = 0; h < D.H; ++h)
                for (int m = 0; m < D.Mmissions; ++m) {
                    const IloBoolVar& v = zz[i][h][m];
                    if (cplex.isExtracted(v) && cplex.getValue(v) > 0.5) kpi.cap[m][h] += 1;
                }
    }
    else {
        // Fallback: infer from executed work
        for (int i = 0; i < D.R; ++i) {
            for (int h = 0; h < D.H; ++h) {
                double B0 = h * D.Delta, B1 = (h + 1) * D.Delta;
                double bestOv = 0.0; int bestM = -1; double bestStart = std::numeric_limits<double>::infinity();
                for (const auto& iv : intervals[i]) {
                    double s = std::get<0>(iv), e = std::get<1>(iv); int m = std::get<2>(iv);
                    double ov = credit_spillover ? overlapLen(s, e, B0, B1) : 0.0;
                    if (ov > bestOv + 1e-9 || (std::abs(ov - bestOv) <= 1e-9 && s < bestStart)) {
                        bestOv = ov; bestM = m; bestStart = s;
                    }
                }
                if (bestM == -1) {
                    double earliest = std::numeric_limits<double>::infinity(); int em = -1;
                    for (const auto& iv : intervals[i]) {
                        double s = std::get<0>(iv); int m = std::get<2>(iv);
                        if (s >= B0 - 1e-9 && s < B1 - 1e-9) { if (s < earliest) { earliest = s; em = m; } }
                    }
                    bestM = em;
                }
                if (bestM >= 0) kpi.cap[bestM][h] += 1;
            }
        }
    }

    // Quota & rolling shortfalls
    double shortQ = 0.0, shortR = 0.0;
    if (!D.q.empty()) {
        for (int m = 0; m < M; ++m) for (int h = 0; h < H; ++h) {
            int req = D.q[m][h];
            int got = kpi.cap[m][h];
            if (req > got) shortQ += (req - got);
        }
    }
    if (!D.Rtot.empty()) {
        for (int m = 0; m < M; ++m) {
            double totMin = 0.0;
            for (int h = 0; h < H; ++h) totMin += Delta * kpi.cap[m][h];
            if (D.Rtot[m] > totMin) shortR += (D.Rtot[m] - totMin);
        }
    }
    kpi.shortQ = shortQ;
    kpi.shortR = shortR;

    // Compose apples-to-apples score
    double score = 0.0;
    for (int m = 0; m < M; ++m) score += D.w[m] * kpi.T[m];
    score += D.lambdaTravel * kpi.travel;
    if (D.etaSwitch > 1e-12 && D.switchPenalty > 1e-12)
        score += D.etaSwitch * D.switchPenalty * static_cast<double>(kpi.switches);
    if (D.rhoQ > 1e-12) score += D.rhoQ * kpi.shortQ;
    if (D.rhoR > 1e-12) score += D.rhoR * kpi.shortR;

    kpi.score = score;
    return kpi;
}
