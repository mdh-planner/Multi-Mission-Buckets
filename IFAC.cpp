// IFAC.cpp
// Multi-Mission Optimization with time-bucket reservations
// Indicator-constraint timing (no big-M).

// Notes:
//  - Indicator constraints: IloIfThen(env, (bin == 1), linear_constraint).
//  - Routing uses standard flow conservation; indicators carry all temporal logic.
//  - Subtours die automatically: along any cycle the chain of ≥ increments
//    forces an infeasible time increase.

#include <ilcplex/ilocplex.h>
#include "definitions.h"
#include "a2a.h"
#include "report.h"
#include "build.h"
#include "lazy_callback.h"


#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <tuple>
#include <random>
#include <numeric>
#include <cassert>
#include <limits>
#include <filesystem>
namespace fs = std::filesystem;

ILOSTLBEGIN

//===============================
// Instance generator
//===============================

static inline double euclid(const Pt& a, const Pt& b) {
	double dx = a.x - b.x, dy = a.y - b.y;
	return std::sqrt(dx * dx + dy * dy);
}

// Evenly spaced mission centers on a circle
static std::vector<Pt> makeMissionCenters(int M, double radius) {
	std::vector<Pt> c(M);
	const double TwoPi = 6.283185307179586;
	for (int m = 0; m < M; ++m) {
		double ang = TwoPi * (double)m / std::max(1, M);
		c[m] = { radius * std::cos(ang), radius * std::sin(ang) };
	}
	return c;
}

// Place depots as requested
static std::vector<Pt> makeDepots(const GeneratorConfig& cfg,
	const std::vector<Pt>& centers,
	std::mt19937& rng) {
	std::vector<Pt> dep(cfg.R);
	std::uniform_real_distribution<double> U(0.0, 1.0);
	const double TwoPi = 6.283185307179586;

	switch (cfg.depot_mode) {
	case GeneratorConfig::DepotMode::NEAR_MISSIONS: {
		for (int i = 0; i < cfg.R; ++i) {
			int m = i % std::max(1, (int)centers.size());
			dep[i] = { centers[m].x * 0.8, centers[m].y * 0.8 };
		}
	} break;
	case GeneratorConfig::DepotMode::CENTRAL: {
		for (int i = 0; i < cfg.R; ++i) dep[i] = { 0.0, 0.0 };
	} break;
	case GeneratorConfig::DepotMode::RANDOM_RING: {
		for (int i = 0; i < cfg.R; ++i) {
			double ang = TwoPi * U(rng);
			dep[i] = { cfg.area_radius * std::cos(ang), cfg.area_radius * std::sin(ang) };
		}
	} break;
	}
	return dep;
}

void generateInstance(ProblemData& D, GeneratorConfig cfg) {
	// Basic sanity + defaults
	if (cfg.M <= 0) cfg.M = 1;
	if (cfg.R <= 0) cfg.R = 1;
	if (cfg.N < cfg.M) cfg.N = cfg.M;                 // ensure at least one task per mission
	if (cfg.flip_bucket < 0) cfg.flip_bucket = std::max(0, cfg.H / 2);
	if (cfg.w.empty()) cfg.w.assign(cfg.M, 1.0);
	if (cfg.Rtot.empty()) cfg.Rtot.assign(cfg.M, 0.0);
	else if ((int)cfg.Rtot.size() != cfg.M) cfg.Rtot.assign(cfg.M, 0.0);
	cfg.surge_robots = std::min(cfg.surge_robots, cfg.R);

	std::mt19937 rng(cfg.seed);
	std::uniform_int_distribution<int> Uservice(cfg.s_min, cfg.s_max);
	std::normal_distribution<double>   N0(0.0, cfg.cluster_sigma);
	std::uniform_real_distribution<double> U01(0.0, 1.0);

	// Fill ProblemData sizes/params
	D.R = cfg.R;
	D.N = cfg.N;
	D.Mmissions = cfg.M;
	D.H = cfg.H;
	D.Delta = cfg.Delta;
	D.eps = cfg.eps;

	D.w = cfg.w;
	D.lambdaTravel = cfg.lambdaTravel;
	D.etaSwitch = cfg.etaSwitch;
	D.switchPenalty = cfg.switchPenalty;

	// Assign tasks to missions ~evenly (shuffle after)
	D.missionOf.assign(D.N, 0);
	for (int j = 0; j < D.N; ++j) D.missionOf[j] = j % D.Mmissions;
	std::shuffle(D.missionOf.begin(), D.missionOf.end(), rng);
	// Ensure every mission has at least one task
	for (int m = 0; m < D.Mmissions && m < D.N; ++m) D.missionOf[m] = m;

	// Build index per mission
	D.tasksOfMission.assign(D.Mmissions, {});
	for (int j = 0; j < D.N; ++j) D.tasksOfMission[D.missionOf[j]].push_back(j);

	// Service times
	D.s.resize(D.N);
	for (int j = 0; j < D.N; ++j) D.s[j] = (double)Uservice(rng);

	//// Geometry: task locations and depots
	//auto centers = makeMissionCenters(D.Mmissions, cfg.area_radius);
	//std::vector<Pt> taskXY(D.N);
	//for (int j = 0; j < D.N; ++j) {
	//    int m = D.missionOf[j];
	//    taskXY[j] = { centers[m].x + N0(rng), centers[m].y + N0(rng) };
	//}
	//auto depotXY = makeDepots(cfg, centers, rng);

	// ---- Geometry: mission centers, tasks around centers, depots ----
	auto centers = makeMissionCenters(D.Mmissions, cfg.area_radius);
	std::vector<Pt> taskXY(D.N);

	for (int j = 0; j < D.N; ++j) {
		int m = D.missionOf[j];
		taskXY[j] = { centers[m].x + N0(rng), centers[m].y + N0(rng) };
	}
	auto depotXY = makeDepots(cfg, centers, rng);

	// keep copies in the ProblemData for later export/plotting
	D.taskXY = taskXY;
	D.depotXY = depotXY;


	// Distance -> time conversion
	auto timeFromDist = [&](double dist, int m_from, int m_to) {
		double t = cfg.base_travel + dist / std::max(1e-6, cfg.speed);
		if (cfg.cross_mission_penalty > 1e-12 && m_from != m_to) t += cfg.cross_mission_penalty;
		return t;
	};

	// Task->Task travel
	D.cTask.assign(D.N, std::vector<double>(D.N, 0.0));
	for (int j = 0; j < D.N; ++j) {
		for (int k = 0; k < D.N; ++k) {
			if (j == k) { D.cTask[j][k] = 0.0; continue; }
			double d = euclid(taskXY[j], taskXY[k]);
			D.cTask[j][k] = timeFromDist(d, D.missionOf[j], D.missionOf[k]);
		}
	}
	// Depot->Task travel
	D.cDep.assign(D.R, std::vector<double>(D.N, 0.0));
	for (int i = 0; i < D.R; ++i) {
		for (int j = 0; j < D.N; ++j) {
			double d = euclid(depotXY[i], taskXY[j]);
			D.cDep[i][j] = timeFromDist(d, -1, D.missionOf[j]);
		}
	}
	// Task->End (return) travel
	D.cEnd.assign(D.R, std::vector<double>(D.N, 0.0));
	for (int i = 0; i < D.R; ++i) {
		for (int j = 0; j < D.N; ++j) {
			double d = euclid(taskXY[j], depotXY[i]);
			D.cEnd[i][j] = timeFromDist(d, D.missionOf[j], -1);
		}
	}

	// Compatibility (random zeros; keep at least one compatible robot per task)
	D.compat.assign(D.R, std::vector<int>(D.N, 1));
	if (cfg.incompat_prob > 0.0) {
		for (int i = 0; i < D.R; ++i) for (int j = 0; j < D.N; ++j) {
			if (U01(rng) < cfg.incompat_prob) D.compat[i][j] = 0;
		}
	}
	if (cfg.ensure_feasible_compat) {
		for (int j = 0; j < D.N; ++j) {
			bool ok = false;
			for (int i = 0; i < D.R; ++i) if (D.compat[i][j]) { ok = true; break; }
			if (!ok) D.compat[j % D.R][j] = 1;
		}
	}

	// Quotas q[m][h]
	D.q.assign(D.Mmissions, std::vector<int>(D.H, 0));
	if (!cfg.manual_q.empty()) {
		// Manual quotas overwrite everything else
		D.q = cfg.manual_q;
		if ((int)D.q.size() != D.Mmissions) D.q.assign(D.Mmissions, std::vector<int>(D.H, 0));
		for (int m = 0; m < D.Mmissions; ++m) {
			if ((int)D.q[m].size() != D.H) D.q[m].assign(D.H, 0);
			for (int h = 0; h < D.H; ++h) D.q[m][h] = std::min(D.q[m][h], D.R);
		}
	}
	else {
		// Flip pattern + optional surge
		int flip = std::min(std::max(0, cfg.flip_bucket), std::max(0, D.H - 1));
		for (int h = 0; h < D.H; ++h) {
			bool early = (h < flip);
			int pri = early ? cfg.early_mission : cfg.late_mission;
			int oth = early ? cfg.late_mission : cfg.early_mission;
			if (pri >= 0 && pri < D.Mmissions) {
				// BUGFIX: late phase now uses cfg.late_robots (was always early_robots)
				D.q[pri][h] = std::min(early ? cfg.early_robots : cfg.late_robots, D.R);
			}
			if (oth >= 0 && oth < D.Mmissions) D.q[oth][h] = std::min(cfg.other_min, D.R);
		}
		if (cfg.use_surge && cfg.surge_bucket >= 0 && cfg.surge_bucket < D.H &&
			cfg.surge_mission >= 0 && cfg.surge_mission < D.Mmissions) {
			D.q[cfg.surge_mission][cfg.surge_bucket] = std::min(cfg.surge_robots, D.R);
			// Optionally zero others in surge bucket:
			for (int m = 0; m < D.Mmissions; ++m) if (m != cfg.surge_mission) D.q[m][cfg.surge_bucket] = 0;
		}
	}

	// Softness & rolling totals
	D.rhoQ = cfg.rhoQ;
	D.Rtot = cfg.Rtot;   // size M, zeros by default
	D.rhoR = cfg.rhoR;
}

// Quota clamping that respects the chosen baseline.
// - Always enforce per-bucket feasibility: Σ_m q[m][h] ≤ R.
// - For baselines where credits come from *work* (not reservations), limit each
//   mission's total ask across buckets by what its tasks could actually supply.
//     * Starts model: cap_m = #tasks in mission m
//     * Presence model: cap_m = Σ_j∈m ceil(s_j / Δ)
//   Also bound by physical max R*H.
//
// For OURS (reservations), only the per-bucket clamp is applied by default.
// If you *want* reservations limited by available work-time, set forcePresenceCap=true.

static int sumRow(const ProblemData& D, int m) {
	int s = 0; for (int h = 0; h < D.H; ++h) s += D.q[m][h]; return s;
}

template <typename T>
static void ensureSize(std::vector<T>& v, int M, const T& val) {
	if ((int)v.size() != M) v.assign(M, val);
}

static void configureBaseline(ProblemData& D, SolveOptions& opt) {
	switch (opt.baseline) {
	case Baseline::OURS:
		opt.useBucketLayer = true;
		opt.useQuotas = true;
		opt.useReservations = true;
		opt.enforceOneMissionPerBucket = true;

		opt.enforceOneStartPerBucket = false; // enforce one task start per bucket
		opt.enforceUpperLimitTaskStartPerBucket = false; // enforce upper limit of number of tasks per bucket defined by opt.K
		opt.K = 2;  // set the max number of starts per (robot, bucket)
		break;

	case Baseline::WEIGHTS_ONLY: // B1
		opt.useBucketLayer = false;
		opt.useQuotas = false;
		if (!opt.alphaWeights.empty()) {
			// If caller provided weights, use them (resize to M if needed)
			ensureSize(opt.alphaWeights, D.Mmissions, 1.0);
			D.w = opt.alphaWeights;
		}
		break;

	case Baseline::MONOLITHIC:   // B5
		opt.useBucketLayer = false;
		opt.useQuotas = false;
		opt.useReservations = false;
		opt.enforceOneMissionPerBucket = false;

		ensureSize(D.w, D.Mmissions, 1.0);
		break;

	case Baseline::DECOUPLED:    // B7
		opt.useBucketLayer = false;
		opt.useQuotas = false;
		if (opt.robotHomeMission.empty()) {
			opt.robotHomeMission.resize(D.R);
			for (int i = 0; i < D.R; ++i) opt.robotHomeMission[i] = i % std::max(1, D.Mmissions);
		}
		// enforce partition via compatibility mask
		for (int i = 0; i < D.R; ++i) for (int j = 0; j < D.N; ++j) {
			if (D.missionOf[j] != opt.robotHomeMission[i]) D.compat[i][j] = 0;
		}
		break;

	case Baseline::MONO_SOFT_BUCKETS:
		opt.useBucketLayer = true;
		opt.useQuotas = true;
		opt.useReservations = false;             // track via z instead of y
		opt.enforceOneMissionPerBucket = false;  // no F1
		opt.enforceOneStartPerBucket = false;    // multiple starts in a bucket allowed
		ensureSize(D.w, D.Mmissions, 1.0);
		break;
	}
}

static void clampQuotasAllBaselines(ProblemData& D,
	const SolveOptions& opt,
	bool creditFromPresence,
	bool forcePresenceCap,
	bool verbose)
{
	if (D.q.empty() || D.Mmissions == 0 || D.H == 0) return;

	// --- 0) Per-bucket feasibility: Σ_m q[m][h] ≤ R ---
	for (int h = 0; h < D.H; ++h) {
		int tot = 0; for (int m = 0; m < D.Mmissions; ++m) tot += D.q[m][h];
		if (tot <= D.R) continue;

		if (verbose) {
			std::cout << endl;
			std::cerr << "[clamp] bucket " << h << ": sum_m q[m][h]=" << tot
				<< " > R=" << D.R << " -> scaling this column\n";

		}

		// Proportional integer scaling of column h
		struct Frac { double frac; int m; };
		double sf = static_cast<double>(D.R) / std::max(1, tot);
		std::vector<int> newCol(D.Mmissions, 0);
		std::vector<Frac> order; order.reserve(D.Mmissions);
		int newSum = 0;

		for (int m = 0; m < D.Mmissions; ++m) {
			double v = sf * D.q[m][h];
			int base = (int)std::floor(v + 1e-12);
			base = std::max(0, std::min(base, D.R));
			newCol[m] = base; newSum += base;
			order.push_back({ v - base, m });
		}
		int rem = D.R - newSum;
		std::sort(order.begin(), order.end(),
			[](const Frac& a, const Frac& b) { return a.frac > b.frac; });
		for (int k = 0; k < rem && k < (int)order.size(); ++k) {
			int m = order[k].m; if (newCol[m] < D.R) newCol[m] += 1;
		}
		for (int m = 0; m < D.Mmissions; ++m) D.q[m][h] = newCol[m];
	}

	// Decide if we apply a per-mission cap in this baseline
	const bool usingReservations = (opt.useBucketLayer && opt.useQuotas && opt.useReservations);
	const bool applyMissionCaps = (!usingReservations) || forcePresenceCap;

	if (!applyMissionCaps) return; // OURS default: mission totals not tied to #tasks

	// --- 1) Per-mission total cap depending on credit model ---
	std::vector<int> cap(D.Mmissions, 0);
	for (int m = 0; m < D.Mmissions; ++m) {
		int c = 0;
		if (!creditFromPresence) {
			// Starts: each task contributes at most once across all buckets
			c = (int)D.tasksOfMission[m].size();
		}
		else {
			// Presence: long tasks can cover multiple buckets
			for (int j : D.tasksOfMission[m]) {
				int buckets = std::max(1, (int)std::ceil(D.s[j] / D.Delta));
				c += buckets;
			}
		}
		// Can never exceed physical R*H
		c = std::min(c, D.R * D.H);
		cap[m] = c;
	}

	for (int m = 0; m < D.Mmissions; ++m) {
		const int orig = sumRow(D, m);
		const int target = std::min(orig, cap[m]);
		if (orig <= target) continue;

		if (verbose) {
			std::cerr << "[clamp] mission " << m
				<< ": sum_h q=" << orig << " > cap=" << cap[m]
				<< " -> scaling this row\n";
		}

		// Proportional integer scaling of row m
		struct Frac { double frac; int h; };
		double sf = (orig > 0) ? (double)target / (double)orig : 0.0;
		std::vector<int> newRow(D.H, 0);
		std::vector<Frac> order; order.reserve(D.H);
		int newSum = 0;

		for (int h = 0; h < D.H; ++h) {
			double v = sf * D.q[m][h];
			int base = (int)std::floor(v + 1e-12);
			base = std::max(0, std::min(base, D.R)); // per-bucket cap
			newRow[h] = base; newSum += base;
			order.push_back({ v - base, h });
		}
		int rem = target - newSum;
		std::sort(order.begin(), order.end(),
			[](const Frac& a, const Frac& b) { return a.frac > b.frac; });
		for (int k = 0; k < rem && k < (int)order.size(); ++k) {
			int h = order[k].h; if (newRow[h] < D.R) newRow[h] += 1;
		}
		for (int h = 0; h < D.H; ++h) D.q[m][h] = newRow[h];
	}
}

// overlap of [a0,a1) and [b0,b1)
static inline double overlapLen(double a0, double a1, double b0, double b1) {
	double L = std::max(a0, b0), R = std::min(a1, b1);
	return (R > L) ? (R - L) : 0.0;
}

// Build intervals per robot of [start,end) with mission id
static std::vector<std::vector<std::tuple<double, double, int>>>
buildRobotIntervals(const ProblemData& D,
	const IloCplex& cplex,
	const IloArray<IloBoolVarArray>& xDep,
	const IloArray< IloArray<IloBoolVarArray> >& xTask,
	const IloArray<IloBoolVarArray>& /*xEnd*/,
	const IloArray<IloNumVarArray>& t)
{
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
		// BUGFIX: use std::get<0> (not std::get[0])
		std::sort(intervals[i].begin(), intervals[i].end(),
			[](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
	}
	return intervals;
}

struct RobotCutState {
	// Where robot i starts the REPLAN from (a "pseudo-depot")
	Pt startXY;            // location at/after tau
	double release;        // earliest time robot i can start new work
	int atTask = -1;       // if finishing a task at tau (in-progress), its id; else -1
};

static double pickCutTime(const ProblemData& D,
	const IloCplex& cplex,
	const IloNumVarArray& Tvars,
	std::mt19937& rng) {
	// upper bound: last mission completion in the solved plan
	double Tmax = 0.0;
	for (int m = 0; m < D.Mmissions; ++m)
		Tmax = std::max(Tmax, cplex.getValue(Tvars[m]));

	if (Tmax <= 0.0) Tmax = D.H * D.Delta;
	std::uniform_real_distribution<double> U(0.0, Tmax);
	return U(rng);
}

struct CutResult {
	std::vector<int> removedTasks;        // tasks started before tau (including in-progress)
	std::vector<RobotCutState> robots;    // per-robot start location and release
};

static CutResult analyzeCutState(
	const ProblemData& D, double tau,
	const IloCplex& cplex,
	const IloArray<IloBoolVarArray>& xDep,
	const IloArray< IloArray<IloBoolVarArray> >& xTask,
	const IloArray<IloNumVarArray>& t
) {
	const int R = D.R, N = D.N;
	CutResult cr;
	cr.robots.assign(R, { /*startXY*/{0,0}, /*release*/tau, /*atTask*/-1 });
	std::vector<char> removed(N, 0);

	// Helper: earliest and last completed task before tau for a robot
	for (int i = 0; i < R; ++i) {
		// default: robot at its depot at time tau
		cr.robots[i].startXY = D.depotXY[i];
		cr.robots[i].release = tau;

		// find first
		int cur = -1;
		for (int j = 0; j < N; ++j)
			if (cplex.getValue(xDep[i][j]) > 0.5) { cur = j; break; }
		if (cur < 0) continue; // robot unused

		int lastDoneTask = -1;
		double lastDoneEnd = 0.0;

		while (true) {
			const double tj = cplex.getValue(t[i][cur]);
			const double ej = tj + D.s[cur];
			const bool started = (tj < tau);
			const bool finished = (ej <= tau + 1e-9);

			if (started) {
				removed[cur] = 1; // task has started -> remove from replanning
			}

			if (finished) {
				// fully completed before tau
				lastDoneTask = cur;
				lastDoneEnd = ej;
			}
			else if (started && ej > tau) {
				// in progress at tau -> must finish at ej, robot only available then
				cr.robots[i].startXY = D.taskXY[cur];
				cr.robots[i].release = ej;
				cr.robots[i].atTask = cur;   // just for reporting
			}

			// move along route
			int nxt = -1;
			for (int k = 0; k < N; ++k)
				if (k != cur && cplex.getValue(xTask[i][cur][k]) > 0.5) { nxt = k; break; }
			if (nxt == -1) break;
			cur = nxt;
		}

		// If not in-progress, set start location to last completed task (if any)
		if (cr.robots[i].atTask < 0 && lastDoneTask >= 0) {
			cr.robots[i].startXY = D.taskXY[lastDoneTask];
			cr.robots[i].release = tau; // available now
		}
		// If robot had no completed or in-progress tasks, it keeps depotXY and release=tau by default
	}

	// Collect removed task ids
	for (int j = 0; j < N; ++j) if (removed[j]) cr.removedTasks.push_back(j);
	return cr;
}

static void addNewTasksToMission(
	ProblemData& D, int mission, int K,
	std::mt19937& rng,
	int s_min = 2, int s_max = 5,          // service time bounds (min)
	double sigma = 120.0                    // spatial spread
) {
	if (K <= 0) return;
	// crude center = mean of existing tasks in mission (fallback to origin)
	Pt ctr{ 0,0 };
	int cnt = 0;
	for (int j = 0; j < D.N; ++j) if (D.missionOf[j] == mission) {
		ctr.x += D.taskXY[j].x; ctr.y += D.taskXY[j].y; ++cnt;
	}
	if (cnt > 0) { ctr.x /= cnt; ctr.y /= cnt; }

	std::normal_distribution<double> N0(0.0, sigma);
	std::uniform_int_distribution<int> Uservice(s_min, s_max);

	// Reserve space
	int oldN = D.N;
	int newN = D.N + K;
	D.s.resize(newN);
	D.missionOf.resize(newN);
	D.taskXY.resize(newN);
	for (auto& row : D.cTask) row.resize(newN, 0.0);
	D.cTask.resize(newN, std::vector<double>(newN, 0.0));
	for (auto& row : D.cDep) row.resize(newN, 0.0);
	for (auto& row : D.cEnd) row.resize(newN, 0.0);
	for (auto& row : D.compat) row.resize(newN, 1);

	// Create tasks
	for (int tId = oldN; tId < newN; ++tId) {
		D.missionOf[tId] = mission;
		D.s[tId] = (double)Uservice(rng);
		D.taskXY[tId] = { ctr.x + N0(rng), ctr.y + N0(rng) };
	}

	// Recompute travel for pairs touching new tasks
	auto timeFromDist = [&](double dist, int m_from, int m_to) {
		// mirror your earlier rule (no extra penalties by default here)
		return dist / std::max(1e-6, 50.0); // or store speed in D if you prefer
	};

	for (int j = 0; j < newN; ++j) {
		for (int k = 0; k < newN; ++k) {
			if (j == k) { D.cTask[j][k] = 0.0; continue; }
			double d = std::hypot(D.taskXY[j].x - D.taskXY[k].x,
				D.taskXY[j].y - D.taskXY[k].y);
			D.cTask[j][k] = timeFromDist(d, D.missionOf[j], D.missionOf[k]);
		}
	}
	for (int i = 0; i < D.R; ++i) {
		for (int j = 0; j < newN; ++j) {
			double d1 = std::hypot(D.depotXY[i].x - D.taskXY[j].x,
				D.depotXY[i].y - D.taskXY[j].y);
			D.cDep[i][j] = timeFromDist(d1, -1, D.missionOf[j]);
			double d2 = std::hypot(D.taskXY[j].x - D.depotXY[i].x,
				D.taskXY[j].y - D.depotXY[i].y);
			D.cEnd[i][j] = timeFromDist(d2, D.missionOf[j], -1);
		}
	}

	// Update index-by-mission list
	D.tasksOfMission.assign(D.Mmissions, {});
	for (int j = 0; j < newN; ++j) D.tasksOfMission[D.missionOf[j]].push_back(j);

	D.N = newN;
}

struct ReplanInstance {
	ProblemData D2;                   // pruned & extended problem
	std::vector<int> mapOld2New;     // old task id -> new id (or -1 if removed)
};

static ReplanInstance buildReplanInstance(
	const ProblemData& D, const CutResult& cr
) {
	ReplanInstance RPI;
	const int N = D.N, R = D.R;

	// Mark keep/remove
	std::vector<char> remove(N, 0);
	for (int j : cr.removedTasks) remove[j] = 1;

	// Build mapping old->new
	RPI.mapOld2New.assign(N, -1);
	int newN = 0;
	for (int j = 0; j < N; ++j) if (!remove[j]) RPI.mapOld2New[j] = newN++;

	// Copy kept tasks into D2
	ProblemData D2 = D;
	D2.N = newN;
	D2.s.assign(newN, 0.0);
	D2.missionOf.assign(newN, 0);
	D2.taskXY.assign(newN, Pt{ 0,0 });
	for (int m = 0; m < D2.Mmissions; ++m) D2.tasksOfMission[m].clear();

	for (int j = 0; j < N; ++j) if (!remove[j]) {
		int jj = RPI.mapOld2New[j];
		D2.s[jj] = D.s[j];
		D2.missionOf[jj] = D.missionOf[j];
		D2.taskXY[jj] = D.taskXY[j];
		D2.tasksOfMission[D2.missionOf[jj]].push_back(jj);
	}

	// Rebuild cTask
	D2.cTask.assign(newN, std::vector<double>(newN, 0.0));
	for (int j = 0; j < N; ++j) if (!remove[j]) {
		int jj = RPI.mapOld2New[j];
		for (int k = 0; k < N; ++k) if (!remove[k]) {
			int kk = RPI.mapOld2New[k];
			D2.cTask[jj][kk] = D.cTask[j][k];
		}
	}

	// Pseudo-depots and release times:
	// set depotXY to robot start positions, and bake release into cDep
	D2.depotXY.resize(R);
	for (int i = 0; i < R; ++i) D2.depotXY[i] = cr.robots[i].startXY;

	D2.cDep.assign(R, std::vector<double>(newN, 0.0));
	D2.cEnd.assign(R, std::vector<double>(newN, 0.0));
	for (int i = 0; i < R; ++i) {
		for (int jj = 0; jj < newN; ++jj) {
			// distance from pseudo-depot to task jj
			const Pt& A = D2.depotXY[i], & B = D2.taskXY[jj];
			double d = std::hypot(A.x - B.x, A.y - B.y);
			double tr = d / std::max(1e-6, 50.0); // or your stored speed rule
			// *** key trick: add release time r_i into depot travel bound ***
			D2.cDep[i][jj] = cr.robots[i].release + tr;

			// end->depot (if you keep end cost) from task to same pseudo-depot
			D2.cEnd[i][jj] = tr; // cost only (no release addition here)
		}
	}

	// Compat: drop removed-task columns
	D2.compat.assign(D.R, std::vector<int>(newN, 1));
	for (int i = 0; i < D.R; ++i)
		for (int j = 0; j < N; ++j) if (!remove[j]) {
			int jj = RPI.mapOld2New[j];
			D2.compat[i][jj] = D.compat[i][j];
		}

	RPI.D2 = std::move(D2);
	return RPI;
}


int main() {


	// -----------------------
	// 1) Build first instance
	// -----------------------
	ProblemData D;
	GeneratorConfig cfg;

	// Example config 
	cfg.R = 4; cfg.N = 20; cfg.M = 3;
	cfg.H = ceil(cfg.N * 0.6);
	//cfg.H = 20;
	cfg.Delta = 10.0; cfg.eps = 1e-3;
	cfg.use_flip = false;  cfg.flip_bucket = 4;
	cfg.use_surge = false; cfg.surge_bucket = 3; cfg.surge_mission = 1; cfg.surge_robots = 3;
	cfg.early_mission = 0; cfg.late_mission = 1;
	cfg.early_robots = 2;  cfg.late_robots = 2; cfg.other_min = 1;
	cfg.incompat_prob = 0;
	cfg.cross_mission_penalty = 0;
	cfg.seed = 12345;

	int optTime = 60;

	for (int type = 0; type < 2; type++) {
		IloEnv env;
		D.reserv = type;
		D.originalN = cfg.N;
		// Manual quotas
		cfg.manual_q.assign(cfg.M, std::vector<int>(cfg.H, 0));
		if (D.reserv) {

			//cfg.manual_q[0][0] = 3;
			cfg.manual_q[1][2] = 3;
			cfg.manual_q[1][3] = 4;
			//cfg.manual_q[1][4] = 2;
			//cfg.manual_q[2][3] = 3;
		}
		generateInstance(D, cfg);

		// Pick baseline and configure
		SolveOptions opt;
		opt.baseline = Baseline::OURS;
		configureBaseline(D, opt);

		// Clamp quotas if bucket layer is on
		if (opt.useBucketLayer && opt.useQuotas) {
			const bool creditFromPresence =
				(opt.baseline == Baseline::MONO_SOFT_BUCKETS) ? false : true;
			const bool forcePresenceCap =
				(opt.baseline == Baseline::OURS) ? false : true;
			clampQuotasAllBaselines(D, opt, creditFromPresence, forcePresenceCap, /*verbose=*/true);
		}

		// -----------------------
		// 2) Build model (plan #1)
		// -----------------------
		auto art = buildModel(D, opt, env);

		IloCplex cplex(art.model);

		// Lazy subtour cuts
		addSubtourLazyCallback(cplex, art.xTask, D.R, D.N);

		// (your CPLEX params)
		//cplex.setParam(IloCplex::Param::RandomSeed, cfg.seed);  // pick any fixed integer
		cplex.setParam(IloCplex::Param::TimeLimit, optTime);
		cplex.setParam(IloCplex::Param::MIP::Display, 2);
		cplex.setParam(IloCplex::Param::Threads, 0);
		cplex.setParam(IloCplex::Param::Parallel, 0);
		cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, 10);
		cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, 50);
		cplex.setParam(IloCplex::Param::MIP::PolishAfter::Time, 10);
		cplex.setParam(IloCplex::Param::Emphasis::MIP, 1);
		cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, 2);
		cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, 2);
		cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, 2);
		cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, 2);
		cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 1);

		if (!cplex.solve()) {
			std::cerr << "CPLEX failed or no solution (plan #1).\n";
			env.end(); return 1;
		}
		std::cout << "Status: " << cplex.getStatus()
			<< "  Objective: " << std::fixed << std::setprecision(3)
			<< cplex.getObjValue() << "\n";

		const bool activeBuckets = hasActiveBucketAccounting(D);
		const bool useY = opt.useBucketLayer && opt.useQuotas && opt.useReservations && activeBuckets;
		const bool useZ = opt.useBucketLayer && opt.useQuotas && !opt.useReservations && activeBuckets;
		const IloArray< IloArray<IloBoolVarArray> >* yPtr = useY ? &art.y : nullptr;
		const IloArray< IloArray<IloBoolVarArray> >* zPtr = useZ ? &art.z : nullptr;

		// ---- Report & export for plan #1
		reportAndExportSolution(
			D, opt, cplex,
			art.xDep, art.xTask, art.xEnd, art.t, art.T,
			yPtr, zPtr,
			/*credit_spillover=*/true,
			/*planID=*/1,
			/*robotRelease=*/nullptr
		);

		// -----------------------
		// 3) Replanning workflow
		// -----------------------
		std::mt19937 rng(1337);

		// Pick a cut time tau based on the solution (you already had this helper)
		//double tau = pickCutTime(D, cplex, art.T, rng) / 2.0;
		double tau = 10.00;
		std::cout << "\n[Replan] Cutting the plan at tau = " << tau << "\n";

		// Analyze executed part & robot states
		auto cut = analyzeCutState(D, tau, cplex, art.xDep, art.xTask, art.t);

		std::vector<double> robotRelease(D.R, 0.0);
		for (int i = 0; i < D.R; ++i) robotRelease[i] = cut.robots[i].release;

		std::cout << "[Replan] Removed (started) tasks: ";
		for (int j : cut.removedTasks) std::cout << j << " ";
		std::cout << "\n";
		for (int i = 0; i < D.R; ++i) {
			std::cout << "  Robot " << i
				<< " release=" << cut.robots[i].release
				<< " at (" << cut.robots[i].startXY.x << "," << cut.robots[i].startXY.y << ")";
			if (cut.robots[i].atTask >= 0) std::cout << " finishing task " << cut.robots[i].atTask;
			std::cout << "\n";
		}


		// Build reduced instance (remaining tasks only, depot moved to pseudo-depots, releases embedded)
		auto RPI = buildReplanInstance(D, cut);

		// Remove all manual quotas for the replan
		// Zero all per-bucket quotas in the replan instance

		//if (!D.reserv) {
		//	RPI.D2.q.assign(RPI.D2.Mmissions, std::vector<int>(RPI.D2.H, 0)); // removes all manual quotas
		//	RPI.D2.Rtot.assign(RPI.D2.Mmissions, 0.0);
		//}
		// After building RPI
		int shift = static_cast<int>(std::ceil(tau / RPI.D2.Delta - 1e-9));
		if (!RPI.D2.q.empty()) {
			for (int m = 0; m < RPI.D2.Mmissions; ++m) {
				std::rotate(RPI.D2.q[m].begin(),
					RPI.D2.q[m].begin() + std::min(shift, RPI.D2.H),
					RPI.D2.q[m].end());
				// Zero out freed tail
				for (int h = RPI.D2.H - shift; h < RPI.D2.H; ++h)
					if (h >= 0) RPI.D2.q[m][h] = 0;
			}
		}


		// (Optional) If you use rolling totals Rtot, scale them to remaining horizon:
		//if (!RPI.D2.Rtot.empty()) {
		//	double frac_remaining = (RPI.D2.H > 0) ? std::max(0.0, (RPI.D2.H - h0) / static_cast<double>(RPI.D2.H)) : 1.0;
		//	for (int m = 0; m < RPI.D2.Mmissions; ++m) {
		//		RPI.D2.Rtot[m] *= frac_remaining; // simple conservative scaling
		//	}
		//}
		int oldN2_beforeAdd = RPI.D2.N;
		// Optionally add new tasks before replanning
		int Knew = 10;
		int missionAdd = 1;
		addNewTasksToMission(RPI.D2, missionAdd, Knew, rng);
		/*Knew = 5;
		missionAdd = 2;
		addNewTasksToMission(RPI.D2, missionAdd, Knew, rng);*/

		// Configure options for replan (copy baseline behavior)
		SolveOptions opt2 = opt;
		configureBaseline(RPI.D2, opt2);
		if (opt2.useBucketLayer && opt2.useQuotas) {
			const bool creditFromPresence2 =
				(opt2.baseline == Baseline::MONO_SOFT_BUCKETS) ? false : true;
			const bool forcePresenceCap2 =
				(opt2.baseline == Baseline::OURS) ? false : true;
			clampQuotasAllBaselines(RPI.D2, opt2, creditFromPresence2, forcePresenceCap2, /*verbose=*/false);
		}

		// -----------------------
		// 4) Build model (plan #2)
		// -----------------------
		auto art2 = buildModel(RPI.D2, opt2, env);

		IloCplex cplex2(art2.model);

		// Lazy subtour cuts
		addSubtourLazyCallback(cplex2, art2.xTask, RPI.D2.R, RPI.D2.N);

		// FEASIBILITY-ORIENTED PROFILE
		cplex2.setParam(IloCplex::Param::TimeLimit, optTime);
		cplex2.setParam(IloCplex::Param::MIP::Display, 2);
		cplex2.setParam(IloCplex::Param::Threads, 12);
		cplex2.setParam(IloCplex::Param::Parallel, 0);
		cplex2.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, 10);
		cplex2.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, 50);
		cplex2.setParam(IloCplex::Param::MIP::PolishAfter::Time, 10);
		cplex2.setParam(IloCplex::Param::Emphasis::MIP, 1);
		cplex2.setParam(IloCplex::Param::MIP::Cuts::Cliques, 2);
		cplex2.setParam(IloCplex::Param::MIP::Cuts::Gomory, 2);
		cplex2.setParam(IloCplex::Param::MIP::Cuts::MIRCut, 2);
		cplex2.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, 2);
		cplex2.setParam(IloCplex::Param::Preprocessing::Presolve, 1);


		if (!cplex2.solve()) {
			std::cerr << "Replan CPLEX failed or no solution (plan #2).\n";
			env.end(); return 1;
		}

		std::cout << "[REPLAN] Status: " << cplex2.getStatus()
			<< "  Obj: " << std::fixed << std::setprecision(3)
			<< cplex2.getObjValue() << "\n";

		const bool activeBuckets2 = hasActiveBucketAccounting(RPI.D2);
		const bool useY2 = opt2.useBucketLayer && opt2.useQuotas && opt2.useReservations && activeBuckets2;
		const bool useZ2 = opt2.useBucketLayer && opt2.useQuotas && !opt2.useReservations && activeBuckets2;
		const IloArray< IloArray<IloBoolVarArray> >* yPtr2 = useY2 ? &art2.y : nullptr;
		const IloArray< IloArray<IloBoolVarArray> >* zPtr2 = useZ2 ? &art2.z : nullptr;

		// ---- Report & export for plan #2
		std::string baselineTag = reportAndExportSolution(
			RPI.D2, opt2, cplex2,
			art2.xDep, art2.xTask, art2.xEnd, art2.t, art2.T,
			yPtr2, zPtr2,
			/*credit_spillover=*/true,
			/*planID=*/2,
			&robotRelease
		);


		// ---- Overwrite existing "2" files with full-horizon (pre-cut + post-cut) ----
		int h0 = static_cast<int>(std::ceil(tau / D.Delta - 1e-9));
		h0 = std::max(0, std::min(h0, D.H));

		// Re-evaluate KPIs so we can read per-bucket capacities for both plans
		auto kpi1 = evaluateApplesToApples(
			D, cplex, art.xDep, art.xTask, art.xEnd, art.t, art.T,
			/*credit_spillover=*/true,
			useY ? &art.y : nullptr,
			useZ ? &art.z : nullptr
		);
		auto kpi2 = evaluateApplesToApples(
			RPI.D2, cplex2, art2.xDep, art2.xTask, art2.xEnd, art2.t, art2.T,
			/*credit_spillover=*/true,
			useY2 ? &art2.y : nullptr,
			useZ2 ? &art2.z : nullptr
		);

			// new tasks in D2 are exactly the indices [oldN2_beforeAdd, RPI.D2.N)
		std::ofstream nf((fs::path(baselineTag) / "sol_new_2.csv").string());
		nf << "task\n";
		for (int jj = oldN2_beforeAdd; jj < RPI.D2.N; ++jj) {
			nf << jj << "\n";
		}


		env.end();
	}
	return 0;
}
