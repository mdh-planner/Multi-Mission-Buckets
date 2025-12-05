#pragma once
#include <vector>

struct Pt { double x, y; };

//===============================
// Problem data
//===============================
struct ProblemData {
    int R = 0;            // robots
    int N = 0;            // tasks
    int Mmissions = 0;    // missions
    bool reserv = false;  // if we have manual reservations or not
    int originalN = 0;    // original task number, before any replans.
    // Time-bucket settings
    int H = 0;            // number of buckets
    double Delta = 10.0;  // bucket length (minutes)
    double eps = 1e-3;    // tighten bucket-upper-bound in (F5)

    // Task / mission data
    std::vector<double> s;                 // service time s_j
    std::vector<int>    missionOf;         // mu(j)
    std::vector<std::vector<double>> cTask;// travel j->k
    std::vector<std::vector<double>> cDep; // depot_i->j
    std::vector<std::vector<double>> cEnd; // j->end_i (optional in cost)

    std::vector<std::vector<int>> compat;  // compat[i][j] in {0,1}

    // Objective parameters
    std::vector<double> w;                 // mission weights
    double lambdaTravel = 0.0;             // weight on travel
    double etaSwitch = 0.0;                // weight on mission switch penalty
    double switchPenalty = 0.0;            // psi when mu(j)!=mu(k)

    // Quotas (per-bucket) and rolling totals
    std::vector<std::vector<int>> q;       // q[m][h] >= 0
    double rhoQ = 0;                     // soft quota penalty if 0, quotas are hard.
    std::vector<double> Rtot;              // total reserved minutes per mission over horizon
    double rhoR = 0.0;                     // soft total penalty

    // Convenience: tasks indexed per mission
    std::vector<std::vector<int>> tasksOfMission;

    // xy for plotting / route export

    std::vector<Pt> taskXY;            // task j -> (x,y)
    std::vector<Pt> depotXY;           // robot i depot -> (x,y)

};

//===============================
// Instance generator
//===============================

struct GeneratorConfig {
    // Sizes
    int R = 3;                 // robots
    int N = 10;                // tasks
    int M = 2;                 // missions
    int H = 6;                 // buckets
    double Delta = 10.0;       // minutes per bucket
    double eps = 1e-3;

    // Service-time generation (minutes)
    int s_min = 2;
    int s_max = 5;             // inclusive upper bound

    // Geometry -> time via speed
    double area_radius = 100.0;   // mission centers on a circle
    double cluster_sigma = 50.0;  // spread of tasks around center
    double speed = 50.0;          // meters per minute
    double base_travel = 0.0;      // constant offset added to all legs

    // Depots placement
    enum class DepotMode { NEAR_MISSIONS, CENTRAL, RANDOM_RING };
    DepotMode depot_mode = DepotMode::NEAR_MISSIONS;

    // Extra penalty if traveling across missions (0 disables)
    double cross_mission_penalty = 0.0;

    // Compatibility
    double incompat_prob = 0.0;    // P(robot i incompatible with task j)
    bool   ensure_feasible_compat = true; // ensure each task has ≥1 compatible robot

  // Objective weights (all ≥ 0) controlling trade-offs in the cost function:

    double lambdaTravel = 0.0005;   // Weight on total travel time.
                                  // Higher -> prioritize shorter routes; lower -> allow longer travel
                                  // if it improves other goals (e.g., earlier mission completion).

    double etaSwitch = 1.0;    // Weight on penalizing mission switches between consecutive tasks.
                                  // Higher -> discourage switching missions along a robot's route.

    double switchPenalty = 100;     // Base penalty applied per mission switch.
                                  // Effective per-switch cost is (etaSwitch * switchPenalty).
                                  // Keep this as a clean "unit cost" and tune etaSwitch for scaling.

    std::vector<double> w;         // size M; if empty => all 1.0

    // Quotas (per-bucket minima). If manual_q present, it overrides flip/surge.
    // manual_q[m][h] in {0..R}
    std::vector<std::vector<int>> manual_q;

    // Priority flip configuration (used only if manual_q empty)
    bool use_flip = true;
    int  flip_bucket = -1;         // default: H/2 (rounded down) if < 0
    int  early_mission = 0;        // priority before flip
    int  late_mission = 1;         // priority after flip
    int  early_robots = 2;         // min robots early for early_mission
    int  late_robots = 2;          // min robots late for late_mission
    int  other_min = 1;            // min robots for the non-priority mission in each phase

    // Surge (used only if manual_q empty)
    bool use_surge = false;
    int  surge_bucket = -1;        // specific bucket gets surge
    int  surge_mission = 0;        // mission that gets surge
    int  surge_robots = 3;         // required robots in surge bucket

    // Softness of quotas & rolling totals
    double rhoQ = 100.0;            // penalty per robot shortfall per bucket (0 => hard)
    std::vector<double> Rtot;      // per-mission minutes reserved across horizon; if empty => none
    double rhoR = 0.5;             // penalty per minute shortfall (0 => hard)

    // RNG
    unsigned seed = 42;
};

//===============================
// Baselines / options
//===============================
enum class Baseline { OURS, WEIGHTS_ONLY, MONOLITHIC, DECOUPLED, MONO_SOFT_BUCKETS };


struct SolveOptions {
    Baseline baseline = Baseline::OURS;
    bool useBucketLayer = false;          // create bucket layer
    bool useQuotas = false;               // (F6)-(F7)
    bool useReservations = false;         // OURS=true; MONO_SOFT_BUCKETS=false
    bool enforceOneMissionPerBucket = true; // F1
    bool enforceOneStartPerBucket = true;   // F2
    bool enforceUpperLimitTaskStartPerBucket = false; // enforce upper limit of number of tasks per bucket defined by opt.K
    int K = 2;  // set the max number of starts per (robot, bucket)
    std::vector<int> robotHomeMission;      // used by DECOUPLED
    std::vector<double> alphaWeights;       // used by WEIGHTS_ONLY
};

void generateInstance(ProblemData& D, GeneratorConfig cfg);
#pragma once
