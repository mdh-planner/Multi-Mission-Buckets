#pragma once
#include <ilcplex/ilocplex.h>
#include "definitions.h"

// Everything needed later for solving + reporting.
struct BuildArtifacts {
    IloModel model;

    // core vars
    IloArray<IloBoolVarArray> xDep;
    IloArray<IloBoolVarArray> xEnd;
    IloArray<IloArray<IloBoolVarArray>> xTask;
    IloArray<IloNumVarArray> t;
    IloNumVarArray T;
    IloNumVarArray rEnd;   // NEW: robot completion times (return to depot)

    // bucket layer (may be empty, depending on options)
    IloArray<IloArray<IloBoolVarArray>> y;   // reservations (OURS)
    IloArray<IloArray<IloBoolVarArray>> b;   // start-bucket selection
    IloArray<IloArray<IloBoolVarArray>> z;   // soft-mono credits
    IloArray<IloNumVarArray> sQ;             // per-bucket quota slack
    IloNumVarArray sR;                        // rolling totals slack

    // convenience flags for main()
    bool useY = false;
    bool useZ = false;

    // CPLEX handle is *not* created here; you do it in main().
    BuildArtifacts(IloEnv& env, int R, int N, int M, int H)
        : model(env),
        xDep(env, R), xEnd(env, R), xTask(env, R),
        t(env, R), T(env, M, 0.0, IloInfinity, ILOFLOAT),
        y(env, R), b(env, R), z(env, R),
        sQ(env, M), sR(env, M, 0.0, IloInfinity, ILOFLOAT) {}
};

// Build the full optimization model (vars + constraints + objective),
// but DO NOT solve. Returns all variables/model so main can:
//   IloCplex cplex(art.model); cplex.solve(); reportAndExportSolution(...);
BuildArtifacts buildModel(const ProblemData& D, const SolveOptions& opt, IloEnv& env);

