#ifndef A2A_H
#define A2A_H

#include <vector>
#include <ilcplex/ilocplex.h>
#include "definitions.h"

// KPI container used by evaluateApplesToApples
struct EvalKPI {
    double rawObj = 0.0;
    double score = 0.0;
    double travel = 0.0;
    int    switches = 0;
    double shortQ = 0.0;
    double shortR = 0.0;
    std::vector<double> T;
    std::vector<std::vector<int>> cap;
};

// True if q/Rtot has any non-zero requirement
bool hasActiveBucketAccounting(const ProblemData& D);

// Compute apples-to-apples KPIs given a solved model
EvalKPI evaluateApplesToApples(
    const ProblemData& D,
    const IloCplex& cplex,
    const IloArray<IloBoolVarArray>& xDep,
    const IloArray< IloArray<IloBoolVarArray> >& xTask,
    const IloArray<IloBoolVarArray>& xEnd,
    const IloArray<IloNumVarArray>& t,
    const IloNumVarArray& Tvars,
    bool credit_spillover,
    const IloArray< IloArray<IloBoolVarArray> >* y_ptr, // may be nullptr
    const IloArray< IloArray<IloBoolVarArray> >* z_ptr  // may be nullptr
);

#endif // A2A_H
