#pragma once
#ifndef REPORT_H
#define REPORT_H

#include <ilcplex/ilocplex.h>
#include "definitions.h"

// Forward declare CPLEX array types we use (already included via ilocplex.h)
const std::string reportAndExportSolution(
    const ProblemData& D,
    const SolveOptions& opt,
    IloCplex& cplex,
    const IloArray<IloBoolVarArray>& xDep,
    const IloArray< IloArray<IloBoolVarArray> >& xTask,
    const IloArray<IloBoolVarArray>& xEnd,
    const IloArray<IloNumVarArray>& t,
    const IloNumVarArray& T,
    const IloArray< IloArray<IloBoolVarArray> >* yPtr, // pass nullptr if not used
    const IloArray< IloArray<IloBoolVarArray> >* zPtr, // pass nullptr if not used
    bool credit_spillover = true,
    int id=1,
    const std::vector<double>* robotRelease = nullptr  //
);

#endif // REPORT_H
