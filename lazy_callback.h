// subtour_lazy.h
#pragma once
#include <ilcplex/ilocplex.h>
#include <vector>
#include <functional>

ILOSTLBEGIN

class SubtourLazyCallback : public IloCplex::LazyConstraintCallbackI {
    const IloArray<IloArray<IloBoolVarArray>>& xTask;
    int R, N;

public:
    SubtourLazyCallback(IloEnv env,
        const IloArray<IloArray<IloBoolVarArray>>& xTask_,
        int R_, int N_)
        : IloCplex::LazyConstraintCallbackI(env),
        xTask(xTask_), R(R_), N(N_) {}

    void main() override {
        IloEnv env = getEnv();

        // 1. Read integer solution
        std::vector<std::vector<std::vector<int>>> arcs(R,
            std::vector<std::vector<int>>(N));

        for (int i = 0; i < R; ++i) {
            for (int j = 0; j < N; ++j) {
                IloNumArray vals(env, N);
                getValues(vals, xTask[i][j]);   // integer values 0/1
                for (int k = 0; k < N; ++k) {
                    if (vals[k] > 0.5) {        // part of the route
                        if (j != k)
                            arcs[i][j].push_back(k);
                    }
                }
                vals.end();
            }
        }

        // 2. Detect subtours (violations)

        std::vector<char> visited(N);
        std::vector<char> inStack(N);
        std::vector<int> stack;

        for (int i = 0; i < R; ++i) {
            std::fill(visited.begin(), visited.end(), 0);
            std::fill(inStack.begin(), inStack.end(), 0);
            stack.clear();

            std::function<bool(int)> dfs = [&](int u) -> bool {
                visited[u] = 1;
                inStack[u] = 1;
                stack.push_back(u);

                for (int v : arcs[i][u]) {
                    if (!visited[v]) {
                        if (dfs(v)) return true;
                    }
                    else if (inStack[v]) {
                        // we found a cycle S
                        std::vector<int> S;
                        auto it = std::find(stack.begin(), stack.end(), v);
                        for (auto it2 = it; it2 != stack.end(); ++it2)
                            S.push_back(*it2);

                        if (S.size() >= 2 && S.size() < N) {
                            IloExpr lhs(env);
                            for (int r = 0; r < R; ++r)
                                for (int a : S)
                                    for (int b : S)
                                        if (a != b)
                                            lhs += xTask[r][a][b];

                            add(lhs <= (int)S.size() - 1);
                            lhs.end();
                            return true; // stop callback
                        }
                    }
                }

                inStack[u] = 0;
                stack.pop_back();
                return false;
            };

            for (int j = 0; j < N; ++j) {
                if (!visited[j]) {
                    if (dfs(j)) return;   // added cut
                }
            }
        }
    }

    IloCplex::CallbackI* duplicateCallback() const override {
        return (new(getEnv()) SubtourLazyCallback(*this));
    }
};

inline IloCplex::Callback addSubtourLazyCallback(
    IloCplex& cplex,
    const IloArray<IloArray<IloBoolVarArray>>& xTask,
    int R, int N
) {
    return cplex.use(new(cplex.getEnv())
        SubtourLazyCallback(cplex.getEnv(), xTask, R, N));
}
