#pragma once
#include "symbolic.hpp"
#include <vector>
#include <map>
#include <string>
#include <memory>

namespace lamina {

struct PiecewiseSolution {
    struct Case {
        std::shared_ptr<SymbolicExpr> condition;
        std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solutions;
    };
    std::vector<Case> cases;
};

class LAMINA_API ParametricSolver {
public:

    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
    solve_system(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    static PiecewiseSolution solve_system_piecewise(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

private:

    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
    solve_linear_parametric(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
    solve_polynomial_parametric(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    static bool is_linear_in_unknowns(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns);
};

}
