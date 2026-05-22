#pragma once
#include "interval.hpp"
#include "symbolic.hpp"
#include <vector>
#include <memory>
#include <string>
#include <utility>

namespace lamina {

struct SignChartEntry {
    Interval interval;
    int sign;
};

struct LAMINA_API PiecewiseIntervalResult {
    struct Case {
        std::shared_ptr<SymbolicExpr> condition;
        IntervalUnion solution;
    };
    std::vector<Case> cases;

    bool is_empty() const { return cases.empty(); }

    bool is_single() const { return cases.size() == 1; }

    IntervalUnion single_solution() const {
        if (cases.empty()) return IntervalUnion::empty();
        return cases[0].solution;
    }
};

class LAMINA_API InequalitySolver {
public:

    static IntervalUnion solve_inequality(
        const std::shared_ptr<SymbolicExpr>& expr,
        InequalityType type,
        const std::string& variable);

    static IntervalUnion solve_inequalities(
        const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                     InequalityType>>& inequalities,
        const std::string& variable);

    static IntervalUnion solve_rational_inequality(
        const std::shared_ptr<SymbolicExpr>& numerator,
        const std::shared_ptr<SymbolicExpr>& denominator,
        InequalityType type,
        const std::string& variable);

    static PiecewiseIntervalResult solve_parametric_inequality(
        const std::shared_ptr<SymbolicExpr>& expr,
        InequalityType type,
        const std::string& variable,
        const std::vector<std::string>& parameters);

private:

    static std::vector<SignChartEntry> build_sign_chart(
        const std::shared_ptr<SymbolicExpr>& poly,
        const std::string& variable,
        const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
        const std::vector<int>& multiplicities);

    static IntervalUnion select_intervals(
        const std::vector<SignChartEntry>& chart,
        InequalityType type,
        const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
        const std::vector<int>& multiplicities);

    static IntervalUnion build_parametric_solution(
        const std::vector<std::shared_ptr<SymbolicExpr>>& symbolic_roots,
        const std::vector<int>& multiplicities,
        int leading_sign,
        InequalityType type);
};

}
