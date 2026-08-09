#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "result.hpp"

class SymbolicExpr;

namespace lamina {

using ExprPtr = std::shared_ptr<SymbolicExpr>;

enum class VerificationStatus {
    NotChecked,
    Verified,
    Refuted,
    Inconclusive
};

template <typename T>
struct Conditional {
    T value;
    std::vector<ExprPtr> conditions;
    VerificationStatus verification = VerificationStatus::NotChecked;
};

struct FiniteSolution {
    ExprPtr value;
    std::size_t multiplicity = 1;
    std::vector<ExprPtr> conditions;
};

struct IntervalSolution {
    ExprPtr lower;
    ExprPtr upper;
    bool lower_closed = true;
    bool upper_closed = true;
    std::vector<ExprPtr> conditions;
};

struct ConditionSet {
    std::string variable;
    ExprPtr predicate;
    std::vector<ExprPtr> conditions;
};

class SolutionSet {
public:
    enum class Kind {
        Empty,
        Universal,
        Finite,
        Intervals,
        Conditional,
        Inconclusive
    };

    static SolutionSet empty() { return SolutionSet(Kind::Empty); }
    static SolutionSet universal() { return SolutionSet(Kind::Universal); }

    static SolutionSet finite(std::vector<FiniteSolution> solutions) {
        SolutionSet set(Kind::Finite);
        set.finite_solutions_ = std::move(solutions);
        return set;
    }

    static SolutionSet intervals(std::vector<IntervalSolution> intervals) {
        SolutionSet set(Kind::Intervals);
        set.intervals_ = std::move(intervals);
        return set;
    }

    static SolutionSet conditional(ConditionSet condition_set) {
        SolutionSet set(Kind::Conditional);
        set.condition_set_ = std::move(condition_set);
        return set;
    }

    static SolutionSet inconclusive(std::string reason) {
        SolutionSet set(Kind::Inconclusive);
        set.reason_ = std::move(reason);
        return set;
    }

    Kind kind() const noexcept { return kind_; }
    const std::vector<FiniteSolution>& finite_solutions() const noexcept {
        return finite_solutions_;
    }
    const std::vector<IntervalSolution>& interval_solutions() const noexcept {
        return intervals_;
    }
    const ConditionSet& condition_set() const noexcept { return condition_set_; }
    const std::string& reason() const noexcept { return reason_; }

private:
    explicit SolutionSet(Kind kind) : kind_(kind) {}

    Kind kind_;
    std::vector<FiniteSolution> finite_solutions_;
    std::vector<IntervalSolution> intervals_;
    ConditionSet condition_set_;
    std::string reason_;
};

struct TransformResult {
    Conditional<ExprPtr> expression;
    std::vector<ExprPtr> convergence_domain;
    std::vector<ExprPtr> roc;
    VerificationStatus verification = VerificationStatus::NotChecked;
};

struct IntegralResult {
    Conditional<ExprPtr> result;
    std::vector<ExprPtr> domain;
    std::vector<ExprPtr> singularities;
    VerificationStatus verification = VerificationStatus::NotChecked;
};

} // namespace lamina
