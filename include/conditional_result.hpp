#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "proof_outcome.hpp"
#include "result.hpp"

namespace LMCAS {

class SymbolicExpr;



using ExprPtr = std::shared_ptr<SymbolicExpr>;

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

struct ParametricSolution {
    ExprPtr value;
    std::vector<std::string> integer_parameters;
    std::vector<ExprPtr> conditions;
};

struct ConditionSet {
    std::string variable;
    ExprPtr predicate;
    std::vector<ExprPtr> conditions;
};

struct EmptySolutions {};
struct UniversalSolutions {};
struct FiniteSolutions {
    std::vector<FiniteSolution> values;
};
struct IntervalSolutions {
    std::vector<IntervalSolution> values;
};
struct ConditionalSolutions {
    ConditionSet value;
};
struct ParametricSolutions {
    std::vector<ParametricSolution> values;
};

using SolutionSet = std::variant<
    EmptySolutions,
    UniversalSolutions,
    FiniteSolutions,
    IntervalSolutions,
    ConditionalSolutions,
    ParametricSolutions>;
using SolveResult = Result<SolutionSet>;

struct EvaluatedTransform {
    ExprPtr expression;
    std::vector<ExprPtr> conditions;
    std::vector<ExprPtr> roc;
};
using TransformEngineResult = Result<Verified<EvaluatedTransform>>;

struct ClosedFormIntegral {
    ExprPtr expression;
    std::vector<ExprPtr> conditions;
};

struct UnevaluatedIntegral {
    ExprPtr integral;
};

using IntegralOutcome = std::variant<
    Verified<ClosedFormIntegral>,
    UnevaluatedIntegral>;
using IntegralResult = Result<IntegralOutcome>;

} // namespace LMCAS
