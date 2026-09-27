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

/**
 * @brief conditions 是在操作上下文假设下须同时满足的条件。
 * 空条件保留原表达式定义域且不增加约束；整数参数遍历全体整数。
 */
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

/**
 * @brief 成功的 SolveResult 表示操作定义域内的完整解集。
 * 完备性未决时返回外层 Inconclusive 错误；EmptySolutions 表示已证明的空解集。
 */
using SolutionSet = std::variant<
    EmptySolutions,
    UniversalSolutions,
    FiniteSolutions,
    IntervalSolutions,
    ConditionalSolutions,
    ParametricSolutions>;
using SolveResult = Result<SolutionSet>;

/**
 * @brief 受检变换成功时返回在 conditions 和收敛域 roc 条件下有效的已求值表达式。
 * 支持性未决或收敛条件无法表示时，返回外层 Inconclusive 错误。
 */
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

/**
 * @brief 成功保留未求值积分。
 * 闭式、收敛性及初等原函数不存在性的证明需另行提供。
 */
struct UnevaluatedIntegral {
    ExprPtr integral;
};

using IntegralOutcome = std::variant<
    Verified<ClosedFormIntegral>,
    UnevaluatedIntegral>;
using IntegralResult = Result<IntegralOutcome>;

} // namespace LMCAS
