#include "internal/calculus_utils_support.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS::calculus_utils_detail {

Result<void> calculus_utils_validate_expr(const std::shared_ptr<SymbolicExpr>& expr,
                                          const std::string& var,
                                          ComputationContext& context,
                                          const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return step;
    }
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "expression cannot be null", operation);
    }
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "variable name cannot be empty", operation);
    }
    return Result<void>::success();
}

Result<void> calculus_utils_validate_two_exprs(
    const std::shared_ptr<SymbolicExpr>& first,
    const std::shared_ptr<SymbolicExpr>& second,
    const std::string& var,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return step;
    }
    if (!first || !LMCAS::detail::node(first) || !second || !LMCAS::detail::node(second)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "expression cannot be null", operation);
    }
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "variable name cannot be empty", operation);
    }
    return Result<void>::success();
}

Result<void> calculus_utils_validate_expr_target(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& target,
    ComputationContext& context,
    const std::string& operation)
{
    auto input = calculus_utils_validate_expr(expr, var, context, operation);
    if (!input) {
        return input;
    }
    if (!target || !LMCAS::detail::node(target)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "target expression cannot be null",
                                     operation);
    }
    return Result<void>::success();
}

Result<void> calculus_utils_validate_expr_bounds(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation)
{
    auto input = calculus_utils_validate_expr(expr, var, context, operation);
    if (!input) {
        return input;
    }
    if (!a || !LMCAS::detail::node(a) || !b || !LMCAS::detail::node(b)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "bounds cannot be null", operation);
    }
    return Result<void>::success();
}
static bool is_infinity_node(const std::shared_ptr<const SymbolicNode>& node)
{
    auto function = std::dynamic_pointer_cast<const FunctionNode>(node);
    return function && function->type() == FunctionNode::FuncType::Infinity;
}

/**
 * @internal
 * @brief 判断表达式是否包含无穷大节点.
 */
bool calculus_utils_is_infinity(const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }

    if (auto f = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return f->type() == FunctionNode::FuncType::Infinity;
    }
    if (auto m = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        for (auto& op : m->operands()) {
            if (is_infinity_node(op)) {
                return true;
            }
        }
    }
    return false;
}

/**
 * @internal
 * @brief 判断两个表达式是否结构相等.
 */
bool calculus_utils_expr_equal(const std::shared_ptr<SymbolicExpr>& a,
                                      const std::shared_ptr<SymbolicExpr>& b)
{
    const bool a_exists = a && LMCAS::detail::node(a);
    const bool b_exists = b && LMCAS::detail::node(b);
    if (!a_exists || !b_exists) {
        return a_exists == b_exists;
    }
    return LMCAS::detail::node(a)->equals(*LMCAS::detail::node(b));
}
std::shared_ptr<SymbolicExpr> calculus_utils_make_abs(
    const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return nullptr;
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(expr)}));
}

}

namespace LMCAS {

std::shared_ptr<SymbolicExpr> log_differentiate(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    if (!f || !LMCAS::detail::node(f)) {
        return nullptr;
    }

    /// 计算 ln(f)
    auto ln_f = SymbolicExpr::ln(f);

    /// 对 ln(f) 求导 - 化简器会自动应用对数规则
    /// (ln(a*b) = ln(a)+ln(b), ln(a^n) = n*ln(a))
    auto ln_f_simplified = ln_f->simplify();
    auto d_ln_f = ln_f_simplified->differentiate(var);

    /// 结果 = f * d/dx[ln(f)]
    auto result = SymbolicExpr::multiply(f, d_ln_f);
    return result->simplify();
}

std::shared_ptr<SymbolicExpr> differential(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    if (!f || !LMCAS::detail::node(f)) {
        return nullptr;
    }

    /// 微分 df = f'(var) * dx,返回系数 f'(var)
    return f->differentiate(var);
}

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::string>> total_differential(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars)
{
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::string>> result;
    if (!f || !LMCAS::detail::node(f)) {
        return result;
    }

    result.reserve(vars.size());
    for (const auto& v : vars) {
        auto partial = f->differentiate(v);
        result.emplace_back(partial, v);
    }
    return result;
}

}
