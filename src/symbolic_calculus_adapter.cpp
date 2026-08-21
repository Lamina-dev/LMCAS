#include "assumption_context.hpp"
#include "integration.hpp"
#include "symbolic_ast.hpp"
#include "visitors/limit_visitor.hpp"

std::shared_ptr<SymbolicExpr> SymbolicExpr::integral(
    std::shared_ptr<SymbolicExpr> operand,
    const std::string& variable_name) {
    return operand ? operand->integrate(variable_name) : nullptr;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::limit_func(
    std::shared_ptr<SymbolicExpr> operand,
    const std::string& variable_name,
    std::shared_ptr<SymbolicExpr> target) {
    return operand ? operand->limit(variable_name, std::move(target)) : nullptr;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::limit(
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& point,
    const std::string& direction) const {
    if (!impl_->root || !point) return nullptr;

    LimitVisitor visitor(variable, lamina::detail::node(point), direction);
    impl_->root->accept(visitor);
    return visitor.get_result()
        ? lamina::detail::make_expression_ptr(visitor.get_result())
        : nullptr;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::integrate(
    const std::string& variable) const {
    if (!impl_->root) return nullptr;

    lamina::Integrator integrator;
    auto result = lamina::detail::make_expression_ptr(
        integrator.integrate(*this, variable));
    return result->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::series(
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& point,
    int order,
    const lamina::AssumptionContext*) const {
    if (!impl_->root || !point) return nullptr;
    if (order < 0) return SymbolicExpr::number(0);

    auto variable_expr = SymbolicExpr::variable(variable);
    auto delta = SymbolicExpr::add(
        variable_expr,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), point));
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    auto derivative = lamina::detail::make_expression_ptr(impl_->root->clone());

    for (int n = 0; n <= order; ++n) {
        if (n > 0) {
            derivative = derivative->differentiate(variable);
            if (!derivative) break;
            derivative = derivative->simplify();
        }

        auto coefficient = derivative->substitute(variable, point);
        if (!coefficient) break;
        coefficient = coefficient->simplify();

        auto term = coefficient;
        if (n > 0) {
            term = SymbolicExpr::multiply(
                term,
                SymbolicExpr::power(delta, SymbolicExpr::number(n)));
        }
        if (n > 1) {
            auto factorial = BigInt::factorial(static_cast<unsigned int>(n));
            term = SymbolicExpr::multiply(
                term,
                SymbolicExpr::number(Rational(BigInt(1), factorial)));
        }
        terms.push_back(lamina::detail::node(term));
    }

    if (terms.empty()) return SymbolicExpr::number(0);
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<AddNode>(terms))->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::make_integral(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable) {
    if (!expression) return nullptr;
    std::vector<std::shared_ptr<const SymbolicNode>> arguments = {
        lamina::detail::node(expression),
        lamina::detail::make_node<VariableNode>(variable),
    };
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Calculus_Integral, arguments));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::make_limit(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& point) {
    if (!expression || !point) return nullptr;
    std::vector<std::shared_ptr<const SymbolicNode>> arguments = {
        lamina::detail::node(expression),
        lamina::detail::make_node<VariableNode>(variable),
        lamina::detail::node(point),
    };
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Limit, arguments));
}
