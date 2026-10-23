#include "internal/rewrite_budget.hpp"
#include <string>
#include <limits>
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/expression_transform.hpp"
#include "internal/visitors/print_visitor.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include "internal/visitors/expand_visitor.hpp"
#include "matcher.hpp"
#include "numeric_evaluation.hpp"
#include "root_of_identity.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"

namespace LMCAS {


bool SymbolicExpr::is_number() const {
    return impl_->root->is_number();
}

bool SymbolicExpr::get_number_value_is_zero() const {
    return impl_->root->is_zero();
}

bool SymbolicExpr::is_big_int() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    return number_node &&
        std::holds_alternative<BigInt>(number_node->value());
}

bool SymbolicExpr::is_rational() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    return number_node &&
        std::holds_alternative<Rational>(number_node->value());
}

bool SymbolicExpr::is_int() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (!number_node) return false;
    if (std::holds_alternative<BigInt>(number_node->value())) return true;
    if (std::holds_alternative<Rational>(number_node->value())) {
        return std::get<Rational>(number_node->value()).is_integer();
    }

    const lmmc_real_t value = std::get<lmmc_real_t>(number_node->value());
    const lmmc_real_t rounded = std::round(value);
    int equal = 0;
    lmmc_double_nearly_equal_tol(value, rounded, 1e-12, 1e-12, &equal);
    return equal != 0;
}

std::variant<int, BigInt, Rational> SymbolicExpr::get_number() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (!number_node) {
        throw std::runtime_error("Expression is not a number");
    }
    if (std::holds_alternative<BigInt>(number_node->value())) {
        return std::get<BigInt>(number_node->value());
    }
    if (std::holds_alternative<Rational>(number_node->value())) {
        return std::get<Rational>(number_node->value());
    }
    return Rational::from_double(std::get<lmmc_real_t>(number_node->value()));
}

int SymbolicExpr::get_int() const {
    if (!is_int()) {
        throw std::runtime_error("Expression is not an integer");
    }
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (std::holds_alternative<BigInt>(number_node->value())) {
        auto value = std::get<BigInt>(number_node->value()).try_to_int64();
        if (!value || *value < std::numeric_limits<int>::min() ||
            *value > std::numeric_limits<int>::max()) {
            throw std::out_of_range("integer value does not fit in int");
        }
        return static_cast<int>(*value);
    }
    if (std::holds_alternative<Rational>(number_node->value())) {
        auto value = std::get<Rational>(number_node->value()).to_bigint().try_to_int64();
        if (!value || *value < std::numeric_limits<int>::min() ||
            *value > std::numeric_limits<int>::max()) {
            throw std::out_of_range("integer value does not fit in int");
        }
        return static_cast<int>(*value);
    }
    const auto value = std::get<lmmc_real_t>(number_node->value());
    if (value < static_cast<lmmc_real_t>(std::numeric_limits<int>::min()) ||
        value > static_cast<lmmc_real_t>(std::numeric_limits<int>::max())) {
        throw std::out_of_range("integer value does not fit in int");
    }
    return static_cast<int>(value);
}

BigInt SymbolicExpr::get_big_int() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (!number_node ||
        !std::holds_alternative<BigInt>(number_node->value())) {
        throw std::runtime_error("Expression is not a BigInt");
    }
    return std::get<BigInt>(number_node->value());
}

Rational SymbolicExpr::get_rational() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (!number_node ||
        !std::holds_alternative<Rational>(number_node->value())) {
        throw std::runtime_error("Expression is not a Rational");
    }
    return std::get<Rational>(number_node->value());
}

Rational SymbolicExpr::convert_rational() const {
    if (is_rational()) return get_rational();
    if (is_big_int()) return Rational(get_big_int());
    return Rational(0);
}


int SymbolicExpr::compare(const std::shared_ptr<SymbolicExpr>& other) const {
    if (!impl_->root || !LMCAS::detail::node(other)) return 0;

    return impl_->root->compare(*LMCAS::detail::node(other));
}


std::string SymbolicExpr::to_string() const {
    PrintVisitor printer;
    if (impl_->root) {
        impl_->root->accept(printer);
        return printer.get_result();
    }
    return "null";
}

lmmc_real_t SymbolicExpr::to_numeric() const {
    auto evaluated = LMCAS::evaluate_numeric(*this);
    if (!evaluated) {
        throw std::runtime_error("numeric evaluation failed: " + evaluated.error().message);
    }
    return static_cast<lmmc_real_t>(evaluated.value().value);
}

bool SymbolicExpr::is_zero() const {
    if (!impl_->root) return false;
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(impl_->root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value())) {
            return std::get<lmmc_real_t>(num->value()) == 0.0;
        }
        if (const auto* integer = std::get_if<BigInt>(&num->value())) {
            return integer->is_zero();
        }
        if (const auto* rational = std::get_if<Rational>(&num->value())) {
            return rational->is_zero();
        }
    }
    return false;
}

bool SymbolicExpr::is_one() const {
    if (!impl_->root) return false;
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(impl_->root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value())) {
            return std::get<lmmc_real_t>(num->value()) == 1.0;
        }
        if (const auto* integer = std::get_if<BigInt>(&num->value())) {
            return *integer == BigInt(1);
        }
        if (const auto* rational = std::get_if<Rational>(&num->value())) {
            return *rational == Rational(1);
        }
    }
    return false;
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_trig() const {
    ComputationContext context;
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
                                 context.limits().max_ast_nodes, "simplify_trig");
    auto res = detail::make_expression_ptr(detail::normalize_node(impl_->root, budget));
    if (!LMCAS::detail::node(res)) return nullptr;

    static const LMCAS::RewriteEngine engine = [] {
        LMCAS::RewriteEngine configured;
        using namespace LMCAS;
        auto x_val = wildcard("x");
        auto x = LMCAS::detail::make_expression_ptr(x_val);

        auto sinx = SymbolicExpr::sin(x);
        auto cosx = SymbolicExpr::cos(x);
        auto n2 = SymbolicExpr::number(2);
        auto sin2 = SymbolicExpr::power(sinx, n2);
        auto cos2 = SymbolicExpr::power(cosx, n2);

        auto pat1 = SymbolicExpr::add(sin2, cos2);
        configured.add_rule(Rule(*pat1, *SymbolicExpr::number(1), {"x"}));

        auto pat2 = SymbolicExpr::add(cos2, sin2);
        configured.add_rule(Rule(*pat2, *SymbolicExpr::number(1), {"x"}));

        auto two_x = SymbolicExpr::multiply(n2, x);
        auto sin2x = SymbolicExpr::sin(two_x);
        auto two_sin_cos = SymbolicExpr::multiply(n2,
            SymbolicExpr::multiply(sinx, cosx));
        configured.add_rule(Rule(*sin2x, *two_sin_cos, {"x"}));

        auto cos2x = SymbolicExpr::cos(two_x);
        auto cos2_sub_sin2 = SymbolicExpr::add(cos2,
            SymbolicExpr::multiply(sin2, SymbolicExpr::number(-1)));
        configured.add_rule(Rule(*cos2x, *cos2_sub_sin2, {"x"}));
        return configured;
    }();

    auto rewritten = engine.apply_checked(*res, context);
    if (!rewritten) return res;
    auto simplified = std::move(rewritten.value());
    auto result_ptr = LMCAS::detail::make_expression_ptr(simplified);
    return detail::make_expression_ptr(detail::normalize_node(detail::node(result_ptr), budget));
}

namespace detail {

SymbolicNodePtr normalize_node(const SymbolicNodePtr& node, RewriteBudget& budget,
                               const FactsQuery& facts, Domain domain) {
    budget.consume();
    if (!node) {
        throw CasError{CasErrc::InvalidArgument, "simplify argument cannot be null",
                       "LMCAS.simplify"};
    }
    budget.measure(node);
    NormalizationVisitor visitor(budget, facts, domain);
    node->accept(visitor);
    auto result = visitor.get_result();
    if (!result) {
        throw CasError{CasErrc::InternalInvariant, "simplify returned null",
                       "LMCAS.simplify"};
    }
    budget.measure(result);
    return result;
}

SymbolicNodePtr expand_node(const SymbolicNodePtr& node, RewriteBudget& budget,
                            const FactsQuery& facts, Domain domain) {
    budget.consume();
    if (!node) {
        throw CasError{CasErrc::InvalidArgument, "expand argument cannot be null",
                       "LMCAS.expand"};
    }
    budget.measure(node);
    ExpandVisitor visitor(budget);
    node->accept(visitor);
    return normalize_node(visitor.get_result(), budget, facts, domain);
}

SymbolicNodePtr differentiate_node(const SymbolicNodePtr& node,
                                   const std::string& variable, RewriteBudget& budget,
                                   const FactsQuery& facts, Domain domain) {
    budget.consume();
    if (!node || variable.empty()) {
        throw CasError{CasErrc::InvalidArgument, "invalid differentiation argument",
                       "LMCAS.differentiate"};
    }
    budget.measure(node);
    DifferentiationVisitor visitor(variable, budget);
    node->accept(visitor);
    return normalize_node(visitor.get_result(), budget, facts, domain);
}

SymbolicNodePtr substitute_node(const SymbolicNodePtr& node, const std::string& variable,
                                const SymbolicNodePtr& replacement, RewriteBudget& budget,
                                const FactsQuery& facts, Domain domain) {
    budget.consume();
    if (!node || !replacement || variable.empty()) {
        throw CasError{CasErrc::InvalidArgument, "invalid substitution argument",
                       "LMCAS.substitute"};
    }
    budget.measure(node);
    budget.measure(replacement);
    auto result = substitute_free(node, variable, replacement, &budget);
    budget.measure(result);
    return normalize_node(result, budget, facts, domain);
}

}


} // namespace LMCAS
