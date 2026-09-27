#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "poly_utils.hpp"

#include <vector>

namespace LMCAS {
namespace {

using CancellationNode = std::shared_ptr<const SymbolicNode>;
using CancellationExpr = std::shared_ptr<SymbolicExpr>;

static CancellationNode positive_denominator_exponent(const NumberNode& exponent) {
    const auto& value = exponent.value();
    if (std::holds_alternative<BigInt>(value)) {
        const auto& integer = std::get<BigInt>(value);
        if (!integer.is_negative()) { return nullptr; }
        return detail::make_node<NumberNode>(BigInt(0) - integer);
    }
    if (std::holds_alternative<lmmc_real_t>(value)) {
        const auto real = std::get<lmmc_real_t>(value);
        if (!(real < 0)) { return nullptr; }
        return detail::make_node<NumberNode>(-real);
    }
    const auto& rational = std::get<Rational>(value);
    if (!rational.get_numerator().is_negative()) { return nullptr; }
    return detail::make_node<NumberNode>(Rational(
        BigInt(0) - rational.get_numerator(), rational.get_denominator()));
}

static void classify_fraction_factor(const CancellationNode& factor,
    std::vector<CancellationNode>& numerator, std::vector<CancellationNode>& denominator) {
    auto power = std::dynamic_pointer_cast<const PowerNode>(factor);
    if (power) {
        auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
        auto positive = exponent ? positive_denominator_exponent(*exponent) : nullptr;
        if (positive) {
            auto number = std::dynamic_pointer_cast<const NumberNode>(positive);
            if (number && number->is_one()) {
                denominator.push_back(power->base());
            } else {
                denominator.push_back(detail::make_node<PowerNode>(power->base(), positive));
            }
            return;
        }
    }
    numerator.push_back(factor);
}

static void separate_num_den(const CancellationNode& node,
    std::vector<CancellationNode>& numerator, std::vector<CancellationNode>& denominator) {
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : product->operands()) {
            classify_fraction_factor(operand, numerator, denominator);
        }
    } else {
        classify_fraction_factor(node, numerator, denominator);
    }
}

static CancellationExpr build_product(const std::vector<CancellationNode>& factors) {
    if (factors.empty()) { return SymbolicExpr::number(1); }
    if (factors.size() == 1) { return detail::make_expression_ptr(factors[0]); }
    return detail::make_expression_ptr(detail::make_node<MultiplyNode>(factors));
}

struct CancellationFraction {
    CancellationExpr numerator;
    CancellationExpr denominator;
};

static bool collect_fraction_terms(const AddNode& addition,
                                   std::vector<CancellationFraction>& terms) {
    bool has_denominator = false;
    for (const auto& term : addition.operands()) {
        std::vector<CancellationNode> numerator, denominator;
        separate_num_den(term, numerator, denominator);
        auto num_expr = build_product(numerator)->simplify();
        auto den_expr = build_product(denominator)->simplify();
        if (!den_expr->is_one()) { has_denominator = true; }
        terms.push_back({num_expr, den_expr});
    }
    return has_denominator;
}

static bool same_denominators(const std::vector<CancellationFraction>& terms) {
    auto first = terms[0].denominator->to_string();
    for (size_t i = 1; i < terms.size(); ++i) {
        if (terms[i].denominator->to_string() != first) { return false; }
    }
    return true;
}

static CancellationExpr distinct_denominator_product(
    const std::vector<CancellationFraction>& terms) {
    auto total_den = SymbolicExpr::number(1);
    std::vector<CancellationExpr> unique_dens;
    for (const auto& term : terms) {
        if (term.denominator->is_one()) { continue; }
        bool found = false;
        for (const auto& denominator : unique_dens) {
            if (denominator->to_string() == term.denominator->to_string()) {
                found = true;
                break;
            }
        }
        if (!found) { unique_dens.push_back(term.denominator); }
    }
    for (const auto& denominator : unique_dens) {
        total_den = SymbolicExpr::multiply(total_den, denominator)->simplify();
    }
    return total_den;
}

static CancellationFraction combine_fraction_terms(
    const std::vector<CancellationFraction>& terms) {
    const bool same = same_denominators(terms);
    auto denominator = same ? terms[0].denominator : distinct_denominator_product(terms);
    std::vector<CancellationNode> operands;
    for (const auto& term : terms) {
        auto numerator = term.numerator;
        if (!same) {
            const auto& total_node = detail::node(denominator);
            const auto& term_node = detail::node(term.denominator);
            auto factor = total_node && term_node &&
                    total_node->compare(*term_node) == 0
                ? SymbolicExpr::number(1)
                : SymbolicExpr::divide(
                      denominator, term.denominator)->simplify();
            numerator = SymbolicExpr::multiply(
                term.numerator, factor)->simplify();
        }
        if (numerator && detail::node(numerator)) {
            operands.push_back(detail::node(numerator));
        }
    }
    if (operands.empty()) { return {SymbolicExpr::number(0), nullptr}; }
    auto numerator = operands.size() == 1 ? detail::make_expression_ptr(operands[0])
        : detail::make_expression_ptr(detail::make_node<AddNode>(operands));
    numerator = same ? numerator->simplify() : numerator->expand()->simplify();
    return {numerator, denominator};
}


template <typename Coefficient>
static void cancel_polynomial_gcd(const Polynomial<Coefficient>& numerator,
    const Polynomial<Coefficient>& denominator, CancellationExpr& cur_num,
    CancellationExpr& cur_den) {
    if (denominator.is_zero()) return;
    auto gcd = Polynomial<Coefficient>::gcd(numerator, denominator);
    if (gcd.degree() >= 1) {
        auto [q_num, r_num] = numerator.div_mod(gcd);
        auto [q_den, r_den] = denominator.div_mod(gcd);
        if (!r_num.is_zero() || !r_den.is_zero()) return;
        cur_num = poly_to_symbolic(q_num)->simplify();
        cur_den = poly_to_symbolic(q_den)->simplify();
    }
}

static CancellationExpr cancel_variable(CancellationExpr& cur_num,
    CancellationExpr& cur_den, const std::string& var, const CancellationExpr& simp) {
    auto num_expanded = cur_num->expand()->simplify();
    auto den_expanded = cur_den->expand()->simplify();
    auto rational_num = symbolic_to_poly<Rational>(num_expanded, var);
    auto rational_den = symbolic_to_poly<Rational>(den_expanded, var);
    if (rational_num && rational_den) {
        if (rational_den.value().is_zero()) { return simp; }
        if (rational_num.value().is_zero()) { return SymbolicExpr::number(0); }
        cancel_polynomial_gcd(rational_num.value(), rational_den.value(), cur_num, cur_den);
        return nullptr;
    }
    auto symbolic_num = symbolic_to_poly<SymbolicPolyCoeff>(num_expanded, var);
    auto symbolic_den = symbolic_to_poly<SymbolicPolyCoeff>(den_expanded, var);
    if (!symbolic_num || !symbolic_den) { return nullptr; }
    if (symbolic_den.value().is_zero()) { return simp; }
    if (symbolic_num.value().is_zero()) { return SymbolicExpr::number(0); }
    cancel_polynomial_gcd(symbolic_num.value(), symbolic_den.value(), cur_num, cur_den);
    return nullptr;
}

static CancellationExpr finish_cancellation(const CancellationExpr& cur_num,
                                             const CancellationExpr& cur_den) {
    if (cur_den->is_one()) { return cur_num; }
    if (detail::node(cur_den)) {
        auto diff = SymbolicExpr::add(detail::make_expression_ptr(detail::node(cur_den)),
            SymbolicExpr::number(-1))->simplify();
        if (diff->is_zero()) {
            return SymbolicExpr::multiply(SymbolicExpr::number(-1), cur_num)->simplify();
        }
    }
    return SymbolicExpr::divide(cur_num, cur_den)->simplify();
}

static CancellationExpr reduce_fraction(CancellationFraction fraction,
    const CancellationExpr& simp, bool from_sum) {
    auto rational_variables = free_variables(detail::node(fraction.numerator));
    const auto denominator_variables = free_variables(detail::node(fraction.denominator));
    rational_variables.insert(denominator_variables.begin(), denominator_variables.end());
    if (!from_sum && rational_variables.empty()) {
        return SymbolicExpr::divide(fraction.numerator, fraction.denominator)->simplify();
    }
    for (const auto& var : rational_variables) {
        try {
            auto terminal = cancel_variable(
                fraction.numerator, fraction.denominator, var, simp);
            if (terminal) { return terminal; }
        } catch (const std::invalid_argument&) {
            continue;
        } catch (const std::out_of_range&) {
            continue;
        } catch (const std::runtime_error&) {
            continue;
        }
    }
    return finish_cancellation(fraction.numerator, fraction.denominator);
}

}

std::shared_ptr<SymbolicExpr> SymbolicExpr::cancel() const {
    auto simp = simplify();
    if (!simp || !detail::node(simp)) { return simp; }
    if (auto addition = std::dynamic_pointer_cast<const AddNode>(detail::node(simp))) {
        std::vector<CancellationFraction> terms;
        if (collect_fraction_terms(*addition, terms)) {
            auto combined = combine_fraction_terms(terms);
            if (!combined.denominator) { return combined.numerator; }
            combined.denominator = combined.denominator->expand()->simplify();
            return reduce_fraction(std::move(combined), simp, true);
        }
    }
    std::vector<CancellationNode> num_factors;
    std::vector<CancellationNode> den_factors;
    separate_num_den(detail::node(simp), num_factors, den_factors);
    if (den_factors.empty()) { return simp; }
    auto numerator = build_product(num_factors)->simplify();
    auto denominator = build_product(den_factors)->simplify();
    if (denominator->is_one()) { return numerator; }
    return reduce_fraction({numerator, denominator}, simp, false);
}

}
