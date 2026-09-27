#include "irrational.hpp"
#include "symbolic.hpp"
#include <utility>

namespace LMCAS {
namespace {

std::shared_ptr<SymbolicExpr> symbolic_named_power(const std::string& name,
                                                  const BigInt& exponent) {
    auto base = SymbolicExpr::variable(name);
    return exponent == BigInt(1) ? base :
        SymbolicExpr::power(base, SymbolicExpr::number(exponent));
}

}

std::shared_ptr<SymbolicExpr> Irrational::Basis::to_symbolic() const {
    std::shared_ptr<SymbolicExpr> result;
    auto append = [&](std::shared_ptr<SymbolicExpr> factor) {
        result = result ? SymbolicExpr::multiply(result, std::move(factor)) : std::move(factor);
    };
    if (!e_power.is_zero()) append(symbolic_named_power("e", e_power));
    if (!pi_power.is_zero()) append(symbolic_named_power("π", pi_power));
    if (radicand != BigInt(1)) append(SymbolicExpr::sqrt(SymbolicExpr::number(radicand)));
    return result ? result : SymbolicExpr::number(1);
}

std::shared_ptr<SymbolicExpr> Irrational::terms_to_symbolic(const Terms& terms) {
    std::shared_ptr<SymbolicExpr> result;
    for (const auto& term : terms) {
        auto expression = term.first.is_constant() ? SymbolicExpr::number(term.second) :
                          term.first.to_symbolic();
        if (!term.first.is_constant() && term.second != Rational(1)) {
            expression = SymbolicExpr::multiply(SymbolicExpr::number(term.second), expression);
        }
        result = result ? SymbolicExpr::add(result, expression) : expression;
    }
    return result ? result : SymbolicExpr::number(0);
}

std::shared_ptr<SymbolicExpr> Irrational::to_symbolic() const {
    const auto rational = as_rational_checked();
    if (rational) return SymbolicExpr::number(rational.value());
    auto numerator = terms_to_symbolic(numerator_);
    if (denominator_.size() == 1 && denominator_.begin()->first.is_constant() &&
        denominator_.begin()->second == Rational(1)) return numerator;
    return SymbolicExpr::divide(numerator, terms_to_symbolic(denominator_));
}

}
