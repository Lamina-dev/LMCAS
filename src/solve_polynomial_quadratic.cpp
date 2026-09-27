#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace LMCAS {
using namespace polynomial_solver_detail;
namespace {
struct NormalizedQuadratic { double a, b, c, original_c; };

std::optional<NormalizedQuadratic> normalize_quadratic(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c) {
    auto av = finite_numeric_value(a);
    auto bv = finite_numeric_value(b);
    auto cv = finite_numeric_value(c);
    if (!is_purely_numeric(a) || !is_purely_numeric(b) || !is_purely_numeric(c)) {
        return std::nullopt;
    }
    if (!av || !bv || !cv || *av == 0.0) {
        return std::nullopt;
    }
    const double scale = std::max({std::abs(*av), std::abs(*bv), std::abs(*cv)});
    const double an = *av / scale, bn = *bv / scale, cn = *cv / scale;
    if ((*av != 0.0 && an == 0.0) || (*bv != 0.0 && bn == 0.0) ||
        (*cv != 0.0 && cn == 0.0)) {
        return std::nullopt;
    }
    return NormalizedQuadratic{an, bn, cn, *cv};
}

std::optional<std::pair<double, double>> stable_quadratic_pair(
    const NormalizedQuadratic& coefficients, double sqrt_discriminant) {
    const auto& [an, bn, cn, cv] = coefficients;
    double first, second;
    if (sqrt_discriminant == 0.0) {
        first = -bn / (2.0 * an);
        second = first;
    } else {
        const double q = -0.5 * (bn + std::copysign(sqrt_discriminant, bn));
        first = q / an;
        second = cn / q;
    }
    if (std::isfinite(first) && std::isfinite(second) &&
        (cv == 0.0 || (first != 0.0 && second != 0.0))) {
        return std::pair<double, double>{first, second};
    }
    return std::nullopt;
}

std::vector<std::shared_ptr<SymbolicExpr>> numeric_biquadratic_roots(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c) {
    auto normalized = normalize_quadratic(a, b, c);
    if (!normalized) {
        return {};
    }
    const auto& [an, bn, cn, cv] = *normalized;
    const double discriminant = std::fma(-4.0 * an, cn, bn * bn);
    if (!std::isfinite(discriminant) || discriminant < 0.0) {
        return {};
    }
    auto pair = stable_quadratic_pair(*normalized, std::sqrt(discriminant));
    if (!pair) {
        return {};
    }
    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    for (double value : {pair->first, pair->second}) {
        if (value >= 0.0) {
            const double positive = std::sqrt(value);
            roots.push_back(SymbolicExpr::number(positive));
            roots.push_back(SymbolicExpr::number(-positive));
        } else {
            auto expression = SymbolicExpr::number(value);
            roots.push_back(SymbolicExpr::sqrt(expression)->simplify());
            roots.push_back(negate(SymbolicExpr::sqrt(expression))->simplify());
        }
    }
    return roots;
}
}

namespace polynomial_solver_detail {
std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c) {
    auto normalized = normalize_quadratic(a, b, c);
    if (normalized) {
        const auto& [an, bn, cn, cv] = *normalized;
        const double discriminant = std::fma(-4.0 * an, cn, bn * bn);
        if (std::isfinite(discriminant) && discriminant >= 0.0) {
            const double square_root = std::sqrt(discriminant);
            if (square_root == 0.0) {
                return {SymbolicExpr::number(-bn / (2.0 * an))};
            }
            auto pair = stable_quadratic_pair(*normalized, square_root);
            if (pair) {
                return {SymbolicExpr::number(pair->first),
                                  SymbolicExpr::number(pair->second)};
            }
        }
    }
    auto b2 = SymbolicExpr::power(b, num(2));
    auto four_ac = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, c));
    auto delta = sub(b2, four_ac)->simplify();

    auto neg_b = negate(b);
    auto two_a = SymbolicExpr::multiply(num(2), a);
    if (delta && delta->is_zero()) {
        return { SymbolicExpr::divide(neg_b, two_a)->simplify() };
    }

    auto sqrt_delta = SymbolicExpr::sqrt(delta);
    auto x1 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, sqrt_delta), two_a)->simplify();
    auto x2 = SymbolicExpr::divide(sub(neg_b, sqrt_delta), two_a)->simplify();

    return { x1, x2 };
}
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_biquadratic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::string&) {
    auto a_simplified = a->simplify();
    if (a_simplified->get_number_value_is_zero()) {
        auto b_simplified = b->simplify();
        if (!b_simplified->get_number_value_is_zero()) {
            return solve_quadratic_internal(b, num(0), c);
        }

        auto c_simplified = c->simplify();
        if (c_simplified->get_number_value_is_zero()) {
            throw std::invalid_argument(
                "solve_biquadratic: identically zero equation has no finite root list");
        }
        if (is_purely_numeric(c_simplified)) {
            return {};
        }
        throw std::invalid_argument(
            "solve_biquadratic: equation is independent of the requested variable");
    }

    auto numeric = numeric_biquadratic_roots(a, b, c);
    if (!numeric.empty()) {
        return numeric;
    }
    auto u_roots = solve_quadratic_internal(a, b, c);

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& u : u_roots) {
        auto pos_root = SymbolicExpr::sqrt(u)->simplify();
        auto neg_root = negate(SymbolicExpr::sqrt(u))->simplify();
        results.push_back(pos_root);
        results.push_back(neg_root);
    }

    return results;
}

}
