#include "internal/complex_quadratic.hpp"

#include "assumption_context.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/rewrite_budget.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "residual_verification.hpp"

#include <new>
#include <optional>
#include <utility>

namespace LMCAS::detail {
namespace {

constexpr const char* kOperation = "solve_complex_quadratic";
using Node = std::shared_ptr<const SymbolicNode>;
using Coefficients = std::array<std::shared_ptr<SymbolicExpr>, 3>;

Node integer(int value) {
    return SymbolicFactory::create_number(BigInt(value));
}

Node add(const Node& left, const Node& right) {
    return SymbolicFactory::create_add({left, right});
}

Node multiply(const Node& left, const Node& right) {
    return SymbolicFactory::create_multiply({left, right});
}

Node negate(const Node& value) {
    return multiply(integer(-1), value);
}

Node divide(const Node& numerator, const Node& denominator) {
    return multiply(numerator, SymbolicFactory::create_power(denominator, integer(-1)));
}

ExpressionResult normalize(const Node& value, ComputationContext& context) {
    auto step = context.consume_steps(1, kOperation);
    if (!step) {
        return ExpressionResult::failure(step.error());
    }
    RewriteBudget budget(context, context.limits().max_recursion_depth,
                         context.limits().max_ast_nodes, kOperation);
    budget.measure(value);
    std::optional<AssumptionFacts> facts;
    if (context.assumptions()) facts.emplace(*context.assumptions());
    NormalizationVisitor visitor(context, facts ? static_cast<const FactsQuery&>(*facts)
                                       : no_facts(), Domain::Real, &budget);
    value->accept(visitor);
    auto result = visitor.get_result();
    if (!result) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
            "quadratic normalization returned null", kOperation);
    }
    budget.measure(result);
    return make_expression_ptr(std::move(result));
}

Result<Coefficients> normalize_coefficients(
    const Coefficients& input, ComputationContext& context) {
    for (const auto& coefficient : input) {
        auto step = context.consume_steps(1, kOperation);
        if (!step) {
            return Result<Coefficients>::failure(step.error());
        }
        if (!coefficient || !node(coefficient)) {
            return Result<Coefficients>::failure(CasErrc::InvalidArgument,
                "quadratic coefficients must not be null", kOperation);
        }
    }
    Coefficients result;
    for (std::size_t i = 0; i < input.size(); ++i) {
        auto simplified = normalize(node(input[i]), context);
        if (!simplified) {
            return Result<Coefficients>::failure(simplified.error());
        }
        result[i] = std::move(simplified.value());
    }
    return result;
}

Result<void> prove_coefficients(const Coefficients& coefficients,
                                const AssumptionContext& facts,
                                ComputationContext& context) {
    auto step = context.consume_steps(4, kOperation);
    if (!step) {
        return step;
    }
    if (node(coefficients[0])->is_zero()) {
        return Result<void>::failure(CasErrc::DomainError,
            "quadratic leading coefficient is zero", kOperation);
    }
    AssumptionFacts adapter(facts);
    for (const auto& coefficient : coefficients) {
        auto real = query_real_value(node(coefficient), adapter, context);
        if (!real) {
            return Result<void>::failure(real.error());
        }
        if (real.value() != Tribool::True) {
            return Result<void>::failure(CasErrc::UnsupportedExpression,
                "quadratic coefficients must be provably real", kOperation);
        }
    }
    auto nonzero = facts.is_nonzero_checked(*coefficients[0], context);
    if (!nonzero) {
        return Result<void>::failure(nonzero.error());
    }
    if (nonzero.value() == Tribool::False) {
        return Result<void>::failure(CasErrc::DomainError,
            "quadratic leading coefficient is zero", kOperation);
    }
    if (nonzero.value() != Tribool::True) {
        return Result<void>::failure(CasErrc::Inconclusive,
            "quadratic leading coefficient is not provably nonzero", kOperation);
    }
    return Result<void>::success();
}

Result<bool> negative_discriminant(const SymbolicExpr& delta,
                                   const AssumptionContext& facts,
                                   ComputationContext& context) {
    auto step = context.consume_steps(2, kOperation);
    if (!step) {
        return Result<bool>::failure(step.error());
    }
    auto nonnegative = facts.is_nonnegative_checked(delta, context);
    if (!nonnegative) {
        return Result<bool>::failure(nonnegative.error());
    }
    if (nonnegative.value() == Tribool::True) {
        return false;
    }
    auto negative = facts.is_negative_checked(delta, context);
    if (!negative) {
        return Result<bool>::failure(negative.error());
    }
    if (negative.value() == Tribool::True) {
        return true;
    }
    return Result<bool>::failure(CasErrc::Inconclusive,
        "quadratic discriminant sign is unknown", kOperation);
}

Result<QuadraticCartesianRoot> normalized_root(
    const Node& real, const Node& imag, ComputationContext& context) {
    auto real_part = normalize(real, context);
    if (!real_part) {
        return Result<QuadraticCartesianRoot>::failure(real_part.error());
    }
    auto imag_part = normalize(imag, context);
    if (!imag_part) {
        return Result<QuadraticCartesianRoot>::failure(imag_part.error());
    }
    return QuadraticCartesianRoot{std::move(real_part.value()), std::move(imag_part.value())};
}

Result<CartesianQuadraticRoots> construct_roots(
    const Coefficients& coefficients, const Node& delta, bool negative,
    ComputationContext& context) {
    const auto denominator = multiply(integer(2), node(coefficients[0]));
    const auto minus_b = negate(node(coefficients[1]));
    const auto radical = SymbolicFactory::create_power(
        negative ? negate(delta) : delta, SymbolicFactory::create_number(Rational(1, 2)));
    const auto real = divide(minus_b, denominator);
    const auto imag = divide(radical, denominator);
    auto first = normalized_root(
        negative ? real : add(real, imag),
        negative ? imag : integer(0), context);
    if (!first) {
        return Result<CartesianQuadraticRoots>::failure(first.error());
    }
    auto second = normalized_root(
        negative ? real : add(real, negate(imag)),
        negative ? negate(imag) : integer(0), context);
    if (!second) {
        return Result<CartesianQuadraticRoots>::failure(second.error());
    }
    return CartesianQuadraticRoots{std::move(first.value()), std::move(second.value())};
}

struct CartesianNodes {
    Node real;
    Node imag;
};

CartesianNodes cartesian_add(const CartesianNodes& left, const CartesianNodes& right) {
    return {add(left.real, right.real), add(left.imag, right.imag)};
}

CartesianNodes cartesian_multiply(const CartesianNodes& left, const CartesianNodes& right) {
    return {add(multiply(left.real, right.real), negate(multiply(left.imag, right.imag))),
            add(multiply(left.real, right.imag), multiply(left.imag, right.real))};
}

Result<void> prove_zero(const Node& residual, ComputationContext& context) {
    auto simplified = normalize(residual, context);
    if (!simplified) {
        return Result<void>::failure(simplified.error());
    }
    auto proof = check_zero_residual(simplified.value(), context);
    if (!proof) {
        return Result<void>::failure(proof.error());
    }
    if (!std::holds_alternative<ProvedZeroResidual>(proof.value())) {
        return Result<void>::failure(CasErrc::Inconclusive,
            "quadratic Cartesian residual could not be proved zero", kOperation);
    }
    return Result<void>::success();
}

Result<void> verify_root(const Coefficients& coefficients,
                         const QuadraticCartesianRoot& root,
                         ComputationContext& context) {
    auto step = context.consume_steps(16, kOperation);
    if (!step) {
        return step;
    }
    const auto zero = integer(0);
    const CartesianNodes a{node(coefficients[0]), zero};
    const CartesianNodes b{node(coefficients[1]), zero};
    const CartesianNodes c{node(coefficients[2]), zero};
    const CartesianNodes z{node(root.real), node(root.imag)};
    const auto residual = cartesian_add(
        cartesian_multiply(cartesian_add(cartesian_multiply(a, z), b), z), c);
    auto real = prove_zero(residual.real, context);
    if (!real) {
        return real;
    }
    return prove_zero(residual.imag, context);
}

Result<CartesianQuadraticRoots> solve_real_quadratic(
    const Coefficients& coefficients, const AssumptionContext& facts,
    ComputationContext& context) {
    auto valid = prove_coefficients(coefficients, facts, context);
    if (!valid) {
        return Result<CartesianQuadraticRoots>::failure(valid.error());
    }
    const auto b = node(coefficients[1]);
    auto delta = normalize(add(multiply(b, b), negate(multiply(integer(4),
        multiply(node(coefficients[0]), node(coefficients[2]))))), context);
    if (!delta) {
        return Result<CartesianQuadraticRoots>::failure(delta.error());
    }
    auto negative = negative_discriminant(*delta.value(), facts, context);
    if (!negative) {
        return Result<CartesianQuadraticRoots>::failure(negative.error());
    }
    auto roots = construct_roots(coefficients, node(delta.value()), negative.value(), context);
    if (!roots) {
        return roots;
    }
    for (const auto& root : roots.value()) {
        auto verified = verify_root(coefficients, root, context);
        if (!verified) {
            return Result<CartesianQuadraticRoots>::failure(verified.error());
        }
    }
    return roots;
}

}

Result<CartesianQuadraticRoots> real_quadratic_roots_checked(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    ComputationContext& context) {
    try {
        auto coefficients = normalize_coefficients({a, b, c}, context);
        if (!coefficients) {
            return Result<CartesianQuadraticRoots>::failure(coefficients.error());
        }
        if (context.assumptions()) {
            return solve_real_quadratic(coefficients.value(), *context.assumptions(), context);
        }
        const AssumptionContext facts;
        return solve_real_quadratic(coefficients.value(), facts, context);
    } catch (const CasError& error) {
        return Result<CartesianQuadraticRoots>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<CartesianQuadraticRoots>::failure(CasErrc::ResourceLimit,
            "quadratic construction allocation failed", kOperation);
    } catch (const std::exception& error) {
        return Result<CartesianQuadraticRoots>::failure(CasErrc::InternalInvariant,
            error.what(), kOperation);
    }
}

}
