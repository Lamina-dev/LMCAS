#include "internal/interval_endpoint.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/exact_root.hpp"
#include <cmath>

namespace LMCAS::detail {

Result<Rational> exact_double_rational(double value,
                                       ComputationContext& context,
                                       const std::string& operation) {
    if (!std::isfinite(value)) {
        return Result<Rational>::failure(
            CasErrc::NumericFailure,
            "finite interval endpoint evaluated to NaN or infinity",
            operation);
    }
    if (value == 0.0) {
        return Result<Rational>::success(Rational(0));
    }

    int exponent = 0;
    std::frexp(value, &exponent);
    const int binary_exponent = exponent - 53;
    const std::size_t required_bits = binary_exponent < 0
        ? static_cast<std::size_t>(-binary_exponent) + 1
        : static_cast<std::size_t>(binary_exponent) + 54;
    auto bit_budget = context.require_integer_bits(required_bits, operation);
    if (!bit_budget) {
        return Result<Rational>::failure(bit_budget.error());
    }

    return Result<Rational>::success(Rational::from_double(value));
}

namespace {


std::optional<ComparableEndpoint> parse_quadratic_surd(
    const std::shared_ptr<const SymbolicNode>& node);

std::optional<ComparableEndpoint> combine_surd_sum(
    const ComparableEndpoint& left, const ComparableEndpoint& right) {
        if (left.infinity != 0 || right.infinity != 0) {
            return std::nullopt;
        }
        if (left.radical_coefficient != Rational(0) &&
            right.radical_coefficient != Rational(0) &&
            left.radicand != right.radicand) {
            return std::nullopt;
        }
        const Rational radicand = left.radical_coefficient != Rational(0)
            ? left.radicand : right.radicand;
        return ComparableEndpoint{
            0,
            left.rational + right.rational,
            left.radical_coefficient + right.radical_coefficient,
            radicand};
}

std::optional<ComparableEndpoint> combine_surd_product(
    const ComparableEndpoint& left, const ComparableEndpoint& right) {
        if (left.infinity != 0 || right.infinity != 0) {
            return std::nullopt;
        }
        const bool left_radical = left.radical_coefficient != Rational(0);
        const bool right_radical = right.radical_coefficient != Rational(0);
        if (left_radical && right_radical && left.radicand != right.radicand) {
            return std::nullopt;
        }
        const Rational radicand = left_radical ? left.radicand : right.radicand;
        const Rational rational = left.rational * right.rational +
            left.radical_coefficient * right.radical_coefficient * radicand;
        const Rational coefficient =
            left.rational * right.radical_coefficient +
            left.radical_coefficient * right.rational;
        return ComparableEndpoint{0, rational, coefficient, radicand};
}

std::optional<ComparableEndpoint> parse_surd_sum(const AddNode& add) {
        ComparableEndpoint result{0, Rational(0), Rational(0), Rational(0)};
        for (const auto& operand : add.operands()) {
            auto parsed = parse_quadratic_surd(operand);
            if (!parsed) {
                return std::nullopt;
            }
            auto combined = combine_surd_sum(result, *parsed);
            if (!combined) {
                return std::nullopt;
            }
            result = std::move(*combined);
        }
        return result;
}

std::optional<ComparableEndpoint> parse_surd_product(const MultiplyNode& multiply) {
        ComparableEndpoint result{0, Rational(1), Rational(0), Rational(0)};
        for (const auto& operand : multiply.operands()) {
            auto parsed = parse_quadratic_surd(operand);
            if (!parsed) {
                return std::nullopt;
            }
            auto combined = combine_surd_product(result, *parsed);
            if (!combined) {
                return std::nullopt;
            }
            result = std::move(*combined);
        }
        return result;
}

std::optional<ComparableEndpoint> parse_surd_power(const PowerNode& power) {
        auto exponent = exact_rational_value(power.exponent());
        if (!exponent) {
            return std::nullopt;
        }
        if (*exponent == Rational(1, 2)) {
            auto radicand = exact_rational_value(power.base());
            if (!radicand || *radicand <= Rational(0)) {
                return std::nullopt;
            }
            return ComparableEndpoint{0, Rational(0), Rational(1), *radicand};
        }
        auto base = parse_quadratic_surd(power.base());
        if (!base || *exponent != Rational(-1)) {
            return std::nullopt;
        }
        const Rational norm = base->rational * base->rational -
            base->radical_coefficient * base->radical_coefficient * base->radicand;
        if (norm == Rational(0)) {
            return std::nullopt;
        }
        return ComparableEndpoint{
            0,
            base->rational / norm,
            (Rational(0) - base->radical_coefficient) / norm,
            base->radicand};
}

std::optional<ComparableEndpoint> parse_quadratic_surd(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) {
        return std::nullopt;
    }
    if (auto number = exact_rational_value(node)) {
        return ComparableEndpoint{0, *number, Rational(0), Rational(0)};
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() != FunctionNode::FuncType::Sqrt ||
            function->arguments().size() != 1) {
            return std::nullopt;
        }
        auto radicand = exact_rational_value(function->arguments()[0]);
        if (!radicand || *radicand <= Rational(0)) {
            return std::nullopt;
        }
        return ComparableEndpoint{0, Rational(0), Rational(1), *radicand};
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return parse_surd_sum(*add);
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return parse_surd_product(*multiply);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return parse_surd_power(*power);
    }
    return std::nullopt;
}

Result<ComparableEndpoint> comparable_number(
    const NumberNode& number, ComputationContext& context,
    const std::string& operation) {
        if (std::holds_alternative<BigInt>(number.value())) {
            return Result<ComparableEndpoint>::success(ComparableEndpoint{
                0, Rational(std::get<BigInt>(number.value())),
                Rational(0), Rational(0)});
        }
        if (std::holds_alternative<Rational>(number.value())) {
            return Result<ComparableEndpoint>::success(ComparableEndpoint{
                0, std::get<Rational>(number.value()), Rational(0), Rational(0)});
        }
        auto rational = exact_double_rational(
            std::get<lmmc_real_t>(number.value()), context, operation);
        if (!rational) {
            return Result<ComparableEndpoint>::failure(rational.error());
        }
        return Result<ComparableEndpoint>::success(
            ComparableEndpoint{0, std::move(rational.value()),
                               Rational(0), Rational(0)});
}

Result<ComparableEndpoint> comparable_surd(
    ComparableEndpoint surd, ComputationContext& context) {
        if (surd.radical_coefficient != Rational(0)) {
            const Rational a = surd.rational;
            const Rational b = surd.radical_coefficient;
            const Rational d = surd.radicand;
            auto algebraic = LMCAS::detail::make_exact_real_algebraic(
                Polynomial<Rational>({
                    a * a - b * b * d,
                    Rational(-2) * a,
                    Rational(1)
                }),
                b > Rational(0) ? 1 : 0, 1, context);
            if (!algebraic) {
                return Result<ComparableEndpoint>::failure(algebraic.error());
            }
            surd.algebraic = std::move(algebraic.value());
        }
        return Result<ComparableEndpoint>::success(std::move(surd));
}

Result<ComparableEndpoint> comparable_root(
    const RootOfNode& root, ComputationContext& context,
    const std::string& operation) {
        auto isolated = LMCAS::detail::isolate_exact_root(
            root.exact_id(), context, operation);
        if (!isolated) {
            return Result<ComparableEndpoint>::failure(isolated.error());
        }
        auto* real = std::get_if<LMCAS::detail::RealIsolation>(
            &isolated.value());
        if (!real) {
            return Result<ComparableEndpoint>::failure(
                CasErrc::DomainError,
                "interval endpoint RootOf is not real",
                operation);
        }
        ComparableEndpoint endpoint_value;
        endpoint_value.algebraic = std::move(real->value);
        return Result<ComparableEndpoint>::success(std::move(endpoint_value));
}

bool outside_real_function_domain(const FunctionNode& function) {
    if (function.arguments().size() != 1) {
        return false;
    }
    auto argument = exact_rational_value(function.arguments()[0]);
    if (!argument) {
        return false;
    }
    if (function.type() == FunctionNode::FuncType::Ln) {
        return *argument <= Rational(0);
    }
    return function.type() == FunctionNode::FuncType::Sqrt &&
        *argument < Rational(0);
}

Result<ComparableEndpoint> unsupported_endpoint(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& operation) {
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return Result<ComparableEndpoint>::failure(
            CasErrc::UnboundSymbol,
            "interval endpoint contains unbound symbol '" + variable->name() + "'",
            operation);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (outside_real_function_domain(*function)) {
            return Result<ComparableEndpoint>::failure(
                CasErrc::DomainError,
                "interval endpoint is outside the real function domain",
                operation);
        }
    }
    return Result<ComparableEndpoint>::failure(
        CasErrc::Inconclusive,
        "exact ordering of a symbolic interval endpoint is not proven",
        operation);
}

Result<ComparableEndpoint> comparable_finite_endpoint(
    const Endpoint& endpoint, ComputationContext& context,
    const std::string& operation) {
    if (!endpoint.value || !LMCAS::detail::node(endpoint.value)) {
        return Result<ComparableEndpoint>::failure(
            CasErrc::InvalidArgument,
            "finite interval endpoint must contain an expression",
            operation);
    }

    auto simplified = endpoint.value->simplify();
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return Result<ComparableEndpoint>::failure(
            CasErrc::InternalInvariant,
            "interval endpoint simplification produced a null expression",
            operation);
    }

    const auto& node = LMCAS::detail::node(simplified);
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return comparable_number(*number, context, operation);
    }
    if (auto surd = parse_quadratic_surd(node)) {
        return comparable_surd(std::move(*surd), context);
    }
    if (auto root = std::dynamic_pointer_cast<const RootOfNode>(node)) {
        return comparable_root(*root, context, operation);
    }
    return unsupported_endpoint(node, operation);
}

}

Result<ComparableEndpoint> comparable_endpoint(
    const Endpoint& endpoint,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return Result<ComparableEndpoint>::failure(step.error());
    }

    if (endpoint.is_neg_infinity && endpoint.is_pos_infinity) {
        return Result<ComparableEndpoint>::failure(
            CasErrc::InvalidArgument,
            "an endpoint cannot be both negative and positive infinity",
            operation);
    }
    if (endpoint.is_neg_infinity || endpoint.is_pos_infinity) {
        if (endpoint.value) {
            return Result<ComparableEndpoint>::failure(
                CasErrc::InvalidArgument,
                "infinite endpoints cannot also contain a finite value",
                operation);
        }
        if (!endpoint.is_open) {
            return Result<ComparableEndpoint>::failure(
                CasErrc::InvalidArgument,
                "infinite interval endpoints must be open",
                operation);
        }
        return Result<ComparableEndpoint>::success(
            ComparableEndpoint{endpoint.is_neg_infinity ? -1 : 1,
                               Rational(0), Rational(0), Rational(0)});
    }
    return comparable_finite_endpoint(endpoint, context, operation);
}

}
