#include "lsr_expr.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <utility>

#include "complex_analysis.hpp"
#include "poly_utils.hpp"
#include "symbolic_ast.hpp"

namespace lamina::lsr {
namespace {

constexpr const char* kSymOperation = "lsr.sym";
constexpr const char* kIntegerOperation = "lsr.integer";
constexpr const char* kRationalOperation = "lsr.rational";
constexpr const char* kApproxOperation = "lsr.approx_real";
constexpr const char* kImaginaryOperation = "lsr.imaginary_unit";
constexpr const char* kComplexOperation = "lsr.complex";
constexpr const char* kRealOperation = "lsr.real";
constexpr const char* kImagOperation = "lsr.imag";
constexpr const char* kConjOperation = "lsr.conj";
constexpr const char* kAbsOperation = "lsr.abs";
constexpr const char* kEquivalentOperation = "lsr.equivalent_core";
constexpr const char* kEquivalentProfileOperation = "lsr.equivalent_core.profile";
constexpr const char* kExprSetOperation = "lsr.expr_set";
constexpr const char* kSolveExprSetOperation = "lsr.solve_expr_set";
constexpr const char* kEvalComplexOperation = "lsr.eval_complex";

ExprResult expression_failure(CasErrc code, std::string message,
                              const char* operation) {
    return ExprResult::failure(code, std::move(message), operation);
}

ExprSetResult expr_set_failure(CasErrc code, std::string message,
                               const char* operation) {
    return ExprSetResult::failure(code, std::move(message), operation);
}

Result<ApproxComplex> complex_failure(CasErrc code, std::string message,
                                      const char* operation) {
    return Result<ApproxComplex>::failure(code, std::move(message), operation);
}

Result<ApproxComplex> eval_complex_failure(const CasError& error) {
    return complex_failure(error.code, error.message, kEvalComplexOperation);
}

ExprResult expr_from_complex_result(const ComplexExprResult& result,
                                    const char* operation) {
    if (!result) {
        return ExprResult::failure(result.error());
    }
    if (!result.value() || !lamina::detail::node(result.value())) {
        return expression_failure(CasErrc::InternalInvariant,
                                  "complex expression result is null",
                                  operation);
    }
    return ExprResult::success(result.value());
}

ApproxReal approx_part(double value) {
    ApproxReal part;
    part.value = value;
    part.absolute_error = std::numeric_limits<double>::epsilon() *
                          std::max(1.0, std::abs(value)) * 4.0;
    part.status = NumericStatus::Finite;
    return part;
}

ApproxComplex approx_complex(double real, double imag) {
    return ApproxComplex{approx_part(real), approx_part(imag)};
}

Result<ApproxComplex> checked_complex(double real, double imag,
                                      const char* operation) {
    if (!std::isfinite(real) || !std::isfinite(imag)) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex evaluation produced a non-finite component",
                               operation);
    }
    return Result<ApproxComplex>::success(approx_complex(real, imag));
}

Result<ApproxComplex> real_to_complex(const Result<ApproxReal>& real) {
    if (!real) {
        return eval_complex_failure(real.error());
    }
    if (!real.value().is_finite()) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex evaluation requires finite real components",
                               kEvalComplexOperation);
    }
    return Result<ApproxComplex>::success(
        ApproxComplex{real.value(), approx_part(0.0)});
}

Result<ApproxComplex> add_complex(const ApproxComplex& lhs,
                                  const ApproxComplex& rhs) {
    auto result = checked_complex(lhs.real.value + rhs.real.value,
                                  lhs.imag.value + rhs.imag.value,
                                  kEvalComplexOperation);
    if (!result) return result;
    result.value().real.absolute_error +=
        lhs.real.absolute_error + rhs.real.absolute_error;
    result.value().imag.absolute_error +=
        lhs.imag.absolute_error + rhs.imag.absolute_error;
    return result;
}

Result<ApproxComplex> multiply_complex(const ApproxComplex& lhs,
                                       const ApproxComplex& rhs) {
    const double a = lhs.real.value;
    const double b = lhs.imag.value;
    const double c = rhs.real.value;
    const double d = rhs.imag.value;
    auto result = checked_complex(a * c - b * d, a * d + b * c,
                                  kEvalComplexOperation);
    if (!result) return result;
    result.value().real.absolute_error +=
        std::abs(c) * lhs.real.absolute_error +
        std::abs(a) * rhs.real.absolute_error +
        std::abs(d) * lhs.imag.absolute_error +
        std::abs(b) * rhs.imag.absolute_error;
    result.value().imag.absolute_error +=
        std::abs(d) * lhs.real.absolute_error +
        std::abs(a) * rhs.imag.absolute_error +
        std::abs(c) * lhs.imag.absolute_error +
        std::abs(b) * rhs.real.absolute_error;
    return result;
}

Result<ApproxComplex> divide_complex(const ApproxComplex& lhs,
                                     const ApproxComplex& rhs) {
    const double c = rhs.real.value;
    const double d = rhs.imag.value;
    const double denom = c * c + d * d;
    if (denom == 0.0 || !std::isfinite(denom)) {
        return complex_failure(CasErrc::DomainError,
                               "complex division by zero or overflow",
                               kEvalComplexOperation);
    }
    return checked_complex((lhs.real.value * c + lhs.imag.value * d) / denom,
                           (lhs.imag.value * c - lhs.real.value * d) / denom,
                           kEvalComplexOperation);
}

bool is_integer_double(double value) {
    return std::isfinite(value) && std::floor(value) == value;
}

bool is_reserved_symbol_name(const std::string& name) {
    return name == "i" || name == "I" || name == "pi" || name == "π" ||
           name == "e";
}

bool is_imaginary_unit_name(const std::string& name) {
    return name == "i" || name == "I";
}

Result<void> validate_eqv_options(const EqvOptions& options) {
    if (options.budget.max_rewrite_steps == 0 ||
        options.budget.max_rewrite_depth == 0 ||
        options.budget.max_node_growth_factor == 0) {
        return Result<void>::failure(
            CasErrc::ResourceLimit,
            "equivalence rewrite budget exhausted before normalization",
            kEquivalentOperation);
    }
    return Result<void>::success();
}

bool exact_integer_node(const std::shared_ptr<const SymbolicNode>& node,
                        int expected) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return false;
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()) == BigInt(expected);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()) == Rational(expected);
    }
    return false;
}

std::optional<int> exact_small_integer_node(
    const std::shared_ptr<const SymbolicNode>& node,
    int min_value,
    int max_value) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;

    BigInt value;
    if (std::holds_alternative<BigInt>(number->value())) {
        value = std::get<BigInt>(number->value());
    } else if (std::holds_alternative<Rational>(number->value())) {
        const Rational& rational = std::get<Rational>(number->value());
        if (!rational.is_integer()) return std::nullopt;
        value = rational.to_BigInt();
    } else {
        return std::nullopt;
    }

    if (value < BigInt(min_value) || value > BigInt(max_value)) {
        return std::nullopt;
    }
    return value.to_int();
}

bool trig_square_argument(const std::shared_ptr<const SymbolicNode>& node,
                          FunctionNode::FuncType type,
                          std::shared_ptr<const SymbolicNode>& argument) {
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power || !exact_integer_node(power->exponent(), 2)) return false;
    auto function = std::dynamic_pointer_cast<const FunctionNode>(power->base());
    if (!function || function->type() != type ||
        function->arguments().size() != 1) {
        return false;
    }
    argument = function->arguments()[0];
    return true;
}

ExprPtr rewrite_trig_basic_identity(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return nullptr;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> rewritten;
        rewritten.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            auto child = rewrite_trig_basic_identity(operand);
            rewritten.push_back(lamina::detail::node(child));
        }

        std::vector<bool> used(rewritten.size(), false);
        std::vector<std::shared_ptr<const SymbolicNode>> result_nodes;
        for (std::size_t i = 0; i < rewritten.size(); ++i) {
            if (used[i]) continue;
            std::shared_ptr<const SymbolicNode> sin_arg;
            std::shared_ptr<const SymbolicNode> cos_arg;
            const bool is_sin_square = trig_square_argument(
                rewritten[i], FunctionNode::FuncType::Sin, sin_arg);
            const bool is_cos_square = trig_square_argument(
                rewritten[i], FunctionNode::FuncType::Cos, cos_arg);
            bool matched = false;
            for (std::size_t j = i + 1; j < rewritten.size(); ++j) {
                if (used[j]) continue;
                std::shared_ptr<const SymbolicNode> other_arg;
                if (is_sin_square &&
                    trig_square_argument(rewritten[j],
                                         FunctionNode::FuncType::Cos,
                                         other_arg) &&
                    sin_arg->equals(*other_arg)) {
                    matched = true;
                } else if (is_cos_square &&
                           trig_square_argument(rewritten[j],
                                                FunctionNode::FuncType::Sin,
                                                other_arg) &&
                           cos_arg->equals(*other_arg)) {
                    matched = true;
                }
                if (matched) {
                    used[i] = true;
                    used[j] = true;
                    result_nodes.push_back(
                        lamina::detail::node(SymbolicExpr::number(1)));
                    break;
                }
            }
            if (!matched && !used[i]) {
                result_nodes.push_back(rewritten[i]);
            }
        }

        if (result_nodes.empty()) {
            return SymbolicExpr::number(0);
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_add(std::move(result_nodes)))->simplify();
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(multiply->operands().size());
        for (const auto& operand : multiply->operands()) {
            auto child = rewrite_trig_basic_identity(operand);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_multiply(std::move(operands)))->simplify();
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = rewrite_trig_basic_identity(power->base());
        auto exponent = rewrite_trig_basic_identity(power->exponent());
        return SymbolicExpr::power(base, exponent)->simplify();
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(function->arguments().size());
        for (const auto& argument : function->arguments()) {
            auto child = rewrite_trig_basic_identity(argument);
            args.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(function->type(),
                                                    std::move(args)))->simplify();
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real_part = rewrite_trig_basic_identity(complex_node->real());
        auto imag_part = rewrite_trig_basic_identity(complex_node->imag());
        auto value = complex(real_part, imag_part);
        if (!value) throw std::runtime_error(value.error().message);
        return value.value()->simplify();
    }

    return lamina::detail::make_expression_ptr(node);
}

ExprPtr rewrite_exp_log_basic_identity(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return nullptr;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            auto child = rewrite_exp_log_basic_identity(operand);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_add(std::move(operands)))->simplify();
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(multiply->operands().size());
        for (const auto& operand : multiply->operands()) {
            auto child = rewrite_exp_log_basic_identity(operand);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_multiply(std::move(operands)))->simplify();
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = rewrite_exp_log_basic_identity(power->base());
        auto exponent = rewrite_exp_log_basic_identity(power->exponent());
        return SymbolicExpr::power(base, exponent)->simplify();
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(function->arguments().size());
        for (const auto& argument : function->arguments()) {
            auto child = rewrite_exp_log_basic_identity(argument);
            args.push_back(lamina::detail::node(child));
        }

        if (function->type() == FunctionNode::FuncType::Exp &&
            args.size() == 1 && exact_integer_node(args[0], 0)) {
            return SymbolicExpr::number(1);
        }
        if (function->type() == FunctionNode::FuncType::Ln &&
            args.size() == 1 && exact_integer_node(args[0], 1)) {
            return SymbolicExpr::number(0);
        }

        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(function->type(),
                                                    std::move(args)))->simplify();
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real_part = rewrite_exp_log_basic_identity(complex_node->real());
        auto imag_part = rewrite_exp_log_basic_identity(complex_node->imag());
        auto value = complex(real_part, imag_part);
        if (!value) throw std::runtime_error(value.error().message);
        return value.value()->simplify();
    }

    return lamina::detail::make_expression_ptr(node);
}

void collect_variable_names(
    const std::shared_ptr<const SymbolicNode>& node,
    std::set<std::string>& variables) {
    if (!node) return;
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        variables.insert(variable->name());
        return;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            collect_variable_names(operand, variables);
        }
        return;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            collect_variable_names(operand, variables);
        }
        return;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        collect_variable_names(power->base(), variables);
        collect_variable_names(power->exponent(), variables);
        return;
    }
    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        collect_variable_names(complex_node->real(), variables);
        collect_variable_names(complex_node->imag(), variables);
    }
}

Result<std::optional<bool>> prove_rational_polynomial_equivalence(
    const ExprPtr& difference,
    ComputationContext& context,
    const EqvOptions& options) {
    if (!difference) {
        return Result<std::optional<bool>>::failure(
            CasErrc::InternalInvariant,
            "equivalence difference is null",
            kEquivalentOperation);
    }

    if (options.budget.max_rewrite_steps < 4) {
        return Result<std::optional<bool>>::failure(
            CasErrc::ResourceLimit,
            "equivalence rewrite budget exhausted before polynomial normalization",
            kEquivalentOperation);
    }

    auto step = context.consume_steps(4, kEquivalentOperation);
    if (!step) return Result<std::optional<bool>>::failure(step.error());

    auto expanded = difference->expand();
    if (!expanded || !lamina::detail::node(expanded)) {
        return Result<std::optional<bool>>::failure(
            CasErrc::InternalInvariant,
            "equivalence expansion returned null",
            kEquivalentOperation);
    }

    std::set<std::string> variables;
    collect_variable_names(lamina::detail::node(expanded), variables);
    if (variables.size() > 1) {
        return Result<std::optional<bool>>::success(std::nullopt);
    }
    const std::string variable = variables.empty() ? "x" : *variables.begin();

    auto recognized = recognize_rational_polynomial(*expanded, variable, context);
    if (!recognized) {
        return Result<std::optional<bool>>::failure(recognized.error());
    }
    if (!recognized.value()) {
        return Result<std::optional<bool>>::success(std::nullopt);
    }
    return Result<std::optional<bool>>::success(recognized.value()->is_zero());
}

Rational polynomial_coeff_or_zero(const Polynomial<Rational>& polynomial,
                                  std::size_t degree) {
    return degree < polynomial.coeffs.size()
        ? polynomial.coeffs[degree]
        : Rational(0);
}

ExprPtr rational_expression(const Rational& value) {
    return SymbolicExpr::number(value);
}

bool exact_rational_sqrt(const Rational& value, Rational& root) {
    if (value < Rational(0)) return false;
    const BigInt numerator_root = value.get_numerator().sqrt();
    const BigInt denominator_root = value.get_denominator().sqrt();
    if (numerator_root * numerator_root != value.get_numerator() ||
        denominator_root * denominator_root != value.get_denominator()) {
        return false;
    }
    root = Rational(numerator_root, denominator_root);
    return true;
}

ExprPtr sqrt_rational_expression(const Rational& value) {
    Rational root;
    if (exact_rational_sqrt(value, root)) {
        return rational_expression(root);
    }
    return SymbolicExpr::sqrt(rational_expression(value))->simplify();
}

ExprResult verified_lsr_complex(ExprPtr real_part, ExprPtr imag_part) {
    auto value = complex(std::move(real_part), std::move(imag_part));
    if (!value) return value;
    auto simplified = value.value()->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
        return expression_failure(CasErrc::InternalInvariant,
                                  "complex root simplification returned null",
                                  kSolveExprSetOperation);
    }
    return ExprResult::success(std::move(simplified));
}

Result<std::optional<ExprSet>> try_lsr_closed_form_rational_poly_roots(
    const ExprPtr& equation,
    const std::string& variable,
    ComputationContext& context) {
    if (!equation) {
        return Result<std::optional<ExprSet>>::failure(
            CasErrc::InvalidArgument, "equation cannot be null",
            kSolveExprSetOperation);
    }

    auto recognized = recognize_rational_polynomial(*equation, variable, context);
    if (!recognized) {
        return Result<std::optional<ExprSet>>::failure(recognized.error());
    }
    if (!recognized.value()) {
        return Result<std::optional<ExprSet>>::success(std::nullopt);
    }

    const Polynomial<Rational>& polynomial = *recognized.value();
    const int degree = polynomial.degree();
    if (degree < 1 || degree > 2) {
        return Result<std::optional<ExprSet>>::success(std::nullopt);
    }

    std::vector<ExprPtr> roots;
    if (degree == 1) {
        const Rational b = polynomial_coeff_or_zero(polynomial, 1);
        if (b.is_zero()) {
            return Result<std::optional<ExprSet>>::success(std::nullopt);
        }
        const Rational c = polynomial_coeff_or_zero(polynomial, 0);
        roots.push_back(rational_expression((-c) / b)->simplify());
    } else {
        const Rational a = polynomial_coeff_or_zero(polynomial, 2);
        const Rational b = polynomial_coeff_or_zero(polynomial, 1);
        const Rational c = polynomial_coeff_or_zero(polynomial, 0);
        if (a.is_zero()) {
            return Result<std::optional<ExprSet>>::success(std::nullopt);
        }

        const Rational two_a = Rational(2) * a;
        const Rational discriminant = b * b - Rational(4) * a * c;
        const Rational real_component = (-b) / two_a;

        if (discriminant < Rational(0)) {
            const Rational positive_discriminant = -discriminant;
            const Rational positive_denominator = two_a.abs();
            auto imag_magnitude = SymbolicExpr::divide(
                sqrt_rational_expression(positive_discriminant),
                rational_expression(positive_denominator))->simplify();
            auto positive = verified_lsr_complex(
                rational_expression(real_component), imag_magnitude);
            if (!positive) {
                return Result<std::optional<ExprSet>>::failure(positive.error());
            }
            auto negative_imag = SymbolicExpr::multiply(
                SymbolicExpr::number(-1), imag_magnitude)->simplify();
            auto negative = verified_lsr_complex(
                rational_expression(real_component), negative_imag);
            if (!negative) {
                return Result<std::optional<ExprSet>>::failure(negative.error());
            }
            roots.push_back(negative.value());
            roots.push_back(positive.value());
        } else {
            auto sqrt_discriminant = sqrt_rational_expression(discriminant);
            auto numerator_left = SymbolicExpr::add(
                rational_expression(-b),
                SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                       sqrt_discriminant))->simplify();
            roots.push_back(SymbolicExpr::divide(
                numerator_left, rational_expression(two_a))->simplify());
            if (!discriminant.is_zero()) {
                auto numerator_right = SymbolicExpr::add(
                    rational_expression(-b), sqrt_discriminant)->simplify();
                roots.push_back(SymbolicExpr::divide(
                    numerator_right, rational_expression(two_a))->simplify());
            }
        }
    }

    auto set = ExprSet::make(std::move(roots));
    if (!set) {
        return Result<std::optional<ExprSet>>::failure(set.error());
    }
    return Result<std::optional<ExprSet>>::success(std::move(set.value()));
}

ExprPtr canonicalize_lsr_complex_product(const SymbolicExpr& expression) {
    const auto& node = lamina::detail::node(expression);

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = exact_small_integer_node(power->exponent(), 0, 16);
        if (exponent) {
            auto canonical_base = canonicalize_lsr_complex_product(
                *lamina::detail::make_expression_ptr(power->base()));
            if (std::dynamic_pointer_cast<const ComplexNode>(
                    lamina::detail::node(canonical_base))) {
                auto result = SymbolicExpr::number(1);
                for (int i = 0; i < *exponent; ++i) {
                    result = canonicalize_lsr_complex_product(
                        *SymbolicExpr::multiply(result, canonical_base));
                }
                return result->simplify();
            }
        }
        return lamina::detail::make_expression_ptr(node);
    }

    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) {
        return lamina::detail::make_expression_ptr(node);
    }

    bool saw_complex = false;
    auto real = SymbolicExpr::number(1);
    auto imag = SymbolicExpr::number(0);

    for (const auto& operand : multiply->operands()) {
        auto canonical_operand = canonicalize_lsr_complex_product(
            *lamina::detail::make_expression_ptr(operand));
        const auto& operand_node = lamina::detail::node(canonical_operand);
        ExprPtr factor_real;
        ExprPtr factor_imag;
        if (auto complex_operand = std::dynamic_pointer_cast<const ComplexNode>(
                operand_node)) {
            saw_complex = true;
            factor_real = lamina::detail::make_expression_ptr(complex_operand->real());
            factor_imag = lamina::detail::make_expression_ptr(complex_operand->imag());
        } else {
            factor_real = canonical_operand;
            factor_imag = SymbolicExpr::number(0);
        }

        auto ac = SymbolicExpr::multiply(real, factor_real);
        auto bd = SymbolicExpr::multiply(imag, factor_imag);
        auto ad = SymbolicExpr::multiply(real, factor_imag);
        auto bc = SymbolicExpr::multiply(imag, factor_real);
        auto next_real = SymbolicExpr::add(
            ac, SymbolicExpr::multiply(SymbolicExpr::number(-1), bd))->simplify();
        auto next_imag = SymbolicExpr::add(ad, bc)->simplify();
        real = std::move(next_real);
        imag = std::move(next_imag);
    }

    if (!saw_complex) {
        return lamina::detail::make_expression_ptr(node);
    }
    auto result = complex(real, imag);
    if (!result) {
        throw std::runtime_error(result.error().message);
    }
    return result.value()->simplify();
}

Result<ApproxComplex> evaluate_complex_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const NumericBindings& bindings,
    ComputationContext& context) {
    auto step = context.consume_steps(1, kEvalComplexOperation);
    if (!step) return Result<ApproxComplex>::failure(step.error());
    if (!node) {
        return complex_failure(CasErrc::InvalidArgument,
                               "expression contains a null node",
                               kEvalComplexOperation);
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real = evaluate_numeric(
            *lamina::detail::make_expression_ptr(complex_node->real()),
            bindings, context);
        if (!real) return eval_complex_failure(real.error());
        auto imag = evaluate_numeric(
            *lamina::detail::make_expression_ptr(complex_node->imag()),
            bindings, context);
        if (!imag) return eval_complex_failure(imag.error());
        if (!real.value().is_finite() || !imag.value().is_finite()) {
            return complex_failure(CasErrc::NumericFailure,
                                   "complex components must be finite",
                                   kEvalComplexOperation);
        }
        return Result<ApproxComplex>::success(
            ApproxComplex{real.value(), imag.value()});
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        auto sum = Result<ApproxComplex>::success(approx_complex(0.0, 0.0));
        for (const auto& operand : add->operands()) {
            auto term = evaluate_complex_node(operand, bindings, context);
            if (!term) return term;
            sum = add_complex(sum.value(), term.value());
            if (!sum) return sum;
        }
        return sum;
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        auto product = Result<ApproxComplex>::success(approx_complex(1.0, 0.0));
        for (const auto& operand : multiply->operands()) {
            auto factor = evaluate_complex_node(operand, bindings, context);
            if (!factor) return factor;
            product = multiply_complex(product.value(), factor.value());
            if (!product) return product;
        }
        return product;
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = evaluate_complex_node(power->base(), bindings, context);
        if (!base) return base;
        auto exponent = evaluate_numeric(
            *lamina::detail::make_expression_ptr(power->exponent()),
            bindings, context);
        if (!exponent) return eval_complex_failure(exponent.error());
        const double exponent_value = exponent.value().value;
        if (!is_integer_double(exponent_value) ||
            std::abs(exponent_value) > 64.0) {
            return complex_failure(CasErrc::UnsupportedExpression,
                                   "complex evaluation only supports integer powers with |n| <= 64",
                                   kEvalComplexOperation);
        }
        int exponent_int = static_cast<int>(exponent_value);
        auto result = Result<ApproxComplex>::success(approx_complex(1.0, 0.0));
        for (int i = 0; i < std::abs(exponent_int); ++i) {
            result = multiply_complex(result.value(), base.value());
            if (!result) return result;
        }
        if (exponent_int < 0) {
            result = divide_complex(approx_complex(1.0, 0.0), result.value());
        }
        return result;
    }

    return real_to_complex(
        evaluate_numeric(*lamina::detail::make_expression_ptr(node),
                         bindings, context));
}

} // namespace

ExprSet::ExprSet(std::vector<ExprPtr> elements)
    : elements_(std::move(elements)) {}

Result<ExprSet> ExprSet::make(std::vector<ExprPtr> elements) {
    std::vector<ExprPtr> unique;
    unique.reserve(elements.size());
    for (auto& element : elements) {
        if (!element) {
            return expr_set_failure(CasErrc::InvalidArgument,
                                    "set<Expr> elements cannot be null",
                                    kExprSetOperation);
        }
        bool duplicate = false;
        for (const auto& existing : unique) {
            if (structurally_equal(*existing, *element)) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            unique.push_back(std::move(element));
        }
    }
    return Result<ExprSet>::success(ExprSet(std::move(unique)));
}

bool ExprSet::contains(const SymbolicExpr& expression) const {
    for (const auto& element : elements_) {
        if (element && structurally_equal(*element, expression)) {
            return true;
        }
    }
    return false;
}

bool ExprSet::subset_of(const ExprSet& other) const {
    for (const auto& element : elements_) {
        if (!element || !other.contains(*element)) {
            return false;
        }
    }
    return true;
}

ExprSet ExprSet::set_union(const ExprSet& other) const {
    std::vector<ExprPtr> result = elements_;
    for (const auto& element : other.elements_) {
        if (element && !contains(*element)) {
            result.push_back(element);
        }
    }
    return ExprSet(std::move(result));
}

ExprSet ExprSet::intersection(const ExprSet& other) const {
    std::vector<ExprPtr> result;
    for (const auto& element : elements_) {
        if (element && other.contains(*element)) {
            result.push_back(element);
        }
    }
    return ExprSet(std::move(result));
}

ExprSet ExprSet::difference(const ExprSet& other) const {
    std::vector<ExprPtr> result;
    for (const auto& element : elements_) {
        if (element && !other.contains(*element)) {
            result.push_back(element);
        }
    }
    return ExprSet(std::move(result));
}

ExprSet ExprSet::symmetric_difference(const ExprSet& other) const {
    auto left_only = difference(other);
    auto right_only = other.difference(*this);
    return left_only.set_union(right_only);
}

ExprResult sym(const std::string& name) {
    if (name.empty()) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "symbol name cannot be empty", kSymOperation);
    }
    if (is_reserved_symbol_name(name)) {
        if (is_imaginary_unit_name(name)) {
            return expression_failure(CasErrc::InvalidArgument,
                                      "imaginary unit symbol is reserved",
                                      kSymOperation);
        }
        return expression_failure(CasErrc::InvalidArgument,
                                  "reserved mathematical constants cannot be shadowed",
                                  kSymOperation);
    }
    try {
        auto expression = SymbolicExpr::variable(name);
        if (!expression) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "symbol factory returned null", kSymOperation);
        }
        return ExprResult::success(std::move(expression));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "symbol allocation failed", kSymOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kSymOperation);
    }
}

ExprResult integer(long long value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "integer allocation failed", kIntegerOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kIntegerOperation);
    }
}

ExprResult integer(const BigInt& value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "integer allocation failed", kIntegerOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kIntegerOperation);
    }
}

ExprResult rational(const Rational& value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "rational allocation failed", kRationalOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kRationalOperation);
    }
}

ExprResult approx_real(double value) {
    if (!std::isfinite(value)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "approximate real must be finite", kApproxOperation);
    }
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "approximate real allocation failed", kApproxOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kApproxOperation);
    }
}

ExprResult imaginary_unit() {
    auto zero = integer(0);
    if (!zero) return ExprResult::failure(zero.error());
    auto one = integer(1);
    if (!one) return ExprResult::failure(one.error());
    return complex(zero.value(), one.value());
}

ExprResult complex(ExprPtr real, ExprPtr imag) {
    if (!real || !imag) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "complex expression parts cannot be null",
                                  kComplexOperation);
    }
    try {
        auto node = SymbolicFactory::create_complex(lamina::detail::node(real),
                                                    lamina::detail::node(imag));
        return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "complex expression allocation failed",
                                  kComplexOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kComplexOperation);
    }
}

ExprResult real(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(real_part_checked(expression, context),
                                    kRealOperation);
}

ExprResult real(const ExprPtr& expression) {
    ComputationContext context;
    return real(expression, context);
}

ExprResult imag(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(imag_part_checked(expression, context),
                                    kImagOperation);
}

ExprResult imag(const ExprPtr& expression) {
    ComputationContext context;
    return imag(expression, context);
}

ExprResult conj(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(conjugate_checked(expression, context),
                                    kConjOperation);
}

ExprResult conj(const ExprPtr& expression) {
    ComputationContext context;
    return conj(expression, context);
}

ExprResult abs(const ExprPtr& expression, ComputationContext& context) {
    auto step = context.consume_steps(1, kAbsOperation);
    if (!step) return ExprResult::failure(step.error());
    auto re = real(expression, context);
    if (!re) return re;
    auto im = imag(expression, context);
    if (!im) return im;
    try {
        auto re_squared = SymbolicExpr::power(re.value(), SymbolicExpr::number(2));
        auto im_squared = SymbolicExpr::power(im.value(), SymbolicExpr::number(2));
        auto sum = SymbolicExpr::add(re_squared, im_squared);
        auto result = SymbolicExpr::sqrt(sum)->simplify();
        if (!result || !lamina::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "complex absolute value construction failed",
                                      kAbsOperation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "complex absolute value allocation failed",
                                  kAbsOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kAbsOperation);
    }
}

ExprResult abs(const ExprPtr& expression) {
    ComputationContext context;
    return abs(expression, context);
}

Result<ApproxReal> evalf(const SymbolicExpr& expression,
                         const NumericBindings& bindings,
                         ComputationContext& context) {
    return evaluate_numeric(expression, bindings, context);
}

Result<ApproxReal> evalf(const SymbolicExpr& expression,
                         const NumericBindings& bindings) {
    ComputationContext context;
    return evalf(expression, bindings, context);
}

Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                   const NumericBindings& bindings,
                                   ComputationContext& context) {
    try {
        if (!lamina::detail::node(expression)) {
            return complex_failure(CasErrc::InvalidArgument,
                                   "cannot evaluate an empty expression as complex",
                                   kEvalComplexOperation);
        }
        return evaluate_complex_node(lamina::detail::node(expression),
                                     bindings, context);
    } catch (const std::bad_alloc&) {
        return complex_failure(CasErrc::ResourceLimit,
                               "complex evaluation allocation failed",
                               kEvalComplexOperation);
    } catch (const std::exception& error) {
        return complex_failure(CasErrc::UnsupportedExpression, error.what(),
                               kEvalComplexOperation);
    }
}

Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                   const NumericBindings& bindings) {
    ComputationContext context;
    return eval_complex(expression, bindings, context);
}

SolveResult solve_set(const ExprPtr& equation,
                      const std::string& variable,
                      ComputationContext& context,
                      const SolveOptions& options) {
    if (!equation) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                    "equation cannot be null", "lsr.solve_set");
    }
    if (variable.empty()) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                    "solve variable cannot be empty",
                                    "lsr.solve_set");
    }
    return solve_dispatch_checked(equation, variable, context, options);
}

SolveResult solve_set(const ExprPtr& equation,
                      const std::string& variable,
                      const SolveOptions& options) {
    ComputationContext context;
    return solve_set(equation, variable, context, options);
}

ExprSetResult expr_set(std::vector<ExprPtr> elements) {
    try {
        return ExprSet::make(std::move(elements));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(CasErrc::ResourceLimit,
                                "set<Expr> allocation failed",
                                kExprSetOperation);
    } catch (const std::exception& error) {
        return expr_set_failure(CasErrc::InvalidArgument, error.what(),
                                kExprSetOperation);
    }
}

ExprSetResult solve_expr_set(const ExprPtr& equation,
                             const std::string& variable,
                             ComputationContext& context,
                             const SolveOptions& options) {
    auto closed_form = try_lsr_closed_form_rational_poly_roots(
        equation, variable, context);
    if (!closed_form) {
        return ExprSetResult::failure(closed_form.error());
    }
    if (closed_form.value()) {
        return ExprSetResult::success(std::move(*closed_form.value()));
    }

    auto solved = solve_set(equation, variable, context, options);
    if (!solved) {
        return ExprSetResult::failure(solved.error());
    }

    const auto& solution_set = solved.value();
    if (solution_set.kind() == SolutionSet::Kind::Empty) {
        return expr_set({});
    }
    if (solution_set.kind() != SolutionSet::Kind::Finite) {
        std::string reason = solution_set.reason();
        if (reason.empty()) {
            reason = "solution set is not a finite enumerable set<Expr>";
        }
        return expr_set_failure(CasErrc::Inconclusive, std::move(reason),
                                kSolveExprSetOperation);
    }

    std::vector<ExprPtr> elements;
    elements.reserve(solution_set.finite_solutions().size());
    for (const auto& solution : solution_set.finite_solutions()) {
        if (!solution.conditions.empty()) {
            return expr_set_failure(CasErrc::Inconclusive,
                                    "conditional finite solutions cannot be lowered to set<Expr>",
                                    kSolveExprSetOperation);
        }
        elements.push_back(solution.value);
    }
    return expr_set(std::move(elements));
}

ExprSetResult solve_expr_set(const ExprPtr& equation,
                             const std::string& variable,
                             const SolveOptions& options) {
    ComputationContext context;
    return solve_expr_set(equation, variable, context, options);
}

const char* error_name(CasErrc code) noexcept {
    switch (code) {
    case CasErrc::InvalidArgument:
        return "InvalidArgument";
    case CasErrc::ParseError:
        return "ParseError";
    case CasErrc::UnboundSymbol:
        return "UnboundSymbol";
    case CasErrc::DomainError:
        return "DomainError";
    case CasErrc::UnsupportedExpression:
        return "UnsupportedExpression";
    case CasErrc::Inconclusive:
        return "Inconclusive";
    case CasErrc::ResourceLimit:
        return "ResourceLimit";
    case CasErrc::Cancelled:
        return "Cancelled";
    case CasErrc::NumericFailure:
        return "NumericFailure";
    case CasErrc::InternalInvariant:
        return "InternalInvariant";
    }
    return "InternalInvariant";
}

const char* error_name(const CasError& error) noexcept {
    if (error.operation == kSymOperation &&
        error.code == CasErrc::InvalidArgument &&
        error.message.find("imaginary unit") != std::string::npos) {
        return "ImaginaryUnitReserved";
    }
    if (error.operation == kEvalComplexOperation &&
        error.code == CasErrc::UnboundSymbol) {
        return "ComplexEvalUnboundSymbol";
    }
    if (error.operation == kComplexOperation &&
        error.code == CasErrc::InvalidArgument) {
        return "ComplexTypeMismatch";
    }
    if (error.operation == kExprSetOperation &&
        error.code == CasErrc::InvalidArgument) {
        return "SetElementTypeMismatch";
    }
    if (error.operation == kSolveExprSetOperation &&
        error.code == CasErrc::Inconclusive) {
        return "SetResultInconclusive";
    }
    if (error.operation == kEquivalentOperation &&
        error.code == CasErrc::ResourceLimit) {
        return "EqvBudgetExceeded";
    }
    if (error.operation == kEquivalentProfileOperation &&
        error.code == CasErrc::UnsupportedExpression) {
        return "EqvRuleDisabled";
    }
    return error_name(error.code);
}

bool structurally_equal(const SymbolicExpr& lhs, const SymbolicExpr& rhs) {
    const auto& left = lamina::detail::node(lhs);
    const auto& right = lamina::detail::node(rhs);
    if (!left || !right) return left == right;
    return left->equals(*right);
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context,
                             const EqvOptions& options) {
    auto options_valid = validate_eqv_options(options);
    if (!options_valid) return Result<bool>::failure(options_valid.error());

    auto step = context.consume_steps(1, kEquivalentOperation);
    if (!step) return Result<bool>::failure(step.error());
    try {
        auto canonical_lhs = canonicalize_lsr_complex_product(lhs);
        auto canonical_rhs = canonicalize_lsr_complex_product(rhs);
        if (structurally_equal(*canonical_lhs, *canonical_rhs)) {
            return Result<bool>::success(true);
        }
        auto difference = SymbolicExpr::add(
            canonical_lhs,
            SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                   canonical_rhs));
        if (!difference) {
            return Result<bool>::failure(CasErrc::InternalInvariant,
                                         "equivalence difference construction failed",
                                         kEquivalentOperation);
        }
        if (difference->simplify()->is_zero()) {
            return Result<bool>::success(true);
        }

        auto polynomial_proof = prove_rational_polynomial_equivalence(
            difference, context, options);
        if (!polynomial_proof) {
            return Result<bool>::failure(polynomial_proof.error());
        }
        if (polynomial_proof.value()) {
            return Result<bool>::success(*polynomial_proof.value());
        }
        if (options.profile == EqvProfile::TrigBasic) {
            if (options.budget.max_rewrite_steps < 8) {
                return Result<bool>::failure(
                    CasErrc::ResourceLimit,
                    "equivalence rewrite budget exhausted before Trig-Basic normalization",
                    kEquivalentOperation);
            }
            auto trig_step = context.consume_steps(8, kEquivalentOperation);
            if (!trig_step) return Result<bool>::failure(trig_step.error());
            auto trig_lhs = rewrite_trig_basic_identity(lamina::detail::node(lhs));
            auto trig_rhs = rewrite_trig_basic_identity(lamina::detail::node(rhs));
            EqvOptions core_options = options;
            core_options.profile = EqvProfile::Core;
            return equivalent_core(*trig_lhs, *trig_rhs, context,
                                   core_options);
        }
        if (options.profile == EqvProfile::ExpLogBasic) {
            if (options.budget.max_rewrite_steps < 8) {
                return Result<bool>::failure(
                    CasErrc::ResourceLimit,
                    "equivalence rewrite budget exhausted before ExpLog-Basic normalization",
                    kEquivalentOperation);
            }
            auto exp_log_step = context.consume_steps(8, kEquivalentOperation);
            if (!exp_log_step) return Result<bool>::failure(exp_log_step.error());
            auto exp_log_lhs = rewrite_exp_log_basic_identity(
                lamina::detail::node(lhs));
            auto exp_log_rhs = rewrite_exp_log_basic_identity(
                lamina::detail::node(rhs));
            EqvOptions core_options = options;
            core_options.profile = EqvProfile::Core;
            return equivalent_core(*exp_log_lhs, *exp_log_rhs, context,
                                   core_options);
        }
        return Result<bool>::success(false);
    } catch (const std::bad_alloc&) {
        return Result<bool>::failure(CasErrc::ResourceLimit,
                                     "equivalence check allocation failed",
                                     kEquivalentOperation);
    } catch (const std::exception& error) {
        return Result<bool>::failure(CasErrc::Inconclusive, error.what(),
                                     kEquivalentOperation);
    }
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context) {
    return equivalent_core(lhs, rhs, context, EqvOptions{});
}

} // namespace lamina::lsr
