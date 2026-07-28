#include "lsr_expr.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <utility>

#include "symbolic_ast.hpp"

namespace lamina::lsr {
namespace {

constexpr const char* kSymOperation = "lsr.sym";
constexpr const char* kIntegerOperation = "lsr.integer";
constexpr const char* kRationalOperation = "lsr.rational";
constexpr const char* kApproxOperation = "lsr.approx_real";
constexpr const char* kImaginaryOperation = "lsr.imaginary_unit";
constexpr const char* kComplexOperation = "lsr.complex";
constexpr const char* kEquivalentOperation = "lsr.equivalent_core";
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
        return Result<ApproxComplex>::failure(real.error());
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

ExprPtr canonicalize_lsr_complex_product(const SymbolicExpr& expression) {
    const auto& node = lamina::detail::node(expression);
    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) {
        return lamina::detail::make_expression_ptr(node);
    }

    bool saw_complex = false;
    auto real = SymbolicExpr::number(1);
    auto imag = SymbolicExpr::number(0);

    for (const auto& operand : multiply->operands()) {
        ExprPtr factor_real;
        ExprPtr factor_imag;
        if (auto complex_operand = std::dynamic_pointer_cast<const ComplexNode>(operand)) {
            saw_complex = true;
            factor_real = lamina::detail::make_expression_ptr(complex_operand->real());
            factor_imag = lamina::detail::make_expression_ptr(complex_operand->imag());
        } else {
            factor_real = lamina::detail::make_expression_ptr(operand);
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
        if (!real) return Result<ApproxComplex>::failure(real.error());
        auto imag = evaluate_numeric(
            *lamina::detail::make_expression_ptr(complex_node->imag()),
            bindings, context);
        if (!imag) return Result<ApproxComplex>::failure(imag.error());
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
        if (!exponent) return Result<ApproxComplex>::failure(exponent.error());
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

bool structurally_equal(const SymbolicExpr& lhs, const SymbolicExpr& rhs) {
    const auto& left = lamina::detail::node(lhs);
    const auto& right = lamina::detail::node(rhs);
    if (!left || !right) return left == right;
    return left->equals(*right);
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context) {
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
        return Result<bool>::success(difference->simplify()->is_zero());
    } catch (const std::bad_alloc&) {
        return Result<bool>::failure(CasErrc::ResourceLimit,
                                     "equivalence check allocation failed",
                                     kEquivalentOperation);
    } catch (const std::exception& error) {
        return Result<bool>::failure(CasErrc::Inconclusive, error.what(),
                                     kEquivalentOperation);
    }
}

} // namespace lamina::lsr
