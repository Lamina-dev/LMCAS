#define _USE_MATH_DEFINES
#include "../include/symbolic_complex.hpp"
#include "../include/symbolic.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <variant>
#include <vector>

namespace lamina {

namespace {

constexpr const char* kComplexAddOperation = "complex_add";
constexpr const char* kComplexSubOperation = "complex_sub";
constexpr const char* kComplexMulOperation = "complex_mul";
constexpr const char* kComplexDivOperation = "complex_div";
constexpr const char* kComplexConjOperation = "complex_conj";
constexpr const char* kComplexAbsOperation = "complex_abs";
constexpr const char* kComplexArgOperation = "complex_arg";
constexpr const char* kComplexPolarOperation = "complex_polar_form";
constexpr const char* kComplexNthRootOperation = "solve_complex_nth_root";
constexpr const char* kComplexQuadraticOperation = "solve_complex_quadratic";
constexpr const char* kComplexLocusCircleOperation = "complex_locus_circle";
constexpr const char* kComplexLocusBisectorOperation =
    "complex_locus_perpendicular_bisector";

Result<void> validate_complex_symbolic(const ComplexSymbolic& z,
                                       const char* name,
                                       ComputationContext& context,
                                       const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!z.real || !lamina::detail::node(z.real)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, std::string(name) + ".real must not be null", operation);
    }
    if (!z.imag || !lamina::detail::node(z.imag)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, std::string(name) + ".imag must not be null", operation);
    }
    return Result<void>::success();
}

Result<void> validate_expr(const std::shared_ptr<SymbolicExpr>& expr,
                           const char* name,
                           ComputationContext& context,
                           const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!expr || !lamina::detail::node(expr)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, std::string(name) + " must not be null", operation);
    }
    return Result<void>::success();
}

bool is_explicit_approx_real(const std::shared_ptr<SymbolicExpr>& expr, double& out) {
    if (!expr || !lamina::detail::node(expr)) return false;
    auto number = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr));
    if (!number || !std::holds_alternative<lmmc_real_t>(number->value())) return false;
    out = static_cast<double>(std::get<lmmc_real_t>(number->value()));
    return std::isfinite(out);
}

template <typename T, typename F>
Result<T> checked_construct(F&& f, const std::string& operation) {
    try {
        return Result<T>::success(f());
    } catch (const detail::ResultPropagation& propagation) {
        return Result<T>::failure(propagation.error());
    } catch (const std::bad_alloc&) {
        return Result<T>::failure(CasErrc::ResourceLimit,
                                  "symbolic complex allocation failed", operation);
    } catch (const std::exception& ex) {
        return Result<T>::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
}

} // namespace

ComplexSymbolicResult complex_add_checked(const ComplexSymbolic& a, const ComplexSymbolic& b,
                                          ComputationContext& context) {
    auto va = validate_complex_symbolic(a, "a", context, kComplexAddOperation);
    if (!va) return ComplexSymbolicResult::failure(va.error());
    auto vb = validate_complex_symbolic(b, "b", context, kComplexAddOperation);
    if (!vb) return ComplexSymbolicResult::failure(vb.error());
    auto budget = context.consume_steps(4, kComplexAddOperation);
    if (!budget) return ComplexSymbolicResult::failure(budget.error());
    return checked_construct<ComplexSymbolic>([&]() {
        return make_complex(SymbolicExpr::add(a.real, b.real), SymbolicExpr::add(a.imag, b.imag));
    }, kComplexAddOperation);
}
ComplexSymbolicResult complex_add_checked(const ComplexSymbolic& a, const ComplexSymbolic& b) {
    ComputationContext context;
    return complex_add_checked(a, b, context);
}

ComplexSymbolicResult complex_sub_checked(const ComplexSymbolic& a, const ComplexSymbolic& b,
                                          ComputationContext& context) {
    auto va = validate_complex_symbolic(a, "a", context, kComplexSubOperation);
    if (!va) return ComplexSymbolicResult::failure(va.error());
    auto vb = validate_complex_symbolic(b, "b", context, kComplexSubOperation);
    if (!vb) return ComplexSymbolicResult::failure(vb.error());
    auto budget = context.consume_steps(6, kComplexSubOperation);
    if (!budget) return ComplexSymbolicResult::failure(budget.error());
    return checked_construct<ComplexSymbolic>([&]() {
        return make_complex(
            SymbolicExpr::add(a.real, SymbolicExpr::multiply(SymbolicExpr::number(-1), b.real)),
            SymbolicExpr::add(a.imag, SymbolicExpr::multiply(SymbolicExpr::number(-1), b.imag)));
    }, kComplexSubOperation);
}
ComplexSymbolicResult complex_sub_checked(const ComplexSymbolic& a, const ComplexSymbolic& b) {
    ComputationContext context;
    return complex_sub_checked(a, b, context);
}

ComplexSymbolicResult complex_mul_checked(const ComplexSymbolic& a, const ComplexSymbolic& b,
                                          ComputationContext& context) {
    auto va = validate_complex_symbolic(a, "a", context, kComplexMulOperation);
    if (!va) return ComplexSymbolicResult::failure(va.error());
    auto vb = validate_complex_symbolic(b, "b", context, kComplexMulOperation);
    if (!vb) return ComplexSymbolicResult::failure(vb.error());
    auto budget = context.consume_steps(10, kComplexMulOperation);
    if (!budget) return ComplexSymbolicResult::failure(budget.error());
    return checked_construct<ComplexSymbolic>([&]() {
        auto real = SymbolicExpr::add(
            SymbolicExpr::multiply(a.real, b.real),
            SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a.imag, b.imag)));
        auto imag = SymbolicExpr::add(SymbolicExpr::multiply(a.real, b.imag),
                                      SymbolicExpr::multiply(a.imag, b.real));
        return make_complex(real, imag);
    }, kComplexMulOperation);
}
ComplexSymbolicResult complex_mul_checked(const ComplexSymbolic& a, const ComplexSymbolic& b) {
    ComputationContext context;
    return complex_mul_checked(a, b, context);
}

ComplexSymbolicResult complex_div_checked(const ComplexSymbolic& a, const ComplexSymbolic& b,
                                          ComputationContext& context) {
    auto va = validate_complex_symbolic(a, "a", context, kComplexDivOperation);
    if (!va) return ComplexSymbolicResult::failure(va.error());
    auto vb = validate_complex_symbolic(b, "b", context, kComplexDivOperation);
    if (!vb) return ComplexSymbolicResult::failure(vb.error());
    if (lamina::detail::node(b.real)->is_zero() && lamina::detail::node(b.imag)->is_zero()) {
        return ComplexSymbolicResult::failure(
            CasErrc::DomainError, "complex division by zero", kComplexDivOperation);
    }
    auto budget = context.consume_steps(16, kComplexDivOperation);
    if (!budget) return ComplexSymbolicResult::failure(budget.error());
    return checked_construct<ComplexSymbolic>([&]() {
        auto denom = SymbolicExpr::add(SymbolicExpr::power(b.real, SymbolicExpr::number(2)),
                                       SymbolicExpr::power(b.imag, SymbolicExpr::number(2)));
        auto real = SymbolicExpr::divide(
            SymbolicExpr::add(SymbolicExpr::multiply(a.real, b.real), SymbolicExpr::multiply(a.imag, b.imag)),
            denom);
        auto imag = SymbolicExpr::divide(
            SymbolicExpr::add(SymbolicExpr::multiply(a.imag, b.real),
                              SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                                     SymbolicExpr::multiply(a.real, b.imag))),
            denom);
        return make_complex(real, imag);
    }, kComplexDivOperation);
}
ComplexSymbolicResult complex_div_checked(const ComplexSymbolic& a, const ComplexSymbolic& b) {
    ComputationContext context;
    return complex_div_checked(a, b, context);
}

ComplexSymbolicResult complex_conj_checked(const ComplexSymbolic& z, ComputationContext& context) {
    auto vz = validate_complex_symbolic(z, "z", context, kComplexConjOperation);
    if (!vz) return ComplexSymbolicResult::failure(vz.error());
    auto budget = context.consume_steps(3, kComplexConjOperation);
    if (!budget) return ComplexSymbolicResult::failure(budget.error());
    return checked_construct<ComplexSymbolic>([&]() {
        return make_complex(z.real, SymbolicExpr::multiply(SymbolicExpr::number(-1), z.imag));
    }, kComplexConjOperation);
}
ComplexSymbolicResult complex_conj_checked(const ComplexSymbolic& z) {
    ComputationContext context;
    return complex_conj_checked(z, context);
}

ExpressionResult complex_abs_checked(const ComplexSymbolic& z, ComputationContext& context) {
    auto vz = validate_complex_symbolic(z, "z", context, kComplexAbsOperation);
    if (!vz) return ExpressionResult::failure(vz.error());
    auto budget = context.consume_steps(8, kComplexAbsOperation);
    if (!budget) return ExpressionResult::failure(budget.error());
    return checked_construct<std::shared_ptr<SymbolicExpr>>([&]() {
        return SymbolicExpr::sqrt(SymbolicExpr::add(
            SymbolicExpr::power(z.real, SymbolicExpr::number(2)),
            SymbolicExpr::power(z.imag, SymbolicExpr::number(2))));
    }, kComplexAbsOperation);
}
ExpressionResult complex_abs_checked(const ComplexSymbolic& z) {
    ComputationContext context;
    return complex_abs_checked(z, context);
}

ExpressionResult complex_arg_checked(const ComplexSymbolic& z, ComputationContext& context) {
    auto vz = validate_complex_symbolic(z, "z", context, kComplexArgOperation);
    if (!vz) return ExpressionResult::failure(vz.error());
    auto budget = context.consume_steps(4, kComplexArgOperation);
    if (!budget) return ExpressionResult::failure(budget.error());
    return checked_construct<std::shared_ptr<SymbolicExpr>>([&]() {
        return SymbolicExpr::atan2(z.imag, z.real);
    }, kComplexArgOperation);
}
ExpressionResult complex_arg_checked(const ComplexSymbolic& z) {
    ComputationContext context;
    return complex_arg_checked(z, context);
}

ComplexSymbolicResult complex_exp_form_checked(std::shared_ptr<SymbolicExpr> r,
                                               std::shared_ptr<SymbolicExpr> theta,
                                               ComputationContext& context) {
    auto vr = validate_expr(r, "r", context, kComplexPolarOperation);
    if (!vr) return ComplexSymbolicResult::failure(vr.error());
    auto vt = validate_expr(theta, "theta", context, kComplexPolarOperation);
    if (!vt) return ComplexSymbolicResult::failure(vt.error());
    auto budget = context.consume_steps(6, kComplexPolarOperation);
    if (!budget) return ComplexSymbolicResult::failure(budget.error());
    return checked_construct<ComplexSymbolic>([&]() {
        auto real = SymbolicExpr::multiply(r, SymbolicExpr::cos(theta));
        auto imag = SymbolicExpr::multiply(r, SymbolicExpr::sin(theta));
        return make_complex(real, imag);
    }, kComplexPolarOperation);
}
ComplexSymbolicResult complex_exp_form_checked(std::shared_ptr<SymbolicExpr> r,
                                               std::shared_ptr<SymbolicExpr> theta) {
    ComputationContext context;
    return complex_exp_form_checked(std::move(r), std::move(theta), context);
}

ComplexSymbolicResult complex_trig_form_checked(std::shared_ptr<SymbolicExpr> r,
                                                std::shared_ptr<SymbolicExpr> theta,
                                                ComputationContext& context) {
    return complex_exp_form_checked(std::move(r), std::move(theta), context);
}
ComplexSymbolicResult complex_trig_form_checked(std::shared_ptr<SymbolicExpr> r,
                                                std::shared_ptr<SymbolicExpr> theta) {
    ComputationContext context;
    return complex_trig_form_checked(std::move(r), std::move(theta), context);
}

ComplexRootsResult solve_complex_nth_root_checked(std::shared_ptr<SymbolicExpr> c, int n,
                                                  ComputationContext& context) {
    auto vc = validate_expr(c, "c", context, kComplexNthRootOperation);
    if (!vc) return ComplexRootsResult::failure(vc.error());
    if (n <= 0) {
        return ComplexRootsResult::failure(
            CasErrc::InvalidArgument, "root degree must be positive", kComplexNthRootOperation);
    }
    auto budget = context.consume_steps(static_cast<std::size_t>(n) * 8 + 8,
                                        kComplexNthRootOperation);
    if (!budget) return ComplexRootsResult::failure(budget.error());

    double val = 0.0;
    if (!is_explicit_approx_real(c, val)) {
        return ComplexRootsResult::failure(
            CasErrc::Inconclusive,
            "checked complex nth roots currently require an explicit approximate real input",
            kComplexNthRootOperation);
    }

    return checked_construct<std::vector<ComplexSymbolic>>([&]() {
        double r = std::abs(val);
        double theta = std::arg(std::complex<double>(val, 0.0));
        std::vector<ComplexSymbolic> roots;
        roots.reserve(static_cast<std::size_t>(n));
        for (int k = 0; k < n; ++k) {
            double root_r = std::pow(r, 1.0 / n);
            double root_theta = (theta + 2 * LMMC_CONST_PI * k) / n;
            roots.push_back(lamina::detail::propagate_result(complex_exp_form_checked(
                SymbolicExpr::number(root_r),
                SymbolicExpr::number(root_theta), context)));
        }
        return roots;
    }, kComplexNthRootOperation);
}
ComplexRootsResult solve_complex_nth_root_checked(std::shared_ptr<SymbolicExpr> c, int n) {
    ComputationContext context;
    return solve_complex_nth_root_checked(std::move(c), n, context);
}

ComplexRootsResult solve_complex_quadratic_checked(
    std::shared_ptr<SymbolicExpr> a, std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c, ComputationContext& context) {
    auto va = validate_expr(a, "a", context, kComplexQuadraticOperation);
    if (!va) return ComplexRootsResult::failure(va.error());
    auto vb = validate_expr(b, "b", context, kComplexQuadraticOperation);
    if (!vb) return ComplexRootsResult::failure(vb.error());
    auto vc = validate_expr(c, "c", context, kComplexQuadraticOperation);
    if (!vc) return ComplexRootsResult::failure(vc.error());
    return checked_construct<std::vector<ComplexSymbolic>>([&]() {
        auto four_ac = SymbolicExpr::multiply(
            SymbolicExpr::number(4), SymbolicExpr::multiply(a, c));
        auto delta = SymbolicExpr::add(
            SymbolicExpr::power(b, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(-1), four_ac));
        auto sqrt_delta = SymbolicExpr::sqrt(delta);
        auto denominator = SymbolicExpr::multiply(SymbolicExpr::number(2), a);
        auto real_denominator =
            make_complex(denominator, SymbolicExpr::number(0));
        auto first_numerator = make_complex(
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(-1), b),
                sqrt_delta),
            SymbolicExpr::number(0));
        auto second_numerator = make_complex(
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(-1), b),
                SymbolicExpr::multiply(SymbolicExpr::number(-1), sqrt_delta)),
            SymbolicExpr::number(0));
        auto first = detail::propagate_result(
            complex_div_checked(first_numerator, real_denominator, context));
        auto second = detail::propagate_result(
            complex_div_checked(second_numerator, real_denominator, context));
        return std::vector<ComplexSymbolic>{std::move(first), std::move(second)};
    }, kComplexQuadraticOperation);
}

ComplexRootsResult solve_complex_quadratic_checked(
    std::shared_ptr<SymbolicExpr> a, std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c) {
    ComputationContext context;
    return solve_complex_quadratic_checked(
        std::move(a), std::move(b), std::move(c), context);
}

ExpressionResult complex_locus_circle_checked(
    const ComplexSymbolic& a, std::shared_ptr<SymbolicExpr> r,
    const std::string& z_var, ComputationContext& context) {
    auto va = validate_complex_symbolic(
        a, "a", context, kComplexLocusCircleOperation);
    if (!va) return ExpressionResult::failure(va.error());
    auto vr = validate_expr(r, "r", context, kComplexLocusCircleOperation);
    if (!vr) return ExpressionResult::failure(vr.error());
    if (z_var.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument, "complex variable must not be empty",
            kComplexLocusCircleOperation);
    }
    return checked_construct<std::shared_ptr<SymbolicExpr>>([&]() {
        auto z = SymbolicExpr::variable(z_var);
        auto difference = detail::propagate_result(complex_sub_checked(
            make_complex(z, SymbolicExpr::number(0)), a, context));
        auto magnitude =
            detail::propagate_result(complex_abs_checked(difference, context));
        return SymbolicExpr::eq(magnitude, r);
    }, kComplexLocusCircleOperation);
}

ExpressionResult complex_locus_circle_checked(
    const ComplexSymbolic& a, std::shared_ptr<SymbolicExpr> r,
    const std::string& z_var) {
    ComputationContext context;
    return complex_locus_circle_checked(a, std::move(r), z_var, context);
}

ExpressionResult complex_locus_perpendicular_bisector_checked(
    const ComplexSymbolic& a, const ComplexSymbolic& b,
    const std::string& z_var, ComputationContext& context) {
    auto va = validate_complex_symbolic(
        a, "a", context, kComplexLocusBisectorOperation);
    if (!va) return ExpressionResult::failure(va.error());
    auto vb = validate_complex_symbolic(
        b, "b", context, kComplexLocusBisectorOperation);
    if (!vb) return ExpressionResult::failure(vb.error());
    if (z_var.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument, "complex variable must not be empty",
            kComplexLocusBisectorOperation);
    }
    return checked_construct<std::shared_ptr<SymbolicExpr>>([&]() {
        auto z = SymbolicExpr::variable(z_var);
        auto z_as_complex = make_complex(z, SymbolicExpr::number(0));
        auto from_a = detail::propagate_result(
            complex_sub_checked(z_as_complex, a, context));
        auto from_b = detail::propagate_result(
            complex_sub_checked(z_as_complex, b, context));
        auto distance_a =
            detail::propagate_result(complex_abs_checked(from_a, context));
        auto distance_b =
            detail::propagate_result(complex_abs_checked(from_b, context));
        return SymbolicExpr::eq(distance_a, distance_b);
    }, kComplexLocusBisectorOperation);
}

ExpressionResult complex_locus_perpendicular_bisector_checked(
    const ComplexSymbolic& a, const ComplexSymbolic& b,
    const std::string& z_var) {
    ComputationContext context;
    return complex_locus_perpendicular_bisector_checked(
        a, b, z_var, context);
}

}
