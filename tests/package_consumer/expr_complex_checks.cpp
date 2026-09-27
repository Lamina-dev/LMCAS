#include "assumption_context.hpp"
#include "expr.hpp"
#include "value.hpp"
#include "symbolic.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

using namespace LMCAS;

static int check_imaginary_product() {
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto i_squared = SymbolicExpr::multiply(i.value(), i.value());
    LMCAS::ComputationContext lsr_context;
    auto i_rule = LMCAS::equivalent_core(
        *i_squared, *SymbolicExpr::number(-1), lsr_context);
    if (!i_rule || !i_rule.value()) {
        std::cerr << "failed to prove LMCAS i*i == -1\n";
        return 11;
    }
    return 0;
}

static int check_imaginary_power() {
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto i_power_two = SymbolicExpr::power(i.value(), SymbolicExpr::number(2));
    LMCAS::ComputationContext lsr_power_context;
    auto i_power_rule = LMCAS::equivalent_core(
        *i_power_two, *SymbolicExpr::number(-1), lsr_power_context);
    if (!i_power_rule || !i_power_rule.value()) {
        std::cerr << "failed to prove LMCAS i^2 == -1\n";
        return 11;
    }
    return 0;
}

static int check_ordinary_i_identity() {
    auto legacy_i = SymbolicExpr::variable("i");
    auto legacy_i_squared = SymbolicExpr::multiply(legacy_i, legacy_i);
    LMCAS::ComputationContext legacy_i_context;
    auto legacy_i_rule = LMCAS::equivalent_core(
        *legacy_i_squared, *SymbolicExpr::number(-1), legacy_i_context);
    if (!legacy_i_rule || legacy_i_rule.value()) {
        std::cerr << "ordinary Expr i was treated as the imaginary unit\n";
        return 11;
    }
    return 0;
}

static int check_null_components() {
    auto null_real_complex =
        LMCAS::complex(nullptr, SymbolicExpr::number(1));
    auto null_imag_complex =
        LMCAS::complex(SymbolicExpr::number(0), nullptr);
    if (null_real_complex || null_imag_complex ||
        std::string(LMCAS::error_name(null_real_complex.error())) !=
            "ComplexTypeMismatch" ||
        std::string(LMCAS::error_name(null_imag_complex.error())) !=
            "ComplexTypeMismatch") {
        std::cerr << "failed to expose LMCAS complex type diagnostics\n";
        return 10;
    }
    return 0;
}

static int check_complex_lowering() {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0),
                                       SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    auto lowered_complex = LMCAS::eval_complex(*three_plus_four_i);
    if (!lowered_complex || !lowered_complex.value().is_finite() ||
        lowered_complex.value().real.value != 3.0 ||
        lowered_complex.value().imag.value != 4.0) {
        std::cerr << "failed to explicitly lower LMCAS Expr to complex\n";
        return 13;
    }
    return 0;
}

static int check_ordinary_i_lowering() {
    auto legacy_i = SymbolicExpr::variable("i");
    auto lowered_ordinary_i = LMCAS::eval_complex(*legacy_i);
    if (lowered_ordinary_i ||
        lowered_ordinary_i.error().code != LMCAS::CasErrc::UnboundSymbol) {
        std::cerr << "ordinary lowercase i complex lowering returned ";
        if (lowered_ordinary_i) {
            std::cerr << "a value\n";
        } else {
            std::cerr << LMCAS::error_name(lowered_ordinary_i.error())
                      << "\n";
        }
        return 13;
    }
    return 0;
}

static int check_complex_arithmetic() {
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto ordinary_multiply_complex = SymbolicExpr::add(
        SymbolicExpr::number(3),
        SymbolicExpr::multiply(SymbolicExpr::number(4), i.value()));
    auto lowered_ordinary_multiply =
        LMCAS::eval_complex(*ordinary_multiply_complex);
    if (!lowered_ordinary_multiply ||
        !lowered_ordinary_multiply.value().is_finite() ||
        lowered_ordinary_multiply.value().real.value != 3.0 ||
        lowered_ordinary_multiply.value().imag.value != 4.0) {
        std::cerr << "failed to lower LMCAS 3 + 4 * i form to complex\n";
        return 13;
    }
    return 0;
}

static int check_zero_reciprocal() {
    auto zero_inverse_complex = LMCAS::eval_complex(
        *SymbolicExpr::power(SymbolicExpr::number(0), SymbolicExpr::number(-1)));
    if (zero_inverse_complex ||
        std::string(LMCAS::error_name(zero_inverse_complex.error())) !=
            "DomainError") {
        std::cerr << "failed to reject LMCAS complex reciprocal of zero\n";
        return 13;
    }
    return 0;
}

static int check_nonfinite_exponent() {
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto nonfinite_complex_power = LMCAS::eval_complex(
        *SymbolicExpr::power(i.value(), SymbolicExpr::infinity()));
    if (nonfinite_complex_power ||
        std::string(LMCAS::error_name(
            nonfinite_complex_power.error())) != "NumericFailure") {
        std::cerr << "failed to reject LMCAS complex non-finite exponent\n";
        return 13;
    }
    return 0;
}

static int check_complex_budget() {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0),
                                       SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    LMCAS::ResourceLimits exhausted_complex_limits;
    exhausted_complex_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_complex_context(exhausted_complex_limits);
    auto exhausted_complex = LMCAS::eval_complex(
        *three_plus_four_i, {}, exhausted_complex_context);
    if (exhausted_complex ||
        std::string(LMCAS::error_name(exhausted_complex.error())) !=
            "ResourceLimit") {
        std::cerr << "failed to expose LMCAS complex resource limits\n";
        return 13;
    }
    return 0;
}

static int check_complex_components() {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0),
                                       SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    auto lsr_real = LMCAS::real(three_plus_four_i);
    auto lsr_imag = LMCAS::imag(three_plus_four_i);
    auto lsr_conj = LMCAS::conj(three_plus_four_i);
    auto lsr_abs = LMCAS::abs(three_plus_four_i);
    auto expected_conj = LMCAS::complex(SymbolicExpr::number(3),
                                              SymbolicExpr::number(-4));
    if (!lsr_real || !lsr_imag || !lsr_conj || !lsr_abs || !expected_conj ||
        !LMCAS::structurally_equal(*lsr_real.value(),
                                         *SymbolicExpr::number(3)) ||
        !LMCAS::structurally_equal(*lsr_imag.value(),
                                         *SymbolicExpr::number(4)) ||
        !LMCAS::structurally_equal(*lsr_conj.value(),
                                         *expected_conj.value())) {
        std::cerr << "failed to call LMCAS complex part facade\n";
        return 14;
    }
    return 0;
}

static int check_null_complex_functions() {
    auto lsr_null_real = LMCAS::real(nullptr);
    auto lsr_null_imag = LMCAS::imag(nullptr);
    auto lsr_null_conj = LMCAS::conj(nullptr);
    auto lsr_null_abs = LMCAS::abs(nullptr);
    if (lsr_null_real || lsr_null_imag || lsr_null_conj || lsr_null_abs ||
        std::string(LMCAS::error_name(lsr_null_real.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_null_imag.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_null_conj.error())) !=
            "InvalidArgument" ||
        std::string(LMCAS::error_name(lsr_null_abs.error())) !=
            "InvalidArgument") {
        std::cerr << "failed to reject null LMCAS complex facade inputs\n";
        return 14;
    }
    return 0;
}

static int check_real_promotion() {
    auto lsr_real_value = SymbolicExpr::number(-5);
    auto lsr_real_value_real = LMCAS::real(lsr_real_value);
    auto lsr_real_value_imag = LMCAS::imag(lsr_real_value);
    auto lsr_real_value_conj = LMCAS::conj(lsr_real_value);
    auto lsr_real_value_abs = LMCAS::abs(lsr_real_value);
    auto lsr_real_value_abs_eval =
        lsr_real_value_abs ? LMCAS::evalf(*lsr_real_value_abs.value())
                           : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                 LMCAS::CasErrc::InternalInvariant,
                                 "abs(-5) construction failed", "consumer");
    if (!lsr_real_value_real || !lsr_real_value_imag ||
        !lsr_real_value_conj || !lsr_real_value_abs_eval) {
        std::cerr << "failed to promote real values through LMCAS complex facade\n";
        return 14;
    }
    if (!LMCAS::structurally_equal(*lsr_real_value_real.value(),
                                         *lsr_real_value) ||
        !LMCAS::structurally_equal(*lsr_real_value_imag.value(),
                                         *SymbolicExpr::number(0)) ||
        !LMCAS::structurally_equal(*lsr_real_value_conj.value(),
                                         *lsr_real_value)) {
        std::cerr << "failed to promote real values through LMCAS complex facade\n";
        return 14;
    }
    if (lsr_real_value_abs_eval.value().value < 4.999999999999 ||
        lsr_real_value_abs_eval.value().value > 5.000000000001) {
        std::cerr << "failed to promote real values through LMCAS complex facade\n";
        return 14;
    }
    return 0;
}

static int check_shifted_complex_roots() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_shifted_complex_roots = LMCAS::solve_expr_set(
        SymbolicExpr::add(
            SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                              SymbolicExpr::multiply(SymbolicExpr::number(2), x)),
            SymbolicExpr::number(2)),
        "x");
    if (!lsr_shifted_complex_roots ||
        lsr_shifted_complex_roots.value().size() != 2) {
        std::cerr << "failed to lower shifted LMCAS complex roots to set<Expr>\n";
        return 12;
    }
    bool saw_negative_one_plus_i = false;
    bool saw_negative_one_minus_i = false;
    for (const auto& root : lsr_shifted_complex_roots.value().elements()) {
        auto lowered_root = LMCAS::eval_complex(*root);
        if (!lowered_root || !lowered_root.value().is_finite()) {
            std::cerr << "failed to evaluate shifted LMCAS complex root\n";
            return 12;
        }
        if (lowered_root.value().real.value == -1.0 &&
            lowered_root.value().imag.value == 1.0) {
            saw_negative_one_plus_i = true;
        }
        if (lowered_root.value().real.value == -1.0 &&
            lowered_root.value().imag.value == -1.0) {
            saw_negative_one_minus_i = true;
        }
    }
    if (!saw_negative_one_plus_i || !saw_negative_one_minus_i) {
        std::cerr << "failed to preserve shifted LMCAS complex root components\n";
        return 12;
    }
    return 0;
}

int run_expr_complex_consumer_checks() {
    const auto checks = {
        check_imaginary_product,
        check_imaginary_power,
        check_ordinary_i_identity,
        check_null_components,
        check_complex_lowering,
        check_ordinary_i_lowering,
        check_complex_arithmetic,
        check_zero_reciprocal,
        check_nonfinite_exponent,
        check_complex_budget,
        check_complex_components,
        check_null_complex_functions,
        check_real_promotion,
        check_shifted_complex_roots,
    };
    for (const auto check : checks) {
        if (const int status = check(); status != 0) {
            return status;
        }
    }
    return 0;
}
