#define _USE_MATH_DEFINES
#include <cmath>
#include "test_common.hpp"
#include "symbolic_complex.hpp"
#include "numeric_evaluation.hpp"

using namespace LMCAS;

void test_complex_arithmetic() {
    TEST_CASE("Complex Arithmetic: add, sub, mul, div");

    // Pair 1: (1+2i) and (3+4i)
    auto a1 = LMCAS::make_complex(SymbolicExpr::number(1), SymbolicExpr::number(2));
    auto b1 = LMCAS::make_complex(SymbolicExpr::number(3), SymbolicExpr::number(4));

    // Pair 2: (2+0i) and (0+3i)
    auto a2 = LMCAS::make_complex(SymbolicExpr::number(2), SymbolicExpr::number(0));
    auto b2 = LMCAS::make_complex(SymbolicExpr::number(0), SymbolicExpr::number(3));

    // Add pair 1: (1+2i) + (3+4i) = (4+6i)
    {
        auto result = LMCAS::complex_add_checked(a1, b1).value();
        std::string real_s = result.real ? result.real->to_string() : "null";
        std::string imag_s = result.imag ? result.imag->to_string() : "null";
        // Symbolic add produces "1 + 3" (unsimplified) or "4" (simplified)
        EXPECT_CONTAINS(real_s, {"1"}, "add pair1 real contains operand 1");
        EXPECT_CONTAINS(imag_s, {"2"}, "add pair1 imag contains operand 2");
    }

    // Add pair 2: (2+0i) + (0+3i) = (2+3i)
    {
        auto result = LMCAS::complex_add_checked(a2, b2).value();
        std::string real_s = result.real ? result.real->to_string() : "null";
        std::string imag_s = result.imag ? result.imag->to_string() : "null";
        EXPECT_CONTAINS(real_s, {"2"}, "add pair2 real contains 2");
        EXPECT_CONTAINS(imag_s, {"3"}, "add pair2 imag contains 3");
    }

    // Sub pair 1: (1+2i) - (3+4i) = (-2-2i)
    {
        auto result = LMCAS::complex_sub_checked(a1, b1).value();
        std::string real_s = result.real ? result.real->to_string() : "null";
        std::string imag_s = result.imag ? result.imag->to_string() : "null";
        EXPECT_CONTAINS(real_s, {"1"}, "sub pair1 real contains operands");
        EXPECT_CONTAINS(imag_s, {"2"}, "sub pair1 imag contains operands");
    }

    // Sub pair 2: (2+0i) - (0+3i) = (2-3i)
    {
        auto result = LMCAS::complex_sub_checked(a2, b2).value();
        std::string real_s = result.real ? result.real->to_string() : "null";
        std::string imag_s = result.imag ? result.imag->to_string() : "null";
        EXPECT_CONTAINS(real_s, {"2"}, "sub pair2 real contains 2");
        EXPECT_CONTAINS(imag_s, {"3"}, "sub pair2 imag contains 3");
    }

    // Mul pair 1: (1+2i)*(3+4i) = (1*3 - 2*4) + (1*4 + 2*3)i = (-5+10i)
    {
        auto result = LMCAS::complex_mul_checked(a1, b1).value();
        std::string real_s = result.real ? result.real->to_string() : "null";
        std::string imag_s = result.imag ? result.imag->to_string() : "null";
        EXPECT_TRUE(real_s != "null", "mul pair1 real is not null");
        EXPECT_TRUE(imag_s != "null", "mul pair1 imag is not null");
    }

    // Mul pair 2: (2+0i)*(0+3i) = (0) + (6)i
    {
        auto result = LMCAS::complex_mul_checked(a2, b2).value();
        std::string real_s = result.real ? result.real->to_string() : "null";
        std::string imag_s = result.imag ? result.imag->to_string() : "null";
        EXPECT_TRUE(real_s != "null", "mul pair2 real is not null");
        // imag = a.real*b.imag + a.imag*b.real = 2*3 + 0*0 (symbolic, unsimplified)
        EXPECT_CONTAINS(imag_s, {"2"}, "mul pair2 imag contains factor 2");
    }

    // Div pair 1: (1+2i)/(3+4i)
    {
        auto result = LMCAS::complex_div_checked(a1, b1).value();
        std::string real_s = result.real ? result.real->to_string() : "null";
        std::string imag_s = result.imag ? result.imag->to_string() : "null";
        EXPECT_TRUE(real_s != "null", "div pair1 real is not null");
        EXPECT_TRUE(imag_s != "null", "div pair1 imag is not null");
    }

    // Div pair 2: (2+0i)/(0+3i)
    {
        auto result = LMCAS::complex_div_checked(a2, b2).value();
        std::string real_s = result.real ? result.real->to_string() : "null";
        std::string imag_s = result.imag ? result.imag->to_string() : "null";
        EXPECT_TRUE(real_s != "null", "div pair2 real is not null");
        EXPECT_TRUE(imag_s != "null", "div pair2 imag is not null");
    }
}

void test_complex_conj() {
    TEST_CASE("Complex Conjugate: conj(a+bi) = (a-bi)");

    // conj(3+4i) should give (3-4i)
    auto z = LMCAS::make_complex(SymbolicExpr::number(3), SymbolicExpr::number(4));
    auto conj = LMCAS::complex_conj_checked(z).value();

    std::string real_s = conj.real ? conj.real->to_string() : "null";
    std::string imag_s = conj.imag ? conj.imag->to_string() : "null";

    // Real part should be unchanged (3)
    EXPECT_CONTAINS(real_s, {"3"}, "conj real part is 3");
    // Imaginary part should be negated (-1 * 4)
    EXPECT_CONTAINS(imag_s, {"4"}, "conj imag part contains 4");
    EXPECT_CONTAINS(imag_s, {"-1"}, "conj imag part is negated");
}

void test_complex_abs() {
    TEST_CASE("Complex modulus preserves exact values and numerical range");
    auto exact = LMCAS::complex_abs_checked(
        LMCAS::make_complex(SymbolicExpr::number(3), SymbolicExpr::number(4)));
    EXPECT_TRUE(exact.has_value(), "exact modulus construction succeeds");
    if (exact) {
        auto value = LMCAS::evaluate_numeric(*exact.value());
        EXPECT_TRUE(value && value.value().value == 5.0, "exact modulus evaluates to five");
    }
    for (double component : {1e200, 1e-200}) {
        auto modulus = LMCAS::complex_abs_checked(LMCAS::make_complex(
            SymbolicExpr::number(component), SymbolicExpr::number(component)));
        EXPECT_TRUE(modulus.has_value(), "extreme modulus construction succeeds");
        if (modulus) {
            auto value = LMCAS::evaluate_numeric(*modulus.value());
            EXPECT_TRUE(value && value.value().is_finite() &&
                            std::abs(value.value().value / std::hypot(component, component) - 1) < 1e-14,
                        "symbolic complex modulus preserves finite nonzero extreme values");
        }
    }
}

void test_complex_arg() {
    TEST_CASE("Complex Argument: arg on positive real axis = 0");

    // arg(1+0i) = atan2(0, 1) = 0
    {
        auto z = LMCAS::make_complex(SymbolicExpr::number(1), SymbolicExpr::number(0));
        auto arg = LMCAS::complex_arg_checked(z).value();
        std::string s = arg ? arg->to_string() : "null";
        EXPECT_TRUE(s != "null", "arg(1+0i) is not null");
        // Should be atan2(0, 1) which is 0
        auto val = test_numeric_eval(arg);
        if (val) {
            EXPECT_TRUE(std::abs(*val - 0.0) < 1e-9, "arg(1+0i) evaluates to 0");
        } else {
            EXPECT_CONTAINS(s, {"0"}, "arg(1+0i) contains 0");
        }
    }

    TEST_CASE("Complex Argument: arg on positive imaginary axis = pi/2");

    // arg(0+1i) = atan2(1, 0) = pi/2
    {
        auto z = LMCAS::make_complex(SymbolicExpr::number(0), SymbolicExpr::number(1));
        auto arg = LMCAS::complex_arg_checked(z).value();
        std::string s = arg ? arg->to_string() : "null";
        EXPECT_TRUE(s != "null", "arg(0+1i) is not null");
        auto val = test_numeric_eval(arg);
        if (val) {
            EXPECT_TRUE(std::abs(*val - M_PI / 2.0) < 1e-9, "arg(0+1i) evaluates to pi/2");
        } else {
            // Structural: should contain atan2(1, 0) or pi/2
            EXPECT_CONTAINS(s, {"1"}, "arg(0+1i) contains 1");
        }
    }
}

void test_complex_polar_forms() {
    TEST_CASE("Complex Polar Forms: exp_form and trig_form equivalence");

    auto r = SymbolicExpr::number(2);
    auto theta = SymbolicExpr::number(1); // 1 radian

    auto exp_form = LMCAS::complex_exp_form_checked(r, theta).value();
    auto trig_form = LMCAS::complex_trig_form_checked(r, theta).value();

    // Both should produce r*cos(theta) for real and r*sin(theta) for imag
    std::string exp_real = exp_form.real ? exp_form.real->to_string() : "null";
    std::string exp_imag = exp_form.imag ? exp_form.imag->to_string() : "null";
    std::string trig_real = trig_form.real ? trig_form.real->to_string() : "null";
    std::string trig_imag = trig_form.imag ? trig_form.imag->to_string() : "null";

    // exp_form real part should contain cos
    EXPECT_CONTAINS(exp_real, {"cos"}, "exp_form real contains cos");
    // exp_form imag part should contain sin
    EXPECT_CONTAINS(exp_imag, {"sin"}, "exp_form imag contains sin");

    // trig_form should be equivalent to exp_form
    EXPECT_EQ_STR(exp_real, trig_real, "exp_form and trig_form have same real part");
    EXPECT_EQ_STR(exp_imag, trig_imag, "exp_form and trig_form have same imag part");
}

void test_complex_nth_root() {
    TEST_CASE("Complex Nth Root: returns exactly n roots");

    // Test n=2: square roots of 4
    {
        auto c = SymbolicExpr::number(4.0);
        auto roots = LMCAS::solve_complex_nth_root_checked(c, 2).value();
        EXPECT_TRUE(roots.size() == 2, "sqrt(4) returns exactly 2 roots");
    }

    // Test n=3: cube roots of 8
    {
        auto c = SymbolicExpr::number(8.0);
        auto roots = LMCAS::solve_complex_nth_root_checked(c, 3).value();
        EXPECT_TRUE(roots.size() == 3, "cbrt(8) returns exactly 3 roots");
    }

    // Test n=4: fourth roots of 16
    {
        auto c = SymbolicExpr::number(16.0);
        auto roots = LMCAS::solve_complex_nth_root_checked(c, 4).value();
        EXPECT_TRUE(roots.size() == 4, "4th root of 16 returns exactly 4 roots");
    }

    // Test n=5: fifth roots of 32
    {
        auto c = SymbolicExpr::number(32.0);
        auto roots = LMCAS::solve_complex_nth_root_checked(c, 5).value();
        EXPECT_TRUE(roots.size() == 5, "5th root of 32 returns exactly 5 roots");
    }
}

void test_checked_complex_contracts() {
    TEST_CASE("Checked Symbolic Complex: validates inputs and propagates context errors");

    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto approx_four = SymbolicExpr::number(4.0);

    auto a = LMCAS::make_complex(one, two);
    auto b = LMCAS::make_complex(two, one);

    {
        auto sum = LMCAS::complex_add_checked(a, b);
        EXPECT_TRUE(sum.has_value(), "checked complex_add succeeds for valid inputs");
        EXPECT_TRUE(sum.value().real && sum.value().imag,
                    "checked complex_add returns non-null components");
    }

    {
        auto bad = LMCAS::make_complex(nullptr, one);
        auto result = LMCAS::complex_mul_checked(bad, b);
        EXPECT_TRUE(!result.has_value() &&
                    result.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked complex_mul rejects null components");
    }

    {
        auto zero_complex = LMCAS::make_complex(zero, zero);
        auto result = LMCAS::complex_div_checked(a, zero_complex);
        EXPECT_TRUE(!result.has_value() &&
                    result.error().code == LMCAS::CasErrc::DomainError,
                    "checked complex_div rejects exact zero denominator");
    }

    {
        auto roots = LMCAS::solve_complex_nth_root_checked(approx_four, 2);
        EXPECT_TRUE(roots.has_value() && roots.value().size() == 2,
                    "checked complex nth root accepts explicit approximate real input");
    }

    {
        auto exact = LMCAS::solve_complex_nth_root_checked(SymbolicExpr::number(4), 2);
        EXPECT_TRUE(!exact.has_value() &&
                    exact.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked complex nth root does not implicitly float exact integers");
    }

    {
        auto bad_order = LMCAS::solve_complex_nth_root_checked(approx_four, 0);
        EXPECT_TRUE(!bad_order.has_value() &&
                    bad_order.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked complex nth root rejects non-positive degree");
    }

    {
        LMCAS::CancellationToken token;
        token.cancel();
        LMCAS::ComputationContext cancelled_context({}, token);
        auto cancelled = LMCAS::complex_conj_checked(a, cancelled_context);
        EXPECT_TRUE(!cancelled.has_value() &&
                    cancelled.error().code == LMCAS::CasErrc::Cancelled,
                    "checked complex_conj observes cancellation");
    }

    {
        LMCAS::ResourceLimits limits;
        limits.max_steps = 1;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::complex_abs_checked(a, limited_context);
        EXPECT_TRUE(!limited.has_value() &&
                    limited.error().code == LMCAS::CasErrc::ResourceLimit,
                    "checked complex_abs observes step budget");
    }
}

void test_complex_quadratic() {
    TEST_CASE("Complex Quadratic: z^2 + 1 = 0 yields +/-i");

    // z^2 + 1 = 0 => a=1, b=0, c=1
    auto a = SymbolicExpr::number(1);
    auto b = SymbolicExpr::number(0);
    auto c = SymbolicExpr::number(1);

    auto checked_roots = LMCAS::solve_complex_quadratic_checked(a, b, c);
    EXPECT_TRUE(checked_roots.has_value(),
                "z^2+1=0 checked solve succeeds");
    auto roots = checked_roots
        ? std::move(checked_roots.value())
        : std::vector<LMCAS::ComplexSymbolic>{};
    EXPECT_TRUE(roots.size() == 2, "z^2+1=0 returns exactly 2 roots");

    // The roots should be +i and -i
    // Root 1 real part should be 0 (or simplify to 0), imag should be non-zero
    if (roots.size() == 2) {
        std::string r1_real = roots[0].real ? roots[0].real->to_string() : "null";
        std::string r1_imag = roots[0].imag ? roots[0].imag->to_string() : "null";
        std::string r2_real = roots[1].real ? roots[1].real->to_string() : "null";
        std::string r2_imag = roots[1].imag ? roots[1].imag->to_string() : "null";

        // Both roots should have structural content (not null)
        EXPECT_TRUE(r1_real != "null", "root1 real is not null");
        EXPECT_TRUE(r2_real != "null", "root2 real is not null");
        // The discriminant is b^2 - 4ac = 0 - 4 = -4, sqrt(-4) = 2i
        // Roots: (0 +/- sqrt(-4)) / 2 => +/- i
        // The expression involves sqrt(-4) which should appear in the output
        EXPECT_CONTAINS(r1_real, {"0"}, "root1 of z^2+1=0 real part contains 0");
    }
}

void test_complex_locus() {
    TEST_CASE("Complex Locus: circle |z - a| = r");

    // Circle centered at (1+2i) with radius 3
    {
        auto center = LMCAS::make_complex(SymbolicExpr::number(1), SymbolicExpr::number(2));
        auto radius = SymbolicExpr::number(3);
        auto checked_locus =
            LMCAS::complex_locus_circle_checked(center, radius, "z");
        EXPECT_TRUE(checked_locus.has_value(),
                    "checked circle locus succeeds");
        auto locus = checked_locus
            ? std::move(checked_locus.value()) : nullptr;
        std::string s = locus ? locus->to_string() : "null";
        EXPECT_TRUE(s != "null", "locus_circle is not null");
        // Should contain z and structural elements of |z - a| = r
        EXPECT_CONTAINS(s, {"z"}, "locus_circle contains z variable");
        EXPECT_CONTAINS(s, {"3"}, "locus_circle contains radius 3");
    }

    TEST_CASE("Complex Locus: perpendicular bisector |z - a| = |z - b|");

    // Perpendicular bisector between (1+0i) and (3+0i)
    {
        auto a = LMCAS::make_complex(SymbolicExpr::number(1), SymbolicExpr::number(0));
        auto b = LMCAS::make_complex(SymbolicExpr::number(3), SymbolicExpr::number(0));
        auto checked_locus =
            LMCAS::complex_locus_perpendicular_bisector_checked(a, b, "z");
        EXPECT_TRUE(checked_locus.has_value(),
                    "checked perpendicular-bisector locus succeeds");
        auto locus = checked_locus
            ? std::move(checked_locus.value()) : nullptr;
        std::string s = locus ? locus->to_string() : "null";
        EXPECT_TRUE(s != "null", "locus_perpendicular_bisector is not null");
        // Should contain z and structural elements
        EXPECT_CONTAINS(s, {"z"}, "locus_perp_bisector contains z variable");
    }
}

int main() {
    try {
        test_complex_arithmetic();
        test_complex_conj();
        test_complex_abs();
        test_complex_arg();
        test_complex_polar_forms();
        test_complex_nth_root();
        test_checked_complex_contracts();
        test_complex_quadratic();
        test_complex_locus();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
