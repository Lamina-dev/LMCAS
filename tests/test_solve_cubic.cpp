#include "test_common.hpp"
#include "solve_polynomial.hpp"
#include "residual_verification.hpp"
#include "root_of_identity.hpp"
#include "expr.hpp"
#include <cmath>
#include <random>
#include <sstream>
#include <algorithm>
#include <array>
#include <string>
#include <stdexcept>

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

static double eval_cubic_at(double a, double b, double c, double d, double x) {
    return a * x * x * x + b * x * x + c * x + d;
}

static void expect_exact_cubic_residual(
    const std::shared_ptr<SymbolicExpr> &root,
    const std::shared_ptr<SymbolicExpr> &a,
    const std::shared_ptr<SymbolicExpr> &b,
    const std::shared_ptr<SymbolicExpr> &c,
    const std::shared_ptr<SymbolicExpr> &d) {
    auto residual = SymbolicExpr::add(
        SymbolicExpr::multiply(a, SymbolicExpr::power(root, num(3))),
        SymbolicExpr::add(
            SymbolicExpr::multiply(b, SymbolicExpr::power(root, num(2))),
            SymbolicExpr::add(SymbolicExpr::multiply(c, root), d)));
    ComputationContext context;
    auto proof = check_zero_residual(residual, context);
    EXPECT_TRUE((proof && std::holds_alternative<ProvedZeroResidual>(proof.value()))) << "every exact cubic root has a proved zero residual";
}

static void expect_rootof_roots(
    const std::vector<std::shared_ptr<SymbolicExpr>> &roots,
    const Polynomial<Rational> &polynomial) {
    auto expression = poly_to_symbolic(polynomial);
    for (int index = 0; index < polynomial.degree(); ++index) {
        auto expected = make_rootof_checked(
            expression, polynomial.variable_name, static_cast<std::size_t>(index));
        EXPECT_TRUE((expected && std::any_of(
                                     roots.begin(), roots.end(), [&](const auto &root) {
                                         return root->compare(expected.value()) == 0;
                                     })))
            << "every exact algebraic identity is present in the solution";
    }
}

TEST(SolveCubic, CubicRationalExtractionAndQuadraticRemainder) {
    auto roots = LMCAS::solve_cubic(num(1), num(-1), num(-2), num(2), "x");
    ASSERT_EQ(roots.size(), 3u) << "(x-1)(x^2-2) has three exact roots";
    auto sqrt_two = SymbolicExpr::sqrt(num(2));
    const std::array<std::shared_ptr<SymbolicExpr>, 3> expected{
        num(1), sqrt_two, SymbolicExpr::multiply(num(-1), sqrt_two)};
    std::array<bool, 3> matched{};
    for (const auto &root : roots) {
        ASSERT_TRUE(root);
        expect_exact_cubic_residual(root, num(1), num(-1), num(-2), num(2));
        bool found = false;
        for (std::size_t index = 0; index < expected.size(); ++index) {
            if (!matched[index] && test_proved_equivalent(root, expected[index])) {
                matched[index] = true;
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "every root matches a distinct expected root";
    }
}

TEST(SolveCubic, CubicIrreducibleComplexRoots) {
    auto roots = LMCAS::solve_cubic(num(1), num(0), num(-2), num(-5), "x");
    EXPECT_TRUE((roots.size() == 3)) << "irreducible cubic has all three exact roots";
    expect_rootof_roots(roots, Polynomial<Rational>(
                                   {Rational(-5), Rational(-2), Rational(0), Rational(1)}, "x"));
    for (const auto &root : roots) {
        expect_exact_cubic_residual(root, num(1), num(0), num(-2), num(-5));
    }
}

TEST(SolveCubic, CubicExactCubeRoots) {
    auto roots = LMCAS::solve_cubic(num(1), num(0), num(0), num(-2), "x");
    EXPECT_TRUE((roots.size() == 3)) << "x^3-2 retains its complex conjugate pair";
    expect_rootof_roots(roots, Polynomial<Rational>(
                                   {Rational(-2), Rational(0), Rational(0), Rational(1)}, "x"));
    for (const auto &root : roots) {
        expect_exact_cubic_residual(root, num(1), num(0), num(0), num(-2));
    }
}

TEST(SolveCubic, CubicTripleRootX13X33x23x10) {
    auto roots = LMCAS::solve_cubic(num(1), num(-3), num(3), num(-1), "x");
    EXPECT_TRUE((roots.size() == 3)) << "Triple root: should return exactly 3 roots";

    for (const auto &root : roots) {
        EXPECT_TRUE((root->is_one())) << "each occurrence of the triple root is exactly one";
        expect_exact_cubic_residual(root, num(1), num(-3), num(3), num(-1));
    }
}

TEST(SolveCubic, CubicDiscriminantClassificationIsScaleRelative) {
    auto roots = LMCAS::solve_cubic(
        num(1),
        num(0),
        num(0),
        SymbolicExpr::number(-1e-200),
        "x");
    EXPECT_TRUE((roots.size() == 3)) << "tiny nonzero constant still has three cubic roots";
    if (roots.size() == 3) {
        const double expected = std::cbrt(1e-200);
        const double real_root = test_numeric_value(roots[0]);
        EXPECT_TRUE((std::abs(real_root - expected) <= expected * 1e-12)) << "tiny nonzero constant must not be classified as a triple zero root";
    }
}

TEST(SolveCubic, CubicCommonCoefficientScaleDoesNotOverflowDepression) {
    const double scale = 1e200;
    auto roots = LMCAS::solve_cubic(
        SymbolicExpr::number(scale),
        SymbolicExpr::number(-6.0 * scale),
        SymbolicExpr::number(11.0 * scale),
        SymbolicExpr::number(-6.0 * scale),
        "x");
    EXPECT_TRUE((roots.size() == 3)) << "scaled (x-1)(x-2)(x-3) returns three roots";
    if (roots.size() == 3) {
        std::vector<double> values{
            test_numeric_value(roots[0]),
            test_numeric_value(roots[1]),
            test_numeric_value(roots[2])};
        std::sort(values.begin(), values.end());
        EXPECT_TRUE((std::isfinite(values[0]) && std::isfinite(values[1]) &&
                     std::isfinite(values[2]) &&
                     std::abs(values[0] - 1.0) < 1e-12 &&
                     std::abs(values[1] - 2.0) < 1e-12 &&
                     std::abs(values[2] - 3.0) < 1e-12))
            << "common coefficient scaling preserves roots 1, 2, and 3";
    }
}

TEST(SolveCubic, CubicD0P0X33x2X12X20) {
    auto roots = LMCAS::solve_cubic(num(1), num(0), num(-3), num(2), "x");
    EXPECT_TRUE((roots.size() == 3)) << "D=0 p!=0: should return exactly 3 roots";

    if (roots.size() == 3) {
        double r1 = test_numeric_value(roots[0]);
        double r2 = test_numeric_value(roots[1]);
        double r3 = test_numeric_value(roots[2]);

        std::vector<double> sorted_roots = {r1, r2, r3};
        std::sort(sorted_roots.begin(), sorted_roots.end());

        EXPECT_TRUE((std::abs(sorted_roots[0] - (-2.0)) < 1e-10)) << "D=0 p!=0: smallest root = -2";
        EXPECT_TRUE((std::abs(sorted_roots[1] - 1.0) < 1e-10)) << "D=0 p!=0: middle root = 1";
        EXPECT_TRUE((std::abs(sorted_roots[2] - 1.0) < 1e-10)) << "D=0 p!=0: largest root = 1";

        for (double r : sorted_roots) {
            double res = eval_cubic_at(1.0, 0.0, -3.0, 2.0, r);
            EXPECT_TRUE((std::abs(res) < 1e-10)) << "D=0 p!=0: root satisfies equation";
        }
    }
}

TEST(SolveCubic, CubicD0CasusIrreducibilisX33x10) {
    auto roots = LMCAS::solve_cubic(num(1), num(0), num(-3), num(1), "x");
    EXPECT_TRUE((roots.size() == 3)) << "D<0: should return exactly 3 roots";
    expect_rootof_roots(roots, Polynomial<Rational>(
                                   {Rational(1), Rational(-3), Rational(0), Rational(1)}, "x"));

    for (const auto &root : roots) {
        expect_exact_cubic_residual(root, num(1), num(0), num(-3), num(1));
    }
}

TEST(SolveCubic, CubicSymbolicCoefficientsX3AX2X10) {
    auto sym_a = SymbolicExpr::variable("a");
    auto roots = LMCAS::solve_cubic(num(1), sym_a, num(1), num(1), "x");
    EXPECT_TRUE((roots.size() == 3)) << "Symbolic: should return exactly 3 roots";

    if (roots.size() == 3) {

        bool found_a = false;
        for (const auto &root : roots) {
            if (root->compare(root->substitute("a", num(0))) != 0) {
                found_a = true;
                break;
            }
        }
        EXPECT_TRUE((found_a)) << "Symbolic: roots contain variable 'a'";
    }
}

TEST(SolveCubic, CubicExactHugeCoefficientAvoidsUnsafeNumericUnderflow) {
    auto leading = SymbolicExpr::number(
        BigInt("1" + std::string(400, '0')));
    auto roots = LMCAS::solve_cubic(leading, num(0), num(0), num(-1), "x");
    EXPECT_TRUE((roots.size() == 3)) << "huge exact cubic still returns every root";
    expect_rootof_roots(roots, Polynomial<Rational>(
                                   {Rational(-1), Rational(0), Rational(0),
                                    Rational(BigInt("1" + std::string(400, '0')))},
                                   "x"));
    for (const auto &root : roots) {
        expect_exact_cubic_residual(root, leading, num(0), num(0), num(-1));
    }
}

TEST(SolveCubic, CubicA0DelegationToQuadratic0X32X24X20) {
    auto roots = LMCAS::solve_cubic(num(0), num(2), num(-4), num(2), "x");
    EXPECT_TRUE((roots.size() >= 1 && roots.size() <= 2)) << "a=0 delegation: returns 1-2 roots (quadratic)";

    if (roots.size() >= 1) {
        double r1 = test_numeric_value(roots[0]);
        EXPECT_TRUE((std::abs(r1 - 1.0) < 1e-10)) << "a=0 delegation: root = 1";
    }
}

TEST(SolveCubic, CubicA0DelegationToLinear0X30X23X60) {
    auto roots = LMCAS::solve_cubic(num(0), num(0), num(3), num(-6), "x");
    EXPECT_TRUE((roots.size() == 1)) << "a=0 b=0 delegation: returns 1 root (linear)";

    if (roots.size() == 1) {
        double r1 = test_numeric_value(roots[0]);
        EXPECT_TRUE((std::abs(r1 - 2.0) < 1e-10)) << "a=0 b=0 delegation: root = 2";
    }
}

TEST(SolveCubic, CubicDegenerationToAConstantEquation) {
    auto no_roots = LMCAS::solve_cubic(
        num(0), num(0), num(0), num(7), "x");
    EXPECT_TRUE((no_roots.empty())) << "nonzero constant cubic degeneration has no roots";

    bool rejected_indeterminate = false;
    try {
        (void)LMCAS::solve_cubic(
            num(0), num(0), num(0), num(0), "x");
    } catch (const std::invalid_argument &) {
        rejected_indeterminate = true;
    }
    EXPECT_TRUE((rejected_indeterminate)) << "identically zero cubic is rejected as indeterminate";
}

TEST(SolveCubic, CubicRootVerification) {
    const double RESIDUAL_TOL = 1e-10;
    const int NUM_CUBIC_TRIALS = 50;
    int cubic_verify_pass_count = 0;
    int cubic_verify_total_roots = 0;

    std::mt19937 rng_cubic(123);
    std::uniform_int_distribution<int> coeff_dist(-10, 10);

    for (int trial = 0; trial < NUM_CUBIC_TRIALS; ++trial) {
        int a_val = coeff_dist(rng_cubic);

        while (a_val == 0)
            a_val = coeff_dist(rng_cubic);
        int b_val = coeff_dist(rng_cubic);
        int c_val = coeff_dist(rng_cubic);
        int d_val = coeff_dist(rng_cubic);

        auto roots = LMCAS::solve_cubic(
            SymbolicExpr::number(static_cast<double>(a_val)),
            SymbolicExpr::number(static_cast<double>(b_val)),
            SymbolicExpr::number(static_cast<double>(c_val)),
            SymbolicExpr::number(static_cast<double>(d_val)), "x");

        if (roots.size() != 3) {
            std::ostringstream msg;
            msg << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
                << ", c=" << c_val << ", d=" << d_val
                << "): expected 3 roots, got " << roots.size();
            ADD_FAILURE() << msg.str();
            continue;
        }

        bool trial_ok = true;
        for (size_t i = 0; i < roots.size(); ++i) {
            double r = test_numeric_value(roots[i]);

            if (std::isnan(r) || std::isinf(r)) {
                continue;
            }
            cubic_verify_total_roots++;
            double residual = eval_cubic_at((double)a_val, (double)b_val,
                                            (double)c_val, (double)d_val, r);
            if (std::abs(residual) >= RESIDUAL_TOL) {
                trial_ok = false;
                std::ostringstream msg;
                msg << "Trial " << trial << " root " << i
                    << " (a=" << a_val << ", b=" << b_val
                    << ", c=" << c_val << ", d=" << d_val
                    << "): |f(r)| = " << std::abs(residual) << " >= 1e-10"
                    << " (r = " << r << ")";
                ADD_FAILURE() << msg.str();
            }
        }
        if (trial_ok) {
            cubic_verify_pass_count++;
        }
    }

    {
        std::ostringstream msg;
        msg << "Cubic root verification: " << cubic_verify_pass_count
            << "/" << NUM_CUBIC_TRIALS << " trials passed ("
            << cubic_verify_total_roots << " real roots verified)";
        EXPECT_TRUE((cubic_verify_pass_count == NUM_CUBIC_TRIALS)) << msg.str();
    }
}

static bool cubic_roots_have_nonfinite_value(double r1, double r2, double r3) {
    return std::isnan(r1) || std::isnan(r2) || std::isnan(r3) ||
           std::isinf(r1) || std::isinf(r2) || std::isinf(r3);
}

TEST(SolveCubic, VietaSFormulasForCubics) {
    const double TOLERANCE = 1e-8;
    const int NUM_TRIALS = 60;
    int vieta_pass_count = 0;

    std::mt19937 rng(42);
    std::uniform_int_distribution<int> root_dist(-5, 5);

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int r1_int = root_dist(rng);
        int r2_int = root_dist(rng);
        int r3_int = root_dist(rng);

        int b_val = -(r1_int + r2_int + r3_int);
        int c_val = r1_int * r2_int + r1_int * r3_int + r2_int * r3_int;
        int d_val = -(r1_int * r2_int * r3_int);

        auto roots = LMCAS::solve_cubic(num(1), num(b_val), num(c_val), num(d_val), "x");

        if (roots.size() != 3) {
            std::ostringstream msg;
            msg << "Trial " << trial << " (b=" << b_val << ", c=" << c_val
                << ", d=" << d_val << "): expected 3 roots, got " << roots.size();
            ADD_FAILURE() << msg.str();
            continue;
        }

        double r1 = test_numeric_value(roots[0]);
        double r2 = test_numeric_value(roots[1]);
        double r3 = test_numeric_value(roots[2]);

        if (cubic_roots_have_nonfinite_value(r1, r2, r3)) {
            std::ostringstream msg;
            msg << "Trial " << trial << " (b=" << b_val << ", c=" << c_val
                << ", d=" << d_val << "): root evaluation produced NaN/Inf";
            ADD_FAILURE() << msg.str();
            continue;
        }

        double sum_roots = r1 + r2 + r3;
        double expected_sum = -(double)b_val;
        bool sum_ok = std::abs(sum_roots - expected_sum) < TOLERANCE;

        double sum_products = r1 * r2 + r1 * r3 + r2 * r3;
        double expected_products = (double)c_val;
        bool products_ok = std::abs(sum_products - expected_products) < TOLERANCE;

        double product_roots = r1 * r2 * r3;
        double expected_product = -(double)d_val;
        bool product_ok = std::abs(product_roots - expected_product) < TOLERANCE;

        if (!sum_ok || !products_ok || !product_ok) {

        } else {
            vieta_pass_count++;
        }
    }

    {
        std::ostringstream msg;
        msg << "Vieta's formulas: " << vieta_pass_count << "/" << NUM_TRIALS << " trials passed";
        EXPECT_TRUE((vieta_pass_count == NUM_TRIALS)) << msg.str();
    }
}

TEST(SolveCubic, FormalResidualCommonDomains) {
    for (const char *source : {
             "ln(x)+x/x-1-ln(x)",
             "exp(ln(x))-x*exp(ln(x))/x",
             "1/(x^2-1)+1/(2*(x+1))-1/(2*(x-1))",
             "1/(sqrt(2)^2)-1/2"}) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << "formal residual parses";
        if (!expression)
            continue;
        ResourceLimits limits;
        limits.max_steps = 10000;
        ComputationContext context(limits);
        auto proof = check_zero_residual(expression.value(), context);
        EXPECT_TRUE((proof.has_value())) << (proof ? "formal residual proof completes within budget" : proof.error().operation + ": " + proof.error().message);
        EXPECT_TRUE((proof && std::holds_alternative<ProvedZeroResidual>(proof.value()))) << "common-domain rational identity is certified";
    }
    auto nonidentity = parse_expr("sqrt(2)+1");
    ComputationContext context;
    auto proof = check_zero_residual(nonidentity.value(), context);
    EXPECT_TRUE((proof && !std::holds_alternative<ProvedZeroResidual>(proof.value()))) << "a nonzero algebraic remainder cannot certify an identity";
    auto oversized = SymbolicExpr::variable("x");
    for (int exponent : {1000, 1000, 1000, 3})
        oversized = SymbolicExpr::power(oversized, SymbolicExpr::number(exponent));
    ComputationContext bounded;
    auto rejected = check_zero_residual(oversized, bounded);
    EXPECT_TRUE((!rejected && rejected.error().code == CasErrc::ResourceLimit)) << "nested powers cannot overflow the formal monomial degree representation";
}
