#include "test_common.hpp"
#include "root_of_utils.hpp"
#include "solve_polynomial.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

using namespace LMCAS;

using LMCAS::Polynomial;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

static double eval_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) return std::nan("");
    return expr->to_numeric();
}

int main() {

    TEST_CASE("RootOf - Checked evaluation preserves failure semantics");
    {
        auto x = SymbolicExpr::variable("x");
        auto exact_poly = SymbolicExpr::add(SymbolicExpr::power(x, num(2)), num(-4));

        LMCAS::ComputationContext valid_context;
        auto valid = LMCAS::rootof_evaluate_checked(
            SymbolicExpr::root_of(exact_poly, "x", 0), valid_context);
        EXPECT_TRUE(valid && std::abs(valid.value() + 2.0) < 1e-10,
                    "checked RootOf verifies an exact real root");

        const BigInt largest_finite_integer =
            (BigInt(1) << 1024) - (BigInt(1) << 971);
        Polynomial<Rational> endpoint_poly(
            {Rational(-largest_finite_integer), Rational(1)}, "x");
        auto endpoint_root = LMCAS::rootof_evaluate_checked(
            SymbolicExpr::root_of(
                LMCAS::poly_to_symbolic(endpoint_poly), "x", 0));
        const std::string endpoint_message =
            endpoint_root
                ? "checked RootOf represents the largest finite double root"
                : "largest finite RootOf failed: " +
                      endpoint_root.error().message;
        EXPECT_TRUE(
            endpoint_root &&
                endpoint_root.value() == std::numeric_limits<double>::max(),
            endpoint_message);

        LMCAS::ComputationContext invalid_index_context;
        auto invalid_index = LMCAS::make_rootof_checked(
            exact_poly, "x", 2, invalid_index_context);
        EXPECT_TRUE(!invalid_index &&
                        invalid_index.error().code == LMCAS::CasErrc::InvalidArgument,
                    "out-of-range RootOf construction is InvalidArgument");

        auto complex_poly = SymbolicExpr::add(SymbolicExpr::power(x, num(2)), num(1));
        LMCAS::ComputationContext complex_context;
        auto complex_expression = SymbolicExpr::root_of(complex_poly, "x", 0);
        auto complex_root = LMCAS::rootof_evaluate_checked(
            complex_expression, complex_context);
        EXPECT_TRUE(!complex_root &&
                        complex_root.error().code == LMCAS::CasErrc::DomainError,
                    "real RootOf wrapper rejects a non-real selected root");
        auto complex_value = LMCAS::rootof_evaluate_complex_checked(
            complex_expression);
        EXPECT_TRUE(complex_value &&
                        std::abs(complex_value.value().real.value) < 1e-9 &&
                        std::abs(complex_value.value().imag.value + 1.0) < 1e-8,
                    "complex RootOf evaluates the first ordered root as -i");
        auto conjugate_value = LMCAS::rootof_evaluate_complex_checked(
            SymbolicExpr::root_of(complex_poly, "x", 1));
        EXPECT_TRUE(conjugate_value &&
                        std::abs(conjugate_value.value().imag.value - 1.0) < 1e-8,
                    "complex RootOf evaluates the conjugate root as +i");

        auto parametric_poly = SymbolicExpr::add(
            SymbolicExpr::power(x, num(2)), SymbolicExpr::variable("a"));
        LMCAS::ComputationContext parametric_context;
        auto parametric = LMCAS::make_rootof_checked(
            parametric_poly, "x", 0, parametric_context);
        EXPECT_TRUE(!parametric &&
                        parametric.error().code == LMCAS::CasErrc::Inconclusive,
                    "parametric RootOf construction is explicitly inconclusive");

        auto approximate_poly = SymbolicExpr::add(
            SymbolicExpr::power(x, num(2)), SymbolicExpr::number(-4.0));
        LMCAS::ComputationContext approximate_context;
        auto approximate = LMCAS::make_rootof_checked(
            approximate_poly, "x", 0, approximate_context);
        EXPECT_TRUE(!approximate &&
                        approximate.error().code == LMCAS::CasErrc::Inconclusive,
                    "ApproxReal coefficients are not silently exactified");

        LMCAS::CancellationToken cancellation;
        cancellation.cancel();
        LMCAS::ComputationContext cancelled_context({}, cancellation);
        auto cancelled = LMCAS::rootof_evaluate_checked(
            SymbolicExpr::root_of(exact_poly, "x", 0), cancelled_context);
        EXPECT_TRUE(!cancelled &&
                        cancelled.error().code == LMCAS::CasErrc::Cancelled,
                    "checked RootOf observes cancellation");

        LMCAS::ResourceLimits limits;
        limits.max_steps = 1;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::rootof_evaluate_checked(
            SymbolicExpr::root_of(exact_poly, "x", 0), limited_context);
        EXPECT_TRUE(!limited &&
                        limited.error().code == LMCAS::CasErrc::ResourceLimit,
                    "checked RootOf enforces traversal budgets");
    }

    TEST_CASE("RootOf - Out-of-range indices are rejected");
    {
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(SymbolicExpr::variable("x"), num(3)),
                SymbolicExpr::variable("x")),
            num(1));
        for (const std::size_t index : {std::size_t(3), std::size_t(5),
                                        std::size_t(100)}) {
            auto result = LMCAS::make_rootof_checked(
                poly_expr, "x", index);
            EXPECT_TRUE(!result &&
                            result.error().code ==
                                LMCAS::CasErrc::InvalidArgument,
                        "out-of-range RootOf index is rejected");
        }
    }

    TEST_CASE("RootOf - Negative indices are rejected at construction");
    {
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(SymbolicExpr::variable("x"), num(3)),
                SymbolicExpr::multiply(num(-2), SymbolicExpr::variable("x"))),
            num(1));

        bool negative_one_rejected = false;
        try {
            (void)SymbolicExpr::root_of(poly_expr, "x", -1);
        } catch (const std::invalid_argument&) {
            negative_one_rejected = true;
        }
        EXPECT_TRUE(negative_one_rejected,
                    "RootOf construction rejects index -1");

        bool negative_ten_rejected = false;
        try {
            (void)SymbolicExpr::root_of(poly_expr, "x", -10);
        } catch (const std::invalid_argument&) {
            negative_ten_rejected = true;
        }
        EXPECT_TRUE(negative_ten_rejected,
                    "RootOf construction rejects index -10");
    }

    TEST_CASE("RootOf - canonical identity ignores scale, repetition, and dummy name");
    {
        auto x = SymbolicExpr::variable("x");
        auto base = SymbolicExpr::add(
            SymbolicExpr::power(x, num(2)), num(-2));
        auto scaled = SymbolicExpr::add(
            SymbolicExpr::multiply(
                num(2), SymbolicExpr::power(x, num(2))),
            num(-4));
        auto y = SymbolicExpr::variable("y");
        auto renamed = SymbolicExpr::add(
            SymbolicExpr::power(y, num(2)), num(-2));
        auto repeated = SymbolicExpr::power(base, num(2));

        auto canonical = SymbolicExpr::root_of(base, "x", 0);
        auto scaled_root = SymbolicExpr::root_of(scaled, "x", 0);
        auto renamed_root = SymbolicExpr::root_of(renamed, "y", 0);
        auto repeated_root = SymbolicExpr::root_of(repeated, "x", 0);
        EXPECT_TRUE(
            LMCAS::detail::node(canonical)->equals(
                *LMCAS::detail::node(scaled_root)) &&
            LMCAS::detail::node(canonical)->equals(
                *LMCAS::detail::node(renamed_root)) &&
            LMCAS::detail::node(canonical)->equals(
                *LMCAS::detail::node(repeated_root)),
            "canonical RootOf identities are structurally equal");
        EXPECT_TRUE(
            LMCAS::detail::node(canonical)->hash() ==
                LMCAS::detail::node(repeated_root)->hash(),
            "canonical RootOf identities share a hash");
    }


    TEST_CASE("RootOf - Parametric coefficients are Inconclusive");
    {
        auto x = SymbolicExpr::variable("x");
        auto a = SymbolicExpr::variable("a");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(3)),
                SymbolicExpr::multiply(a, x)),
            num(1));
        auto result = LMCAS::make_rootof_checked(
            poly_expr, "x", 0);
        EXPECT_TRUE(!result &&
                        result.error().code == LMCAS::CasErrc::Inconclusive,
                    "parametric RootOf construction is Inconclusive");
    }

    TEST_CASE("RootOf - Multiple parametric coefficients are Inconclusive");
    {
        auto x = SymbolicExpr::variable("x");
        auto b = SymbolicExpr::variable("b");
        auto c = SymbolicExpr::variable("c");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(2)),
                SymbolicExpr::multiply(b, x)),
            c);
        auto result = LMCAS::make_rootof_checked(
            poly_expr, "x", 0);
        EXPECT_TRUE(!result &&
                        result.error().code == LMCAS::CasErrc::Inconclusive,
                    "multiple parametric coefficients are Inconclusive");
    }

    TEST_CASE("RootOf - Simplify degree-2 polynomial to closed-form");
    {

        auto x = SymbolicExpr::variable("x");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::power(x, num(2)),
            num(-4));

        auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
        auto simplified_k0 = LMCAS::rootof_simplify(rootof_k0);

        std::string s0 = simplified_k0->to_string();
        EXPECT_TRUE(s0.find("RootOf") == std::string::npos,
            "rootof_simplify(degree-2, k=0) returns non-RootOf expression");

        double val0 = eval_numeric(simplified_k0);
        EXPECT_TRUE(!std::isnan(val0) && std::abs(val0 - (-2.0)) < 1e-10,
            "rootof_simplify(x^2-4, x, 0) = -2");

        auto rootof_k1 = SymbolicExpr::root_of(poly_expr, "x", 1);
        auto simplified_k1 = LMCAS::rootof_simplify(rootof_k1);

        std::string s1 = simplified_k1->to_string();
        EXPECT_TRUE(s1.find("RootOf") == std::string::npos,
            "rootof_simplify(degree-2, k=1) returns non-RootOf expression");

        double val1 = eval_numeric(simplified_k1);
        EXPECT_TRUE(!std::isnan(val1) && std::abs(val1 - 2.0) < 1e-10,
            "rootof_simplify(x^2-4, x, 1) = 2");
    }

    TEST_CASE("RootOf - Non-polynomial construction is Inconclusive");
    {
        auto x = SymbolicExpr::variable("x");
        auto unsupported = SymbolicExpr::add(
            SymbolicExpr::power(x, num(2)), SymbolicExpr::sin(x));
        auto root = LMCAS::make_rootof_checked(
            unsupported, "x", 0);
        EXPECT_TRUE(!root &&
                        root.error().code == LMCAS::CasErrc::Inconclusive,
                    "non-polynomial RootOf input is rejected");
    }

    TEST_CASE("RootOf - Higher-degree simplification preserves canonical identity");
    {
        auto x = SymbolicExpr::variable("x");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::power(x, num(3)),
                    SymbolicExpr::multiply(
                        num(-6), SymbolicExpr::power(x, num(2)))),
                SymbolicExpr::multiply(num(11), x)),
            num(-6));
        for (int index = 0; index < 3; ++index) {
            auto root = SymbolicExpr::root_of(poly_expr, "x", index);
            auto simplified = LMCAS::rootof_simplify(root);
            EXPECT_TRUE(
                LMCAS::detail::node(simplified)->equals(
                    *LMCAS::detail::node(root)),
                "cubic RootOf simplification preserves identity");
            auto value = LMCAS::rootof_evaluate_checked(root);
            EXPECT_TRUE(value &&
                            std::abs(value.value() -
                                     static_cast<double>(index + 1)) < 1e-8,
                        "cubic RootOf exact order is stable");
        }
    }

    TEST_CASE("RootOf - Degree-four identities retain exact order");
    {
        auto x = SymbolicExpr::variable("x");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(4)),
                SymbolicExpr::multiply(
                    num(-5), SymbolicExpr::power(x, num(2)))),
            num(4));
        const double expected[] = {-2.0, -1.0, 1.0, 2.0};
        for (int index = 0; index < 4; ++index) {
            auto root = SymbolicExpr::root_of(poly_expr, "x", index);
            auto simplified = LMCAS::rootof_simplify(root);
            EXPECT_TRUE(
                LMCAS::detail::node(simplified)->equals(
                    *LMCAS::detail::node(root)),
                "quartic RootOf simplification preserves identity");
            auto value = LMCAS::rootof_evaluate_checked(root);
            EXPECT_TRUE(value &&
                            std::abs(value.value() - expected[index]) < 1e-8,
                        "quartic RootOf exact order is stable");
        }
    }

    TEST_CASE("RootOf - Valid index on numeric polynomial evaluates correctly");
    {

        auto x = SymbolicExpr::variable("x");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::power(x, num(3)),
                    SymbolicExpr::multiply(num(-6), SymbolicExpr::power(x, num(2)))),
                SymbolicExpr::multiply(num(11), x)),
            num(-6));

        auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
        auto result0 = LMCAS::rootof_evaluate_checked(rootof_k0);
        EXPECT_TRUE(result0.has_value(),
            "checked RootOf evaluation with valid k=0 returns a value");
        if (result0.has_value()) {
            EXPECT_TRUE(std::abs(result0.value() - 1.0) < 1e-10,
                "checked RootOf evaluation of cubic k=0 is the smallest root");
        }

        auto rootof_k2 = SymbolicExpr::root_of(poly_expr, "x", 2);
        auto result2 = LMCAS::rootof_evaluate_checked(rootof_k2);
        EXPECT_TRUE(result2.has_value(),
            "checked RootOf evaluation with valid k=2 returns a value");
        if (result2.has_value()) {
            EXPECT_TRUE(std::abs(result2.value() - 3.0) < 1e-10,
                "checked RootOf evaluation of cubic k=2 is the largest root");
        }
    }


    TEST_CASE("RootOf - Certified higher-degree complex ordering");
    {
        Polynomial<Rational> cubic(
            {Rational(-1), Rational(0), Rational(0), Rational(1)}, "x");
        auto cubic_expr = LMCAS::poly_to_symbolic(cubic);
        const double cubic_expected[][2] = {
            {1.0, 0.0},
            {-0.5, -std::sqrt(3.0) / 2.0},
            {-0.5, std::sqrt(3.0) / 2.0}};
        for (int index = 0; index < 3; ++index) {
            auto root = SymbolicExpr::root_of(cubic_expr, "x", index);
            auto value = LMCAS::rootof_evaluate_complex_checked(root);
            EXPECT_TRUE(
                value &&
                    std::abs(value.value().real.value -
                             cubic_expected[index][0]) < 1e-9 &&
                    std::abs(value.value().imag.value -
                             cubic_expected[index][1]) < 1e-9,
                "x^3-1 RootOf index has certified real-first complex order");
            if (index != 0) {
                LMCAS::ComputationContext real_context;
                auto real = LMCAS::rootof_evaluate_checked(
                    root, real_context);
                EXPECT_TRUE(
                    !real &&
                        real.error().code == LMCAS::CasErrc::DomainError,
                    "real evaluation rejects a selected cubic non-real root");
            }
        }

        Polynomial<Rational> shifted(
            {Rational(10), Rational(-6), Rational(3),
             Rational(0), Rational(1)}, "x");
        auto shifted_expr = LMCAS::poly_to_symbolic(shifted);
        const double shifted_expected[][2] = {
            {-1.0, -2.0}, {-1.0, 2.0},
            {1.0, -1.0}, {1.0, 1.0}};
        for (int index = 0; index < 4; ++index) {
            auto value = LMCAS::rootof_evaluate_complex_checked(
                SymbolicExpr::root_of(shifted_expr, "x", index));
            EXPECT_TRUE(
                value &&
                    std::abs(value.value().real.value -
                             shifted_expected[index][0]) < 1e-9 &&
                    std::abs(value.value().imag.value -
                             shifted_expected[index][1]) < 1e-9,
                "shifted quartic roots use exact real/imaginary lexicographic order");
        }
    }


    TEST_CASE("RootOf - Clustered complex rectangles remain distinct");
    {
        const Rational small(1, 1048576);
        Polynomial<Rational> polynomial =
            Polynomial<Rational>({small, Rational(0), Rational(1)}, "x") *
            Polynomial<Rational>({Rational(4) * small,
                                  Rational(0), Rational(1)}, "x");
        auto expression = LMCAS::poly_to_symbolic(polynomial);
        const double expected[] = {
            -1.0 / 512.0, -1.0 / 1024.0,
            1.0 / 1024.0, 1.0 / 512.0};
        for (int index = 0; index < 4; ++index) {
            auto value = LMCAS::rootof_evaluate_complex_checked(
                SymbolicExpr::root_of(expression, "x", index));
            EXPECT_TRUE(
                value &&
                    std::abs(value.value().real.value) < 1e-12 &&
                    std::abs(value.value().imag.value - expected[index]) <
                        1e-10 &&
                    value.value().imag.absolute_error > 0.0,
                "clustered imaginary roots retain distinct certified enclosures");
        }
    }

    TEST_CASE("RootOf - Generic quintic complex isolation and resources");
    {
        Polynomial<Rational> polynomial(
            {Rational(1), Rational(1), Rational(0),
             Rational(0), Rational(0), Rational(1)}, "x");
        auto expression = LMCAS::poly_to_symbolic(polynomial);
        std::vector<LMCAS::ApproxComplex> values;
        for (int index = 0; index < 5; ++index) {
            auto value = LMCAS::rootof_evaluate_complex_checked(
                SymbolicExpr::root_of(expression, "x", index));
            EXPECT_TRUE(
                value && std::isfinite(value.value().real.value) &&
                    std::isfinite(value.value().imag.value) &&
                    value.value().real.absolute_error >= 0.0 &&
                    value.value().imag.absolute_error >= 0.0,
                "generic quintic RootOf has a finite certified complex enclosure");
            if (value) values.push_back(value.value());
        }
        EXPECT_TRUE(
            values.size() == 5 &&
                std::abs(values[0].imag.value) < 1e-12 &&
                std::abs(values[1].real.value - values[2].real.value) <
                    1e-10 &&
                std::abs(values[1].imag.value + values[2].imag.value) <
                    1e-10 &&
                std::abs(values[3].real.value - values[4].real.value) <
                    1e-10 &&
                std::abs(values[3].imag.value + values[4].imag.value) <
                    1e-10,
            "generic quintic roots are one real root plus certified conjugate pairs");

        auto selected = SymbolicExpr::root_of(expression, "x", 1);
        LMCAS::CancellationToken cancellation;
        cancellation.cancel();
        LMCAS::ComputationContext cancelled({}, cancellation);
        auto cancelled_value =
            LMCAS::rootof_evaluate_complex_checked(selected, cancelled);
        EXPECT_TRUE(
            !cancelled_value &&
                cancelled_value.error().code == LMCAS::CasErrc::Cancelled,
            "higher-degree complex isolation propagates cancellation");

        LMCAS::ResourceLimits limits;
        limits.max_steps = 8;
        LMCAS::ComputationContext limited(limits);
        auto limited_value =
            LMCAS::rootof_evaluate_complex_checked(selected, limited);
        EXPECT_TRUE(
            !limited_value &&
                limited_value.error().code == LMCAS::CasErrc::ResourceLimit,
            "higher-degree complex isolation enforces step budgets");

        auto repeated_expression = SymbolicExpr::power(expression, num(2));
        for (int index = 0; index < 5; ++index) {
            auto canonical = SymbolicExpr::root_of(expression, "x", index);
            auto repeated = SymbolicExpr::root_of(
                repeated_expression, "x", index);
            EXPECT_TRUE(
                LMCAS::detail::node(canonical)->equals(
                    *LMCAS::detail::node(repeated)),
                "square-free canonicalization preserves quintic complex identity");
        }
    }
    return TEST_REPORT();
}
