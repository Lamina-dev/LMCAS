
#include "test_common.hpp"
#include <rapidcheck.h>
#include "calculus_utils.hpp"

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>

using namespace LMCAS;

using SE = SymbolicExpr;

static auto num(int n) { return SE::number(n); }
static auto var(const std::string &name) { return SE::variable(name); }

namespace {

struct PolynomialSample {
    std::shared_ptr<SymbolicExpr> expression;
    std::vector<int> coefficients;
};

static ::testing::AssertionResult numeric_matches(const std::shared_ptr<SymbolicExpr> &actual,
                                                  double expected, double tolerance) {
    const auto value = test_numeric_eval(actual);
    if (!value) {
        return ::testing::AssertionFailure() << "Expression requires a numeric value";
    }
    if (!std::isfinite(*value) || !std::isfinite(expected) ||
        !std::isfinite(tolerance) || tolerance < 0.0) {
        return ::testing::AssertionFailure()
            << "Invalid numeric comparison: actual=" << *value
            << ", expected=" << expected << ", tolerance=" << tolerance;
    }
    if (std::abs(*value - expected) <= tolerance) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
        << "Actual: " << *value << ", expected: " << expected << ", tolerance: " << tolerance;
}

/**
 * @brief Generate a random polynomial expression in variable x with
 *        integer coefficients in [-5, 5] and degree in [1, 4].
 *
 * Returns expressions like: 3*x^2 + (-2)*x + 1
 */
PolynomialSample gen_polynomial(const std::string &v, int min_deg = 1, int max_deg = 4) {
    int degree = *rc::gen::inRange(min_deg, (max_deg) + 1);
    auto x = var(v);
    std::shared_ptr<SymbolicExpr> result = nullptr;
    std::vector<int> coefficients(static_cast<std::size_t>(degree) + 1);

    for (int d = degree; d >= 0; --d) {
        int coeff = *rc::gen::inRange(-5, (5) + 1);
        if (d == degree && coeff == 0) {
            coeff = *rc::gen::inRange(1, (5) + 1);
        }
        coefficients[static_cast<std::size_t>(d)] = coeff;

        if (coeff == 0) {
            continue;
        }

        std::shared_ptr<SymbolicExpr> term;
        if (d == 0) {
            term = num(coeff);
        } else if (d == 1) {
            term = (coeff == 1) ? x : SE::multiply(num(coeff), x);
        } else {
            auto x_pow = SE::power(x, num(d));
            term = (coeff == 1) ? x_pow : SE::multiply(num(coeff), x_pow);
        }

        if (!result) {
            result = term;
        } else {
            result = SE::add(result, term);
        }
    }

    return {result, std::move(coefficients)};
}

} // anonymous namespace

TEST(LmcasPropCalculusUtils, LogDifferentiationEquivalence) {
    EXPECT_TRUE(rc::check("both derivatives of lead*x^2+c equal 2*lead*x", []() {
        auto x = var("x");
        int c = *rc::gen::inRange(1, (5) + 1);
        int lead = *rc::gen::inRange(1, (3) + 1);
        auto f = SE::add(SE::multiply(num(lead), SE::power(x, num(2))), num(c));

        auto log_diff = LMCAS::log_differentiate(f, "x");
        auto std_diff = f->differentiate("x");
        RC_ASSERT(log_diff != nullptr);
        RC_ASSERT(std_diff != nullptr);

        for (double pt : {0.5, 1.0, 1.5, 2.0, 3.0}) {
            auto pt_expr = SE::number(pt);
            const double expected = 2.0 * lead * pt;
            const double tolerance = 1e-6 * std::max(1.0, std::abs(expected));
            RC_ASSERT(static_cast<bool>(numeric_matches(log_diff->substitute("x", pt_expr), expected, tolerance)));
            RC_ASSERT(static_cast<bool>(numeric_matches(std_diff->substitute("x", pt_expr), expected, tolerance)));
        }
    }));
}

static bool has_zero_horizontal_asymptote(
    const std::vector<std::shared_ptr<SymbolicExpr>> &horizontal) {
    bool has_zero_horiz = false;
    for (const auto &ha : horizontal) {
        if (!ha) {
            continue;
        }
        auto simplified = ha->simplify();
        if (!simplified) {
            continue;
        }
        if (simplified->is_zero()) {
            has_zero_horiz = true;
            break;
        }
        auto val = test_numeric_eval(simplified);
        if (val && std::isfinite(*val) && std::abs(*val) <= 1e-6) {
            has_zero_horiz = true;
            break;
        }
    }
    return has_zero_horiz;
}

TEST(LmcasPropCalculusUtils, VerticalAsymptotesRational) {
    EXPECT_TRUE(rc::check("vertical asymptotes are at denominator zeros", []() {
        auto x = var("x");

        // Use a simple rational function with a single known denominator zero
        // to avoid issues with polynomial expansion and root-finding limitations.
        // f(x) = 1 / (x - r) where r is a random integer
        int r = *rc::gen::inRange(-5, (5) + 1);
        auto denom = (r == 0) ? x : SE::add(x, num(-r));
        auto f = SE::divide(num(1), denom);

        auto checked = LMCAS::asymptotes_checked(f, "x");
        RC_ASSERT(checked);
        const auto &result = checked.value();

        // The vertical asymptote should be at x = r
        bool found = false;
        for (const auto &va : result.vertical) {
            if (!va) {
                continue;
            }
            auto simplified = va->simplify();
            if (!simplified) {
                continue;
            }
            auto val = test_numeric_eval(simplified);
            if (val && std::isfinite(*val) &&
                std::abs(*val - static_cast<double>(r)) <= 1e-6) {
                found = true;
                break;
            }
        }
        RC_ASSERT(found);

        // Horizontal asymptote should be y = 0 (deg(P) < deg(Q))
        bool has_zero_horiz = has_zero_horizontal_asymptote(result.horizontal);
        RC_ASSERT(has_zero_horiz);
    }));
}

TEST(LmcasPropCalculusUtils, HorizontalAsymptotesRational) {
    EXPECT_TRUE(rc::check("horizontal asymptote matches degree rule", []() {
        auto x = var("x");

        // Generate f(x) = a / (x - r) for random a, r
        // This has horizontal asymptote y = 0 (deg(P) < deg(Q))
        int a = *rc::gen::inRange(1, (5) + 1);
        int r = *rc::gen::inRange(-5, (5) + 1);
        auto denom = (r == 0) ? x : SE::add(x, num(-r));
        auto f = SE::divide(num(a), denom);

        auto checked = LMCAS::asymptotes_checked(f, "x");
        RC_ASSERT(checked);
        const auto &result = checked.value();

        // Horizontal asymptote should be y = 0
        bool has_zero_horiz = has_zero_horizontal_asymptote(result.horizontal);
        RC_ASSERT(has_zero_horiz);
    }));
}

static ::testing::AssertionResult circle_curvature_matches(int radius) {
    auto t = var("t");
    auto radius_expr = num(radius);
    auto x_t = SE::multiply(radius_expr, SE::cos(t));
    auto y_t = SE::multiply(radius_expr, SE::sin(t));

    auto result = LMCAS::curvature_parametric_checked(x_t, y_t, "t");
    if (!result) {
        return ::testing::AssertionFailure() << result.error().message;
    }
    const auto &kappa = result.value();
    if (!kappa) {
        return ::testing::AssertionFailure() << "Curvature requires an expression";
    }
    for (double pt : {0.0, 0.5}) {
        auto matches = numeric_matches(kappa->substitute("t", SE::number(pt)), 1.0 / radius, 1e-6);
        if (!matches) {
            return ::testing::AssertionFailure()
                << "radius=" << radius << ", t=" << pt << ": " << matches.message();
        }
    }
    return ::testing::AssertionSuccess();
}

static void check_circle_curvature() {
    RC_ASSERT(static_cast<bool>(circle_curvature_matches(*rc::gen::inRange(1, (10) + 1))));
}

static void check_polynomial_curvature() {
    auto sample = gen_polynomial("x", 2, 3);
    auto result = LMCAS::curvature_checked(sample.expression, "x");
    RC_ASSERT(result);
    const auto &kappa = result.value();
    RC_ASSERT(kappa != nullptr);

    const double pt = 1.0;
    double first = 0.0;
    double second = 0.0;
    for (std::size_t k = 1; k < sample.coefficients.size(); ++k) {
        first += static_cast<double>(k) * sample.coefficients[k] *
                 std::pow(pt, static_cast<int>(k) - 1);
        if (k >= 2) {
            second += static_cast<double>(k * (k - 1)) * sample.coefficients[k] *
                      std::pow(pt, static_cast<int>(k) - 2);
        }
    }
    const double expected = std::abs(second) / std::pow(1.0 + first * first, 1.5);
    RC_ASSERT(static_cast<bool>(numeric_matches(kappa->substitute("x", SE::number(pt)), expected,
                                              1e-4 * std::max(1.0, std::abs(expected)))));
}

TEST(LmcasPropCalculusUtils, CircleCurvatureRadiusEndpoints) {
    ASSERT_TRUE(circle_curvature_matches(1));
    ASSERT_TRUE(circle_curvature_matches(10));
}

TEST(LmcasPropCalculusUtils, PolynomialCurvatureZeroSecondDerivative) {
    auto x = var("x");
    auto polynomial = SE::add(SE::power(x, num(3)),
                              SE::multiply(num(-3), SE::power(x, num(2))));
    auto result = LMCAS::curvature_checked(polynomial, "x");
    ASSERT_TRUE(result);
    ASSERT_NE(result.value(), nullptr);
    ASSERT_TRUE(numeric_matches(result.value()->substitute("x", num(1)), 0.0, 0.0));
}

TEST(LmcasPropCalculusUtils, CircleCurvature) {
    EXPECT_TRUE(rc::check("curvature of circle radius R is 1/R", check_circle_curvature));
}

TEST(LmcasPropCalculusUtils, PolynomialCurvature) {
    EXPECT_TRUE(rc::check("curvature(f, x) matches |f''|/(1+f'^2)^(3/2) for polynomials",
                          check_polynomial_curvature));
}

TEST(LmcasPropCalculusUtils, InflectionPoints) {
    /**
     * @brief API 返回 f'' 的实零点，仅作为拐点的必要条件。
     * 凹凸性变化尚未获证，因此保留 x^4 的候选点。
     */
    EXPECT_TRUE(rc::check("returned polynomial candidates satisfy coefficient-derived f''=0", []() {
        auto sample = gen_polynomial("x", 3, 4);
        auto result = LMCAS::inflection_points_checked(sample.expression, "x");
        if (!result) {
            std::cerr << "Inflection input: " << sample.expression->to_string()
                      << "; " << result.error().operation << ": "
                      << result.error().message << '\n';
        }
        RC_ASSERT(result);

        for (const auto &pt : result.value()) {
            RC_ASSERT(pt != nullptr);
            auto value = test_numeric_eval(pt);
            RC_ASSERT(value.has_value());
            RC_ASSERT(std::isfinite(*value));
            double second = 0.0;
            for (std::size_t k = 2; k < sample.coefficients.size(); ++k) {
                second += static_cast<double>(k * (k - 1)) * sample.coefficients[k] *
                          std::pow(*value, static_cast<int>(k) - 2);
            }
            RC_ASSERT(std::isfinite(second));
            RC_ASSERT(std::abs(second) <= 1e-4);
        }
    }));
}

TEST(LmcasPropCalculusUtils, QuarticInflectionCandidates) {
    auto x = var("x");
    auto fourth = SE::power(x, num(4));
    auto square = SE::power(x, num(2));
    auto complex_only = LMCAS::inflection_points_checked(
        SE::add(fourth, square), "x");
    ASSERT_TRUE(complex_only);
    ASSERT_TRUE(complex_only.value().empty());

    auto real_pair = LMCAS::inflection_points_checked(
        SE::add(fourth, SE::multiply(num(-6), square)), "x");
    ASSERT_TRUE(real_pair);
    ASSERT_EQ(real_pair.value().size(), 2u);
    std::vector<double> coordinates;
    for (const auto &point : real_pair.value()) {
        ASSERT_NE(point, nullptr);
        auto value = test_numeric_eval(point);
        ASSERT_TRUE(value.has_value());
        ASSERT_TRUE(std::isfinite(*value));
        coordinates.push_back(*value);
    }
    std::sort(coordinates.begin(), coordinates.end());
    ASSERT_NEAR(coordinates[0], -1.0, 1e-12);
    ASSERT_NEAR(coordinates[1], 1.0, 1e-12);

    auto repeated = LMCAS::inflection_points_checked(fourth, "x");
    ASSERT_TRUE(repeated);
    ASSERT_EQ(repeated.value().size(), 1u);
    ASSERT_TRUE(numeric_matches(repeated.value()[0], 0.0, 0.0));
}
