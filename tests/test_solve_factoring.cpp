#include "test_common.hpp"
#include "solve_polynomial.hpp"
#include "root_of_utils.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <algorithm>
#include <map>
#include <set>

using namespace LMCAS;

static Polynomial<Rational> poly_from_roots(const std::vector<Rational> &roots, const std::string &var = "x") {
    Polynomial<Rational> result({Rational(1)}, var);
    for (const auto &r : roots) {
        Polynomial<Rational> factor({-r, Rational(1)}, var);
        result = result * factor;
    }
    return result;
}

static Polynomial<Rational> reconstruct_square_free_factors(
    const std::vector<std::pair<Polynomial<Rational>, int>> &factors,
    const std::string &variable = "x") {
    Polynomial<Rational> reconstructed({Rational(1)}, variable);
    for (const auto &[factor, multiplicity] : factors) {
        for (int i = 0; i < multiplicity; ++i) {
            reconstructed = reconstructed * factor;
        }
    }
    return reconstructed;
}

static Polynomial<SymbolicPolyCoeff> to_symbolic_poly(const Polynomial<Rational> &rat_poly) {
    std::vector<SymbolicPolyCoeff> sym_coeffs;
    sym_coeffs.reserve(rat_poly.coeffs.size());
    for (const auto &c : rat_poly.coeffs) {
        sym_coeffs.push_back(SymbolicPolyCoeff(SymbolicExpr::number(c)));
    }
    return Polynomial<SymbolicPolyCoeff>(sym_coeffs, rat_poly.variable_name);
}

static int count_occurrences(const std::vector<Rational> &vec, const Rational &val) {
    int count = 0;
    for (const auto &v : vec) {
        if (v == val)
            count++;
    }
    return count;
}

static double eval_poly_at(const Polynomial<Rational> &poly, double x) {
    double result = 0.0;
    double x_pow = 1.0;
    for (size_t i = 0; i < poly.coeffs.size(); ++i) {
        result += poly.coeffs[i].to_double() * x_pow;
        x_pow *= x;
    }
    return result;
}

static bool is_rootof(const std::shared_ptr<SymbolicExpr> &expr) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    return std::dynamic_pointer_cast<const RootOfNode>(
               LMCAS::detail::node(expr)) != nullptr;
}

TEST(SolveFactoring, ZeroConstantFactorOutX) {
    auto p = poly_from_roots({Rational(0), Rational(0), Rational(0), Rational(2), Rational(-1)});

    auto roots = find_rational_roots(p);
    ASSERT_EQ(roots.size(), 5U);
    EXPECT_EQ(count_occurrences(roots, Rational(0)), 3);
    EXPECT_EQ(count_occurrences(roots, Rational(2)), 1);
    EXPECT_EQ(count_occurrences(roots, Rational(-1)), 1);

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    ASSERT_EQ(all_roots.size(), 5U);
    for (const auto &root : all_roots) {
        double value = test_numeric_value(root);
        ASSERT_TRUE(std::isfinite(value));
        EXPECT_NEAR(eval_poly_at(p, value), 0.0, 1e-8);
    }
}

TEST(SolveFactoring, ZeroConstantSingleX) {
    Polynomial<Rational> p({Rational(0), Rational(1), Rational(0), Rational(1)}, "x");

    auto roots = find_rational_roots(p);
    EXPECT_TRUE((roots.size() == 1)) << "find_rational_roots: only root 0 found";
    EXPECT_TRUE((roots[0] == Rational(0))) << "find_rational_roots: root is 0";
}

TEST(SolveFactoring, RepeatedRationalRoots) {
    auto p = poly_from_roots({Rational(1), Rational(1), Rational(1), Rational(-2), Rational(-2)});

    auto roots = find_rational_roots(p);
    EXPECT_TRUE((count_occurrences(roots, Rational(1)) == 3)) << "find_rational_roots: root 1 has multiplicity 3";

    EXPECT_TRUE((roots.size() == 5)) << "find_rational_roots: returns all 5 rational roots of (x-1)^3(x+2)^2";
    EXPECT_TRUE((count_occurrences(roots, Rational(-2)) == 2)) << "find_rational_roots: root -2 has multiplicity 2";

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    ASSERT_EQ(all_roots.size(), 5U);

    int count_1 = 0;
    int count_neg2 = 0;
    for (const auto &root : all_roots) {
        double value = test_numeric_value(root);
        ASSERT_TRUE(std::isfinite(value));
        EXPECT_NEAR(eval_poly_at(p, value), 0.0, 1e-8);
        if (std::abs(value - 1.0) < 1e-8) {
            ++count_1;
        }
        if (std::abs(value + 2.0) < 1e-8) {
            ++count_neg2;
        }
    }
    EXPECT_EQ(count_1, 3);
    EXPECT_EQ(count_neg2, 2);
}

TEST(SolveFactoring, RepeatedRationalRootsFractional) {
    auto p = poly_from_roots({Rational(1, 2), Rational(1, 2), Rational(3)});

    auto roots = find_rational_roots(p);
    EXPECT_TRUE((count_occurrences(roots, Rational(1, 2)) == 2)) << "find_rational_roots: root 1/2 has multiplicity 2";
    EXPECT_TRUE((roots.size() == 3)) << "find_rational_roots: returns all 3 rational roots of (x-1/2)^2(x-3)";

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    ASSERT_EQ(all_roots.size(), 3U);
    int half_count = 0;
    int three_count = 0;
    for (const auto &root : all_roots) {
        double value = test_numeric_value(root);
        ASSERT_TRUE(std::isfinite(value));
        EXPECT_NEAR(eval_poly_at(p, value), 0.0, 1e-8);
        half_count += std::abs(value - 0.5) < 1e-8;
        three_count += std::abs(value - 3.0) < 1e-8;
    }
    EXPECT_EQ(half_count, 2);
    EXPECT_EQ(three_count, 1);
}

TEST(SolveFactoring, FullyReducibleDegree6) {
    auto p = poly_from_roots({Rational(1), Rational(2), Rational(3),
                              Rational(-1), Rational(-2), Rational(-3)});

    auto roots = find_rational_roots(p);
    ASSERT_EQ(roots.size(), 6U);
    for (int root : {-3, -2, -1, 1, 2, 3}) {
        EXPECT_EQ(
            count_occurrences(roots, Rational(root)), 1);
    }

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    ASSERT_EQ(all_roots.size(), 6U);
    for (const auto &root : all_roots) {
        double value = test_numeric_value(root);
        ASSERT_TRUE(std::isfinite(value));
        EXPECT_NEAR(eval_poly_at(p, value), 0.0, 1e-8);
    }
}

TEST(SolveFactoring, FullyReducibleDegree6WithFractions) {
    auto p = poly_from_roots({Rational(1, 2), Rational(1, 3), Rational(-1),
                              Rational(2), Rational(-5), Rational(7)});

    auto rational_roots = find_rational_roots(p);
    ASSERT_EQ(rational_roots.size(), 6U);
    for (const auto &expected : {
             Rational(1, 2), Rational(1, 3), Rational(-1),
             Rational(2), Rational(-5), Rational(7)}) {
        EXPECT_EQ(
            count_occurrences(rational_roots, expected), 1);
    }
    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    ASSERT_EQ(all_roots.size(), 6U);

    for (const auto &root : all_roots) {
        double value = test_numeric_value(root);
        ASSERT_TRUE(std::isfinite(value));
        EXPECT_NEAR(eval_poly_at(p, value), 0.0, 1e-6);
    }
}

TEST(SolveFactoring, PartiallyReducibleLinearPlusIrreducible) {
    Polynomial<Rational> linear({Rational(-1), Rational(1)}, "x");
    Polynomial<Rational> quintic({Rational(1), Rational(1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
    auto p = linear * quintic;

    auto roots = find_rational_roots(p);
    EXPECT_TRUE((count_occurrences(roots, Rational(1)) == 1)) << "find_rational_roots: finds rational root 1";

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    ASSERT_EQ(all_roots.size(), 6U);

    int rootof_count = 0;
    int numeric_count = 0;
    for (const auto &root : all_roots) {
        if (is_rootof(root)) {
            ++rootof_count;
            continue;
        }
        double value = test_numeric_value(root);
        ASSERT_TRUE(std::isfinite(value));
        EXPECT_NEAR(eval_poly_at(p, value), 0.0, 1e-8);
        EXPECT_NEAR(value, 1.0, 1e-8);
        ++numeric_count;
    }
    EXPECT_EQ(rootof_count, 5);
    EXPECT_EQ(numeric_count, 1);
}

TEST(SolveFactoring, PartiallyReducibleTwoLinearPlusIrreducible) {
    Polynomial<Rational> lin1({Rational(-2), Rational(1)}, "x");
    Polynomial<Rational> lin2({Rational(3), Rational(1)}, "x");

    Polynomial<Rational> quintic({Rational(1), Rational(0), Rational(0), Rational(0), Rational(1), Rational(1)}, "x");
    auto p = lin1 * lin2 * quintic;

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    ASSERT_EQ(all_roots.size(), 7U);

    int rational_verified = 0;
    for (const auto &root : all_roots) {
        if (is_rootof(root)) {
            continue;
        }
        double value = test_numeric_value(root);
        ASSERT_TRUE(std::isfinite(value));
        EXPECT_NEAR(eval_poly_at(p, value), 0.0, 1e-8);
        EXPECT_TRUE(std::abs(value - 2.0) < 1e-8 ||
                    std::abs(value + 3.0) < 1e-8);
        ++rational_verified;
    }
    EXPECT_EQ(rational_verified, 2);
}

TEST(SolveFactoring, SquareFreeAlready) {
    auto polynomial =
        poly_from_roots({Rational(1), Rational(2), Rational(3)});

    auto factors = square_free_factorization(polynomial);

    ASSERT_EQ(factors.size(), 1U);
    EXPECT_EQ(factors[0].second, 1);
    EXPECT_TRUE(factors[0].first == polynomial);
    EXPECT_TRUE(reconstruct_square_free_factors(factors) == polynomial);
}

TEST(SolveFactoring, SquareFreeWithRepeated) {
    auto polynomial = poly_from_roots(
        {Rational(1), Rational(1), Rational(2),
         Rational(2), Rational(2)});

    auto factors = square_free_factorization(polynomial);

    ASSERT_EQ(factors.size(), 2U);
    std::set<int> multiplicities;
    for (const auto &[factor, multiplicity] : factors) {
        EXPECT_EQ(factor.degree(), 1);
        multiplicities.insert(multiplicity);
    }
    EXPECT_EQ(multiplicities, (std::set<int>{2, 3}));
    EXPECT_TRUE(reconstruct_square_free_factors(factors) == polynomial);
}

TEST(SolveFactoring, SquareFreeHighMultiplicity) {
    auto polynomial = poly_from_roots(
        {Rational(1), Rational(1), Rational(1), Rational(1)});

    auto factors = square_free_factorization(polynomial);

    ASSERT_EQ(factors.size(), 1U);
    EXPECT_EQ(factors[0].first.degree(), 1);
    EXPECT_EQ(factors[0].second, 4);
    EXPECT_TRUE(reconstruct_square_free_factors(factors) == polynomial);
}

TEST(SolveFactoring, SquareFreeIrreducible) {
    Polynomial<Rational> polynomial(
        {Rational(1), Rational(0), Rational(1)}, "x");

    auto factors = square_free_factorization(polynomial);

    ASSERT_EQ(factors.size(), 1U);
    EXPECT_EQ(factors[0].second, 1);
    EXPECT_TRUE(factors[0].first == polynomial);
}
