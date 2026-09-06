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

static Polynomial<Rational> poly_from_roots(const std::vector<Rational>& roots, const std::string& var = "x") {
    Polynomial<Rational> result({Rational(1)}, var);
    for (const auto& r : roots) {
        Polynomial<Rational> factor({-r, Rational(1)}, var);
        result = result * factor;
    }
    return result;
}

static Polynomial<SymbolicPolyCoeff> to_symbolic_poly(const Polynomial<Rational>& rat_poly) {
    std::vector<SymbolicPolyCoeff> sym_coeffs;
    sym_coeffs.reserve(rat_poly.coeffs.size());
    for (const auto& c : rat_poly.coeffs) {
        sym_coeffs.push_back(SymbolicPolyCoeff(SymbolicExpr::number(c)));
    }
    return Polynomial<SymbolicPolyCoeff>(sym_coeffs, rat_poly.variable_name);
}

static int count_occurrences(const std::vector<Rational>& vec, const Rational& val) {
    int count = 0;
    for (const auto& v : vec) {
        if (v == val) count++;
    }
    return count;
}

static double eval_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) return 0.0;

    if (auto n = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
        if (std::holds_alternative<lmmc_real_t>(n->value())) return std::get<lmmc_real_t>(n->value());
        if (std::holds_alternative<BigInt>(n->value())) return std::get<BigInt>(n->value()).to_double();
        if (std::holds_alternative<Rational>(n->value())) return std::get<Rational>(n->value()).to_double();
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        double result = 0.0;
        for (auto& op : add->operands()) {
            result += eval_numeric(LMCAS::detail::make_expression_ptr(op));
        }
        return result;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        double result = 1.0;
        for (auto& op : mul->operands()) {
            result *= eval_numeric(LMCAS::detail::make_expression_ptr(op));
        }
        return result;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        double base = eval_numeric(LMCAS::detail::make_expression_ptr(pow->base()));
        double exp = eval_numeric(LMCAS::detail::make_expression_ptr(pow->exponent()));
        if (base < 0.0 && std::abs(exp - std::round(exp)) > 1e-15) {
            double denom = std::round(1.0 / exp);
            if (std::abs(exp * denom - 1.0) < 1e-12 && ((int)denom % 2 == 1)) {
                return -std::pow(-base, exp);
            }
            return std::nan("");
        }
        return std::pow(base, exp);
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        if (func->arguments().size() == 1) {
            double arg = eval_numeric(LMCAS::detail::make_expression_ptr(func->arguments()[0]));
            switch (func->type()) {
                case FunctionNode::FuncType::Sin: return std::sin(arg);
                case FunctionNode::FuncType::Cos: return std::cos(arg);
                case FunctionNode::FuncType::Tan: return std::tan(arg);
                case FunctionNode::FuncType::Exp: return std::exp(arg);
                case FunctionNode::FuncType::Ln: return std::log(arg);
                case FunctionNode::FuncType::Sqrt:
                    if (arg < 0.0) return std::nan("");
                    return std::sqrt(arg);
                case FunctionNode::FuncType::Abs: return std::abs(arg);
                case FunctionNode::FuncType::ArcCos: return std::acos(arg);
                case FunctionNode::FuncType::ArcSin: return std::asin(arg);
                case FunctionNode::FuncType::ArcTan: return std::atan(arg);
                default: break;
            }
        }
    }

    if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
        return std::nan("");
    }

    return std::nan("");
}

static double eval_poly_at(const Polynomial<Rational>& poly, double x) {
    double result = 0.0;
    double x_pow = 1.0;
    for (size_t i = 0; i < poly.coeffs.size(); ++i) {
        result += poly.coeffs[i].to_double() * x_pow;
        x_pow *= x;
    }
    return result;
}

static bool is_rootof(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) return false;
    return std::dynamic_pointer_cast<const RootOfNode>(
               LMCAS::detail::node(expr)) != nullptr;
}

void test_zero_constant_factor_out_x() {
    TEST_CASE("Zero constant term: x^3(x-2)(x+1) - factor out x^3");

    auto p = poly_from_roots({Rational(0), Rational(0), Rational(0), Rational(2), Rational(-1)});

    auto roots = find_rational_roots(p);
    EXPECT_TRUE(count_occurrences(roots, Rational(0)) == 3,
        "find_rational_roots: x^3 factor gives root 0 with multiplicity 3");

    EXPECT_TRUE(roots.size() >= 4,
        "find_rational_roots: finds at least 4 roots (3 zeros + at least 1 from quotient)");

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    EXPECT_TRUE(all_roots.size() == 5,
        "solve_by_factoring: returns all 5 roots for x^3(x-2)(x+1)");

    int verified = 0;
    for (const auto& root : all_roots) {
        double r_val = eval_numeric(root);
        if (!std::isnan(r_val)) {
            double residual = eval_poly_at(p, r_val);
            if (std::abs(residual) < 1e-8) {
                verified++;
            }
        }
    }
    EXPECT_TRUE(verified == 5,
        "solve_by_factoring: all 5 roots satisfy the polynomial");
}

void test_zero_constant_single_x() {
    TEST_CASE("Zero constant term: x(x^2+1) - factor out single x");

    Polynomial<Rational> p({Rational(0), Rational(1), Rational(0), Rational(1)}, "x");

    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 1, "find_rational_roots: only root 0 found");
    EXPECT_TRUE(roots[0] == Rational(0), "find_rational_roots: root is 0");
}

void test_repeated_rational_roots() {
    TEST_CASE("Repeated rational roots: (x-1)^3(x+2)^2");

    auto p = poly_from_roots({Rational(1), Rational(1), Rational(1), Rational(-2), Rational(-2)});

    auto roots = find_rational_roots(p);
    EXPECT_TRUE(count_occurrences(roots, Rational(1)) == 3,
        "find_rational_roots: root 1 has multiplicity 3");

    EXPECT_TRUE(roots.size() == 5,
        "find_rational_roots: returns all 5 rational roots of (x-1)^3(x+2)^2");
    EXPECT_TRUE(count_occurrences(roots, Rational(-2)) == 2,
        "find_rational_roots: root -2 has multiplicity 2");

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    EXPECT_TRUE(all_roots.size() == 5,
        "solve_by_factoring: returns all 5 roots for (x-1)^3(x+2)^2");

    int count_1 = 0, count_neg2 = 0;
    for (const auto& root : all_roots) {
        double r_val = eval_numeric(root);
        if (!std::isnan(r_val)) {
            if (std::abs(r_val - 1.0) < 1e-8) count_1++;
            if (std::abs(r_val + 2.0) < 1e-8) count_neg2++;
        }
    }
    EXPECT_TRUE(count_1 == 3,
        "solve_by_factoring: root 1 appears 3 times");
    EXPECT_TRUE(count_neg2 == 2,
        "solve_by_factoring: root -2 appears 2 times");
}

void test_repeated_rational_roots_fractional() {
    TEST_CASE("Repeated rational roots: (x-1/2)^2(x-3)");

    auto p = poly_from_roots({Rational(1, 2), Rational(1, 2), Rational(3)});

    auto roots = find_rational_roots(p);
    EXPECT_TRUE(count_occurrences(roots, Rational(1, 2)) == 2,
        "find_rational_roots: root 1/2 has multiplicity 2");
    EXPECT_TRUE(roots.size() == 3,
        "find_rational_roots: returns all 3 rational roots of (x-1/2)^2(x-3)");

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    EXPECT_TRUE(all_roots.size() == 3,
        "solve_by_factoring: returns all 3 roots for (x-1/2)^2(x-3)");
}

void test_fully_reducible_degree6() {
    TEST_CASE("Fully reducible degree-6: (x-1)(x-2)(x-3)(x+1)(x+2)(x+3)");

    auto p = poly_from_roots({Rational(1), Rational(2), Rational(3),
                              Rational(-1), Rational(-2), Rational(-3)});

    auto roots = find_rational_roots(p);

    EXPECT_TRUE(roots.size() >= 2,
        "find_rational_roots: finds at least 2 roots before stopping at degree 4 quotient");

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    EXPECT_TRUE(all_roots.size() == 6,
        "solve_by_factoring: returns exactly 6 roots for degree-6 fully reducible");

    int verified = 0;
    for (const auto& root : all_roots) {
        double r_val = eval_numeric(root);
        if (!std::isnan(r_val)) {
            double residual = eval_poly_at(p, r_val);
            if (std::abs(residual) < 1e-8) {
                verified++;
            }
        }
    }
    EXPECT_TRUE(verified == 6,
        "solve_by_factoring: all 6 roots satisfy the polynomial numerically");
}

void test_fully_reducible_degree6_with_fractions() {
    TEST_CASE("Fully reducible degree-6 with fractional roots: (x-1/2)(x-1/3)(x+1)(x+2)(x-5)(x+7)");

    auto p = poly_from_roots({Rational(1, 2), Rational(1, 3), Rational(-1),
                              Rational(2), Rational(-5), Rational(7)});

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    EXPECT_TRUE(all_roots.size() == 6,
        "solve_by_factoring: returns exactly 6 roots for degree-6 with fractions");

    int verified = 0;
    for (const auto& root : all_roots) {
        double r_val = eval_numeric(root);
        if (!std::isnan(r_val)) {
            double residual = eval_poly_at(p, r_val);
            if (std::abs(residual) < 1e-6) {
                verified++;
            }
        }
    }
    EXPECT_TRUE(verified == 6,
        "solve_by_factoring: all 6 fractional roots satisfy the polynomial");
}

void test_partially_reducible_linear_plus_irreducible() {
    TEST_CASE("Partially reducible: (x-1)(x^5+x+1) - one rational root + irreducible quintic");

    Polynomial<Rational> linear({Rational(-1), Rational(1)}, "x");
    Polynomial<Rational> quintic({Rational(1), Rational(1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
    auto p = linear * quintic;

    auto roots = find_rational_roots(p);
    EXPECT_TRUE(count_occurrences(roots, Rational(1)) == 1,
        "find_rational_roots: finds rational root 1");

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    EXPECT_TRUE(all_roots.size() == 6,
        "solve_by_factoring: returns 6 roots (1 rational + 5 RootOf)");

    int rootof_count = 0;
    int numeric_count = 0;
    for (const auto& root : all_roots) {
        if (is_rootof(root)) {
            rootof_count++;
        } else {
            double r_val = eval_numeric(root);
            if (!std::isnan(r_val) && std::abs(eval_poly_at(p, r_val)) < 1e-8) {
                numeric_count++;
            }
        }
    }
    EXPECT_TRUE(rootof_count == 5,
        "solve_by_factoring: 5 roots are RootOf (irreducible quintic)");
    EXPECT_TRUE(numeric_count >= 1,
        "solve_by_factoring: at least 1 root is a numeric value (the rational root)");
}

void test_partially_reducible_two_linear_plus_irreducible() {
    TEST_CASE("Partially reducible: (x-2)(x+3)(x^5+x^4+1) - two rational roots + irreducible quintic");

    Polynomial<Rational> lin1({Rational(-2), Rational(1)}, "x");
    Polynomial<Rational> lin2({Rational(3), Rational(1)}, "x");

    Polynomial<Rational> quintic({Rational(1), Rational(0), Rational(0), Rational(0), Rational(1), Rational(1)}, "x");
    auto p = lin1 * lin2 * quintic;

    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");
    EXPECT_TRUE(all_roots.size() == 7,
        "solve_by_factoring: returns 7 roots (2 rational + 5 RootOf)");

    int rational_verified = 0;
    for (const auto& root : all_roots) {
        if (!is_rootof(root)) {
            double r_val = eval_numeric(root);
            if (!std::isnan(r_val)) {
                if (std::abs(r_val - 2.0) < 1e-8 || std::abs(r_val + 3.0) < 1e-8) {
                    rational_verified++;
                }
            }
        }
    }
    EXPECT_TRUE(rational_verified >= 2,
        "solve_by_factoring: both rational roots (2 and -3) are found");
}

void test_square_free_already() {
    TEST_CASE("Square-free factorization: polynomial already square-free");

    auto p = poly_from_roots({Rational(1), Rational(2), Rational(3)});

    auto factors = square_free_factorization(p);

    EXPECT_TRUE(factors.size() == 1,
        "square_free_factorization: single factor for square-free poly");
    if (!factors.empty()) {
        EXPECT_TRUE(factors[0].second == 1,
            "square_free_factorization: multiplicity is 1");
        EXPECT_TRUE(factors[0].first.degree() == 3,
            "square_free_factorization: factor has degree 3");
    }
}

void test_square_free_with_repeated() {
    TEST_CASE("Square-free factorization: (x-1)^2(x-2)^3");

    auto p = poly_from_roots({Rational(1), Rational(1), Rational(2), Rational(2), Rational(2)});

    auto factors = square_free_factorization(p);

    EXPECT_TRUE(factors.size() >= 1,
        "square_free_factorization: returns at least one factor pair");

    int total_degree = 0;
    for (const auto& [factor, mult] : factors) {
        total_degree += factor.degree() * mult;
    }
    EXPECT_TRUE(total_degree == 5,
        "square_free_factorization: total degree (sum of factor.degree * mult) equals 5");

    std::set<int> mults;
    for (const auto& [factor, mult] : factors) {
        mults.insert(mult);
    }
    EXPECT_TRUE(mults.size() == factors.size(),
        "square_free_factorization: all multiplicities are distinct");
}

void test_square_free_high_multiplicity() {
    TEST_CASE("Square-free factorization: (x-1)^4");

    auto p = poly_from_roots({Rational(1), Rational(1), Rational(1), Rational(1)});

    auto factors = square_free_factorization(p);
    EXPECT_TRUE(!factors.empty(),
        "square_free_factorization: returns factors for (x-1)^4");

    int total_degree = 0;
    for (const auto& [factor, mult] : factors) {
        total_degree += factor.degree() * mult;
    }
    EXPECT_TRUE(total_degree == 4,
        "square_free_factorization: total degree is 4 for (x-1)^4");
}

void test_square_free_irreducible() {
    TEST_CASE("Square-free factorization: x^2+1 (irreducible, already square-free)");

    Polynomial<Rational> p({Rational(1), Rational(0), Rational(1)}, "x");

    auto factors = square_free_factorization(p);
    EXPECT_TRUE(factors.size() == 1,
        "square_free_factorization: single factor for irreducible x^2+1");
    if (!factors.empty()) {
        EXPECT_TRUE(factors[0].second == 1,
            "square_free_factorization: multiplicity is 1");
        EXPECT_TRUE(factors[0].first.degree() == 2,
            "square_free_factorization: factor has degree 2");
    }
}

void test_solve_by_factoring_zero_constant() {
    TEST_CASE("solve_by_factoring: x^2(x-1)(x-2)(x-3) - zero constant, degree 5");

    auto p = poly_from_roots({Rational(0), Rational(0), Rational(1), Rational(2), Rational(3)});
    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");

    EXPECT_TRUE(all_roots.size() == 5,
        "solve_by_factoring: returns 5 roots for x^2(x-1)(x-2)(x-3)");

    int verified = 0;
    for (const auto& root : all_roots) {
        double r_val = eval_numeric(root);
        if (!std::isnan(r_val)) {
            double residual = eval_poly_at(p, r_val);
            if (std::abs(residual) < 1e-8) {
                verified++;
            }
        }
    }
    EXPECT_TRUE(verified == 5,
        "solve_by_factoring: all 5 roots satisfy the polynomial");
}

void test_solve_by_factoring_repeated_roots() {
    TEST_CASE("solve_by_factoring: (x-1)^3(x+2)^2 - repeated roots, degree 5");

    auto p = poly_from_roots({Rational(1), Rational(1), Rational(1), Rational(-2), Rational(-2)});
    auto sym_poly = to_symbolic_poly(p);
    auto all_roots = solve_by_factoring(sym_poly, "x");

    EXPECT_TRUE(all_roots.size() == 5,
        "solve_by_factoring: returns 5 roots (counting multiplicity)");

    int count_1 = 0, count_neg2 = 0;
    for (const auto& root : all_roots) {
        double r_val = eval_numeric(root);
        if (!std::isnan(r_val)) {
            if (std::abs(r_val - 1.0) < 1e-8) count_1++;
            if (std::abs(r_val + 2.0) < 1e-8) count_neg2++;
        }
    }
    EXPECT_TRUE(count_1 == 3,
        "solve_by_factoring: root 1 appears 3 times");
    EXPECT_TRUE(count_neg2 == 2,
        "solve_by_factoring: root -2 appears 2 times");
}

int main() {
    std::cout << "=== Preprocessing Edge Cases Tests ===\n\n";

    test_zero_constant_factor_out_x();
    test_zero_constant_single_x();

    test_repeated_rational_roots();
    test_repeated_rational_roots_fractional();

    test_fully_reducible_degree6();
    test_fully_reducible_degree6_with_fractions();

    test_partially_reducible_linear_plus_irreducible();
    test_partially_reducible_two_linear_plus_irreducible();

    test_square_free_already();
    test_square_free_with_repeated();
    test_square_free_high_multiplicity();
    test_square_free_irreducible();

    test_solve_by_factoring_zero_constant();
    test_solve_by_factoring_repeated_roots();

    return TEST_REPORT();
}
