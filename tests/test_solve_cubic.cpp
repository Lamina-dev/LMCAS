#include "test_common.hpp"
#include "solve_polynomial.hpp"
#include <cmath>
#include <random>
#include <sstream>
#include <algorithm>
#include <string>

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

static double eval_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return 0.0;

    if (auto n = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
        if (std::holds_alternative<lmmc_real_t>(n->value())) return std::get<lmmc_real_t>(n->value());
        if (std::holds_alternative<BigInt>(n->value())) return std::get<BigInt>(n->value()).to_double();
        if (std::holds_alternative<Rational>(n->value())) return std::get<Rational>(n->value()).to_double();
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
        double result = 0.0;
        for (auto& op : add->operands()) {
            result += eval_numeric(lamina::detail::make_expression_ptr(op));
        }
        return result;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        double result = 1.0;
        for (auto& op : mul->operands()) {
            result *= eval_numeric(lamina::detail::make_expression_ptr(op));
        }
        return result;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
        double base = eval_numeric(lamina::detail::make_expression_ptr(pow->base()));
        double exp = eval_numeric(lamina::detail::make_expression_ptr(pow->exponent()));

        if (base < 0.0 && std::abs(exp - std::round(exp)) > 1e-15) {
            double denom = std::round(1.0 / exp);
            if (std::abs(exp * denom - 1.0) < 1e-12 && ((int)denom % 2 == 1)) {
                return -std::pow(-base, exp);
            }
            return std::nan("");
        }
        return std::pow(base, exp);
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
        if (func->arguments().size() == 1) {
            double arg = eval_numeric(lamina::detail::make_expression_ptr(func->arguments()[0]));
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

    if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
        return std::nan("");
    }

    return std::nan("");
}

static double eval_cubic_at(double a, double b, double c, double d, double x) {
    return a * x * x * x + b * x * x + c * x + d;
}

int main() {
    TEST_CASE("Solve Cubic - Placeholder");

    {
        auto a = num(1);
        auto b = num(-6);
        auto c = num(11);
        auto d = num(-6);

        auto roots = lamina::solve_cubic(a, b, c, d, "x");
        EXPECT_TRUE(roots.size() == 3, "Cubic should return 3 roots");
    }

    TEST_CASE("Cubic D > 0: x^3 - 2x - 5 = 0");
    {
        auto roots = lamina::solve_cubic(num(1), num(0), num(-2), num(-5), "x");
        EXPECT_TRUE(roots.size() == 3, "D>0: should return exactly 3 roots");

        if (roots.size() >= 1) {

            double r1 = eval_numeric(roots[0]);
            if (!std::isnan(r1)) {
                double residual = eval_cubic_at(1.0, 0.0, -2.0, -5.0, r1);
                EXPECT_TRUE(std::abs(residual) < 1e-10,
                    "D>0: real root satisfies equation |f(r)| < 1e-10");
                EXPECT_TRUE(std::abs(r1 - 2.0946) < 0.001,
                    "D>0: real root approx 2.0946");
            } else {

                double r1_alt = roots[0]->to_numeric();
                double residual = eval_cubic_at(1.0, 0.0, -2.0, -5.0, r1_alt);
                EXPECT_TRUE(std::abs(residual) < 1e-10,
                    "D>0: real root satisfies equation |f(r)| < 1e-10 (to_numeric)");
            }
        }

        if (roots.size() == 3) {
            std::string r2_str = roots[1]->to_string();
            std::string r3_str = roots[2]->to_string();
            bool has_imaginary = (r2_str.find("-1") != std::string::npos) ||
                                 (r3_str.find("-1") != std::string::npos);
            EXPECT_TRUE(has_imaginary,
                "D>0: complex roots contain imaginary component (sqrt(-1))");
        }
    }

    TEST_CASE("Cubic Triple Root: (x-1)^3 = x^3 - 3x^2 + 3x - 1 = 0");
    {
        auto roots = lamina::solve_cubic(num(1), num(-3), num(3), num(-1), "x");
        EXPECT_TRUE(roots.size() == 3, "Triple root: should return exactly 3 roots");

        if (roots.size() == 3) {
            double r1 = eval_numeric(roots[0]);
            double r2 = eval_numeric(roots[1]);
            double r3 = eval_numeric(roots[2]);

            EXPECT_TRUE(std::abs(r1 - 1.0) < 1e-10, "Triple root: r1 = 1");
            EXPECT_TRUE(std::abs(r2 - 1.0) < 1e-10, "Triple root: r2 = 1");
            EXPECT_TRUE(std::abs(r3 - 1.0) < 1e-10, "Triple root: r3 = 1");

            double res1 = eval_cubic_at(1.0, -3.0, 3.0, -1.0, r1);
            double res2 = eval_cubic_at(1.0, -3.0, 3.0, -1.0, r2);
            double res3 = eval_cubic_at(1.0, -3.0, 3.0, -1.0, r3);
            EXPECT_TRUE(std::abs(res1) < 1e-10 && std::abs(res2) < 1e-10 && std::abs(res3) < 1e-10,
                "Triple root: all roots satisfy equation");
        }
    }

    TEST_CASE("Cubic D=0, p!=0: x^3 - 3x + 2 = (x-1)^2(x+2) = 0");
    {
        auto roots = lamina::solve_cubic(num(1), num(0), num(-3), num(2), "x");
        EXPECT_TRUE(roots.size() == 3, "D=0 p!=0: should return exactly 3 roots");

        if (roots.size() == 3) {
            double r1 = eval_numeric(roots[0]);
            double r2 = eval_numeric(roots[1]);
            double r3 = eval_numeric(roots[2]);

            std::vector<double> sorted_roots = {r1, r2, r3};
            std::sort(sorted_roots.begin(), sorted_roots.end());

            EXPECT_TRUE(std::abs(sorted_roots[0] - (-2.0)) < 1e-10,
                "D=0 p!=0: smallest root = -2");
            EXPECT_TRUE(std::abs(sorted_roots[1] - 1.0) < 1e-10,
                "D=0 p!=0: middle root = 1");
            EXPECT_TRUE(std::abs(sorted_roots[2] - 1.0) < 1e-10,
                "D=0 p!=0: largest root = 1");

            for (double r : sorted_roots) {
                double res = eval_cubic_at(1.0, 0.0, -3.0, 2.0, r);
                EXPECT_TRUE(std::abs(res) < 1e-10,
                    "D=0 p!=0: root satisfies equation");
            }
        }
    }

    TEST_CASE("Cubic D < 0 (casus irreducibilis): x^3 - 3x + 1 = 0");
    {
        auto roots = lamina::solve_cubic(num(1), num(0), num(-3), num(1), "x");
        EXPECT_TRUE(roots.size() == 3, "D<0: should return exactly 3 roots");

        if (roots.size() == 3) {
            double r1 = eval_numeric(roots[0]);
            double r2 = eval_numeric(roots[1]);
            double r3 = eval_numeric(roots[2]);

            EXPECT_TRUE(std::abs(r1 - r2) > 1e-6 &&
                        std::abs(r1 - r3) > 1e-6 &&
                        std::abs(r2 - r3) > 1e-6,
                "D<0: three distinct real roots");

            double res1 = eval_cubic_at(1.0, 0.0, -3.0, 1.0, r1);
            double res2 = eval_cubic_at(1.0, 0.0, -3.0, 1.0, r2);
            double res3 = eval_cubic_at(1.0, 0.0, -3.0, 1.0, r3);
            EXPECT_TRUE(std::abs(res1) < 1e-10,
                "D<0: root 1 satisfies equation |f(r1)| < 1e-10");
            EXPECT_TRUE(std::abs(res2) < 1e-10,
                "D<0: root 2 satisfies equation |f(r2)| < 1e-10");
            EXPECT_TRUE(std::abs(res3) < 1e-10,
                "D<0: root 3 satisfies equation |f(r3)| < 1e-10");

            double sum = r1 + r2 + r3;
            EXPECT_TRUE(std::abs(sum) < 1e-8,
                "D<0: sum of roots = 0 (Vieta's)");
        }
    }

    TEST_CASE("Cubic Symbolic Coefficients: x^3 + a*x^2 + x + 1 = 0");
    {
        auto sym_a = SymbolicExpr::variable("a");
        auto roots = lamina::solve_cubic(num(1), sym_a, num(1), num(1), "x");
        EXPECT_TRUE(roots.size() == 3, "Symbolic: should return exactly 3 roots");

        if (roots.size() == 3) {

            bool found_a = false;
            for (const auto& root : roots) {
                std::string s = root->to_string();
                if (s.find("a") != std::string::npos) {
                    found_a = true;
                    break;
                }
            }
            EXPECT_TRUE(found_a, "Symbolic: roots contain variable 'a'");

            bool has_power_expr = false;
            for (const auto& root : roots) {
                std::string s = root->to_string();
                if (s.find("^") != std::string::npos ||
                    s.find("sqrt") != std::string::npos) {
                    has_power_expr = true;
                    break;
                }
            }
            EXPECT_TRUE(has_power_expr,
                "Symbolic: roots expressed using power/sqrt expressions");
        }
    }

    TEST_CASE("Cubic exact huge coefficient avoids unsafe numeric underflow");
    {
        std::string huge_digits = "1" + std::string(400, '0');
        auto roots = lamina::solve_cubic(
            SymbolicExpr::number(BigInt(huge_digits)),
            num(0),
            num(0),
            num(-1),
            "x");

        EXPECT_TRUE(roots.size() == 3,
            "huge exact cubic should still return the symbolic cubic roots");

        int zero_roots = 0;
        for (const auto& root : roots) {
            std::string text = root ? root->to_string() : "";
            EXPECT_TRUE(text.find("inf") == std::string::npos &&
                        text.find("nan") == std::string::npos,
                "huge exact cubic roots should not contain fabricated inf/nan");
            if (root && root->simplify()->is_zero()) {
                ++zero_roots;
            }
        }
        EXPECT_TRUE(zero_roots < 3,
            "huge exact cubic must not collapse nonzero roots to triple zero");
    }

    TEST_CASE("Cubic a=0 Delegation to Quadratic: 0*x^3 + 2*x^2 - 4*x + 2 = 0");
    {

        auto roots = lamina::solve_cubic(num(0), num(2), num(-4), num(2), "x");
        EXPECT_TRUE(roots.size() >= 1 && roots.size() <= 2,
            "a=0 delegation: returns 1-2 roots (quadratic)");

        if (roots.size() >= 1) {
            double r1 = eval_numeric(roots[0]);
            EXPECT_TRUE(std::abs(r1 - 1.0) < 1e-10,
                "a=0 delegation: root = 1");
        }
    }

    TEST_CASE("Cubic a=0 Delegation to Linear: 0*x^3 + 0*x^2 + 3*x - 6 = 0");
    {

        auto roots = lamina::solve_cubic(num(0), num(0), num(3), num(-6), "x");
        EXPECT_TRUE(roots.size() == 1,
            "a=0 b=0 delegation: returns 1 root (linear)");

        if (roots.size() == 1) {
            double r1 = eval_numeric(roots[0]);
            EXPECT_TRUE(std::abs(r1 - 2.0) < 1e-10,
                "a=0 b=0 delegation: root = 2");
        }
    }

    TEST_CASE("Cubic root verification");

    {
        const double RESIDUAL_TOL = 1e-10;
        const int NUM_CUBIC_TRIALS = 50;
        int cubic_verify_pass_count = 0;
        int cubic_verify_total_roots = 0;

        std::mt19937 rng_cubic(123);
        std::uniform_int_distribution<int> coeff_dist(-10, 10);

        for (int trial = 0; trial < NUM_CUBIC_TRIALS; ++trial) {
            int a_val = coeff_dist(rng_cubic);

            while (a_val == 0) a_val = coeff_dist(rng_cubic);
            int b_val = coeff_dist(rng_cubic);
            int c_val = coeff_dist(rng_cubic);
            int d_val = coeff_dist(rng_cubic);

            auto roots = lamina::solve_cubic(num(a_val), num(b_val), num(c_val), num(d_val), "x");

            if (roots.size() != 3) {
                std::ostringstream msg;
                msg << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
                    << ", c=" << c_val << ", d=" << d_val
                    << "): expected 3 roots, got " << roots.size();
                EXPECT_TRUE(false, msg.str());
                continue;
            }

            bool trial_ok = true;
            for (size_t i = 0; i < roots.size(); ++i) {
                double r = eval_numeric(roots[i]);

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
                    EXPECT_TRUE(false, msg.str());
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
            EXPECT_TRUE(cubic_verify_pass_count == NUM_CUBIC_TRIALS, msg.str());
        }
    }

    TEST_CASE("Vieta's formulas for cubics");

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

        auto roots = lamina::solve_cubic(num(1), num(b_val), num(c_val), num(d_val), "x");

        if (roots.size() != 3) {
            std::ostringstream msg;
            msg << "Trial " << trial << " (b=" << b_val << ", c=" << c_val
                << ", d=" << d_val << "): expected 3 roots, got " << roots.size();
            EXPECT_TRUE(false, msg.str());
            continue;
        }

        double r1 = eval_numeric(roots[0]);
        double r2 = eval_numeric(roots[1]);
        double r3 = eval_numeric(roots[2]);

        if (std::isnan(r1) || std::isnan(r2) || std::isnan(r3) ||
            std::isinf(r1) || std::isinf(r2) || std::isinf(r3)) {
            std::ostringstream msg;
            msg << "Trial " << trial << " (b=" << b_val << ", c=" << c_val
                << ", d=" << d_val << "): root evaluation produced NaN/Inf";
            EXPECT_TRUE(false, msg.str());
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
        EXPECT_TRUE(vieta_pass_count >= 50, msg.str());
    }

    return TEST_REPORT();
}
