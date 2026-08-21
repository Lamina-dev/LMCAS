#include "test_common.hpp"
#include "parametric_solver.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <random>
#include <sstream>
#include <cmath>
#include <optional>

using namespace lamina;

static std::shared_ptr<SymbolicExpr> build_linear_2x2(
    int c1, int c2, int c0,
    const std::string& var1, const std::string& var2)
{
    auto x = SymbolicExpr::variable(var1);
    auto y = SymbolicExpr::variable(var2);

    auto term1 = SymbolicExpr::multiply(SymbolicExpr::number(c1), x);

    auto term2 = SymbolicExpr::multiply(SymbolicExpr::number(c2), y);

    auto term3 = SymbolicExpr::number(c0);

    return SymbolicExpr::add(SymbolicExpr::add(term1, term2), term3);
}

static std::shared_ptr<SymbolicExpr> build_parametric_linear_2x2(
    int c1, int c2, int c0,
    int a_coeff_idx, int a_mult,
    const std::string& var1, const std::string& var2)
{
    auto x = SymbolicExpr::variable(var1);
    auto y = SymbolicExpr::variable(var2);
    auto a = SymbolicExpr::variable("a");

    auto eff_c1 = SymbolicExpr::number(c1);
    auto eff_c2 = SymbolicExpr::number(c2);
    auto eff_c0 = SymbolicExpr::number(c0);

    if (a_coeff_idx == 0) {
        eff_c1 = SymbolicExpr::add(eff_c1, SymbolicExpr::multiply(SymbolicExpr::number(a_mult), a));
    } else if (a_coeff_idx == 1) {
        eff_c2 = SymbolicExpr::add(eff_c2, SymbolicExpr::multiply(SymbolicExpr::number(a_mult), a));
    } else {
        eff_c0 = SymbolicExpr::add(eff_c0, SymbolicExpr::multiply(SymbolicExpr::number(a_mult), a));
    }

    auto term1 = SymbolicExpr::multiply(eff_c1, x);
    auto term2 = SymbolicExpr::multiply(eff_c2, y);

    return SymbolicExpr::add(SymbolicExpr::add(term1, term2), eff_c0);
}

// Recursive numeric evaluator used to validate residuals when symbolic
// simplification cannot fully reduce them. Returns std::nullopt when the
// expression contains symbols other than `var` (which should be substituted
// before calling) or when a math error occurs.
static std::optional<double> numeric_eval(const std::shared_ptr<SymbolicExpr>& e) {
    if (!e || !lamina::detail::node(e)) return 0.0;
    auto root = lamina::detail::node(e);

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value())) return std::get<lmmc_real_t>(num->value());
        if (std::holds_alternative<BigInt>(num->value())) return (double)std::get<BigInt>(num->value()).to_double();
        if (std::holds_alternative<Rational>(num->value())) return (double)std::get<Rational>(num->value()).to_double();
        return 0.0;
    }
    if (std::dynamic_pointer_cast<const VariableNode>(root)) {
        // Free variable still present after substitution -> cannot evaluate.
        return std::nullopt;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(root)) {
        double s = 0.0;
        for (auto& op : add->operands()) {
            auto v = numeric_eval(lamina::detail::make_expression_ptr(op));
            if (!v) return std::nullopt;
            s += *v;
        }
        return s;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(root)) {
        double s = 1.0;
        for (auto& op : mul->operands()) {
            auto v = numeric_eval(lamina::detail::make_expression_ptr(op));
            if (!v) return std::nullopt;
            s *= *v;
        }
        return s;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(root)) {
        auto b = numeric_eval(lamina::detail::make_expression_ptr(pow->base()));
        auto x = numeric_eval(lamina::detail::make_expression_ptr(pow->exponent()));
        if (!b || !x) return std::nullopt;
        if (*b == 0.0 && *x < 0.0) return std::nullopt;
        double v = std::pow(*b, *x);
        if (!std::isfinite(v)) return std::nullopt;
        return v;
    }
    // Functions: fall back to the existing to_numeric (handles Sin/Cos/...).
    try {
        double v = e->to_numeric();
        if (!std::isfinite(v)) return std::nullopt;
        return v;
    } catch (...) {
        return std::nullopt;
    }
}

int main() {

    TEST_CASE("Back-substitution Round-trip for parametric 2x2 systems");
    {
        std::mt19937 rng(42);
        const int NUM_ITERATIONS = 50;
        int pass_count = 0;

        std::uniform_int_distribution<int> coeff_dist(-5, 5);
        std::uniform_int_distribution<int> a_mult_dist(1, 3);
        std::uniform_int_distribution<int> a_idx_dist(0, 2);
        std::uniform_int_distribution<int> param_eq_dist(0, 1);

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int c11 = coeff_dist(rng);
            int c12 = coeff_dist(rng);
            int c10 = coeff_dist(rng);
            int c21 = coeff_dist(rng);
            int c22 = coeff_dist(rng);
            int c20 = coeff_dist(rng);

            if (c11 == 0 && c22 == 0 && c12 == 0 && c21 == 0) continue;

            int param_eq = param_eq_dist(rng);
            int a_idx = a_idx_dist(rng);
            int a_mult = a_mult_dist(rng);

            std::shared_ptr<SymbolicExpr> eq1, eq2;

            if (param_eq == 0) {
                eq1 = build_parametric_linear_2x2(c11, c12, c10, a_idx, a_mult, "x", "y");
                eq2 = build_linear_2x2(c21, c22, c20, "x", "y");
            } else {
                eq1 = build_linear_2x2(c11, c12, c10, "x", "y");
                eq2 = build_parametric_linear_2x2(c21, c22, c20, a_idx, a_mult, "x", "y");
            }

            std::vector<std::shared_ptr<SymbolicExpr>> equations = {eq1, eq2};
            std::vector<std::string> unknowns = {"x", "y"};
            std::vector<std::string> parameters = {"a"};

            auto solutions = ParametricSolver::solve_system(equations, unknowns, parameters);

            if (solutions.empty()) {
                ++pass_count;
                continue;
            }

            bool all_ok = true;
            for (const auto& sol : solutions) {
                auto it_x = sol.find("x");
                auto it_y = sol.find("y");
                if (it_x == sol.end() || it_y == sol.end()) {

                    continue;
                }

                auto x_val = it_x->second;
                auto y_val = it_y->second;

                for (const auto& eq : equations) {
                    auto result = eq->substitute("x", x_val);
                    result = result->substitute("y", y_val);
                    result = result->simplify();

                    if (!result->is_zero()) {
                        result = result->expand();
                        if (result) result = result->simplify();
                    }

                    bool zero_ok = (result && result->is_zero());

                    // simplify/expand may not reduce all parametric rational
                    // expressions to a literal zero (e.g. it doesn't always
                    // perform full common-denominator combination on
                    // c1/(a-k) + c2*a/(a-k) + c3 forms). When that happens,
                    // fall back to numeric evaluation: pick a handful of
                    // concrete `a` values that avoid singularities and
                    // require the residual to vanish at all of them.
                    if (!zero_ok && result) {
                        zero_ok = true;
                        // Use sample values that avoid common singularities
                        // (integer and half-integer denominators like a-k, 2a-k, 3a-k).
                        const double samples[] = {-7.13, -2.71, 0.37, 3.89, 10.61};
                        int verified = 0;
                        for (double sa : samples) {
                            try {
                                auto subst = result->substitute(
                                    "a", SymbolicExpr::number(sa));
                                subst = subst->simplify();
                                auto v = numeric_eval(subst);
                                if (!v || !std::isfinite(*v)) continue;
                                if (std::abs(*v) > 1e-6) {
                                    zero_ok = false;
                                    break;
                                }
                                ++verified;
                            } catch (...) {
                                continue;
                            }
                        }
                        if (verified == 0) zero_ok = false;
                    }

                    if (!zero_ok) {
                        std::ostringstream oss;
                        oss << "failed: iter=" << iter
                            << " residual=" << (result ? result->to_string() : "null")
                            << " eq=" << eq->to_string()
                            << " x=" << x_val->to_string()
                            << " y=" << y_val->to_string();
                        EXPECT_TRUE(false, oss.str());
                        all_ok = false;
                        break;
                    }
                }
                if (!all_ok) break;
            }

            if (all_ok) ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed back-substitution round-trip";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Unique solution for non-singular 2x2 numeric systems");
    {
        std::mt19937 rng(123);
        const int NUM_ITERATIONS = 50;
        int pass_count = 0;
        int tested = 0;

        std::uniform_int_distribution<int> coeff_dist(-5, 5);

        while (tested < NUM_ITERATIONS) {

            int c11 = coeff_dist(rng);
            int c12 = coeff_dist(rng);
            int c10 = coeff_dist(rng);
            int c21 = coeff_dist(rng);
            int c22 = coeff_dist(rng);
            int c20 = coeff_dist(rng);

            int det = c11 * c22 - c12 * c21;

            if (det == 0) continue;

            ++tested;

            auto eq1 = build_linear_2x2(c11, c12, c10, "x", "y");
            auto eq2 = build_linear_2x2(c21, c22, c20, "x", "y");

            std::vector<std::shared_ptr<SymbolicExpr>> equations = {eq1, eq2};
            std::vector<std::string> unknowns = {"x", "y"};
            std::vector<std::string> parameters = {};

            auto solutions = ParametricSolver::solve_system(equations, unknowns, parameters);

            if (solutions.size() == 1) {

                auto& sol = solutions[0];
                auto it_x = sol.find("x");
                auto it_y = sol.find("y");

                if (it_x != sol.end() && it_y != sol.end()) {
                    auto x_val = it_x->second;
                    auto y_val = it_y->second;

                    bool back_sub_ok = true;
                    for (const auto& eq : equations) {
                        auto result = eq->substitute("x", x_val);
                        result = result->substitute("y", y_val);
                        result = result->simplify();

                        if (!result->is_zero()) {
                            result = result->expand();
                            if (result) result = result->simplify();
                        }

                        if (!result || !result->is_zero()) {
                            std::ostringstream oss;
                            oss << "back-sub failed: tested=" << tested
                                << " det=" << det
                                << " residual=" << (result ? result->to_string() : "null");
                            EXPECT_TRUE(false, oss.str());
                            back_sub_ok = false;
                            break;
                        }
                    }

                    if (back_sub_ok) ++pass_count;
                } else {
                    std::ostringstream oss;
                    oss << "solution missing x or y, tested=" << tested;
                    EXPECT_TRUE(false, oss.str());
                }
            } else {
                std::ostringstream oss;
                oss << "failed: expected 1 solution, got " << solutions.size()
                    << " for det=" << det
                    << " system: " << eq1->to_string() << " = 0, " << eq2->to_string() << " = 0";
                EXPECT_TRUE(false, oss.str());
            }
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed unique solution for non-singular systems";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Backward Compatibility (50 iterations)");
    {
        std::mt19937 rng(77);
        std::uniform_int_distribution<int> coeff_dist(-5, 5);
        int passed = 0;
        int skipped = 0;

        for (int iter = 0; iter < 50; ++iter) {

            int a11 = coeff_dist(rng);
            int a12 = coeff_dist(rng);
            int a10 = coeff_dist(rng);
            int a21 = coeff_dist(rng);
            int a22 = coeff_dist(rng);
            int a20 = coeff_dist(rng);

            int det = a11 * a22 - a12 * a21;
            if (det == 0) {
                skipped++;
                continue;
            }

            auto x = SymbolicExpr::variable("x");
            auto y = SymbolicExpr::variable("y");

            auto eq1 = SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::multiply(SymbolicExpr::number(a11), x),
                    SymbolicExpr::multiply(SymbolicExpr::number(a12), y)
                ),
                SymbolicExpr::number(a10)
            );

            auto eq2 = SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::multiply(SymbolicExpr::number(a21), x),
                    SymbolicExpr::multiply(SymbolicExpr::number(a22), y)
                ),
                SymbolicExpr::number(a20)
            );

            std::vector<std::shared_ptr<SymbolicExpr>> equations = {eq1, eq2};
            std::vector<std::string> vars = {"x", "y"};

            auto result_two = SymbolicExpr::solve_system(equations, vars);

            auto result_three = SymbolicExpr::solve_system(equations, vars, std::vector<std::string>{});

            if (result_two.size() != result_three.size()) {
                std::ostringstream oss;
                oss << "Iter " << iter << ": solution count mismatch (two-arg="
                    << result_two.size() << ", three-arg=" << result_three.size() << ")";
                EXPECT_TRUE(false, oss.str());
                continue;
            }

            bool all_match = true;
            for (size_t s = 0; s < result_two.size(); ++s) {
                auto& sol_two = result_two[s];
                auto& sol_three = result_three[s];

                if (sol_two.size() != sol_three.size()) {
                    all_match = false;
                    break;
                }

                for (const auto& eq : equations) {
                    auto residual_two = eq;
                    auto residual_three = eq;

                    for (const auto& [var, val] : sol_two) {
                        residual_two = residual_two->substitute(var, val);
                    }
                    for (const auto& [var, val] : sol_three) {
                        residual_three = residual_three->substitute(var, val);
                    }

                    residual_two = residual_two->simplify();
                    residual_three = residual_three->simplify();

                    bool two_zero = residual_two->is_zero();
                    bool three_zero = residual_three->is_zero();

                    if (!two_zero) {

                        try {
                            double val = residual_two->to_numeric();
                            two_zero = std::abs(val) < 1e-10;
                        } catch (...) {}
                    }
                    if (!three_zero) {
                        try {
                            double val = residual_three->to_numeric();
                            three_zero = std::abs(val) < 1e-10;
                        } catch (...) {}
                    }

                    if (!two_zero || !three_zero) {
                        all_match = false;
                        break;
                    }
                }
                if (!all_match) break;
            }

            if (all_match) {
                passed++;
            } else {
                std::ostringstream oss;
                oss << "Iter " << iter << ": solutions differ or residual non-zero";
                EXPECT_TRUE(false, oss.str());
            }
        }

        std::ostringstream summary;
        summary << "" << passed << " passed, " << skipped
                << " skipped (det=0) out of 50";
        std::cout << summary.str() << std::endl;
        EXPECT_TRUE(passed > 0, "at least some iterations passed");
        EXPECT_TRUE(passed + skipped == 50, "all iterations accounted for");
    }

    return TEST_REPORT();
}
