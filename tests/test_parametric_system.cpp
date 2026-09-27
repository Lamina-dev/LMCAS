#include "test_common.hpp"
#include "parametric_solver.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "expr.hpp"
#include <random>
#include <sstream>
#include <cmath>
#include <optional>

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> build_linear_2x2(
    int c1, int c2, int c0,
    const std::string &var1, const std::string &var2) {
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
    const std::string &var1, const std::string &var2) {
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

/// 递归数值求值器用于验证符号化简后仍保留结构的残差.
/// std::nullopt 表示表达式仍包含待替换符号或数值求值产生定义域诊断.
static std::optional<double> numeric_eval(const std::shared_ptr<SymbolicExpr> &e);

static std::optional<double> evaluate_system_add(const AddNode &add) {

    double s = 0.0;
    for (auto &op : add.operands()) {
        auto v = numeric_eval(LMCAS::detail::make_expression_ptr(op));
        if (!v) {
            return std::nullopt;
        }
        s += *v;
    }
    return s;
}

static std::optional<double> evaluate_system_mul(const MultiplyNode &mul) {

    double s = 1.0;
    for (auto &op : mul.operands()) {
        auto v = numeric_eval(LMCAS::detail::make_expression_ptr(op));
        if (!v) {
            return std::nullopt;
        }
        s *= *v;
    }
    return s;
}

static std::optional<double> evaluate_system_pow(const PowerNode &pow) {

    auto b = numeric_eval(LMCAS::detail::make_expression_ptr(pow.base()));
    auto x = numeric_eval(LMCAS::detail::make_expression_ptr(pow.exponent()));
    if (!b || !x) {
        return std::nullopt;
    }
    if (*b == 0.0 && *x < 0.0) {
        return std::nullopt;
    }
    double v = std::pow(*b, *x);
    if (!std::isfinite(v)) {
        return std::nullopt;
    }
    return v;
}

static std::optional<double> numeric_eval(const std::shared_ptr<SymbolicExpr> &e) {
    if (!e || !LMCAS::detail::node(e)) {
        return std::nullopt;
    }
    auto root = LMCAS::detail::node(e);

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(root)) {
        return test_numeric_number(*num);
    }
    if (std::dynamic_pointer_cast<const VariableNode>(root)) {
        return std::nullopt;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(root)) {
        return evaluate_system_add(*add);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(root)) {
        return evaluate_system_mul(*mul);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(root)) {
        return evaluate_system_pow(*pow);
    }
    // Functions: fall back to the existing to_numeric (handles Sin/Cos/...).
    try {
        double v = e->to_numeric();
        if (!std::isfinite(v)) {
            return std::nullopt;
        }
        return v;
    } catch (...) {
        return std::nullopt;
    }
}

static bool sampled_parametric_residual_is_zero(
    const std::shared_ptr<SymbolicExpr> &result) {

    bool zero_ok = true;
    const double samples[] = {-7.13, -2.71, 0.37, 3.89, 10.61};
    int verified = 0;
    for (double sa : samples) {
        try {
            auto subst = result->substitute(
                "a", SymbolicExpr::number(sa));
            subst = subst->simplify();
            auto v = numeric_eval(subst);
            if (!v || !std::isfinite(*v)) {
                continue;
            }
            if (std::abs(*v) > 1e-6) {
                zero_ok = false;
                break;
            }
            ++verified;
        } catch (...) {
            continue;
        }
    }
    if (verified == 0) {
        zero_ok = false;
    }
    return zero_ok;
}

static bool parametric_equation_residual_is_zero(
    const std::shared_ptr<SymbolicExpr> &eq,
    const std::shared_ptr<SymbolicExpr> &x_val,
    const std::shared_ptr<SymbolicExpr> &y_val, int iter) {

    auto result = eq->substitute("x", x_val);
    result = result->substitute("y", y_val);
    result = result->simplify();

    if (!result->is_zero()) {
        result = result->expand();
        if (result) {
            result = result->simplify();
        }
    }

    bool zero_ok = (result && result->is_zero());

    /**
     * @brief 符号化简未得到零时，在避开奇点的多个 a 值上检查残差均为零。
     * simplify/expand 可能未对 c1/(a-k) + c2*a/(a-k) + c3 等表达式完全通分。
     */
    if (!zero_ok && result) {
        zero_ok = sampled_parametric_residual_is_zero(result);
    }

    if (!zero_ok) {
        std::ostringstream oss;
        oss << "failed: iter=" << iter
            << " residual=" << (result ? result->to_string() : "null")
            << " eq=" << eq->to_string()
            << " x=" << x_val->to_string()
            << " y=" << y_val->to_string();
        ADD_FAILURE() << oss.str();
        return false;
    }
    return true;
}

static bool numeric_system_residuals_are_zero(
    const std::vector<std::shared_ptr<SymbolicExpr>> &equations,
    const std::shared_ptr<SymbolicExpr> &x_val,
    const std::shared_ptr<SymbolicExpr> &y_val, int tested, int det) {
    bool back_sub_ok = true;
    for (const auto &eq : equations) {
        auto result = eq->substitute("x", x_val);
        result = result->substitute("y", y_val);
        result = result->simplify();

        if (!result->is_zero()) {
            result = result->expand();
            if (result) {
                result = result->simplify();
            }
        }

        if (!result || !result->is_zero()) {
            std::ostringstream oss;
            oss << "back-sub failed: tested=" << tested
                << " det=" << det
                << " residual=" << (result ? result->to_string() : "null");
            ADD_FAILURE() << oss.str();
            back_sub_ok = false;
            break;
        }
    }

    return back_sub_ok;
}

static bool parametric_solutions_satisfy_equations(
    const ParametricSolutionList &solutions,
    const std::vector<std::shared_ptr<SymbolicExpr>> &equations, int iter) {
    bool all_ok = true;
    for (const auto &sol : solutions) {
        auto it_x = sol.find("x");
        auto it_y = sol.find("y");
        if (it_x == sol.end() || it_y == sol.end()) {

            continue;
        }

        auto x_val = it_x->second;
        auto y_val = it_y->second;

        for (const auto &eq : equations) {
            if (!parametric_equation_residual_is_zero(eq, x_val, y_val, iter)) {
                all_ok = false;
                break;
            }
        }
        if (!all_ok) {
            break;
        }
    }

    return all_ok;
}

TEST(ParametricSystem, BackSubstitutionRoundTripForParametric2x2Systems) {
    std::mt19937 rng(42);
    const int NUM_ITERATIONS = 50;

    std::uniform_int_distribution<int> coeff_dist(-5, 5);
    std::uniform_int_distribution<int> a_mult_dist(1, 3);
    std::uniform_int_distribution<int> a_idx_dist(0, 2);
    std::uniform_int_distribution<int> param_eq_dist(0, 1);

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);

        int c11 = coeff_dist(rng);
        int c12 = coeff_dist(rng);
        int c10 = coeff_dist(rng);
        int c21 = coeff_dist(rng);
        int c22 = coeff_dist(rng);
        int c20 = coeff_dist(rng);

        EXPECT_FALSE(c11 == 0 && c22 == 0 && c12 == 0 && c21 == 0)
            << "back-substitution round-trip requires a nonzero system";
        if (c11 == 0 && c22 == 0 && c12 == 0 && c21 == 0) {
            continue;
        }

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
            continue;
        }

        EXPECT_TRUE(parametric_solutions_satisfy_equations(
            solutions, equations, iter));
    }
}

TEST(ParametricSystem, UniqueSolutionForNonSingular2x2NumericSystems) {
    std::mt19937 rng(123);
    const int NUM_ITERATIONS = 50;
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

        if (det == 0) {
            continue;
        }

        ++tested;
        SCOPED_TRACE(::testing::Message() << "tested=" << tested << " det=" << det);

        auto eq1 = build_linear_2x2(c11, c12, c10, "x", "y");
        auto eq2 = build_linear_2x2(c21, c22, c20, "x", "y");

        std::vector<std::shared_ptr<SymbolicExpr>> equations = {eq1, eq2};
        std::vector<std::string> unknowns = {"x", "y"};
        std::vector<std::string> parameters = {};

        auto solutions = ParametricSolver::solve_system(equations, unknowns, parameters);

        if (solutions.size() == 1) {

            auto &sol = solutions[0];
            auto it_x = sol.find("x");
            auto it_y = sol.find("y");

            if (it_x != sol.end() && it_y != sol.end()) {
                auto x_val = it_x->second;
                auto y_val = it_y->second;

                EXPECT_TRUE(numeric_system_residuals_are_zero(
                    equations, x_val, y_val, tested, det));
            } else {
                std::ostringstream oss;
                oss << "solution missing x or y, tested=" << tested;
                ADD_FAILURE() << oss.str();
            }
        } else {
            std::ostringstream oss;
            oss << "failed: expected 1 solution, got " << solutions.size()
                << " for det=" << det
                << " system: " << eq1->to_string() << " = 0, " << eq2->to_string() << " = 0";
            ADD_FAILURE() << oss.str();
        }
    }

    EXPECT_EQ(tested, NUM_ITERATIONS);
}

TEST(ParametricSystem, CheckedParametricPolynomialContract) {
    auto invalid = ParametricSolver::solve_polynomial_parametric_checked(
        {}, {"x"}, {});
    EXPECT_TRUE((!invalid &&
                 invalid.error().code == CasErrc::InvalidArgument))
        << "empty parametric polynomial input is invalid";

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext context(limits);
    auto x = SymbolicExpr::variable("x");
    auto limited = ParametricSolver::solve_polynomial_parametric_checked(
        {SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-1))},
        {"x"}, {}, context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == CasErrc::ResourceLimit))
        << "parametric polynomial solve preserves exhausted budget";
    auto parameter = parse_expr("a*x-1");
    auto projection = solve_finite_checked(parameter.value(), "x");
    EXPECT_TRUE((!projection && projection.error().code == CasErrc::Inconclusive)) << "an unspecialized parameter pivot has no unconditional finite projection";
    auto zero_product = parse_expr("a*x");
    auto generic_branch = ParametricSolver::solve_polynomial_parametric_checked(
        {zero_product.value()}, {"x"}, {"a"});
    EXPECT_TRUE((!generic_branch && generic_branch.error().code == CasErrc::Inconclusive)) << "univariate branching preserves a possibly-zero leading factor";
    auto inconsistent = ParametricSolver::solve_polynomial_parametric_checked(
        {SymbolicExpr::power(x, SymbolicExpr::number(2)),
         SymbolicExpr::add(x, SymbolicExpr::number(-1))},
        {"x"}, {});
    EXPECT_TRUE((inconsistent && inconsistent.value().empty())) << "a proved empty polynomial branch cannot become a free unknown";
}

TEST(ParametricSystem, ParameterOnlyConstantTermIsPreservedExactly) {
    auto equation = parse_expr("2*x+a^2+sin(a)");
    auto expected = parse_expr("-(a^2+sin(a))/2");
    ASSERT_TRUE(equation) << equation.error().message;
    ASSERT_TRUE(expected) << expected.error().message;

    auto solutions = ParametricSolver::solve_system(
        {equation.value()}, {"x"}, {"a"});

    ASSERT_EQ(solutions.size(), 1U);
    auto value = solutions.front().find("x");
    ASSERT_NE(value, solutions.front().end());
    EXPECT_TRUE(test_proved_equivalent(value->second, expected.value()));
    auto residual = equation.value()
                        ->substitute("x", value->second)
                        ->simplify();
    ASSERT_NE(residual, nullptr);
    EXPECT_TRUE(residual->is_zero() ||
                test_proved_equivalent(
                    residual, SymbolicExpr::number(0)));
}
