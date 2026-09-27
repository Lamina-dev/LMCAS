#include "test_common.hpp"

#include "value.hpp"
#include "numeric_evaluation.hpp"

using namespace LMCAS;

Value cas_simplify(const std::vector<Value> &args) {
    if (args.empty()) {
        return Value(SymbolicExpr::number(0));
    }
    auto checked = args[0].as_symbolic_checked();
    if (!checked) throw checked.error();
    auto expr = checked.value();
    return Value(expr->simplify());
}

Value cas_differentiate(const std::vector<Value> &args) {
    if (args.size() < 2) {
        return Value(SymbolicExpr::number(0));
    }
    auto checked = args[0].as_symbolic_checked();
    if (!checked) throw checked.error();
    auto expr = checked.value();
    std::string var = args[1].to_string();
    return Value(expr->differentiate(var));
}

Value cas_solve(const std::vector<Value> &args) {
    if (args.size() < 2) {
        return Value();
    }
    auto checked = args[0].as_symbolic_checked();
    if (!checked) throw checked.error();
    auto expr = checked.value();
    std::string var = args[1].to_string();

    auto solutions = LMCAS::solve_finite_checked(expr, var).value();

    std::vector<Value> val_sols;
    for (const auto &s : solutions) {
        val_sols.push_back(Value(s));
    }
    return Value(val_sols);
}

TEST(Main, Simplify11) {
    auto one = SymbolicExpr::number(1);
    auto expr1 = SymbolicExpr::add(one, one);
    std::vector<Value> args = {Value(expr1)};
    Value res = cas_simplify(args);
    EXPECT_EQ((res.to_string()), ("2")) << "1+1 = 2";
}

TEST(Main, Sqrt4) {
    auto four = SymbolicExpr::number(4);
    auto expr = SymbolicExpr::sqrt(four);
    std::vector<Value> args = {Value(expr)};
    Value res = cas_simplify(args);
    EXPECT_EQ((res.to_string()), ("2")) << "sqrt(4) = 2";
}

TEST(Main, Sqrt8) {
    auto eight = SymbolicExpr::number(8);
    auto expr = SymbolicExpr::sqrt(eight);
    std::vector<Value> args = {Value(expr)};
    Value res = cas_simplify(args);

    auto actual = res.as_symbolic_checked();
    ASSERT_TRUE(actual);
    ASSERT_TRUE(detail::node(actual.value()));
    auto two = SymbolicExpr::number(2);
    auto expected = SymbolicExpr::multiply(two, SymbolicExpr::power(
                                                    two, SymbolicExpr::number(Rational(1, 2))))
                        ->simplify();
    EXPECT_TRUE((actual.value()->simplify()->compare(expected) == 0)) << "sqrt(8) = 2*sqrt(2)";
    auto numeric = test_numeric_eval(actual.value());
    ASSERT_TRUE((numeric.has_value())) << "sqrt(8) can be evaluated numerically";
    if (numeric) {
        {
            const double actual_value = (*numeric);
            const double expected_value = (std::sqrt(8.0));
            const double tolerance = (1e-12);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
    }
}

TEST(Main, BigintMulViaCas) {
    auto num = SymbolicExpr::number(BigInt("123456789123456789"));
    auto expr = SymbolicExpr::multiply(num, num);
    std::vector<Value> args = {Value(expr)};
    Value res = cas_simplify(args);
    {
        const std::string actual_text = (res.to_string());
        for (const auto &token : std::vector<std::string>{"152415787"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "BigInt mul result start" << ": missing " << token << " in " << actual_text;
        }
    }
}

TEST(Main, Differentiation) {
    auto x = SymbolicExpr::variable("x");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    std::vector<Value> args = {Value(x2), Value("x")};
    Value res = cas_differentiate(args);
    EXPECT_EQ((res.to_string()), ("2*x")) << "Diff x^2";

    auto sin_x = SymbolicExpr::sin(x);
    args = {Value(sin_x), Value("x")};
    res = cas_differentiate(args);
    EXPECT_EQ((res.to_string()), ("cos(x)")) << "Diff sin(x)";
}

TEST(Main, Solve2x60) {
    auto x = SymbolicExpr::variable("x");

    auto eq = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(-6));
    std::vector<Value> args = {Value(eq), Value("x")};
    Value res = cas_solve(args);
    {
        const std::string actual_text = (res.to_string());
        for (const auto &token : std::vector<std::string>{"3"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "Solution contains 3" << ": missing " << token << " in " << actual_text;
        }
    }
}

TEST(Main, MatrixAddition) {
    auto m1 = SymbolicExpr::matrix({{SymbolicExpr::number(1), SymbolicExpr::number(2)},
                                    {SymbolicExpr::number(3), SymbolicExpr::number(4)}});
    auto m2 = SymbolicExpr::matrix({{SymbolicExpr::number(3), SymbolicExpr::number(4)},
                                    {SymbolicExpr::number(5), SymbolicExpr::number(6)}});
    auto mat_add = SymbolicExpr::add(m1, m2);
    std::vector<Value> args = {Value(mat_add)};
    Value res = cas_simplify(args);
    {
        const std::string actual_text = (res.to_string());
        for (const auto &token : std::vector<std::string>{"4", "6", "8", "10"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "Matrix elements sum" << ": missing " << token << " in " << actual_text;
        }
    }
}

TEST(Main, Logarithm) {
    auto x = SymbolicExpr::variable("x");

    auto x_squared = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto log_expr = SymbolicExpr::log(x_squared, x);
    auto simplified = log_expr->simplify();
    auto legal = evaluate_numeric(*simplified, {{"x", 2}});
    ASSERT_TRUE((legal.has_value())) << "positive nonunit log base is legal";
    if (legal) {
        EXPECT_NEAR(legal.value().value, 2.0, 1e-12) << "log base preserves its value";
    }
    for (double value : {0.0, 1.0, -2.0}) {
        auto invalid = evaluate_numeric(*simplified, {{"x", value}});
        EXPECT_TRUE((!invalid && invalid.error().code == CasErrc::DomainError)) << "logarithm simplification retains the base domain";
    }
}
