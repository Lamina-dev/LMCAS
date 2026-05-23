#include "test_common.hpp"

#include "value.hpp"

Value cas_simplify(const std::vector<Value>& args) {
    if (args.empty()) return Value(SymbolicExpr::number(0));
    auto expr = args[0].as_symbolic();
    if (!expr) return args[0];
    return Value(expr->simplify());
}

Value cas_differentiate(const std::vector<Value>& args) {
    if (args.size() < 2) return Value(SymbolicExpr::number(0));
    auto expr = args[0].as_symbolic();
    std::string var = args[1].to_string();
    return Value(expr->differentiate(var));
}

Value cas_solve(const std::vector<Value>& args) {
    if (args.size() < 2) return Value();
    auto expr = args[0].as_symbolic();
    std::string var = args[1].to_string();

    auto solutions = SymbolicExpr::solve(expr, var);

    std::vector<Value> val_sols;
    for (const auto& s : solutions) {
        val_sols.push_back(Value(s));
    }
    return Value(val_sols);
}

int main() {
    auto x = SymbolicExpr::variable("x");

    TEST_CASE("Simplify 1 + 1");
    {
        auto one = SymbolicExpr::number(1);
        auto expr1 = SymbolicExpr::add(one, one);
        std::vector<Value> args = { Value(expr1) };
        Value res = cas_simplify(args);
        EXPECT_EQ_STR(res.to_string(), "2", "1+1 = 2");
    }

    TEST_CASE("Sqrt(4)");
    {
        auto four = SymbolicExpr::number(4);
        auto expr = SymbolicExpr::sqrt(four);
        std::vector<Value> args = { Value(expr) };
        Value res = cas_simplify(args);
        EXPECT_EQ_STR(res.to_string(), "2", "sqrt(4) = 2");
    }

    TEST_CASE("Sqrt(8)");
    {
        auto eight = SymbolicExpr::number(8);
        auto expr = SymbolicExpr::sqrt(eight);
        std::vector<Value> args = { Value(expr) };
        Value res = cas_simplify(args);

        std::string s = res.to_string();
        std::cout << "[DEBUG] sqrt(8) -> " << s << std::endl;
        bool ok = (s.find("2") != std::string::npos) && (s.find("8") == std::string::npos);
        EXPECT_TRUE(ok, "sqrt(8) = 2*sqrt(2) (checks for 2 and no 8)");
    }

    TEST_CASE("BigInt Mul via CAS");
    {
        auto num = SymbolicExpr::number(BigInt("123456789123456789"));
        auto expr = SymbolicExpr::multiply(num, num);
        std::vector<Value> args = { Value(expr) };
        Value res = cas_simplify(args);
        EXPECT_CONTAINS(res.to_string(), {"152415787"}, "BigInt mul result start");
    }

    TEST_CASE("Differentiation");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        std::vector<Value> args = { Value(x2), Value("x") };
        Value res = cas_differentiate(args);
        EXPECT_EQ_STR(res.to_string(), "2*x", "Diff x^2");

        auto sin_x = SymbolicExpr::sin(x);
        args = { Value(sin_x), Value("x") };
        res = cas_differentiate(args);
        EXPECT_EQ_STR(res.to_string(), "cos(x)", "Diff sin(x)");
    }

    TEST_CASE("Solve 2x - 6 = 0");
    {
        auto eq = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(-6)
        );
        std::vector<Value> args = { Value(eq), Value("x") };
        Value res = cas_solve(args);
        EXPECT_CONTAINS(res.to_string(), {"3"}, "Solution contains 3");
    }

    TEST_CASE("Matrix Addition");
    {
        auto m1 = SymbolicExpr::matrix({
            {SymbolicExpr::number(1), SymbolicExpr::number(2)},
            {SymbolicExpr::number(3), SymbolicExpr::number(4)}
        });
        auto m2 = SymbolicExpr::matrix({
            {SymbolicExpr::number(3), SymbolicExpr::number(4)},
            {SymbolicExpr::number(5), SymbolicExpr::number(6)}
        });
        auto mat_add = SymbolicExpr::add(m1, m2);
        std::vector<Value> args = { Value(mat_add) };
        Value res = cas_simplify(args);
        EXPECT_CONTAINS(res.to_string(), {"4", "6", "8", "10"}, "Matrix elements sum");
    }

    TEST_CASE("Logarithm");
    {
        auto x_squared = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto log_expr = SymbolicExpr::log(x_squared, x);
        std::vector<Value> args = { Value(log_expr) };
        Value res = cas_simplify(args);
        EXPECT_EQ_STR(res.to_string(), "2", "log(x^2, x) = 2");
    }

    return TEST_REPORT();
}
