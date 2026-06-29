/**
 * @file test_z_transform.cpp
 * @brief Z 变换单元测试。
 */
#include "test_common.hpp"
#include "transform_engine.hpp"
#include "symbolic_ast.hpp"

using namespace lamina;

int main() {
    TEST_CASE("Z transform: unit step (constant 1)");
    {
        // Z{1} = z/(z-1)
        auto one = SymbolicExpr::number(1);
        auto result = z_transform(one, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{1} should not be null");
        auto z = SymbolicExpr::variable("z");
        auto expected = SymbolicExpr::divide(z,
            SymbolicExpr::add(z, SymbolicExpr::number(-1)));
        EXPECT_EQ_EXPR(result->simplify(), expected->simplify(),
            "Z{1} = z/(z-1)");
    }

    TEST_CASE("Z transform: exponential sequence a^n");
    {
        // Z{(1/2)^n} = z/(z - 1/2)
        auto half = SymbolicExpr::number(0.5);
        auto n = SymbolicExpr::variable("n");
        auto f_n = SymbolicExpr::power(half, n);
        auto result = z_transform(f_n, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{(1/2)^n} should not be null");
        auto z = SymbolicExpr::variable("z");
        auto expected = SymbolicExpr::divide(z,
            SymbolicExpr::add(z, SymbolicExpr::number(-0.5)));
        EXPECT_EQ_EXPR(result->simplify(), expected->simplify(),
            "Z{(1/2)^n} = z/(z-0.5)");
    }

    TEST_CASE("Z transform: polynomial sequence n");
    {
        // Z{n} = z/(z-1)^2
        auto n = SymbolicExpr::variable("n");
        auto result = z_transform(n, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{n} should not be null");
        auto z = SymbolicExpr::variable("z");
        auto db = SymbolicExpr::add(z, SymbolicExpr::number(-1));
        auto expected = SymbolicExpr::divide(z,
            SymbolicExpr::power(db, SymbolicExpr::number(2)));
        EXPECT_EQ_EXPR(result->simplify(), expected->simplify(),
            "Z{n} = z/(z-1)^2");
    }

    TEST_CASE("Z transform: sinusoidal sequence sin(w*n)");
    {
        // Z{sin(w*n)} = z*sin(w) / (z^2 - 2z*cos(w) + 1)
        auto w = SymbolicExpr::variable("w");
        auto n = SymbolicExpr::variable("n");
        auto f_n = SymbolicExpr::sin(SymbolicExpr::multiply(w, n));
        auto result = z_transform(f_n, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{sin(w*n)} should not be null");
        // Verify it's not an unevaluated node
        EXPECT_FALSE(
            std::dynamic_pointer_cast<TransformNode>(result->root) != nullptr,
            "Z{sin(w*n)} should be evaluated (not unevaluated)");
        // Check the result contains expected components
        auto str = result->simplify()->to_string();
        EXPECT_CONTAINS(str, {"z", "sin"}, "Z{sin(w*n)} contains z and sin");
    }

    TEST_CASE("Z transform: cosine sequence cos(w*n)");
    {
        // Z{cos(w*n)} = z*(z - cos(w)) / (z^2 - 2z*cos(w) + 1)
        auto w = SymbolicExpr::variable("w");
        auto n = SymbolicExpr::variable("n");
        auto f_n = SymbolicExpr::cos(SymbolicExpr::multiply(w, n));
        auto result = z_transform(f_n, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{cos(w*n)} should not be null");
        EXPECT_FALSE(
            std::dynamic_pointer_cast<TransformNode>(result->root) != nullptr,
            "Z{cos(w*n)} should be evaluated");
        auto str = result->simplify()->to_string();
        EXPECT_CONTAINS(str, {"z", "cos"}, "Z{cos(w*n)} contains z and cos");
    }

    TEST_CASE("Z transform: constant sequence c");
    {
        // Z{5} = 5*z/(z-1)
        auto five = SymbolicExpr::number(5);
        auto result = z_transform(five, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{5} should not be null");
        auto z = SymbolicExpr::variable("z");
        auto expected = SymbolicExpr::multiply(five,
            SymbolicExpr::divide(z,
                SymbolicExpr::add(z, SymbolicExpr::number(-1))));
        EXPECT_EQ_EXPR(result->simplify(), expected->simplify(),
            "Z{5} = 5*z/(z-1)");
    }

    TEST_CASE("Z transform: linearity (sum)");
    {
        // Z{1 + n} = Z{1} + Z{n} = z/(z-1) + z/(z-1)^2
        auto one = SymbolicExpr::number(1);
        auto n = SymbolicExpr::variable("n");
        auto f_n = SymbolicExpr::add(one, n);
        auto result = z_transform(f_n, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{1+n} should not be null");
        EXPECT_FALSE(
            std::dynamic_pointer_cast<TransformNode>(result->root) != nullptr,
            "Z{1+n} should be evaluated");
    }

    TEST_CASE("Z transform: linearity (scalar multiple)");
    {
        // Z{3*a^n} = 3*z/(z-a)
        auto three = SymbolicExpr::number(3);
        auto a = SymbolicExpr::variable("a");
        auto n = SymbolicExpr::variable("n");
        auto f_n = SymbolicExpr::multiply(three,
            SymbolicExpr::power(a, n));
        auto result = z_transform(f_n, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{3*a^n} should not be null");
        EXPECT_FALSE(
            std::dynamic_pointer_cast<TransformNode>(result->root) != nullptr,
            "Z{3*a^n} should be evaluated");
    }

    TEST_CASE("Z transform: unevaluated for unknown form");
    {
        // Z{ln(n)} should return unevaluated TransformNode
        auto n = SymbolicExpr::variable("n");
        auto f_n = SymbolicExpr::ln(n);
        auto result = z_transform(f_n, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{ln(n)} should not be null");
        auto tn = std::dynamic_pointer_cast<TransformNode>(result->root);
        EXPECT_TRUE(tn != nullptr, "Z{ln(n)} should be unevaluated TransformNode");
        if (tn) {
            EXPECT_TRUE(
                tn->transform_type == TransformNode::TransformType::ZTransform,
                "Transform type should be ZTransform");
        }
    }

    TEST_CASE("Z transform: n*a^n");
    {
        // Z{n*a^n} = a*z/(z-a)^2
        auto a = SymbolicExpr::variable("a");
        auto n = SymbolicExpr::variable("n");
        auto f_n = SymbolicExpr::multiply(n, SymbolicExpr::power(a, n));
        auto result = z_transform(f_n, "n", "z");
        EXPECT_TRUE(result != nullptr, "Z{n*a^n} should not be null");
        EXPECT_FALSE(
            std::dynamic_pointer_cast<TransformNode>(result->root) != nullptr,
            "Z{n*a^n} should be evaluated");
        auto str = result->simplify()->to_string();
        EXPECT_CONTAINS(str, {"z", "a"}, "Z{n*a^n} contains z and a");
    }

    return TEST_REPORT();
}
