#include "test_common.hpp"
#include <iostream>
#include <limits>
#include <set>
#include <map>
#include "value.hpp"
#include "symbolic.hpp"

using namespace LMCAS;

int main() {
    Value v1(10);
    Value v2(20);
    Value v3(10);
    Value v4 = nullptr;
    Value v5 = nullptr;
    Value v6(3.14);

    EXPECT_TRUE(v1 == v3, "v1 should equal v3");
    EXPECT_TRUE(!(v1 == v2), "v1 should not equal v2");
    EXPECT_TRUE(v1 < v2, "v1 should be less than v2");
    EXPECT_TRUE(!(v2 < v1), "v2 should not be less than v1");
    EXPECT_TRUE(v4 == v5, "null should equal null");
    EXPECT_TRUE(!(v4 < v5), "null should not be less than null");

    EXPECT_TRUE(v1 < v6 || v6 < v1, "Different types should be comparable");

    auto x = SymbolicExpr::variable("x");
    auto expr1 = SymbolicExpr::add(x, SymbolicExpr::number(1));

    auto y = SymbolicExpr::variable("x");
    auto expr2 = SymbolicExpr::add(y, SymbolicExpr::number(1));

    Value sym1(expr1);
    Value sym2(expr2);

    EXPECT_TRUE(sym1 == sym2, "SymbolicExpr AST structure should compare equal");
    EXPECT_TRUE(!(sym1 < sym2) && !(sym2 < sym1), "SymbolicExpr AST structure < should be false for equal trees");

    auto numeric_expr = SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::number(1));
    Value symbolic_number(numeric_expr);
    EXPECT_NEAR(symbolic_number.as_number(), 2.0, 1e-12,
                "Value::as_number evaluates finite symbolic numeric expressions");

    bool symbolic_var_threw = false;
    try {
        (void)Value(SymbolicExpr::variable("z")).as_number();
    } catch (const std::runtime_error&) {
        symbolic_var_threw = true;
    }
    EXPECT_TRUE(symbolic_var_threw,
                "Value::as_number rejects unbound symbolic expressions");

    bool null_symbolic_threw = false;
    try {
        (void)Value(std::shared_ptr<SymbolicExpr>{});
    } catch (const std::invalid_argument&) {
        null_symbolic_threw = true;
    }
    EXPECT_TRUE(null_symbolic_threw,
                "Value rejects null symbolic expressions at construction");
    EXPECT_TRUE(sym1.kind() == Value::Type::Symbolic,
                "Value exposes its type without mutable discriminant access");
    EXPECT_TRUE(
        std::holds_alternative<std::shared_ptr<SymbolicExpr>>(sym1.storage()),
        "Value exposes read-only storage for inspection");

    auto null_number = v4.as_number_checked();
    EXPECT_TRUE(!null_number &&
                    null_number.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked numeric conversion reports incompatible values");
    auto null_symbolic = v4.as_symbolic_checked();
    EXPECT_TRUE(!null_symbolic &&
                    null_symbolic.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked symbolic conversion reports incompatible values");

    Value short_vector(std::vector<Value>{Value(1)});
    Value long_vector(std::vector<Value>{Value(1), Value(2)});
    auto mismatched_sum = short_vector.vector_add_checked(long_vector);
    EXPECT_TRUE(!mismatched_sum &&
                    mismatched_sum.error().code ==
                        LMCAS::CasErrc::DimensionMismatch,
                "checked vector addition reports dimension mismatch");
    Value text_vector(std::vector<Value>{Value("not numeric")});
    auto invalid_dot = short_vector.dot_product_checked(text_vector);
    EXPECT_TRUE(!invalid_dot &&
                    invalid_dot.error().code ==
                        LMCAS::CasErrc::InvalidArgument,
                "checked dot product reports nonnumeric elements");

    bool nan_threw = false;
    try {
        (void)Value(std::numeric_limits<lmmc_real_t>::quiet_NaN());
    } catch (const std::invalid_argument&) {
        nan_threw = true;
    }
    EXPECT_TRUE(nan_threw, "Value rejects NaN at the construction boundary");

    Value infinity(std::numeric_limits<lmmc_real_t>::infinity());
    EXPECT_TRUE(infinity.is_infinity(), "Value preserves infinity as a dedicated type");
    EXPECT_TRUE(infinity.as_symbolic_compatible(),
                "positive infinity reports symbolic compatibility");
    EXPECT_TRUE(infinity.to_string() == "inf", "positive infinity has stable formatting");
    EXPECT_TRUE(infinity.as_symbolic()->to_string() == "inf",
                "positive infinity converts to a symbolic infinity node");

    Value negative_infinity(-std::numeric_limits<lmmc_real_t>::infinity());
    EXPECT_TRUE(negative_infinity.is_infinity(),
                "Value preserves negative infinity as a dedicated type");
    EXPECT_TRUE(negative_infinity.as_symbolic_compatible(),
                "negative infinity reports symbolic compatibility");
    EXPECT_TRUE(negative_infinity.to_string() == "-inf",
                "negative infinity has stable formatting");
    auto negative_infinity_symbolic = negative_infinity.as_symbolic();
    EXPECT_TRUE(negative_infinity_symbolic->compare(SymbolicExpr::infinity(-1)) == 0,
                "negative infinity converts to the canonical symbolic representation");

    Value approximate_half(static_cast<lmmc_real_t>(0.5));
    auto approximate_symbolic = approximate_half.as_symbolic();
    EXPECT_TRUE(approximate_symbolic->compare(SymbolicExpr::number(0.5)) == 0,
                "finite floating-point Values remain approximate symbolic numbers");
    EXPECT_TRUE(approximate_symbolic->compare(
                    SymbolicExpr::number(Rational(1, 2))) != 0,
                "finite floating-point Values are not exactified as rationals");

    bool nested_ragged_threw = false;
    try {
        (void)Value(std::vector<std::vector<Value>>{
            {Value(1), Value(2)},
            {Value(3)}
        });
    } catch (const std::invalid_argument&) {
        nested_ragged_threw = true;
    }
    EXPECT_TRUE(nested_ragged_threw,
                "direct matrix construction rejects rows with different widths");

    bool array_ragged_threw = false;
    try {
        Value row1(std::vector<Value>{Value(1), Value(2)});
        Value row2(std::vector<Value>{Value(3)});
        (void)Value(std::vector<Value>{row1, row2});
    } catch (const std::invalid_argument&) {
        array_ragged_threw = true;
    }
    EXPECT_TRUE(array_ragged_threw,
                "nested-array matrix construction rejects rows with different widths");

    Value matrix_left(std::vector<std::vector<Value>>{
        {Value(1), Value(2), Value(3)},
        {Value(4), Value(5), Value(6)}
    });
    Value matrix_right(std::vector<std::vector<Value>>{
        {Value(1), Value(0)},
        {Value(0), Value(1)},
        {Value(1), Value(1)}
    });
    EXPECT_TRUE(matrix_left.matrix_multiply(matrix_right).to_string() == "[[4, 5], [10, 11]]",
                "rectangular matrix multiplication remains supported");
    auto checked_matrix_product =
        matrix_left.matrix_multiply_checked(matrix_right);
    EXPECT_TRUE(
        checked_matrix_product &&
            checked_matrix_product.value().to_string() ==
                "[[4, 5], [10, 11]]",
        "checked rectangular matrix multiplication returns its product");
    auto invalid_matrix_product =
        matrix_right.matrix_multiply_checked(matrix_right);
    EXPECT_TRUE(
        !invalid_matrix_product &&
            invalid_matrix_product.error().code ==
                LMCAS::CasErrc::DimensionMismatch,
        "checked matrix multiplication reports incompatible dimensions");

    std::set<Value> s;
    s.insert(v1);
    s.insert(v2);
    s.insert(v3);
    s.insert(v4);
    s.insert(v5);
    s.insert(v6);
    s.insert(sym1);
    s.insert(sym2);

    EXPECT_TRUE(s.size() == 5, "Set should have 5 unique elements (10, 20, null, 3.14, x+1)");

    return TEST_REPORT();
}
