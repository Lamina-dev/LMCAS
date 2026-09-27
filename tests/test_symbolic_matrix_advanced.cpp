#include "test_common.hpp"
#include "symbolic_matrix.hpp"
#include "numeric_evaluation.hpp"
#include "residual_verification.hpp"

#include <cmath>
#include <limits>

#include <utility>
#include <vector>
using namespace LMCAS;

TEST(SymbolicMatrixAdvanced, MatrixMultiply) {
    // A * I = A (identity multiplication)
    {
        auto a = SymbolicExpr::variable("a");
        auto b = SymbolicExpr::variable("b");
        auto c = SymbolicExpr::variable("c");
        auto d = SymbolicExpr::variable("d");

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A_data = {
            {a, b},
            {c, d}};
        auto A = SymbolicExpr::matrix(A_data);

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> I_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}};
        auto I = SymbolicExpr::matrix(I_data);

        auto result = LMCAS::matrix_multiply_checked(A, I);
        ASSERT_TRUE(result.has_value());
        EXPECT_TRUE(test_same_expression(result.value(), A));
    }

    // General 2x2 multiplication with numeric matrices
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(2)},
            {SymbolicExpr::number(3), SymbolicExpr::number(4)}};
        auto A = SymbolicExpr::matrix(A_data);

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> B_data = {
            {SymbolicExpr::number(5), SymbolicExpr::number(6)},
            {SymbolicExpr::number(7), SymbolicExpr::number(8)}};
        auto B = SymbolicExpr::matrix(B_data);

        auto result = LMCAS::matrix_multiply_checked(A, B);
        ASSERT_TRUE(result.has_value());
        auto expected = SymbolicExpr::matrix({
            {SymbolicExpr::number(19), SymbolicExpr::number(22)},
            {SymbolicExpr::number(43), SymbolicExpr::number(50)}});
        EXPECT_TRUE(test_same_expression(result.value(), expected));
    }
}

TEST(SymbolicMatrixAdvanced, MatrixMultiplyCallerBudgets) {
    const auto one = SymbolicExpr::number(1);
    const auto two = SymbolicExpr::number(2);
    const auto matrix = SymbolicExpr::matrix({{one, one}, {one, one}});
    ComputationContext context;
    auto product = matrix_multiply_checked(matrix, matrix, context);
    ASSERT_TRUE(product);
    const auto expected = SymbolicExpr::matrix({{two, two}, {two, two}});
    EXPECT_EQ(detail::node(product.value())->compare(*detail::node(expected)), 0);

    for (auto member : {&ResourceLimits::max_steps, &ResourceLimits::max_ast_nodes,
                        &ResourceLimits::max_recursion_depth, &ResourceLimits::max_expansion_terms}) {
        ResourceLimits limits;
        limits.*member = member == &ResourceLimits::max_steps ? 2 : 0;
        ComputationContext limited(limits);
        auto exhausted = matrix_multiply_checked(matrix, matrix, limited);
        ASSERT_FALSE(exhausted);
        EXPECT_EQ(exhausted.error().code, CasErrc::ResourceLimit);
    }

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    auto interrupted = matrix_multiply_checked(matrix, matrix, cancelled);
    ASSERT_FALSE(interrupted);
    EXPECT_EQ(interrupted.error().code, CasErrc::Cancelled);

    ResourceLimits limits;
    limits.max_ast_nodes = 0;
    ComputationContext invalid_shape(limits);
    auto mismatch = matrix_multiply_checked(
        matrix, SymbolicExpr::matrix({{one}, {one}, {one}}), invalid_shape);
    ASSERT_FALSE(mismatch);
    EXPECT_EQ(mismatch.error().code, CasErrc::InvalidArgument);
}

TEST(SymbolicMatrixAdvanced, MatrixMultiplyExpansionDimensions) {
    const auto one = SymbolicExpr::number(1);
    ResourceLimits limits;
    limits.max_expansion_terms = 3;
    ComputationContext output_limited(limits);
    auto square = SymbolicExpr::matrix({{one, one}, {one, one}});
    auto output = matrix_multiply_checked(square, square, output_limited);
    ASSERT_FALSE(output);
    EXPECT_EQ(output.error().code, CasErrc::ResourceLimit);

    ComputationContext inner_limited(limits);
    auto row = SymbolicExpr::matrix({{one, one, one, one}});
    auto column = SymbolicExpr::matrix({{one}, {one}, {one}, {one}});
    auto inner = matrix_multiply_checked(row, column, inner_limited);
    ASSERT_FALSE(inner);
    EXPECT_EQ(inner.error().code, CasErrc::ResourceLimit);
}

static void expect_scaled_matrix_entries(const MatrixNode &matrix, int coefficient,
                                         ComputationContext &context) {
    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t col = 0; col < 3; ++col) {
            const int expected = row == 0 && col == 1 ? 2 * coefficient : 0;
            auto proof = check_zero_residual(SymbolicExpr::add(
                                                 detail::make_expression_ptr(matrix.get(row, col)),
                                                 SymbolicExpr::number(-expected)),
                                             context);
            EXPECT_TRUE((proof && std::holds_alternative<ProvedZeroResidual>(proof.value()))) << "scalar matrix multiplication has its exact entries";
        }
    }
}

static void check_scaled_rectangular_products(const ExprPtr &a, ComputationContext &context) {
    for (int coefficient : {0, -2}) {
        for (bool scalar_first : {false, true}) {
            auto scalar = SymbolicExpr::number(coefficient);
            auto scaled = (scalar_first ? SymbolicExpr::multiply(scalar, a)
                                        : SymbolicExpr::multiply(a, scalar))
                              ->simplify();
            auto result = std::dynamic_pointer_cast<const MatrixNode>(detail::node(scaled));
            EXPECT_TRUE((result && result->rows() == 2 && result->cols() == 3)) << "scalar matrix multiplication preserves rectangular shape";
            if (!result) {
                continue;
            }
            expect_scaled_matrix_entries(*result, coefficient, context);
        }
    }
}

TEST(SymbolicMatrixAdvanced, MatrixProductShapes) {
    ComputationContext context;
    auto zero = SymbolicExpr::number(0);
    auto a = SymbolicExpr::matrix({{zero, SymbolicExpr::number(2), zero}, {zero, zero, zero}});
    auto b = SymbolicExpr::matrix({{zero}, {SymbolicExpr::number(3)}, {zero}});
    auto product = matrix_multiply_checked(a, b);
    ASSERT_TRUE((product.has_value())) << "sparse rectangular matrix product succeeds";
    if (!product) {
        return;
    }
    auto matrix = std::dynamic_pointer_cast<const MatrixNode>(detail::node(product.value()));
    EXPECT_TRUE((matrix && matrix->rows() == 2 && matrix->cols() == 1)) << "checked rectangular product is an explicit matrix";
    if (!matrix) {
        return;
    }
    for (std::size_t row = 0; row < 2; ++row) {
        auto entry = detail::make_expression_ptr(matrix->get(row, 0));
        auto proof = check_zero_residual(SymbolicExpr::add(entry,
                                                           SymbolicExpr::number(row == 0 ? -6 : 0)),
                                         context);
        EXPECT_TRUE((proof && std::holds_alternative<ProvedZeroResidual>(proof.value()))) << "sparse rectangular product has its exact entries";
    }
    check_scaled_rectangular_products(a, context);
}

TEST(SymbolicMatrixAdvanced, MatrixDeterminant) {
    // 2x2 determinant: det([[a,b],[c,d]]) = ad - bc
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(3), SymbolicExpr::number(8)},
            {SymbolicExpr::number(4), SymbolicExpr::number(6)}};
        auto m = SymbolicExpr::matrix(m_data);
        auto det = LMCAS::matrix_determinant_checked(m).value();
        // det = 3*6 - 8*4 = 18 - 32 = -14
        auto simplified = det ? det->simplify() : nullptr;
        EXPECT_TRUE(test_same_expression(simplified, SymbolicExpr::number(-14))) << "Det 2x2: 3*6 - 8*4 = -14";
    }

    // 3x3 determinant with known value
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(6), SymbolicExpr::number(1), SymbolicExpr::number(1)},
            {SymbolicExpr::number(4), SymbolicExpr::number(-2), SymbolicExpr::number(5)},
            {SymbolicExpr::number(2), SymbolicExpr::number(8), SymbolicExpr::number(7)}};
        auto m = SymbolicExpr::matrix(m_data);
        auto det = LMCAS::matrix_determinant_checked(m).value();
        // det = 6*(-2*7 - 5*8) - 1*(4*7 - 5*2) + 1*(4*8 - (-2)*2)
        //     = 6*(-14-40) - 1*(28-10) + 1*(32+4)
        //     = 6*(-54) - 18 + 36
        //     = -324 - 18 + 36 = -306
        auto simplified = det ? det->simplify() : nullptr;
        EXPECT_TRUE(test_same_expression(simplified, SymbolicExpr::number(-306))) << "Det 3x3: known determinant = -306";
    }
}

TEST(SymbolicMatrixAdvanced, MatrixInverse) {
    // 2x2 invertible matrix: verify A * A_inv structure
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(4), SymbolicExpr::number(7)},
            {SymbolicExpr::number(2), SymbolicExpr::number(6)}};
        auto A = SymbolicExpr::matrix(m_data);
        auto A_inv = LMCAS::matrix_inverse_checked(A).value();
        EXPECT_TRUE((A_inv != nullptr)) << "matrix_inverse returns non-null for invertible 2x2";

        // Verify A * A_inv = I
        auto product = LMCAS::matrix_multiply_checked(A, A_inv).value();
        auto simplified = product ? product->simplify() : nullptr;

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}};
        auto identity = SymbolicExpr::matrix(id_data);
        EXPECT_TRUE(test_same_expression(simplified, identity)) << "A * A_inv = Identity for 2x2 invertible matrix";
    }
}

static void expect_numeric_matrix(const ExpressionResult &result,
                                  const double (&expected)[2][2], double tolerance,
                                  const std::string &message) {
    if (!result) {
        std::cerr << message << ": code=" << static_cast<int>(result.error().code)
                  << " operation=" << result.error().operation
                  << " message=" << result.error().message << '\n';
    }
    ASSERT_TRUE((result.has_value())) << message + " construction succeeds";
    if (!result) {
        return;
    }
    auto matrix = std::dynamic_pointer_cast<const MatrixNode>(
        detail::node(result.value()->simplify()));
    EXPECT_TRUE((matrix && matrix->rows() == 2 && matrix->cols() == 2)) << message + " has the required dimensions";
    if (!matrix || matrix->rows() != 2 || matrix->cols() != 2) {
        return;
    }
    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t col = 0; col < 2; ++col) {
            auto entry = evaluate_numeric(*detail::make_expression_ptr(matrix->get(row, col)));
            if (!entry) {
                std::cerr << message << ": code=" << static_cast<int>(entry.error().code)
                          << " operation=" << entry.error().operation
                          << " message=" << entry.error().message << '\n';
            }
            ASSERT_TRUE((entry.has_value())) << message + " entry evaluates";
            if (entry) {
                {
                    const double actual_value = (entry.value().value);
                    const double expected_value = (expected[row][col]);
                    SCOPED_TRACE(message + " entry has its mathematical value");
                    EXPECT_TRUE(std::isfinite(actual_value));
                    EXPECT_GE(tolerance, 0.0);
                    EXPECT_NEAR(actual_value, expected_value, tolerance);
                }
            }
        }
    }
}

TEST(SymbolicMatrixAdvanced, MatrixRotation) {
    const double identity[2][2] = {{1, 0}, {0, 1}};
    expect_numeric_matrix(matrix_rotation_checked(0.0), identity, 0.0,
                          "zero-angle rotation is the identity");
    const double quarter_turn[2][2] = {{0, -1}, {1, 0}};
    expect_numeric_matrix(matrix_rotation_checked(std::acos(0.0)), quarter_turn, 1e-12,
                          "quarter-turn rotation");
}

static void expect_reflection_coordinates(
    const std::shared_ptr<const MatrixNode> &matrix,
    const double (&expected)[2], size_t angle_index) {
    for (size_t row = 0; row < 2; ++row) {
        auto value = test_numeric_eval(
            LMCAS::detail::make_expression_ptr(matrix->get(row, 0)));
        EXPECT_TRUE((value && std::isfinite(*value) &&
                     std::abs(*value - expected[row]) < 1e-12))
            << (angle_index ==
                        0
                    ? "x-axis reflection reverses the normal component"
                    : "diagonal reflection swaps the coordinates");
    }
}

TEST(SymbolicMatrixAdvanced, MatrixReflection) {
    auto vector = SymbolicExpr::matrix({{SymbolicExpr::number(2)}, {SymbolicExpr::number(3)}});
    const double angles[] = {0.0, std::acos(-1.0) / 4.0};
    const double expected[][2] = {{2.0, -3.0}, {3.0, 2.0}};
    for (size_t i = 0; i < 2; ++i) {
        auto reflection = matrix_reflection_checked(angles[i]);
        ASSERT_TRUE((reflection.has_value())) << "reflection construction succeeds";
        if (!reflection) {
            continue;
        }
        auto reflected = matrix_multiply_checked(reflection.value(), vector);
        ASSERT_TRUE((reflected.has_value())) << "reflection applies to a vector";
        if (!reflected) {
            continue;
        }
        auto matrix = std::dynamic_pointer_cast<const MatrixNode>(
            LMCAS::detail::node(reflected.value()->simplify()));
        EXPECT_TRUE((matrix && matrix->rows() == 2 && matrix->cols() == 1)) << "reflection preserves vector dimensions";
        if (!matrix || matrix->rows() != 2 || matrix->cols() != 1) {
            continue;
        }
        expect_reflection_coordinates(matrix, expected[i], i);
    }
}

TEST(SymbolicMatrixAdvanced, MatrixReflectionLargeFiniteAngle) {
    const double angle = 1e308;
    auto reflection = matrix_reflection_checked(angle);
    ASSERT_TRUE(reflection.has_value());
    auto matrix = std::dynamic_pointer_cast<const MatrixNode>(
        detail::node(reflection.value()->simplify()));
    ASSERT_TRUE(matrix);
    ASSERT_EQ(matrix->rows(), 2u);
    ASSERT_EQ(matrix->cols(), 2u);
    double entries[2][2];
    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t col = 0; col < 2; ++col) {
            auto entry = evaluate_numeric(*detail::make_expression_ptr(matrix->get(row, col)));
            ASSERT_TRUE(entry.has_value());
            entries[row][col] = entry.value().value;
            EXPECT_TRUE(std::isfinite(entries[row][col]));
        }
    }
    const double axis[] = {std::cos(angle), std::sin(angle)};
    const double normal[] = {-axis[1], axis[0]};
    for (std::size_t row = 0; row < 2; ++row) {
        EXPECT_NEAR(entries[row][0] * axis[0] + entries[row][1] * axis[1],
                    axis[row], 1e-12);
        EXPECT_NEAR(entries[row][0] * normal[0] + entries[row][1] * normal[1],
                    -normal[row], 1e-12);
    }
    const double identity[2][2] = {{1, 0}, {0, 1}};
    expect_numeric_matrix(matrix_multiply_checked(reflection.value(), reflection.value()),
                          identity, 1e-12, "reflection squared is the identity");
}

TEST(SymbolicMatrixAdvanced, MatrixReflectionCheckedErrors) {
    for (double angle : {std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()}) {
        auto invalid = matrix_reflection_checked(angle);
        ASSERT_FALSE(invalid.has_value());
        EXPECT_EQ(invalid.error().code, CasErrc::InvalidArgument);
        EXPECT_EQ(invalid.error().operation, "matrix_reflection");
    }

    CancellationToken cancellation;
    ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = matrix_reflection_checked(1e308, 2, cancelled_context);
    ASSERT_FALSE(cancelled.has_value());
    EXPECT_EQ(cancelled.error().code, CasErrc::Cancelled);
    EXPECT_EQ(cancelled.error().operation, "matrix_reflection");

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext limited_context(limits);
    auto limited = matrix_reflection_checked(1e308, 2, limited_context);
    ASSERT_FALSE(limited.has_value());
    EXPECT_EQ(limited.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(limited.error().operation, "matrix_reflection");

    auto unsupported = matrix_reflection_checked(1e308, 3);
    ASSERT_FALSE(unsupported.has_value());
    EXPECT_EQ(unsupported.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(unsupported.error().operation, "matrix_reflection");
}

TEST(SymbolicMatrixAdvanced, MatrixScaling) {
    const double scaled[2][2] = {{2, 0}, {0, 3}};
    expect_numeric_matrix(matrix_scaling_checked(2.0, 3.0), scaled, 0.0,
                          "coordinate scaling");
    const double identity[2][2] = {{1, 0}, {0, 1}};
    expect_numeric_matrix(matrix_scaling_checked(1.0, 1.0), identity, 0.0,
                          "unit scaling is the identity");
}

TEST(SymbolicMatrixAdvanced, MatrixTransformsRejectNonfiniteArguments) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();

    for (double theta : {nan, infinity, -infinity}) {
        auto result = matrix_rotation_checked(theta, 2);
        ASSERT_FALSE(result.has_value());
        EXPECT_EQ(result.error().code, CasErrc::InvalidArgument);
        EXPECT_EQ(result.error().operation, "matrix_rotation");
    }
    for (const auto &scale : std::vector<std::pair<double, double>>{
             {nan, 1.0}, {infinity, 1.0}, {-infinity, 1.0},
             {1.0, nan}, {1.0, infinity}, {1.0, -infinity}}) {
        auto result = matrix_scaling_checked(
            scale.first, scale.second, 2);
        ASSERT_FALSE(result.has_value());
        EXPECT_EQ(result.error().code, CasErrc::InvalidArgument);
        EXPECT_EQ(result.error().operation, "matrix_scaling");
    }

    auto bad_rotation_dimension = matrix_rotation_checked(nan, 3);
    ASSERT_FALSE(bad_rotation_dimension.has_value());
    EXPECT_EQ(bad_rotation_dimension.error().code,
              CasErrc::DimensionMismatch);
    auto bad_scaling_dimension = matrix_scaling_checked(nan, 1.0, 3);
    ASSERT_FALSE(bad_scaling_dimension.has_value());
    EXPECT_EQ(bad_scaling_dimension.error().code,
              CasErrc::DimensionMismatch);

    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto cancelled_rotation =
        matrix_rotation_checked(nan, 2, cancelled);
    ASSERT_FALSE(cancelled_rotation.has_value());
    EXPECT_EQ(cancelled_rotation.error().code, CasErrc::Cancelled);
}
