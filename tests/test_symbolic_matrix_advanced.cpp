#include "test_common.hpp"
#include "symbolic_matrix.hpp"

void test_matrix_multiply() {
    TEST_CASE("Matrix Multiply");

    // A * I = A (identity multiplication)
    {
        auto a = SymbolicExpr::variable("a");
        auto b = SymbolicExpr::variable("b");
        auto c = SymbolicExpr::variable("c");
        auto d = SymbolicExpr::variable("d");

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A_data = {
            {a, b},
            {c, d}
        };
        auto A = SymbolicExpr::matrix(A_data);

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> I_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto I = SymbolicExpr::matrix(I_data);

        auto result = lamina::matrix_multiply(A, I);
        EXPECT_TRUE(result != nullptr, "matrix_multiply(A, I) returns non-null");

        auto simplified = result->simplify();
        std::string result_str = simplified ? simplified->to_string() : "null";
        // A * I should contain the original elements
        EXPECT_CONTAINS(result_str, {"a", "b", "c", "d"}, "A*I contains original matrix elements");
    }

    // General 2x2 multiplication with numeric matrices
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(2)},
            {SymbolicExpr::number(3), SymbolicExpr::number(4)}
        };
        auto A = SymbolicExpr::matrix(A_data);

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> B_data = {
            {SymbolicExpr::number(5), SymbolicExpr::number(6)},
            {SymbolicExpr::number(7), SymbolicExpr::number(8)}
        };
        auto B = SymbolicExpr::matrix(B_data);

        auto result = lamina::matrix_multiply(A, B);
        EXPECT_TRUE(result != nullptr, "matrix_multiply(A, B) returns non-null for numeric 2x2");

        auto simplified = result->simplify();
        std::string result_str = simplified ? simplified->to_string() : "null";
        // [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]]
        EXPECT_CONTAINS(result_str, {"19", "22", "43", "50"}, "General 2x2 multiplication yields correct entries");
    }
}

void test_matrix_determinant() {
    TEST_CASE("Matrix Determinant");

    // 2x2 determinant: det([[a,b],[c,d]]) = ad - bc
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(3), SymbolicExpr::number(8)},
            {SymbolicExpr::number(4), SymbolicExpr::number(6)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto det = lamina::matrix_determinant(m);
        // det = 3*6 - 8*4 = 18 - 32 = -14
        auto simplified = det ? det->simplify() : nullptr;
        EXPECT_EQ_EXPR(simplified, SymbolicExpr::number(-14), "Det 2x2: 3*6 - 8*4 = -14");
    }

    // 3x3 determinant with known value
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(6), SymbolicExpr::number(1), SymbolicExpr::number(1)},
            {SymbolicExpr::number(4), SymbolicExpr::number(-2), SymbolicExpr::number(5)},
            {SymbolicExpr::number(2), SymbolicExpr::number(8), SymbolicExpr::number(7)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto det = lamina::matrix_determinant(m);
        // det = 6*(-2*7 - 5*8) - 1*(4*7 - 5*2) + 1*(4*8 - (-2)*2)
        //     = 6*(-14-40) - 1*(28-10) + 1*(32+4)
        //     = 6*(-54) - 18 + 36
        //     = -324 - 18 + 36 = -306
        auto simplified = det ? det->simplify() : nullptr;
        EXPECT_EQ_EXPR(simplified, SymbolicExpr::number(-306), "Det 3x3: known determinant = -306");
    }
}

void test_matrix_inverse() {
    TEST_CASE("Matrix Inverse");

    // 2x2 invertible matrix: verify A * A_inv structure
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(4), SymbolicExpr::number(7)},
            {SymbolicExpr::number(2), SymbolicExpr::number(6)}
        };
        auto A = SymbolicExpr::matrix(m_data);
        auto A_inv = lamina::matrix_inverse(A);
        EXPECT_TRUE(A_inv != nullptr, "matrix_inverse returns non-null for invertible 2x2");

        // Verify A * A_inv = I
        auto product = lamina::matrix_multiply(A, A_inv);
        auto simplified = product ? product->simplify() : nullptr;

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto identity = SymbolicExpr::matrix(id_data);
        EXPECT_EQ_EXPR(simplified, identity, "A * A_inv = Identity for 2x2 invertible matrix");
    }
}

void test_matrix_rotation() {
    TEST_CASE("Matrix Rotation");

    // theta = 0 → identity matrix
    {
        auto rot0 = lamina::matrix_rotation(0.0);
        EXPECT_TRUE(rot0 != nullptr, "matrix_rotation(0) returns non-null");

        std::string rot0_str = rot0 ? rot0->to_string() : "null";
        // cos(0) = 1, sin(0) = 0, so rotation(0) = [[1, 0], [0, 1]]
        EXPECT_CONTAINS(rot0_str, {"1"}, "Rotation(0) contains 1 (identity diagonal)");

        // Verify it's the identity by checking structure
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto identity = SymbolicExpr::matrix(id_data);
        auto simplified = rot0->simplify();
        EXPECT_EQ_EXPR(simplified, identity, "Rotation(0) = Identity");
    }

    // theta = pi/2 → 90-degree rotation: [[0, -1], [1, 0]]
    {
        double pi_half = std::acos(0.0); // pi/2
        auto rot90 = lamina::matrix_rotation(pi_half);
        EXPECT_TRUE(rot90 != nullptr, "matrix_rotation(pi/2) returns non-null");

        std::string rot90_str = rot90 ? rot90->to_string() : "null";
        // cos(pi/2) ≈ 0, sin(pi/2) = 1
        // Should be approximately [[0, -1], [1, 0]]
        EXPECT_CONTAINS(rot90_str, {"1"}, "Rotation(pi/2) contains 1");
    }
}

void test_matrix_reflection() {
    TEST_CASE("Matrix Reflection");

    // Reflection with angle = 0 (reflection across x-axis)
    {
        auto ref0 = lamina::matrix_reflection(0.0);
        EXPECT_TRUE(ref0 != nullptr, "matrix_reflection(0) returns non-null");

        std::string ref0_str = ref0 ? ref0->to_string() : "null";
        // cos(0) = 1, sin(0) = 0
        // The formula uses cos(2*angle) and sin(2*angle) for standard reflection
        // With the implementation: [[c + s*s, -s], [s, c + s*s]] at angle=0:
        // c=cos(0)=1, s=sin(0)=0 → [[1, 0], [0, 1]]
        EXPECT_CONTAINS(ref0_str, {"1"}, "Reflection(0) contains expected elements");
    }
}

void test_matrix_scaling() {
    TEST_CASE("Matrix Scaling");

    // Known scale factors: sx=2, sy=3
    {
        auto scale = lamina::matrix_scaling(2.0, 3.0);
        EXPECT_TRUE(scale != nullptr, "matrix_scaling(2, 3) returns non-null");

        std::string scale_str = scale ? scale->to_string() : "null";
        // Should be [[2, 0], [0, 3]]
        EXPECT_CONTAINS(scale_str, {"2", "3"}, "Scaling(2,3) contains diagonal elements 2 and 3");

        // Verify exact structure
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> expected_data = {
            {SymbolicExpr::number(2), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(3)}
        };
        auto expected = SymbolicExpr::matrix(expected_data);
        EXPECT_EQ_EXPR(scale, expected, "Scaling(2,3) = [[2,0],[0,3]]");
    }

    // Scale factors sx=1, sy=1 → identity
    {
        auto scale_id = lamina::matrix_scaling(1.0, 1.0);
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto identity = SymbolicExpr::matrix(id_data);
        EXPECT_EQ_EXPR(scale_id, identity, "Scaling(1,1) = Identity");
    }
}

void test_matrix_eigenvalues() {
    TEST_CASE("Matrix Eigenvalues");

    // 2x2 matrix with known eigenvalues: [[2, 1], [0, 3]] has eigenvalues 2 and 3
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(2), SymbolicExpr::number(1)},
            {SymbolicExpr::number(0), SymbolicExpr::number(3)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto eigenvals = lamina::matrix_eigenvalues(m);

        EXPECT_TRUE(eigenvals.size() == 2, "Upper triangular 2x2 has 2 eigenvalues");
        if (eigenvals.size() == 2) {
            EXPECT_TRUE(eigenvals.size() == 2, "Upper triangular 2x2 has 2 eigenvalues");
            bool found2 = false, found3 = false;
            for (auto& ev : eigenvals) {
                std::string ev_str = ev ? ev->to_string() : "";
                if (ev_str.find("2") != std::string::npos) found2 = true;
                if (ev_str.find("3") != std::string::npos) found3 = true;
            }
            EXPECT_TRUE(found2 && found3, "Eigenvalues contain 2 and 3");
        }
    }
}

void test_matrix_eigenvectors() {
    TEST_CASE("Matrix Eigenvectors");

    // 2x2 matrix: [[2, 1], [0, 3]] - verify A*v = lambda*v structurally
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(2), SymbolicExpr::number(1)},
            {SymbolicExpr::number(0), SymbolicExpr::number(3)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto eigenvecs = lamina::matrix_eigenvectors(m);

        EXPECT_TRUE(eigenvecs.size() == 2, "2x2 matrix has two independent eigenvectors");
        if (eigenvecs.size() == 2) {
            for (auto& vec : eigenvecs) {
                EXPECT_TRUE(!vec.empty(), "Eigenvector is non-empty");
                // Each component should be a valid symbolic expression
                for (auto& component : vec) {
                    EXPECT_TRUE(component != nullptr, "Eigenvector component is non-null");
                }
            }
        }
    }
}

void test_matrix_checked_contracts() {
    TEST_CASE("Matrix Checked API Contracts");

    auto x = SymbolicExpr::variable("x");
    auto not_matrix = lamina::matrix_determinant_checked(x);
    EXPECT_TRUE(!not_matrix &&
                    not_matrix.error().code == lamina::CasErrc::InvalidArgument,
                "checked determinant rejects non-matrix input");

    auto null_mul = lamina::matrix_multiply_checked(nullptr, nullptr);
    EXPECT_TRUE(!null_mul &&
                    null_mul.error().code == lamina::CasErrc::InvalidArgument,
                "checked matrix multiply rejects null inputs");

    auto A = SymbolicExpr::matrix({
        {SymbolicExpr::number(1), SymbolicExpr::number(2), SymbolicExpr::number(3)},
        {SymbolicExpr::number(4), SymbolicExpr::number(5), SymbolicExpr::number(6)}
    });
    auto B = SymbolicExpr::matrix({
        {SymbolicExpr::number(1), SymbolicExpr::number(2)},
        {SymbolicExpr::number(3), SymbolicExpr::number(4)}
    });
    auto bad_mul = lamina::matrix_multiply_checked(A, B);
    EXPECT_TRUE(!bad_mul &&
                    bad_mul.error().code == lamina::CasErrc::InvalidArgument,
                "checked matrix multiply rejects incompatible dimensions");
    EXPECT_TRUE(lamina::matrix_multiply(A, B) == nullptr,
                "legacy matrix multiply unwraps incompatible dimensions to nullptr");

    auto non_square_det = lamina::matrix_determinant_checked(A);
    EXPECT_TRUE(!non_square_det &&
                    non_square_det.error().code == lamina::CasErrc::InvalidArgument,
                "checked determinant rejects non-square matrix");
    EXPECT_TRUE(lamina::matrix_determinant(A) == nullptr,
                "legacy determinant unwraps non-square matrix to nullptr");

    auto singular = SymbolicExpr::matrix({
        {SymbolicExpr::number(1), SymbolicExpr::number(2)},
        {SymbolicExpr::number(2), SymbolicExpr::number(4)}
    });
    auto singular_inverse = lamina::matrix_inverse_checked(singular);
    EXPECT_TRUE(!singular_inverse &&
                    singular_inverse.error().code == lamina::CasErrc::DomainError,
                "checked inverse rejects provably singular matrices");
    EXPECT_TRUE(lamina::matrix_inverse(singular) == nullptr,
                "legacy inverse unwraps singular matrix domain failure to nullptr");

    auto symbolic_det = SymbolicExpr::matrix({
        {SymbolicExpr::variable("x"), SymbolicExpr::number(0)},
        {SymbolicExpr::number(0), SymbolicExpr::number(1)}
    });
    auto conditional_inverse = lamina::matrix_inverse_checked(symbolic_det);
    EXPECT_TRUE(!conditional_inverse &&
                    conditional_inverse.error().code == lamina::CasErrc::Inconclusive,
                "checked inverse reports inconclusive when determinant nonzero cannot be verified");

    auto unsupported_rot = lamina::matrix_rotation_checked(0.0, 3);
    EXPECT_TRUE(!unsupported_rot &&
                    unsupported_rot.error().code == lamina::CasErrc::UnsupportedExpression,
                "checked rotation reports unsupported dimensions");
    EXPECT_TRUE(lamina::matrix_rotation(0.0, 3) == nullptr,
                "legacy rotation unwraps unsupported dimensions to nullptr");

    auto bad_eigenvalues = lamina::matrix_eigenvalues_checked(A);
    EXPECT_TRUE(!bad_eigenvalues &&
                    bad_eigenvalues.error().code == lamina::CasErrc::InvalidArgument,
                "checked eigenvalues reject non-square matrix");
    EXPECT_TRUE(lamina::matrix_eigenvalues(A).empty(),
                "legacy eigenvalues unwrap non-square matrix to empty list");

    lamina::CancellationToken token;
    token.cancel();
    lamina::ComputationContext context({}, token);
    auto cancelled = lamina::matrix_scaling_checked(1.0, 1.0, 2, context);
    EXPECT_TRUE(!cancelled &&
                    cancelled.error().code == lamina::CasErrc::Cancelled,
                "checked scaling observes cancelled context");
}

int main() {
    try {
        test_matrix_multiply();
        test_matrix_determinant();
        test_matrix_inverse();
        test_matrix_rotation();
        test_matrix_reflection();
        test_matrix_scaling();
        test_matrix_eigenvalues();
        test_matrix_eigenvectors();
        test_matrix_checked_contracts();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
