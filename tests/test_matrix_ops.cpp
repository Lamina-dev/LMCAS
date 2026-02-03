#include "test_common.hpp"

int main() {
    TEST_CASE("Matrix Transpose & Inverse & RREF");

    // 1. Transpose: [[1, 2], [3, 4]] -> [[1, 3], [2, 4]]
    {
         std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m1_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(2)},
            {SymbolicExpr::number(3), SymbolicExpr::number(4)}
        };
        auto m1 = SymbolicExpr::matrix(m1_data);
        auto t1 = SymbolicExpr::transpose(m1);
        
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m1_t_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(3)},
            {SymbolicExpr::number(2), SymbolicExpr::number(4)}
        };
        auto m1_t_expected = SymbolicExpr::matrix(m1_t_data);
        
        EXPECT_EQ_EXPR(t1, m1_t_expected, "Transpose 2x2");
    }

    // 2. Inverse of [[4, 7], [2, 6]]
    // Det = 24 - 14 = 10
    // Adj = [[6, -7], [-2, 4]]
    // Inv = [[0.6, -0.7], [-0.2, 0.4]]
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m2_data = {
            {SymbolicExpr::number(4), SymbolicExpr::number(7)},
            {SymbolicExpr::number(2), SymbolicExpr::number(6)}
        };
        auto m2 = SymbolicExpr::matrix(m2_data);
        auto inv2 = SymbolicExpr::inverse(m2);
        
        // Check M * M^-1 == I
        auto prod = SymbolicExpr::multiply(m2, inv2)->simplify();
        // std::cout << "Prod: " << prod->to_string() << std::endl;
        
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto id = SymbolicExpr::matrix(id_data);
        EXPECT_EQ_EXPR(prod, id, "Inverse * Original == Identity");
    }

    // 3. RREF of [[1, 2, 3], [0, 1, 4], [5, 6, 0]]
    // Should solve to identity if full rank? Or something simpler
    // Let's try simple RREF: [[2, 4], [1, 5]]
    // R1/=2 -> [1, 2]
    // R2 = R2 - 1*R1 -> [0, 3]
    // R2/=3 -> [0, 1]
    // R1 = R1 - 2*R2 -> [1, 0]
    // Result: I
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto id = SymbolicExpr::matrix(id_data);  

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m3_data = {
            {SymbolicExpr::number(2), SymbolicExpr::number(4)},
            {SymbolicExpr::number(1), SymbolicExpr::number(5)}
        };
        auto m3 = SymbolicExpr::matrix(m3_data);
        auto rref3 = SymbolicExpr::rref(m3);
        EXPECT_EQ_EXPR(rref3, id, "RREF of 2x2 full rank matrix");
    }

    return TEST_REPORT();
}
