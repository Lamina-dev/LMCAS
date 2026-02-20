#include "test_common.hpp"

int main() {
    TEST_CASE("Matrix Transpose & Inverse & RREF");

    
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

    
    
    
    
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m2_data = {
            {SymbolicExpr::number(4), SymbolicExpr::number(7)},
            {SymbolicExpr::number(2), SymbolicExpr::number(6)}
        };
        auto m2 = SymbolicExpr::matrix(m2_data);
        auto inv2 = SymbolicExpr::inverse(m2);
        
        
        auto prod = SymbolicExpr::multiply(m2, inv2)->simplify();
        
        
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> id_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto id = SymbolicExpr::matrix(id_data);
        EXPECT_EQ_EXPR(prod, id, "Inverse * Original == Identity");
    }

    
    
    
    
    
    
    
    
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
