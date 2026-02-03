#include "test_common.hpp"

int main() {
    TEST_CASE("Solve System (Direct & Linear)");

    // 1. Linear System 2x2
    // 2x + y = 5
    // x - y = 1 => x=2, y=1
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        
        // 2x + y - 5
        auto eq1 = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::add(y, SymbolicExpr::number(-5))
        );
        
        // x - y - 1
        auto eq2 = SymbolicExpr::add(
            x,
            SymbolicExpr::add(SymbolicExpr::multiply(y, SymbolicExpr::number(-1)), SymbolicExpr::number(-1))
        );
        
        auto solutions = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});
        EXPECT_TRUE(solutions.size() == 1, "Solutions size should be 1");
        if (solutions.size() > 0) {
            EXPECT_EQ_EXPR(solutions[0]["x"], SymbolicExpr::number(2), "x should be 2");
            EXPECT_EQ_EXPR(solutions[0]["y"], SymbolicExpr::number(1), "y should be 1");
        }
    }
    
    // 2. Linear System with parameter
    // x + ay = 0
    // x - ay = 2 => 2x=2->x=1, 2ay=-2->ay=-1->y=-1/a
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto a = SymbolicExpr::variable("a");
        
        auto eq1 = SymbolicExpr::add(x, SymbolicExpr::multiply(a, y));
        auto eq2 = SymbolicExpr::add(x, SymbolicExpr::add(SymbolicExpr::multiply(a, SymbolicExpr::multiply(y, SymbolicExpr::number(-1))), SymbolicExpr::number(-2))); // x - ay - 2 = 0
        
        auto solutions = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});
        //std::cout << "Solving parametric system... solutions found: " << solutions.size() << std::endl;
        
        EXPECT_TRUE(solutions.size() == 1, "Solutions size should be 1");
        
        if (solutions.size() > 0) {
            //std::cout << "x = " << solutions[0]["x"]->to_string() << std::endl;
            //std::cout << "y = " << solutions[0]["y"]->to_string() << std::endl;
            EXPECT_EQ_EXPR(solutions[0]["x"], SymbolicExpr::number(1), "x should be 1");
            
            // y = -1/a
            auto expected_y = SymbolicExpr::power(SymbolicExpr::multiply(a, SymbolicExpr::number(-1)), SymbolicExpr::number(-1))->simplify(); 
            // wait, -1/a = -(a^-1)
            auto expected_y2 = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(a, SymbolicExpr::number(-1)))->simplify();
            EXPECT_EQ_EXPR(solutions[0]["y"], expected_y2, "y should be -a^-1");
        }
    }
    
    // 3. Eigenvalues
    // Matrix [[1, 0], [0, 2]] -> lambda = 1, 2
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(2)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto evals = SymbolicExpr::eigenvalues(m);
        // Should return a Vector expr containing 1 and 2
        
        // Output might be [1, 2] or [2, 1]
        // Currently solve returns vector<shared_ptr> 
        if (evals && evals->type == SymbolicExpr::Type::Vector) {
             //std::cout << "Eigenvalues: " << evals->to_string() << std::endl;
             // Check if it contains 1 and 2
             bool has1 = false, has2 = false;
             for(auto& op : evals->operands) {
                 if (op->to_string() == "1") has1 = true;
                 if (op->to_string() == "2") has2 = true;
             }
             EXPECT_TRUE(has1 && has2, "Eigenvalues should contain 1 and 2");
        } else {
            EXPECT_TRUE(false, "Eigenvalues return type mismatch");
        }
    }

    // 4. Eigenvectors
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(2)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto eigs = SymbolicExpr::eigenvectors(m);
        
        // Should have 2 eigenvalues
        // We might get them in any order
        EXPECT_TRUE(eigs.size() == 2, "Should return 2 eigenvectors");
        
        for(auto& p : eigs) {
            //std::cout << "Lambda: " << p.first->to_string() << " Vectors: ";
            //for(auto& v : p.second) std::cout << v->to_string() << " ";
            //std::cout << std::endl;
            
            //std::cout << "Checking Lambda: " << p.first->to_string() << std::endl;
            if (p.first->to_string() == "1") {
                // Expect [1, 0]
                EXPECT_TRUE(p.second.size() == 1, "size 1");
                auto& v = p.second[0];
                EXPECT_EQ_EXPR(v->operands[0], SymbolicExpr::number(1), "v[0]=1");
                EXPECT_EQ_EXPR(v->operands[1], SymbolicExpr::number(0), "v[1]=0");
            } else if (p.first->to_string() == "2") {
                // Expect [0, 1]
                EXPECT_TRUE(p.second.size() == 1, "size 1");
                auto& v = p.second[0];
                EXPECT_EQ_EXPR(v->operands[0], SymbolicExpr::number(0), "v[0]=0");
                EXPECT_EQ_EXPR(v->operands[1], SymbolicExpr::number(1), "v[1]=1");
            }
        }
    }
    return TEST_REPORT();
}
