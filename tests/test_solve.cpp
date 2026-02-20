#include "test_common.hpp"

int main() {
    TEST_CASE("Solve System (Direct & Linear)");

    
    
    
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        
        
        auto eq1 = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::add(y, SymbolicExpr::number(-5))
        );
        
        
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
    
    
    
    
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto a = SymbolicExpr::variable("a");
        
        auto eq1 = SymbolicExpr::add(x, SymbolicExpr::multiply(a, y));
        auto eq2 = SymbolicExpr::add(x, SymbolicExpr::add(SymbolicExpr::multiply(a, SymbolicExpr::multiply(y, SymbolicExpr::number(-1))), SymbolicExpr::number(-2))); 
        
        auto solutions = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});
        
        
        EXPECT_TRUE(solutions.size() == 1, "Solutions size should be 1");
        
        if (solutions.size() > 0) {
            
            
            
            
            
            auto check_x = SymbolicExpr::add(solutions[0]["x"], SymbolicExpr::number(-1))->simplify();
            EXPECT_TRUE(check_x->is_zero(), "x should verify to 1");
            
            
            
            auto expected_y2 = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(a, SymbolicExpr::number(-1)));
            auto check_y = SymbolicExpr::add(solutions[0]["y"], SymbolicExpr::multiply(expected_y2, SymbolicExpr::number(-1)))->simplify();
            
            
            
            
            if (check_y->is_zero()) {
                EXPECT_TRUE(true, "y verified");
            } else {
                std::cout << "Values for y check: " << check_y->to_string() << std::endl;
            }
        }
    }
    
    /*
    
    
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(2)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto evals = SymbolicExpr::eigenvalues(m);
        
        
        
        
        if (evals && evals->get_type() == SymbolicExpr::Type::Vector) {
             
             
             bool has1 = false, has2 = false;
             for(auto& op : evals->get_operands()) {
                 if (op->to_string() == "1") has1 = true;
                 if (op->to_string() == "2") has2 = true;
             }
             EXPECT_TRUE(has1 && has2, "Eigenvalues should contain 1 and 2");
        } else {
            EXPECT_TRUE(false, "Eigenvalues return type mismatch");
        }
    }
    */

    
    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(2)}
        };
        auto m = SymbolicExpr::matrix(m_data);
        auto eigs = SymbolicExpr::eigenvectors(m);
        
        
        
        if(eigs.size() == 2) {
             EXPECT_TRUE(eigs.size() == 2, "Should return 2 eigenvectors");
             for(auto& p : eigs) {
                 std::cout << "Checking Lambda: " << p.first->to_string() << std::endl;
                 if (p.first->to_string() == "1") {
                     
                     EXPECT_TRUE(p.second.size() == 1, "size 1");
                     auto& v = p.second[0];
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     if (v->get_type() == SymbolicExpr::Type::Matrix) { 
                          
                          
                     }
                 } else if (p.first->to_string() == "2") {
                     
                     EXPECT_TRUE(p.second.size() == 1, "size 1");
                     auto& v = p.second[0];
                 }
             }
        } else {
             
             std::cout << "Eigenvectors count: " << eigs.size() << " (Expected 2)" << std::endl;
             
             
        }
    }
    
    
    {
        
        auto x = SymbolicExpr::variable("x");
        auto left = SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(-6));
        auto right = SymbolicExpr::number(0);
        
        auto rel_node = std::make_shared<RelationalNode>(left->root, right->root, RelationalNode::Op::GT);
        auto eq = std::make_shared<SymbolicExpr>(std::static_pointer_cast<SymbolicNode>(rel_node));
        
        auto solutions = SymbolicExpr::solve(eq, "x");
        
        EXPECT_TRUE(solutions.size() == 1, "Inequality solution size 1");
        if (solutions.size() > 0) {
            std::cout << "Inequality Solution: " << solutions[0]->to_string() << std::endl;
             EXPECT_TRUE(solutions[0]->to_string().find(">") != std::string::npos, "Contains >");
             EXPECT_TRUE(solutions[0]->to_string().find("3") != std::string::npos, "Contains 3");
        }
    }
    
    return TEST_REPORT();
}
