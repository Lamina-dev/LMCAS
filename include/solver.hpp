#pragma once
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <map>

namespace lamina {

class LAMINA_API Solver {
public:
    
    
    
    
    
    static std::map<std::string, SymbolicExpr> solve_linear_system(
        const std::vector<SymbolicExpr>& equations, 
        const std::vector<std::string>& variables);

    
    
    static std::vector<SymbolicExpr> groebner_basis(
        const std::vector<SymbolicExpr>& polynomials,
        const std::vector<std::string>& variables);

    
    
    
    
    static std::vector<std::map<std::string, SymbolicExpr>> solve_polynomial_system(
        const std::vector<SymbolicExpr>& equations,
        const std::vector<std::string>& variables);
};

} 
