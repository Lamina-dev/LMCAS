#pragma once
#include "symbolic.hpp"
#include "matcher.hpp"

namespace lamina {

class LAMINA_API Integrator {
    RewriteEngine engine;

public:
    Integrator();

    
    SymbolicExpr integrate(const SymbolicExpr& expr, const std::string& var_name);
    
    
    SymbolicExpr integrate_def(const SymbolicExpr& expr, const std::string& var_name, 
                              const SymbolicExpr& lower, const SymbolicExpr& upper);
};

} 
