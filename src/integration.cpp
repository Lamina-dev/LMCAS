#include "integration.hpp"

namespace lamina {

Integrator::Integrator() {
    
    
    
    
    
    
    
    
    
    
}

SymbolicExpr Integrator::integrate(const SymbolicExpr& expr, const std::string& var_name) {
    
    
    
    
    
    
    
    
    
    if (!expr.root) return expr;
    
    
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for (auto& op : add->operands) {
            new_ops.push_back(integrate(SymbolicExpr(op), var_name).root);
        }
        return SymbolicExpr(SymbolicFactory::create_add(new_ops));
    }
    
    
    
    
    
    
    auto x = SymbolicExpr(SymbolicFactory::create_variable(var_name));
    auto n = wildcard("n");
    
    
    auto x_pow_n = SymbolicExpr(std::make_shared<PowerNode>(x.root, n.root));
    
    MatchMap bindings;
    std::unordered_set<std::string> w = {"n"};
    
    if (Matcher::match(x_pow_n, expr, w, bindings)) {
        
        
        SymbolicExpr val_n = bindings["n"];
        
        
        
    }
    
    
    std::vector<std::shared_ptr<SymbolicNode>> args;
    args.push_back(expr.root);
    args.push_back(SymbolicFactory::create_variable(var_name));
    return SymbolicExpr(std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args)); 
}

SymbolicExpr Integrator::integrate_def(const SymbolicExpr& expr, const std::string& var_name, 
                                     const SymbolicExpr& lower, const SymbolicExpr& upper) {
    
    auto F = integrate(expr, var_name);
    
    
    
    
    return F; 
}

} 
