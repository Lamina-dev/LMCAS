
#include "symbolic.hpp"

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_ln() const {
    if (operands.empty()) return std::make_shared<SymbolicExpr>(*this);
    auto val = operands[0]->simplify();
    
    // ln(1) = 0
    if (val->is_number() && val->convert_rational() == ::Rational(1)) {
        return SymbolicExpr::number(0);
    }
    
    // ln(e) = 1
    if (val->type == SymbolicExpr::Type::Variable && val->identifier == "e") {
        return SymbolicExpr::number(1);
    }
    
    // ln(e^x) = x
    if (val->type == SymbolicExpr::Type::Power && val->operands.size() == 2) {
        auto base = val->operands[0]; 
        auto exp = val->operands[1];
        if (base->type == SymbolicExpr::Type::Variable && base->identifier == "e") {
            return exp;
        }
    }
    
    // ln(x^n) = n*ln(x)
    if (val->type == SymbolicExpr::Type::Power && val->operands.size() == 2) {
         return SymbolicExpr::multiply(val->operands[1], SymbolicExpr::ln(val->operands[0]))->simplify();
    }
    
    return SymbolicExpr::ln(val);
}
