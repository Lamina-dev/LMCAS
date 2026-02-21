
// Static helper implementation for integral
std::shared_ptr<SymbolicExpr> SymbolicExpr::integral(std::shared_ptr<SymbolicExpr> op, const std::string& var) {
    if (!op) return nullptr;
    return op->integrate(var);
}

// Static helper implementation for limit_func
std::shared_ptr<SymbolicExpr> SymbolicExpr::limit_func(std::shared_ptr<SymbolicExpr> op, const std::string& var, std::shared_ptr<SymbolicExpr> target) {
    if (!op || !target) return nullptr;
    // Call the member function limit()
    // limit takes a double/int/SymbolicExpr?
    // Let's check the signature of limit member function in symbolic.hpp/cpp
    // It seems limit takes (const std::shared_ptr<SymbolicExpr>& target, const std::string& direction)
    return op->limit(target);
}
