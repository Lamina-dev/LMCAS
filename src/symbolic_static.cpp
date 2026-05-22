std::shared_ptr<SymbolicExpr> SymbolicExpr::integral(std::shared_ptr<SymbolicExpr> op, const std::string& var) {
    if (!op) return nullptr;
    return op->integrate(var);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::limit_func(std::shared_ptr<SymbolicExpr> op, const std::string& var, std::shared_ptr<SymbolicExpr> target) {
    if (!op || !target) return nullptr;

    return op->limit(target);
}
