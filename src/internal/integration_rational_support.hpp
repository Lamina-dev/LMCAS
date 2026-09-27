#pragma once
#include "internal/integration_support.hpp"

namespace LMCAS {
inline std::shared_ptr<SymbolicExpr> rd_var_minus(const std::string& var, const Rational& r) {
    auto v = SymbolicExpr::variable(var);
    if (r == Rational(0)) { return v; }
    auto neg_r = SymbolicExpr::number(Rational(0) - r);
    return SymbolicExpr::add(v, neg_r);
}


}
