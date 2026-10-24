#pragma once

#include "symbolic.hpp"
#include "rational_polynomial.hpp"
#include "poly_utils.hpp"
#include "rational.hpp"
#include <vector>
#include <utility>
#include <memory>
#include <string>

namespace LMCAS {

LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_cubic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::string& var);

LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_quartic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e,
    const std::string& var);

LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_biquadratic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::string& var);


/**
 * @brief 通过因式分解求解多项式方程。
 * @param poly 符号系数多项式
 * @param var 求解变量名
 * @return 所有根的符号表达式列表
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> solve_by_factoring(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var);

}
