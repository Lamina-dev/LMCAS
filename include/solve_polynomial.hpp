#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "rational.hpp"
#include <vector>
#include <utility>
#include <memory>
#include <string>

namespace lamina {

// 求解 ax^3 + bx^2 + cx + d = 0，返回 1-3 个符号根
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_cubic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::string& var);

// 求解 ax^4 + bx^3 + cx^2 + dx + e = 0
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_quartic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e,
    const std::string& var);

// 双二次快捷: ax^4 + bx^2 + c = 0
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_biquadratic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::string& var);

// 有理根定理：尝试找到有理根
LAMINA_API std::vector<Rational> find_rational_roots(const Polynomial<Rational>& poly);

// 无平方因子分解：将多项式分解为 f₁·f₂²·f₃³·… 的形式
// 返回 (因子, 重数) 对的向量，每个因子是无平方的
LAMINA_API std::vector<std::pair<Polynomial<Rational>, int>> square_free_factorization(
    const Polynomial<Rational>& poly);

// 因式分解后递归求解
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_by_factoring(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var);

} // namespace lamina
