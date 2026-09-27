#pragma once

#include "polynomial.hpp"
#include "rational.hpp"
#include <utility>
#include <vector>

namespace LMCAS {

/**
 * @brief 查找有理系数多项式的所有有理根。
 * @param poly 有理系数多项式
 * @return 有理根列表，保留重数；零多项式和常数多项式返回空列表。
 * 公共非零有理系数先通过首一化消去；一次余式直接精确求解。
 */
LMCAS_API std::vector<Rational> find_rational_roots(const Polynomial<Rational>& poly);

/**
 * @brief 对有理系数多项式进行无平方因子分解。
 * @param poly 有理系数多项式
 * @return 因子与重数的列表
 */
LMCAS_API std::vector<std::pair<Polynomial<Rational>, int>> square_free_factorization(
    const Polynomial<Rational>& poly);

}
