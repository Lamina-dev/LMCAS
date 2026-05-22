/**
 * @file root_of_utils.hpp
 * @brief RootOf 表达式的构造、化简与数值求值。
 */
#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "lmmc/config.h"
#include <vector>
#include <memory>
#include <string>
#include <optional>

namespace lamina {

/**
 * @brief 对 RootOf 表达式进行数值求值。
 * @param rootof_expr RootOf 符号表达式
 * @return 若求值成功则返回数值结果，否则返回 nullopt
 */
LAMINA_API std::optional<lmmc_real_t> rootof_evaluate(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

/**
 * @brief 化简 RootOf 表达式（如可用根式表示则转换为闭式）。
 * @param rootof_expr RootOf 符号表达式
 * @return 化简后的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

/**
 * @brief 为不可约多项式构造 RootOf 解表达式列表。
 * @param poly 符号系数多项式
 * @param var 求解变量名
 * @return RootOf 表达式列表，每个元素对应多项式的一个根
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> make_rootof_solutions(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var);

}
