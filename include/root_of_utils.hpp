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

// RootOf 数值求值：返回第 k 个根的数值近似
LAMINA_API std::optional<lmmc_real_t> rootof_evaluate(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

// RootOf 化简：检查是否可以用闭合形式表示
LAMINA_API std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

// 生成 RootOf 表达式列表（对不可约多项式）
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> make_rootof_solutions(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var);

} // namespace lamina
