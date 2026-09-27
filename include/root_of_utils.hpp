/**
 * @file root_of_utils.hpp
 * @brief RootOf 恒等关系、求值及代数闭式化简。
 */
#pragma once

#include "root_of_identity.hpp"
#include "poly_utils.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>

namespace LMCAS {


/**
 * @brief 化简 RootOf 表达式（如可用根式表示则转换为闭式）。
 * @param rootof_expr RootOf 符号表达式
 * @return 化简后的符号表达式
 */
LMCAS_API std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);


}
