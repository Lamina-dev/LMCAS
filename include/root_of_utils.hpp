/**
 * @file root_of_utils.hpp
 * @brief RootOf 表达式的构造、化简与数值求值。
 */
#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "numeric_evaluation.hpp"
#include "lmmc/config.h"
#include <vector>
#include <memory>
#include <string>
#include <optional>

namespace LMCAS {

using RootOfEvaluationResult = Result<lmmc_real_t>;
using RootOfComplexEvaluationResult = Result<ApproxComplex>;
using RootOfConstructionResult = Result<std::shared_ptr<SymbolicExpr>>;

LMCAS_API RootOfConstructionResult make_rootof_checked(
    const std::shared_ptr<SymbolicExpr>& polynomial,
    const std::string& variable,
    std::size_t index,
    ComputationContext& context);

LMCAS_API RootOfConstructionResult make_rootof_checked(
    const std::shared_ptr<SymbolicExpr>& polynomial,
    const std::string& variable,
    std::size_t index);


/**
 * The polynomial must be structurally provable as an exact rational
 * polynomial. Real evaluation requires the selected root to be real;
 * selecting a non-real root returns `DomainError`. Malformed expressions
 * and invalid indices return `InvalidArgument`.
 */
LMCAS_API RootOfEvaluationResult rootof_evaluate_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr,
    ComputationContext& context);

LMCAS_API RootOfEvaluationResult rootof_evaluate_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

LMCAS_API RootOfComplexEvaluationResult rootof_evaluate_complex_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr,
    ComputationContext& context);

LMCAS_API RootOfComplexEvaluationResult rootof_evaluate_complex_checked(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

/**
 * @brief 化简 RootOf 表达式（如可用根式表示则转换为闭式）。
 * @param rootof_expr RootOf 符号表达式
 * @return 化简后的符号表达式
 */
LMCAS_API std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

/**
 * @brief 为不可约多项式构造 RootOf 解表达式列表。
 * @param poly 符号系数多项式
 * @param var 求解变量名
 * @return RootOf 表达式列表，每个元素对应多项式的一个根
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> make_rootof_solutions(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var);

}
