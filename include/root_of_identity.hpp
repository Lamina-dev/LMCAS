/**
 * @file root_of_identity.hpp
 * @brief 精确 RootOf 标识构造及受检数值求值。
 */
#pragma once

#include "symbolic.hpp"
#include "computation_context.hpp"
#include "numeric_value.hpp"
#include "polynomial_conversion.hpp"
#include <vector>
#include <cstddef>
#include <memory>
#include <string>

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
 * @brief 对 RootOf 表达式进行实数求值。
 * @pre 多项式须能从结构上证明为精确有理系数多项式，所选根须为实数。
 * @return 选中非实根时返回 DomainError；表达式格式或索引无效时返回 InvalidArgument。
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
 * @brief 为不可约多项式构造 RootOf 解表达式列表。
 * @param poly 符号系数多项式
 * @param var 求解变量名
 * @return RootOf 表达式列表，每个元素对应多项式的一个根
 */
LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> make_rootof_solutions(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var);

}
