/**
 * @file polynomial_conversion.hpp
 * @brief 符号表达式与一元多项式的转换。
 */
#pragma once

#include "computation_context.hpp"
#include "polynomial.hpp"
#include "symbolic.hpp"

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace LMCAS {

/**
 * @brief 将符号表达式转换为一元多项式。
 *
 * 支持 BigInt、Rational 和 SymbolicPolyCoeff 系数。
 * 不支持的子树或无法表示的系数使整个转换返回 UnsupportedExpression；
 * 成功返回的零多项式表示数学上的零。空输入或变量名返回 InvalidArgument。
 * 幂展开次数须为小于 1000 的非负整数；更大次数或超出 int 次数表示范围的
 * 系数存储返回 ResourceLimit。有限 double 按精确 binary64 值转换；
 * 证明源表达式精确性请用 recognize_rational_polynomial()。
 */
template <typename T>
LMCAS_API Result<Polynomial<T>> symbolic_to_poly(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable);

/**
 * @brief 将支持的一元多项式转换为符号表达式。
 *
 * 支持 BigInt、Rational 和 SymbolicPolyCoeff 系数。
 */
template <typename T>
LMCAS_API std::shared_ptr<SymbolicExpr> poly_to_symbolic(
    const Polynomial<T>& polynomial);

using OptionalRationalPolynomial = std::optional<Polynomial<Rational>>;


/**
 * @brief 证明表达式是精确的一元有理系数多项式。
 *
 * 成功且 optional 为空表示表达式超出当前结构支持域。
 * 精确系数仅取自精确数值节点；遍历与展开计入上下文预算。
 */
LMCAS_API Result<OptionalRationalPolynomial> recognize_rational_polynomial(
    const SymbolicExpr& expression,
    const std::string& variable,
    ComputationContext& context);

/** @brief 符号多项式算法使用的表达式系数。 */
struct SymbolicPolyCoeff {
    std::shared_ptr<SymbolicExpr> val;

    SymbolicPolyCoeff() : val(SymbolicExpr::number(0)) {}
    explicit SymbolicPolyCoeff(int value) : val(SymbolicExpr::number(value)) {}
    SymbolicPolyCoeff(std::shared_ptr<SymbolicExpr> value) : val(std::move(value)) {}

    bool operator==(const SymbolicPolyCoeff& other) const {
        if (!val || !other.val) return false;
        if (val == other.val) return true;
        auto difference = SymbolicExpr::add(
            val,
            SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1)));
        return difference->simplify()->is_zero();
    }

    bool operator!=(const SymbolicPolyCoeff& other) const {
        return !(*this == other);
    }

    SymbolicPolyCoeff operator+(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::add(val, other.val));
    }

    SymbolicPolyCoeff operator-(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::add(
            val,
            SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1))));
    }

    SymbolicPolyCoeff operator*(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::multiply(val, other.val));
    }

    SymbolicPolyCoeff operator/(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::divide(val, other.val));
    }

    SymbolicPolyCoeff operator-() const {
        return SymbolicPolyCoeff(
            SymbolicExpr::multiply(val, SymbolicExpr::number(-1)));
    }

    std::string to_string() const {
        return val ? val->to_string() : "0";
    }

    friend SymbolicPolyCoeff abs(const SymbolicPolyCoeff& value) {
        return value;
    }

    bool operator<(const SymbolicPolyCoeff&) const {
        return false;
    }
};

/** @brief 精确提取系数；不支持的值返回错误，保持数值完整。 */
template <typename T>
Result<T> extract_coeff_value(const std::shared_ptr<SymbolicExpr>& coefficient);

template <>
LMCAS_API Result<SymbolicPolyCoeff> extract_coeff_value<SymbolicPolyCoeff>(
    const std::shared_ptr<SymbolicExpr>& coefficient);

template <>
LMCAS_API Result<BigInt> extract_coeff_value<BigInt>(
    const std::shared_ptr<SymbolicExpr>& coefficient);

template <>
LMCAS_API Result<Rational> extract_coeff_value<Rational>(
    const std::shared_ptr<SymbolicExpr>& coefficient);

/** @brief 判断表达式是否含指定变量的自由出现。 */
LMCAS_API bool contains(
    const SymbolicExpr& expression,
    const std::string& variable);


}
