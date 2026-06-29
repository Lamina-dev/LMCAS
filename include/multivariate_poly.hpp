/**
 * @file multivariate_poly.hpp
 * @brief 多元多项式类 MultiPoly，稀疏表示，支持任意数量变量上的环运算。
 *
 * 以 (Monomial, Rational) 对的有序列表表示多元多项式，单项式按指定的
 * MonomialOrder 排序。提供加减乘、求值、度数查询、一元转换等完整接口。
 */
#pragma once

#include "monomial_order.hpp"
#include "rational.hpp"
#include "polynomial.hpp"

#include <map>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#ifdef LAMINA_CORE_EXPORTS
#define LAMINA_API __declspec(dllexport)
#else
#define LAMINA_API __declspec(dllimport)
#endif
#else
#define LAMINA_API
#endif

namespace lamina {

/**
 * @brief 多元多项式类，稀疏表示
 *
 * 以 (Monomial, Rational) 对的有序列表表示多元多项式。
 * 单项式按指定的 MonomialOrder 排序，系数为有理数。
 * 支持加减乘、求值、度数查询、一元转换等完整环运算接口。
 */
class LAMINA_API MultiPoly {
public:
    /// 项类型：(单项式, 系数)
    using Term = std::pair<Monomial, Rational>;

    /// @brief 构造零多项式
    MultiPoly();

    /**
     * @brief 从项列表和变量名列表构造多元多项式
     * @param[in] terms 项列表，每项为 (单项式, 系数) 对
     * @param[in] vars 变量名列表，确定各分量对应的变量
     * @param[in] order 单项式序类型，默认分次逆字典序
     */
    MultiPoly(std::vector<Term> terms, std::vector<std::string> vars,
              MonomialOrderType order = MonomialOrderType::GrevLex);

    /**
     * @brief 从常数构造多元多项式
     * @param[in] constant 常数值
     * @param[in] vars 变量名列表
     */
    explicit MultiPoly(const Rational& constant, const std::vector<std::string>& vars);

    /// --- 算术运算 ---

    /** @brief 多项式加法 */
    MultiPoly operator+(const MultiPoly& other) const;

    /** @brief 多项式减法 */
    MultiPoly operator-(const MultiPoly& other) const;

    /** @brief 多项式乘法 */
    MultiPoly operator*(const MultiPoly& other) const;

    /** @brief 一元取负 */
    MultiPoly operator-() const;

    /** @brief 判等 */
    bool operator==(const MultiPoly& other) const;

    /** @brief 判不等 */
    bool operator!=(const MultiPoly& other) const;

    /**
     * @brief 标量乘法
     * @param[in] scalar 有理数标量
     * @return 各项系数乘以 scalar 后的多项式
     */
    MultiPoly operator*(const Rational& scalar) const;

    /**
     * @brief 精确除法（假设整除）
     * @param[in] divisor 除数多项式
     * @return 商多项式
     * @throw std::runtime_error 除数为零或不整除时抛出
     */
    MultiPoly exact_div(const MultiPoly& divisor) const;

    /// --- 求值 ---

    /**
     * @brief 将指定变量代入有理数值，返回降维后的多项式
     * @param[in] var 待代入的变量名
     * @param[in] val 代入值
     * @return 代入后的多项式（变量数减一或不变）
     */
    MultiPoly eval(const std::string& var, const Rational& val) const;

    /**
     * @brief 将多个变量同时代入有理数值
     * @param[in] substitution 变量名到值的映射
     * @return 代入后的多项式
     */
    MultiPoly eval(const std::map<std::string, Rational>& substitution) const;

    /// --- 度数查询 ---

    /**
     * @brief 计算多项式的全次数
     * @return 所有项中单项式全次数的最大值；零多项式返回 -1
     */
    int total_degree() const;

    /**
     * @brief 计算关于指定变量的次数
     * @param[in] var 变量名
     * @return 该变量在所有项中指数的最大值；零多项式返回 -1
     */
    int degree(const std::string& var) const;

    /**
     * @brief 获取关于指定变量的首项系数
     *
     * 将多项式视为 var 的一元多项式，返回最高次项的系数（为辅助变量的多项式）。
     * @param[in] var 主变量名
     * @return 首项系数多项式
     */
    MultiPoly leading_coeff(const std::string& var) const;

    /// --- 转换 ---

    /**
     * @brief 转换为一元 Polynomial<Rational>
     *
     * 仅当多项式实际只含一个变量时有效。
     * @return 对应的一元多项式
     * @throw std::invalid_argument 含多个有效变量时抛出
     */
    Polynomial<Rational> to_univariate() const;

    /**
     * @brief 从一元多项式构造多元多项式
     * @param[in] poly 一元多项式
     * @param[in] var 对应的变量名
     * @return 等价的 MultiPoly 表示
     */
    static MultiPoly from_univariate(const Polynomial<Rational>& poly,
                                     const std::string& var);

    /// --- 属性查询 ---

    /** @brief 判断是否为零多项式 */
    bool is_zero() const;

    /** @brief 判断是否为常数多项式 */
    bool is_constant() const;

    /** @brief 判断是否实质为一元多项式（仅一个变量出现非零指数） */
    bool is_univariate() const;

    /** @brief 判断是否为齐次多项式（所有项全次数相同） */
    bool is_homogeneous() const;

    /** @brief 获取非零项数 */
    int num_terms() const;

    /** @brief 获取变量个数 */
    int num_vars() const;

    /**
     * @brief 获取变量名列表
     * @return 变量名列表的常引用
     */
    const std::vector<std::string>& variables() const;

    /**
     * @brief 获取项列表
     * @return 有序项列表的常引用
     */
    const std::vector<Term>& terms() const;

    /// --- 工具方法 ---

    /**
     * @brief 提取数值内容（所有系数的 GCD）
     * @return 所有系数绝对值的最大公约数（有理数形式）
     */
    Rational numeric_content() const;

    /**
     * @brief 使多项式本原化（除以数值内容）
     * @return 本原多项式（数值内容为 1）
     */
    MultiPoly make_primitive() const;

    /**
     * @brief 转换为可读字符串表示
     * @return 多项式的字符串形式
     */
    std::string to_string() const;

private:
    std::vector<Term> terms_;           ///< 有序项列表（按 order_ 降序）
    std::vector<std::string> vars_;     ///< 变量名列表
    MonomialOrderType order_;           ///< 单项式序类型

    /**
     * @internal
     * @brief 规范化：合并同类项、去零项、按单项式序排序
     */
    void normalize();
};

} // namespace lamina
