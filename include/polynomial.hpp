/**
 * @file polynomial.hpp
 * @brief 一元多项式模板类 Polynomial<T>,支持四则运算,GCD,求导,求值.
 */
#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include "bigint.hpp"
#include "rational.hpp"

namespace LMCAS {

/**
 * @brief 计算两个系数的最大公约数(泛型实现)
 * @tparam T 系数类型
 * @param a 第一个系数
 * @param b 第二个系数
 * @return gcd(a, b)
 */
template<typename T>
T gcd_coeff_impl(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        return std::abs(std::gcd(a, b));
    }
    else if constexpr (std::is_same_v<T, BigInt>) {
        return BigInt::gcd(a, b);
    }
    else if constexpr (std::is_same_v<T, Rational>) {
        return Rational(1);
    }
    else {
        return T(1);
    }
}

/**
 * @brief 一元多项式模板类
 * @tparam CoeffType 系数类型(支持 BigInt,Rational,int 等)
 *
 * 系数按升幂存储:coeffs[i] 对应 x^i 的系数.
 * 非零操作数必须属于同一命名变量环（包括显式常数），否则算术、
 * 带余除法、伪除法和 GCD 抛出 std::invalid_argument，判等返回 false。
 * 加减乘中的零元素嵌入非零侧变量环；双方为零时保留左侧变量环。
 * 零被除数的商和余数采用除数变量环；零除数始终非法。
 */
template <typename CoeffType>
class Polynomial {
public:

    std::vector<CoeffType> coeffs;   ///< 系数数组,coeffs[i] 为 x^i 的系数
    std::string variable_name;       ///< 变量名

    /**
     * @brief 构造零多项式
     * @param var 变量名,默认 "x"
     */
    Polynomial(const std::string& var = "x") : variable_name(var) {}

    /**
     * @brief 从系数向量构造多项式
     * @param c 系数向量(升幂排列)
     * @param var 变量名
     */
    Polynomial(const std::vector<CoeffType>& c, const std::string& var = "x")
        : coeffs(c), variable_name(var) {
        trim();
    }

    /**
     * @brief 从单个常数构造常数多项式
     * @param c 常数值
     * @param var 变量名
     */
    Polynomial(const CoeffType& c, const std::string& var = "x")
        : coeffs{c}, variable_name(var) {
        trim();
    }

    /** @brief 去除高次零系数 */
    void trim() {
        while (coeffs.size() > 0 && coeffs.back() == CoeffType(0)) {
            coeffs.pop_back();
        }
    }

    /** @brief 判断是否为零多项式 */
    bool is_zero() const {
        return coeffs.empty();
    }

    /**
     * @brief 获取多项式的次数
     * @return 次数;零多项式返回 -1
     */
    int degree() const {
        if (coeffs.empty()) return -1;
        return static_cast<int>(coeffs.size()) - 1;
    }

    /**
     * @brief 获取首项系数
     * @return 最高次项的系数;零多项式返回 0
     */
    CoeffType lead_coeff() const {
        if (coeffs.empty()) return CoeffType(0);
        return coeffs.back();
    }

    /**
     * @brief 多项式加法
     * @param other 加数
     * @return 和多项式
     */
    Polynomial operator+(const Polynomial& other) const {
        require_compatible_variables(other, "Polynomial::operator+");
        const auto& result_var = is_zero() && !other.is_zero()
            ? other.variable_name : variable_name;
        Polynomial res(result_var);
        size_t n = std::max(coeffs.size(), other.coeffs.size());
        res.coeffs.resize(n);

        for (size_t i = 0; i < n; ++i) {
            CoeffType a = (i < coeffs.size()) ? coeffs[i] : CoeffType(0);
            CoeffType b = (i < other.coeffs.size()) ? other.coeffs[i] : CoeffType(0);
            res.coeffs[i] = a + b;
        }
        res.trim();
        return res;
    }

    /**
     * @brief 多项式减法
     * @param other 减数
     * @return 差多项式
     */
    Polynomial operator-(const Polynomial& other) const {
         require_compatible_variables(other, "Polynomial::operator-");
         const auto& result_var = is_zero() && !other.is_zero()
             ? other.variable_name : variable_name;
         Polynomial res(result_var);
         size_t n = std::max(coeffs.size(), other.coeffs.size());
         res.coeffs.resize(n);

         for (size_t i = 0; i < n; ++i) {
             CoeffType a = (i < coeffs.size()) ? coeffs[i] : CoeffType(0);
             CoeffType b = (i < other.coeffs.size()) ? other.coeffs[i] : CoeffType(0);
             res.coeffs[i] = a - b;
         }
         res.trim();
         return res;
    }

    /**
     * @brief 多项式乘法
     * @param other 乘数
     * @return 积多项式
     */
    Polynomial operator*(const Polynomial& other) const {
        require_compatible_variables(other, "Polynomial::operator*");
        const auto& result_var = is_zero() && !other.is_zero()
            ? other.variable_name : variable_name;
        Polynomial res(result_var);
        if (is_zero() || other.is_zero()) { return res; }
        res.coeffs.resize(coeffs.size() + other.coeffs.size() - 1, CoeffType(0));

        for (size_t i = 0; i < coeffs.size(); ++i) {
            for (size_t j = 0; j < other.coeffs.size(); ++j) {
                res.coeffs[i + j] = res.coeffs[i + j] + coeffs[i] * other.coeffs[j];
            }
        }
        res.trim();
        return res;
    }

    /**
     * @brief 多项式判等
     * @param other 比较对象
     * @return 同一命名变量环内系数完全相同返回 true；所有零多项式相等
     */
    bool operator==(const Polynomial& other) const {
         if (degree() != other.degree()) {
             return false;
         }
         if (!variables_compatible(other)) {
             return false;
         }
         for (size_t i = 0; i < coeffs.size(); ++i) {
             if (coeffs[i] != other.coeffs[i]) {
                 return false;
             }
         }
         return true;
    }

    /**
     * @brief 多项式带余除法
     * @param other 除数多项式
     * @return pair(商, 余数)
     * @throw std::runtime_error 除数为零多项式时抛出
     * @throw std::invalid_argument 两个非零操作数的变量名不同时抛出
     */
    std::pair<Polynomial, Polynomial> div_mod(const Polynomial& other) const {
        if (other.is_zero()) throw std::runtime_error("Division by zero polynomial");
        require_compatible_variables(other, "Polynomial::div_mod");
        if (is_zero()) {
            return {Polynomial(other.variable_name), Polynomial(other.variable_name)};
        }

        Polynomial quotient(variable_name);
        Polynomial remainder = *this;

        int deg_rem = remainder.degree();
        int deg_div = other.degree();
        CoeffType lc_div = other.lead_coeff();

        if (deg_rem < deg_div) {
             return {quotient, remainder};
        }

        quotient.coeffs.resize(deg_rem - deg_div + 1, CoeffType(0));

        while (deg_rem >= deg_div && !remainder.is_zero()) {

             CoeffType factor = remainder.lead_coeff() / lc_div;
             int deg_diff = deg_rem - deg_div;

             quotient.coeffs[deg_diff] = factor;

             Polynomial term(variable_name);
             term.coeffs.resize(deg_diff + 1, CoeffType(0));
             term.coeffs[deg_diff] = factor;

             Polynomial subtrahend = other * term;
             remainder = remainder - subtrahend;
             deg_rem = remainder.degree();
        }

        quotient.trim();
        return {quotient, remainder};
    }

    /**
     * @brief 计算两个系数的 GCD(静态方法)
     * @param a 第一个系数
     * @param b 第二个系数
     * @return gcd(a, b)
     */
    static CoeffType gcd_coeff(const CoeffType& a, const CoeffType& b) {
        if constexpr (std::is_same_v<CoeffType, BigInt>) {
            return BigInt::gcd(a, b);
        } else if constexpr (std::is_same_v<CoeffType, Rational>) {
            return CoeffType(1);
        } else {

             if constexpr (std::is_integral_v<CoeffType>) {
                 return std::gcd(a, b);
             }

             return CoeffType(1);
        }
    }

    /**
     * @brief 计算多项式的容度(所有系数的 GCD)
     * @return 容度值
     */
    CoeffType content() const;

    /**
     * @brief 计算多项式的本原部分(除以容度)
     * @return 本原多项式
     */
    Polynomial primitive_part() const;

    /**
     * @brief 通过系数缩放执行伪除法,运算保持在原系数环
     * @param other 除数多项式
     * @return pair(伪商, 伪余数)
     * @throw std::runtime_error 除数为零多项式时抛出
     * @throw std::invalid_argument 两个非零操作数的变量名不同时抛出
     */
    std::pair<Polynomial, Polynomial> pseudo_div_mod(const Polynomial& other) const;

    /**
     * @brief 对多项式求值(Horner 法)
     * @param val 自变量的值
     * @return 多项式在 val 处的值
     */
    CoeffType eval(const CoeffType& val) const {
        if (coeffs.empty()) return CoeffType(0);
        CoeffType res = coeffs.back();
        for (int i = (int)coeffs.size() - 2; i >= 0; --i) {
            res = res * val + coeffs[i];
        }
        return res;
    }

    /**
     * @brief 求导
     * @return 导多项式
     */
    Polynomial differentiate() const {
        if (degree() < 1) return Polynomial(variable_name);
        Polynomial res(variable_name);
        res.coeffs.resize(coeffs.size() - 1);
        for (size_t i = 1; i < coeffs.size(); ++i) {
            CoeffType c = coeffs[i];

            res.coeffs[i-1] = c * CoeffType(i);
        }
        res.trim();
        return res;
    }

    /**
     * @brief 首一化(使首项系数为 1)
     * @return 首一多项式
     */
    Polynomial make_monic() const {
        if (is_zero()) return *this;
        CoeffType lc = lead_coeff();
        if (lc == CoeffType(1)) return *this;

        Polynomial res(variable_name);
        res.coeffs.reserve(coeffs.size());
        for (const auto& c : coeffs) {
            res.coeffs.push_back(c / lc);
        }
        return res;
    }

    /**
     * @brief 伪除法取余数
     * @param other 除数多项式
     * @return 伪余数
     */
    Polynomial pseudo_div_mod_rem(const Polynomial& other) const;

    /**
     * @brief 计算两个多项式的最大公因式
     * @param a 第一个多项式
     * @param b 第二个多项式
     * @return gcd(a, b)
     * @throw std::invalid_argument 两个非零操作数的变量名不同时抛出
     */
    static Polynomial gcd(Polynomial a, Polynomial b);

    /**
     * @brief 计算无平方因子部分
     * @return 去除重因子后的多项式
     */
    Polynomial square_free_part() const;

    /**
     * @brief 转换为字符串表示
     * @return 多项式的字符串形式
     */
    std::string to_string() const {
        if (is_zero()) {
            return "0";
        }
        std::string s;
        for (int i = degree(); i >= 0; --i) {
            if (coeffs[i] == CoeffType(0)) {
                continue;
            }
            append_string_term(s, coeffs[i], i);
        }
        return s.empty() ? "0" : s;
    }

private:
    bool variables_compatible(const Polynomial& other) const {
        return variable_name == other.variable_name || is_zero() || other.is_zero();
    }

    void require_compatible_variables(const Polynomial& other, const char* operation) const {
        if (!variables_compatible(other)) {
            throw std::invalid_argument(
                std::string(operation) + ": variable name mismatch ('" +
                variable_name + "' vs '" + other.variable_name + "')");
        }
    }

    void append_string_term(std::string& text, CoeffType coefficient, int power) const {
        const bool positive = !(coefficient < CoeffType(0));
        if (!text.empty()) {
            text += positive ? " + " : " - ";
        } else if (!positive) {
            text += "-";
        }
        if (!positive) {
            coefficient = coefficient * CoeffType(-1);
        }
        if (!(coefficient == CoeffType(1)) || power == 0) {
            if constexpr (std::is_same_v<CoeffType, BigInt> ||
                          std::is_same_v<CoeffType, Rational>) {
                text += coefficient.to_string();
            } else {
                text += std::to_string(coefficient);
            }
            if (power > 0) {
                text += "*";
            }
        }
        if (power > 0) {
            text += variable_name;
        }
        if (power > 1) {
            text += "^" + std::to_string(power);
        }
    }
};

/**
 * @brief 多项式流输出运算符
 * @tparam T 系数类型
 * @param os 输出流
 * @param p 多项式
 * @return 输出流引用
 */
template<typename T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& p) {
    if (p.is_zero()) {
        return os << "0";
    }
    for (int i = p.degree(); i >= 0; --i) {
        T c = p.coeffs[i];
        if (c == T(0)) {
            continue;
        }

        if (i < p.degree()) {
            if (c > T(0)) {
                os << " + ";
            } else {
                os << " - ";
            }
        } else if (c < T(0)) {
            os << "-";
        }

        T abs_c = (c < T(0)) ? (c * T(-1)) : c;
        if (abs_c != T(1) || i == 0) {
            os << abs_c;
        }

        if (i > 0) {
            os << p.variable_name;
        }
        if (i > 1) {
            os << "^" << i;
        }
    }
    return os;
}

}

#include "polynomial_division.hpp"
