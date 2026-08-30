/**
 * @file irrational.hpp
 * @brief 无理数表示类 Irrational,支持 sqrt,pi,e 及其线性组合.
 */
#pragma once
#define _USE_MATH_DEFINES
#include "symbolic.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#ifndef LMMC_E
#define LMMC_E 2.71828182845904523536
#endif

/**
 * @brief 计算浮点数的平方根(封装 LMMC 接口)
 * @param x 被开方数
 * @return 平方根值
 */
inline lmmc_real_t _irrational_sqrt(lmmc_real_t x) {
    lmmc_real_t res;
    LMMC_REAL_SQRT(&res, &x);
    return res;
}

/** @brief 无理数表示类,支持 sqrtn,pi,e,log 及其线性组合的精确表示与运算 */
class Irrational {
public:

    /**
     * @brief 将无理数转换为符号表达式
     * @return 对应的符号表达式
     */
    std::shared_ptr<SymbolicExpr> to_symbolic() const {
        auto is_zero_tol = [](lmmc_real_t x) -> bool {
            int eq;
            lmmc_double_nearly_equal_tol(x, 0.0, 1e-15, 1e-15, &eq);
            return eq != 0;
        };
        auto is_one_tol = [](lmmc_real_t x) -> bool {
            int eq;
            lmmc_double_nearly_equal_tol(x, 1.0, 1e-15, 1e-15, &eq);
            return eq != 0;
        };
        switch (type) {
            case Type::SQRT: {

                auto sqrtExpr = SymbolicExpr::sqrt(SymbolicExpr::number(static_cast<int>(radicand)));
                if (is_zero_tol(coefficient)) {
                    return SymbolicExpr::number(0);
                } else if (is_one_tol(coefficient)) {
                    return sqrtExpr;
                } else {
                    return SymbolicExpr::multiply(SymbolicExpr::number(::Rational::from_double(coefficient)), sqrtExpr);
                }
            }
            case Type::PI:

                if (is_zero_tol(coefficient)) {
                    return SymbolicExpr::number(0);
                } else if (is_one_tol(coefficient)) {
                    return SymbolicExpr::variable("π");
                } else {
                    return SymbolicExpr::multiply(SymbolicExpr::number(::Rational::from_double(coefficient)), SymbolicExpr::variable("π"));
                }
            case Type::E:
                if (is_zero_tol(coefficient)) {
                    return SymbolicExpr::number(0);
                } else if (is_one_tol(coefficient)) {
                    return SymbolicExpr::variable("e");
                } else {
                    return SymbolicExpr::multiply(SymbolicExpr::number(::Rational::from_double(coefficient)), SymbolicExpr::variable("e"));
                }
            case Type::LOG:

                if (is_zero_tol(coefficient)) {
                    return SymbolicExpr::number(0);
                } else {
                    return SymbolicExpr::multiply(SymbolicExpr::number(::Rational::from_double(coefficient)), SymbolicExpr::variable("log(" + std::to_string(radicand) + ")"));
                }
            case Type::COMPLEX: {
                /// 重建完整线性组合,使 constant_term 与全部基向量共同参与结果.
                std::shared_ptr<SymbolicExpr> result;
                auto add_term = [&](const std::shared_ptr<SymbolicExpr>& term) {
                    if (!result) {
                        result = term;
                    } else {
                        result = SymbolicExpr::add(result, term);
                    }
                };

                if (!is_zero_tol(constant_term)) {
                    add_term(SymbolicExpr::number(::Rational::from_double(constant_term)));
                }

                for (const auto& [key, coeff] : coefficients) {
                    if (is_zero_tol(coeff)) continue;
                    std::shared_ptr<SymbolicExpr> basis;
                    if (key == "pi") {
                        basis = SymbolicExpr::variable("π");
                    } else if (key == "e") {
                        basis = SymbolicExpr::variable("e");
                    } else if (key.substr(0, 4) == "sqrt") {
                        long long n = std::stoll(key.substr(4));
                        basis = SymbolicExpr::sqrt(SymbolicExpr::number(static_cast<int>(n)));
                    } else {
                        /// 将未知基向量保留为变量,维持原始符号信息.
                        basis = SymbolicExpr::variable(key);
                    }
                    std::shared_ptr<SymbolicExpr> term = is_one_tol(coeff)
                        ? basis
                        : SymbolicExpr::multiply(SymbolicExpr::number(::Rational::from_double(coeff)), basis);
                    add_term(term);
                }

                return result ? result : SymbolicExpr::number(0);
            }
            default:
                return SymbolicExpr::number(0);
        }
    }
    /** @brief 无理数的内部类型 */
    enum class Type {
        SQRT,     ///< 平方根形式 coeff * sqrtradicand
        PI,       ///< pi 的倍数
        E,        ///< e 的倍数
        LOG,      ///< 对数形式
        COMPLEX   ///< 多项线性组合
    };

private:
    Type type;

    lmmc_real_t coefficient;
    long long radicand;

    std::map<std::string, lmmc_real_t> coefficients;
    lmmc_real_t constant_term;

    static std::pair<long long, long long> simplify_sqrt(long long n) {
        long long perfect_square = 1;
        long long remainder = n;

        for (long long i = 2; i * i <= n; ++i) {
            while (remainder % (i * i) == 0) {
                perfect_square *= i;
                remainder /= (i * i);
            }
        }
        return {perfect_square, remainder};
    }

public:

    /** @brief 默认构造,初始化为零值的 COMPLEX 类型 */
    Irrational() : type(Type::COMPLEX), coefficient(0), radicand(1), constant_term(0) {}

    /**
     * @brief 构造平方根形式无理数 coeff * sqrtn
     * @param n 被开方数
     * @param coeff 系数,默认为 1.0
     * @return 化简后的无理数对象
     */
    static Irrational sqrt(long long n, lmmc_real_t coeff = 1.0) {
        Irrational result;
        result.type = Type::SQRT;

        auto [perfect, remainder] = simplify_sqrt(n);
        result.coefficient = coeff * perfect;
        result.radicand = remainder;
        result.constant_term = 0;

        return result;
    }

    /**
     * @brief 构造 pi 的倍数
     * @param coeff 系数,默认为 1.0
     * @return coeff * pi
     */
    static Irrational pi(lmmc_real_t coeff = 1.0) {
        Irrational result;
        result.type = Type::PI;
        result.coefficient = coeff;
        result.radicand = 1;
        result.constant_term = 0;
        return result;
    }

    /**
     * @brief 构造 e 的倍数
     * @param coeff 系数,默认为 1.0
     * @return coeff * e
     */
    static Irrational e(lmmc_real_t coeff = 1.0) {
        Irrational result;
        result.type = Type::E;
        result.coefficient = coeff;
        result.radicand = 1;
        result.constant_term = 0;
        return result;
    }

    /**
     * @brief 构造有理常数(退化为有理数的无理数表示)
     * @param value 常数值
     * @return 常数无理数对象
     */
    static Irrational constant(lmmc_real_t value) {
        Irrational result;
        result.type = Type::COMPLEX;
        result.coefficient = 0;
        result.radicand = 1;
        result.constant_term = value;
        return result;
    }

    /** @brief 将当前无理数转换为 COMPLEX 线性组合形式 */
    void to_complex() {
        if (type == Type::COMPLEX) return;

        coefficients.clear();
        constant_term = 0;

        switch (type) {
            case Type::SQRT:
                if (radicand == 1) {
                    constant_term = coefficient;
                } else {
                    coefficients["sqrt" + std::to_string(radicand)] = coefficient;
                }
                break;
            case Type::PI:
                coefficients["pi"] = coefficient;
                break;
            case Type::E:
                coefficients["e"] = coefficient;
                break;
            default:
                break;
        }
        type = Type::COMPLEX;
    }

    Irrational operator+(const Irrational& other) const {
        Irrational result = *this;
        Irrational other_copy = other;

        result.to_complex();
        other_copy.to_complex();

        result.constant_term += other_copy.constant_term;

        for (const auto& [key, coeff]: other_copy.coefficients) {
            result.coefficients[key] += coeff;
        }

        return result;
    }

    Irrational operator-(const Irrational& other) const {
        Irrational result = *this;
        Irrational other_copy = other;

        result.to_complex();
        other_copy.to_complex();

        result.constant_term -= other_copy.constant_term;

        for (const auto& [key, coeff]: other_copy.coefficients) {
            result.coefficients[key] -= coeff;
        }

        return result;
    }

    Irrational operator*(lmmc_real_t scalar) const {
        Irrational result = *this;

        if (type == Type::COMPLEX) {
            result.constant_term *= scalar;
            for (auto& [key, coeff]: result.coefficients) {
                coeff *= scalar;
            }
        } else {
            result.coefficient *= scalar;
        }

        return result;
    }

    Irrational operator*(const Irrational& other) const {

        if (type == Type::COMPLEX && coefficients.empty()) {
            return other * constant_term;
        }
        if (other.type == Type::COMPLEX && other.coefficients.empty()) {
            return *this * other.constant_term;
        }

        if (type == Type::SQRT && other.type == Type::SQRT) {
            return Irrational::sqrt(radicand * other.radicand,
                                    coefficient * other.coefficient);
        }

        return Irrational::constant(to_double() * other.to_double());
    }

    Irrational operator/(const Irrational& other) const {

        if (other.type == Type::COMPLEX && other.coefficients.empty() && other.constant_term != 0) {
            return *this * (1.0 / other.constant_term);
        }

        lmmc_real_t other_val = other.to_double();
        lmmc_real_t abs_other;
        LMMC_REAL_ABS(&abs_other, &other_val);
        if (abs_other < 1e-15) {
            throw std::runtime_error("Irrational: division by zero");
        }
        return Irrational::constant(to_double() / other_val);
    }

    Irrational operator-() const {
        return *this * (-1.0);
    }

    bool operator==(const Irrational& other) const {
        int eq;
        lmmc_double_nearly_equal(to_double(), other.to_double(), &eq);
        return eq != 0;
    }

    bool operator<(const Irrational& other) const {
        return to_double() < other.to_double();
    }

    bool operator<=(const Irrational& other) const {
        return *this < other || *this == other;
    }

    bool operator>(const Irrational& other) const {
        return other < *this;
    }

    bool operator>=(const Irrational& other) const {
        return *this > other || *this == other;
    }

    /**
     * @brief 转换为浮点数近似值
     * @return 双精度浮点近似
     */
    lmmc_real_t to_double() const {
        switch (type) {
            case Type::SQRT:
                if (radicand == 1) {
                    return coefficient;
                }
                {
                    lmmc_real_t rad = radicand;
                    return coefficient * _irrational_sqrt(rad);
                }
            case Type::PI:
                return coefficient * LMMC_PI;
            case Type::E:
                return coefficient * LMMC_E;
            case Type::LOG:
                {
                    lmmc_real_t rad = radicand;
                    lmmc_real_t res;
                    LMMC_REAL_LOG(&res, &rad);
                    return coefficient * res;
                }
            case Type::COMPLEX: {
                lmmc_real_t result = constant_term;
                for (const auto& [key, coeff]: coefficients) {
                    if (key == "pi") {
                        result += coeff * LMMC_PI;
                    } else if (key == "e") {
                        result += coeff * LMMC_E;
                    } else if (key.substr(0, 4) == "sqrt") {
                        long long n = std::stoll(key.substr(4));
                        lmmc_real_t n_real = n;
                        result += coeff * _irrational_sqrt(n_real);
                    }
                }
                return result;
            }
            default:
                return 0.0;
        }
    }

    /**
     * @brief 转换为可读字符串(如 "2sqrt3","pi/2")
     * @return 格式化字符串
     */
    std::string to_string() const {
        auto round_val = [](lmmc_real_t x) -> lmmc_real_t {
            lmmc_real_t res, half = 0.5;
            if (x < 0) half = -0.5;
            LMMC_REAL_ADD(&res, &x, &half);
            long long iptr = (long long)res;
            return (lmmc_real_t)iptr;
        };
        auto abs_val = [](lmmc_real_t x) -> lmmc_real_t {
            lmmc_real_t res;
            LMMC_REAL_ABS(&res, &x);
            return res;
        };
        auto is_zero_tol = [](lmmc_real_t x) -> bool {
            int eq;
            lmmc_double_nearly_equal_tol(x, 0.0, 1e-15, 1e-15, &eq);
            return eq != 0;
        };
        auto is_equal_tol = [](lmmc_real_t a, lmmc_real_t b) -> bool {
            int eq;
            lmmc_double_nearly_equal_tol(a, b, 1e-15, 1e-15, &eq);
            return eq != 0;
        };
        switch (type) {
            case Type::SQRT:
                if (radicand == 1) {

                    if (is_equal_tol(coefficient, round_val(coefficient))) {
                        return std::to_string(static_cast<int>(round_val(coefficient)));
                    }
                    return std::to_string(coefficient);
                }
                if (is_equal_tol(coefficient, 1.0)) {
                    return "√" + std::to_string(radicand);
                }
                if (is_equal_tol(coefficient, -1.0)) {
                    return "-√" + std::to_string(radicand);
                }

                if (is_equal_tol(coefficient, round_val(coefficient))) {
                    return std::to_string(static_cast<int>(round_val(coefficient))) + "√" + std::to_string(radicand);
                } else {

                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(6) << coefficient;
                    std::string temp = oss.str();
                    temp.erase(temp.find_last_not_of('0') + 1);
                    if (temp.back() == '.') temp.pop_back();
                    return temp + "√" + std::to_string(radicand);
                }

            case Type::PI:
                if (is_equal_tol(coefficient, 1.0)) {
                    return "π";
                }
                if (is_equal_tol(coefficient, -1.0)) {
                    return "-π";
                }

                if (is_equal_tol(coefficient, round_val(coefficient))) {
                    return std::to_string(static_cast<int>(round_val(coefficient))) + "π";
                } else {

                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(6) << coefficient;
                    std::string temp = oss.str();
                    temp.erase(temp.find_last_not_of('0') + 1);
                    if (temp.back() == '.') temp.pop_back();
                    return temp + "π";
                }

            case Type::E:
                if (is_equal_tol(coefficient, 1.0)) {
                    return "e";
                }
                if (is_equal_tol(coefficient, -1.0)) {
                    return "-e";
                }

                if (is_equal_tol(coefficient, round_val(coefficient))) {
                    return std::to_string(static_cast<int>(round_val(coefficient))) + "e";
                } else {

                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(6) << coefficient;
                    std::string temp = oss.str();
                    temp.erase(temp.find_last_not_of('0') + 1);
                    if (temp.back() == '.') temp.pop_back();
                    return temp + "e";
                }

            case Type::LOG:
                if (is_equal_tol(coefficient, 1.0)) {
                    return "log(" + std::to_string(radicand) + ")";
                }
                if (is_equal_tol(coefficient, -1.0)) {
                    return "-log(" + std::to_string(radicand) + ")";
                }

                if (is_equal_tol(coefficient, round_val(coefficient))) {
                    return std::to_string(static_cast<int>(round_val(coefficient))) + "log(" + std::to_string(radicand) + ")";
                } else {

                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(6) << coefficient;
                    std::string temp = oss.str();
                    temp.erase(temp.find_last_not_of('0') + 1);
                    if (temp.back() == '.') temp.pop_back();
                    return temp + "log(" + std::to_string(radicand) + ")";
                }

            case Type::COMPLEX: {
                std::string result;
                bool first = true;

                if (!is_zero_tol(constant_term)) {
                    if (is_equal_tol(constant_term, round_val(constant_term))) {

                        result += std::to_string(static_cast<int>(round_val(constant_term)));
                    } else {

                        std::ostringstream oss;
                        oss << std::fixed << std::setprecision(6) << constant_term;
                        std::string temp = oss.str();

                        temp.erase(temp.find_last_not_of('0') + 1);
                        if (temp.back() == '.') temp.pop_back();
                        result += temp;
                    }
                    first = false;
                }

                for (const auto& [key, coeff]: coefficients) {
                    if (is_zero_tol(coeff)) continue;

                    if (!first && coeff > 0) result += " + ";
                    else if (!first && coeff < 0)
                        result += " - ";

                    lmmc_real_t abs_coeff = abs_val(coeff);
                    std::string term;

                    if (key == "pi") {
                        if (is_equal_tol(abs_coeff, 1.0)) {
                            term = "π";
                        } else if (is_equal_tol(abs_coeff, round_val(abs_coeff))) {
                            term = std::to_string(static_cast<int>(round_val(abs_coeff))) + "π";
                        } else {
                            term = std::to_string(abs_coeff) + "π";
                        }
                    } else if (key == "e") {
                        if (is_equal_tol(abs_coeff, 1.0)) {
                            term = "e";
                        } else if (is_equal_tol(abs_coeff, round_val(abs_coeff))) {
                            term = std::to_string(static_cast<int>(round_val(abs_coeff))) + "e";
                        } else {
                            term = std::to_string(abs_coeff) + "e";
                        }
                    } else if (key.substr(0, 4) == "sqrt") {
                        long long n = std::stoll(key.substr(4));
                        if (is_equal_tol(abs_coeff, 1.0)) {
                            term = "√" + std::to_string(n);
                        } else if (is_equal_tol(abs_coeff, round_val(abs_coeff))) {
                            term = std::to_string(static_cast<int>(round_val(abs_coeff))) + "√" + std::to_string(n);
                        } else {
                            term = std::to_string(abs_coeff) + "√" + std::to_string(n);
                        }
                    } else {
                        /// 将未知基向量按原始 key 输出,与 to_symbolic 的 COMPLEX 分支保持一致.
                        if (is_equal_tol(abs_coeff, 1.0)) {
                            term = key;
                        } else if (is_equal_tol(abs_coeff, round_val(abs_coeff))) {
                            term = std::to_string(static_cast<int>(round_val(abs_coeff))) + key;
                        } else {
                            term = std::to_string(abs_coeff) + key;
                        }
                    }

                    if (first && coeff < 0) result += "-";
                    result += term;
                    first = false;
                }

                return result.empty() ? "0" : result;
            }
            default:
                return "0";
        }
    }

    /**
     * @brief 判断是否为零
     * @return 若数值近似为零则返回 true
     */
    bool is_zero() const {
        int eq;
        lmmc_double_nearly_equal_tol(to_double(), 0.0, 1e-15, 1e-15, &eq);
        return eq != 0;
    }

    /**
     * @brief 判断是否可精确表示为有理数
     * @return 若为纯常数项(无无理部分)则返回 true
     */
    bool is_rational() const {
        if (type == Type::COMPLEX) {
            return coefficients.empty();
        }
        return false;
    }

    /** @brief 化简,移除系数为零的项 */
    void simplify() {
        if (type == Type::COMPLEX) {
            auto it = coefficients.begin();
            while (it != coefficients.end()) {
                int eq;
                lmmc_double_nearly_equal_tol(it->second, 0.0, 1e-15, 1e-15, &eq);
                if (eq) {
                    it = coefficients.erase(it);
                } else {
                    ++it;
                }
            }

            int eq_const;
            lmmc_double_nearly_equal_tol(constant_term, 0.0, 1e-15, 1e-15, &eq_const);
            if (coefficients.empty() && eq_const) {
                constant_term = 0.0;
            }
        }
    }

    /**
     * @brief 判断是否为正数
     * @return 若数值大于零则返回 true
     */
    bool is_positive() const {
        int eq;
        lmmc_double_nearly_equal_tol(to_double(), 0.0, 1e-15, 1e-15, &eq);
        return !eq && to_double() > 0;
    }

    /**
     * @brief 判断是否为负数
     * @return 若数值小于零则返回 true
     */
    bool is_negative() const {
        int eq;
        lmmc_double_nearly_equal_tol(to_double(), 0.0, 1e-15, 1e-15, &eq);
        return !eq && to_double() < 0;
    }

    /**
     * @brief 取绝对值
     * @return 绝对值无理数对象
     */
    Irrational abs() const {
        if (is_negative()) {
            return -*this;
        }
        return *this;
    }

    /**
     * @brief 计算整数次幂
     * @param exponent 指数
     * @return 幂运算结果
     */
    Irrational pow(int exponent) const {
        if (exponent == 0) {
            return Irrational::constant(1.0);
        }
        if (exponent == 1) {
            return *this;
        }
        if (exponent == 2 && type == Type::SQRT) {

            return Irrational::constant(coefficient * coefficient * radicand);
        }

        lmmc_real_t res = 1.0;
        lmmc_real_t base = to_double();
        int e = exponent > 0 ? exponent : -exponent;
        while (e > 0) {
            if (e % 2 == 1) LMMC_REAL_MUL(&res, &res, &base);
            LMMC_REAL_MUL(&base, &base, &base);
            e /= 2;
        }
        if (exponent < 0) {
            lmmc_real_t one = 1.0;
            LMMC_REAL_DIV(&res, &one, &res);
        }
        return Irrational::constant(res);
    }

    friend std::ostream& operator<<(std::ostream& os, const Irrational& ir) {
        os << ir.to_string();
        return os;
    }

    /**
     * @brief 获取无理数的内部类型
     * @return 类型枚举值
     */
    Type get_type() const { return type; }
};
