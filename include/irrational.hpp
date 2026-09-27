/** @file irrational.hpp
 * @brief 精确实根式及根式、pi、e 的有理组合。
 */
#pragma once

#include "assumption.hpp"
#include "computation_context.hpp"
#include "rational.hpp"
#include "lmcas_export.hpp"
#include <map>
#include <memory>
#include <iosfwd>
#include <string>
#include <utility>

namespace LMCAS {

class SymbolicExpr;

/**
 * @brief 保留精确系数和无平方因子的被开方数。
 *
 * 商保留精确分子与分母，仅 to_double() 执行舍入；相等性比较使用精确交叉乘积。
 * 排序及 is_positive()/is_negative() 使用近似判断，不能作为精确证明。
 * sign_checked()/abs_checked() 使用认证包围区间；含 pi/e 的组合在 4096 位时仍跨零，
 * 返回 Inconclusive。纯根式持续细化，直至完成认证或耗尽预算。
 * 无上下文的基本运算与 BigInt 一样不限预算；受检便捷入口及 abs() 使用默认 ComputationContext。
 */
class LMCAS_API Irrational {
public:
    enum class Type { SQRT, PI, E, LOG, COMPLEX };

    Irrational() = default;
    /** @brief 负被开方数返回 DomainError；预算耗尽时返回错误，保持精确性。 */
    static Result<Irrational> sqrt_checked(const BigInt& n, lmmc_real_t coefficient,
                                           ComputationContext& context);
    static Result<Irrational> sqrt_checked(const BigInt& n,
                                           lmmc_real_t coefficient = 1.0);
    static Irrational pi(lmmc_real_t coefficient = 1.0);
    static Irrational e(lmmc_real_t coefficient = 1.0);
    static Irrational constant(lmmc_real_t value);
    static Irrational constant(const Rational& value);

    void to_complex() { type_ = Type::COMPLEX; }
    Type get_type() const { return type_; }
    Irrational operator+(const Irrational& other) const;
    Irrational operator-(const Irrational& other) const;
    Irrational operator*(const Irrational& other) const;
    Irrational operator*(const Rational& scalar) const;
    Irrational operator*(lmmc_real_t scalar) const;
    Irrational operator/(const Irrational& other) const;
    Irrational operator-() const;
    Irrational pow(const BigInt& exponent) const;
    bool operator==(const Irrational& other) const;
    bool operator<(const Irrational& other) const;
    bool operator<=(const Irrational& other) const;
    bool operator>(const Irrational& other) const;
    bool operator>=(const Irrational& other) const;
    bool is_zero() const;
    bool is_rational() const;
    Result<Rational> as_rational_checked() const;
    bool is_positive() const;
    bool is_negative() const;
    /** @brief 认证符号，成功值为 Negative、Zero 或 Positive。 */
    Result<Sign> sign_checked(ComputationContext& context) const;
    Result<Sign> sign_checked() const;
    Result<Irrational> abs_checked(ComputationContext& context) const;
    Result<Irrational> abs_checked() const;
    /**
     * @brief 计算经过认证的绝对值。
     * @throws std::domain_error 对应 DomainError。
     * @throws std::length_error 对应 ResourceLimit。
     * @throws std::logic_error 对应 InternalInvariant。
     * @throws std::runtime_error 对应其余受检失败（包括 Inconclusive、Cancelled），附带原因。
     */
    Irrational abs() const;
    void simplify();
    lmmc_real_t to_double() const;
    std::string to_string() const;
    std::shared_ptr<SymbolicExpr> to_symbolic() const;
    friend LMCAS_API std::ostream& operator<<(std::ostream& os, const Irrational& value);

private:
    /** @brief 用任意精度整数保存数学指数与被开方数。 */
    struct Basis {
        BigInt radicand{1};
        BigInt pi_power{0};
        BigInt e_power{0};
        bool operator<(const Basis& other) const;
        bool operator==(const Basis& other) const;
        bool is_constant() const;
        std::string to_string() const;
        lmmc_real_t to_double() const;
        std::shared_ptr<SymbolicExpr> to_symbolic() const;
    };
    using Terms = std::map<Basis, Rational>;
    Type type_ = Type::COMPLEX;
    Terms numerator_;
    Terms denominator_{{Basis{}, Rational(1)}};

    static Terms multiply_terms(const Terms& left, const Terms& right);
    static void add_term(Terms& terms, const Basis& basis, const Rational& coefficient);
    static std::string terms_to_string(const Terms& terms);
    static lmmc_real_t terms_to_double(const Terms& terms);
    static std::shared_ptr<SymbolicExpr> terms_to_symbolic(const Terms& terms);
    struct EnclosureRound;
    static std::pair<Rational, Rational> terms_bounds(
        const Terms& terms, EnclosureRound& round);
    static Result<std::pair<BigInt, BigInt>> extract_square_factors(
        const BigInt& n, ComputationContext& context);
};

} // namespace LMCAS
