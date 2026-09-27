#pragma once

#include "computation_context.hpp"
#include "rational.hpp"
#include <algorithm>
#include <limits>
#include <utility>

namespace LMCAS::detail {

using ExactConstantInterval = std::pair<Rational, Rational>;
Result<ExactConstantInterval> exact_atan_bounds(
    const Rational& value, std::size_t precision_bits, ComputationContext& context);
Result<ExactConstantInterval> exact_pi_bounds(
    std::size_t precision_bits, ComputationContext& context);
Result<ExactConstantInterval> exact_e_bounds(
    std::size_t precision_bits, ComputationContext& context);

/**
 * @brief 常数及根式包围区间共用的内部算术策略。
 * 保守界覆盖未约分 Rational 的交叉乘积及结果；受检入口统一捕获 CasError，仅返回完整结果。
 */
class ExactBoundArithmetic {
public:
    ExactBoundArithmetic(ComputationContext& context, const char* operation)
        : context_(context), operation_(operation) {}

    void require(Result<void> result) const {
        if (!result) throw result.error();
    }
    void step() const { require(context_.consume_steps(1, operation_)); }
    void bits(std::size_t count) const {
        require(context_.require_integer_bits(count, operation_));
    }
    std::size_t growth(std::size_t a, std::size_t b) const {
        if (b > std::numeric_limits<std::size_t>::max() - a)
            throw CasError{CasErrc::ResourceLimit, "integer size overflow", operation_};
        const auto size = a + b;
        bits(size);
        return size;
    }
    void operands(const Rational& a, const Rational& b, bool addition = false) const {
        step();
        const auto a_bits = std::max(a.get_numerator().bit_length(),
                                     a.get_denominator().bit_length());
        const auto b_bits = std::max(b.get_numerator().bit_length(),
                                     b.get_denominator().bit_length());
        const auto size = growth(a_bits, b_bits);
        if (addition) growth(size, 1);
    }
    Rational add(const Rational& a, const Rational& b) const {
        operands(a, b, true); return a + b;
    }
    Rational sub(const Rational& a, const Rational& b) const {
        operands(a, b, true); return a - b;
    }
    Rational mul(const Rational& a, const Rational& b) const {
        operands(a, b); return a * b;
    }
    Rational div(const Rational& a, const Rational& b) const {
        operands(a, b);
        if (b.is_zero())
            throw CasError{CasErrc::DomainError, "zero enclosure divisor", operation_};
        return a / b;
    }
    bool less(const Rational& a, const Rational& b) const {
        operands(a, b); return a < b;
    }
    BigInt dyadic_scale(std::size_t precision) const {
        step();
        growth(precision, 1);
        if (precision > std::numeric_limits<mp_size_t>::max())
            throw CasError{CasErrc::ResourceLimit, "shift size overflow", operation_};
        return BigInt(1) << static_cast<mp_size_t>(precision);
    }
private:
    ComputationContext& context_;
    const char* operation_;
};

}
