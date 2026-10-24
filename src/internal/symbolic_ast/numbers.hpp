/** @file internal/symbolic_ast/numbers.hpp */
#pragma once
#include "internal/symbolic_ast/base.hpp"
#include <stdexcept>

namespace LMCAS {

class NumberNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

    const std::variant<BigInt, Rational> value_;

    explicit NumberNode(const BigInt& v) : value_(v) {}
    explicit NumberNode(const Rational& v) : value_(v) {}
    explicit NumberNode(lmmc_real_t v);
    explicit NumberNode(std::variant<BigInt, Rational> v)
        : value_(std::move(v)) {}

public:
    const std::variant<BigInt, Rational>& value() const noexcept {
        return value_;
    }

    int type_priority() const override { return -10; }

protected:
    std::size_t compute_hash() const override {
        if (std::holds_alternative<Rational>(value_)) {
            return std::get<Rational>(value_).hash();
        }
        return Rational(std::get<BigInt>(value_)).hash();
    }

    int compare_same_type(const SymbolicNode& other) const override;

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {

        if (std::holds_alternative<BigInt>(value_)) {
            return LMCAS::detail::make_node<NumberNode>(std::get<BigInt>(value_));
        }
        return LMCAS::detail::make_node<NumberNode>(std::get<Rational>(value_));
    }

    bool is_number() const override { return true; }
    bool is_zero() const override {

        if (std::holds_alternative<BigInt>(value_)) {
            return std::get<BigInt>(value_) == BigInt(0);
        }
        return std::get<Rational>(value_) == Rational(0);
    }
    bool is_one() const override {
        if (std::holds_alternative<BigInt>(value_)) {
            return std::get<BigInt>(value_) == BigInt(1);
        }
        return std::get<Rational>(value_) == Rational(1);
    }
    bool is_negative_one() const;
    bool is_positive() const override {
        if (std::holds_alternative<BigInt>(value_)) {
            const auto& b = std::get<BigInt>(value_);
            return !b.is_negative() && !(b == BigInt(0));
        }
        return std::get<Rational>(value_) > Rational(0);
    }
};

/**
 * @brief 复数节点，表示 real + imag * i。
 */
class ComplexNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> real_;
    const std::shared_ptr<const SymbolicNode> imag_;

    ComplexNode(std::shared_ptr<const SymbolicNode> r, std::shared_ptr<const SymbolicNode> i)
        : real_(std::move(r)), imag_(std::move(i)) {
        if (!real_ || !imag_) {
            throw std::invalid_argument("ComplexNode real and imaginary parts cannot be null");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& real() const noexcept { return real_; }
    const std::shared_ptr<const SymbolicNode>& imag() const noexcept { return imag_; }

    int type_priority() const override { return -5; } /**< 优先级介于 NumberNode 与 VariableNode 之间。 */

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, real_->hash());
        hash_combine(seed, imag_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const ComplexNode&>(other);
        int cmp_r = real_->compare(*o.real_);
        if (cmp_r != 0) {
            return cmp_r;
        }
        return imag_->compare(*o.imag_);
    }

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<ComplexNode>(real_->clone(), imag_->clone());
    }

    bool is_zero() const override {
        return real_->is_zero() && imag_->is_zero();
    }
};

}
