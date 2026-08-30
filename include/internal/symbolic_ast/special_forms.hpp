/** @file internal/symbolic_ast/special_forms.hpp */
#pragma once

#include "integral.hpp"
#include "../../limit_direction.hpp"
#include "../exact_root_id.hpp"

#include <cstddef>


/** Unevaluated limit syntax. The variable is bound only in the body. */
class LimitNode : public SymbolicNode {
private:
    template <typename Node, typename... Args>
    friend std::shared_ptr<const Node>
        lamina::detail::make_node(Args&&... args);

    const std::shared_ptr<const SymbolicNode> body_;
    const std::string variable_;
    const std::shared_ptr<const SymbolicNode> point_;
    const LimitDirection direction_;

    LimitNode(std::shared_ptr<const SymbolicNode> body,
              std::string variable,
              std::shared_ptr<const SymbolicNode> point,
              LimitDirection direction)
        : body_(std::move(body)), variable_(std::move(variable)),
          point_(std::move(point)), direction_(direction) {
        if (!body_ || !point_) {
            throw std::invalid_argument("limit body and point cannot be null");
        }
        if (variable_.empty()) {
            throw std::invalid_argument("limit variable cannot be empty");
        }
    }

public:
    const auto& body() const noexcept { return body_; }
    const std::string& variable() const noexcept { return variable_; }
    const auto& point() const noexcept { return point_; }
    LimitDirection direction() const noexcept { return direction_; }
    int type_priority() const override { return 109; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        hash_combine(seed, body_->hash());
        hash_combine(seed, std::hash<std::string>{}(variable_));
        hash_combine(seed, point_->hash());
        hash_combine(seed, static_cast<std::size_t>(direction_));
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& limit = static_cast<const LimitNode&>(other);
        int comparison = variable_.compare(limit.variable_);
        if (comparison != 0) return comparison;
        if (direction_ != limit.direction_) {
            return static_cast<int>(direction_) < static_cast<int>(limit.direction_) ? -1 : 1;
        }
        comparison = point_->compare(*limit.point_);
        return comparison != 0 ? comparison : body_->compare(*limit.body_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override {
        lamina::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<LimitNode>(
            body_->clone(), variable_, point_->clone(), direction_);
    }
};

/** Canonical exact algebraic root identity. */
class RootOfNode : public SymbolicNode {
private:
    template <typename Node, typename... Args>
    friend std::shared_ptr<const Node>
        lamina::detail::make_node(Args&&... args);

    const lamina::detail::ExactRootId id_;
    const std::string display_variable_;

    explicit RootOfNode(
        lamina::detail::ExactRootId id,
        std::string display_variable = "_root")
        : id_(std::move(id)),
          display_variable_(std::move(display_variable)) {
        if (id_.polynomial.degree() <= 0 ||
            id_.index >= static_cast<std::size_t>(id_.polynomial.degree()) ||
            display_variable_.empty()) {
            throw std::invalid_argument("RootOf identity is invalid");
        }
    }

public:
    const lamina::detail::ExactRootId& exact_id() const noexcept { return id_; }
    std::size_t index() const noexcept { return id_.index; }
    const std::string& variable() const noexcept {
        return display_variable_;
    }
    std::shared_ptr<const SymbolicNode> polynomial() const {
        const auto variable_node =
            SymbolicFactory::create_variable(variable());
        std::vector<std::shared_ptr<const SymbolicNode>> terms;
        for (std::size_t degree = 0;
             degree < id_.polynomial.coeffs.size(); ++degree) {
            const auto& coefficient = id_.polynomial.coeffs[degree];
            if (coefficient == Rational(0)) continue;
            auto coefficient_node =
                SymbolicFactory::create_number(coefficient);
            if (degree == 0) {
                terms.push_back(coefficient_node);
                continue;
            }
            auto variable_power = degree == 1
                ? variable_node
                : SymbolicFactory::create_power(
                      variable_node,
                      SymbolicFactory::create_number(
                          BigInt(static_cast<long long>(degree))));
            terms.push_back(
                coefficient == Rational(1)
                    ? variable_power
                    : SymbolicFactory::create_multiply(
                          {coefficient_node, variable_power}));
        }
        return SymbolicFactory::create_add(std::move(terms));
    }
    int type_priority() const override { return 110; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        hash_combine(seed, id_.hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        return id_.compare(static_cast<const RootOfNode&>(other).id_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override {
        lamina::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<RootOfNode>(
            id_, display_variable_);
    }
};
