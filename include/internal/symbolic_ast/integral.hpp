/** @file internal/symbolic_ast/integral.hpp */
#pragma once

#include "relations.hpp"

class IntegralNode : public SymbolicNode {
private:
    template <typename Node, typename... Args>
    friend std::shared_ptr<const Node>
        lamina::detail::make_node(Args&&... args);
    const std::shared_ptr<const SymbolicNode> body_;
    const std::string variable_;
    const std::shared_ptr<const SymbolicNode> lower_;
    const std::shared_ptr<const SymbolicNode> upper_;

    IntegralNode(std::shared_ptr<const SymbolicNode> body,
                 std::string variable,
                 std::shared_ptr<const SymbolicNode> lower = nullptr,
                 std::shared_ptr<const SymbolicNode> upper = nullptr)
        : body_(std::move(body)), variable_(std::move(variable)),
          lower_(std::move(lower)), upper_(std::move(upper)) {
        if (!body_) throw std::invalid_argument("integral body cannot be null");
        if (variable_.empty()) throw std::invalid_argument("integral variable cannot be empty");
        if (static_cast<bool>(lower_) != static_cast<bool>(upper_)) {
            throw std::invalid_argument("integral bounds must both be present or absent");
        }
    }

public:
    const auto& body() const noexcept { return body_; }
    const std::string& variable() const noexcept { return variable_; }
    const auto& lower() const noexcept { return lower_; }
    const auto& upper() const noexcept { return upper_; }
    bool is_definite() const noexcept { return lower_ && upper_; }
    int type_priority() const override { return 108; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        hash_combine(seed, body_->hash());
        hash_combine(seed, std::hash<std::string>{}(variable_));
        if (lower_) hash_combine(seed, lower_->hash());
        if (upper_) hash_combine(seed, upper_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& integral = static_cast<const IntegralNode&>(other);
        int comparison = variable_.compare(integral.variable_);
        if (comparison != 0) return comparison;
        comparison = body_->compare(*integral.body_);
        if (comparison != 0) return comparison;
        if (is_definite() != integral.is_definite()) return is_definite() ? 1 : -1;
        if (!is_definite()) return 0;
        comparison = lower_->compare(*integral.lower_);
        return comparison != 0 ? comparison : upper_->compare(*integral.upper_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override {
        lamina::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<IntegralNode>(
            body_->clone(), variable_,
            lower_ ? lower_->clone() : nullptr,
            upper_ ? upper_->clone() : nullptr);
    }
};
