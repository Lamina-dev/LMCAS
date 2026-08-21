/** @file internal/symbolic_ast/arithmetic.hpp */
#pragma once
#include "numbers.hpp"

/**
 * @brief 变量节点，表示一个符号变量。
 */
class VariableNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::string name_;

    explicit VariableNode(std::string n) : name_(std::move(n)) {}

public:
    const std::string& name() const noexcept { return name_; }

    int type_priority() const override { return 10; }

protected:
    std::size_t compute_hash() const override {
        return std::hash<std::string>{}(name_);
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const VariableNode&>(other);
        return name_.compare(o.name_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override { return lamina::detail::make_node<VariableNode>(name_); }
};

/**
 * @brief 加法节点，表示多个操作数的求和。
 *
 * 构造时自动扁平化嵌套的 AddNode，并按规范顺序排序操作数。
 */
class AddNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::vector<std::shared_ptr<const SymbolicNode>> operands_;

    explicit AddNode(std::vector<std::shared_ptr<const SymbolicNode>> ops)
        : operands_([&ops]() {
        if (ops.empty()) {
            throw std::invalid_argument("AddNode requires at least one operand");
        }
        std::vector<std::shared_ptr<const SymbolicNode>> flattened;
        for (const auto& op : ops) {
            if (!op) {
                throw std::invalid_argument("AddNode operand cannot be null");
            }
            if (auto add = std::dynamic_pointer_cast<const AddNode>(op)) {
                flattened.insert(flattened.end(), add->operands_.begin(), add->operands_.end());
            } else {
                flattened.push_back(op);
            }
        }
        std::sort(flattened.begin(), flattened.end(), [](const auto& a, const auto& b) {
            bool a_num = a->is_number();
            bool b_num = b->is_number();
            if (a_num != b_num) return !a_num;
            return a->compare(*b) < 0;
        });
        return flattened;
    }()) {}

public:
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands() const noexcept {
        return operands_;
    }

    int type_priority() const override { return 5; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        for (const auto& op : operands_) {
            hash_combine(seed, op->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const AddNode&>(other);
        if (operands_.size() != o.operands_.size()) {
            return operands_.size() < o.operands_.size() ? -1 : 1;
        }
        for (size_t i = 0; i < operands_.size(); ++i) {
            int cmp = operands_[i]->compare(*o.operands_[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(operands_.size());
        for (const auto& op : operands_) new_ops.push_back(op->clone());
        return lamina::detail::make_node<AddNode>(std::move(new_ops));
    }
};

/**
 * @brief 乘法节点，表示多个操作数的乘积。
 *
 * 构造时自动扁平化嵌套的 MultiplyNode，并按规范顺序排序操作数。
 */
class MultiplyNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::vector<std::shared_ptr<const SymbolicNode>> operands_;

    explicit MultiplyNode(std::vector<std::shared_ptr<const SymbolicNode>> ops)
        : operands_([&ops]() {
        if (ops.empty()) {
            throw std::invalid_argument("MultiplyNode requires at least one operand");
        }
        std::vector<std::shared_ptr<const SymbolicNode>> flattened;
        for (const auto& op : ops) {
            if (!op) {
                throw std::invalid_argument("MultiplyNode operand cannot be null");
            }
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(op)) {
                flattened.insert(flattened.end(), mul->operands_.begin(), mul->operands_.end());
            } else {
                flattened.push_back(op);
            }
        }
        std::sort(flattened.begin(), flattened.end(), [](const auto& a, const auto& b) {
            return a->compare(*b) < 0;
        });
        return flattened;
    }()) {}

public:
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands() const noexcept {
        return operands_;
    }

    int type_priority() const override { return 4; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        for (const auto& op : operands_) {
            hash_combine(seed, op->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const MultiplyNode&>(other);
        if (operands_.size() != o.operands_.size()) {
            return operands_.size() < o.operands_.size() ? -1 : 1;
        }
        for (size_t i = 0; i < operands_.size(); ++i) {
            int cmp = operands_[i]->compare(*o.operands_[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(operands_.size());
        for (const auto& op : operands_) new_ops.push_back(op->clone());
        return lamina::detail::make_node<MultiplyNode>(std::move(new_ops));
    }
};

/**
 * @brief 幂运算节点，表示 base^exponent。
 */
class PowerNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> base_;
    const std::shared_ptr<const SymbolicNode> exponent_;

    PowerNode(std::shared_ptr<const SymbolicNode> b, std::shared_ptr<const SymbolicNode> e)
        : base_(std::move(b)), exponent_(std::move(e)) {
        if (!base_ || !exponent_) {
            throw std::invalid_argument("PowerNode base and exponent cannot be null");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& base() const noexcept { return base_; }
    const std::shared_ptr<const SymbolicNode>& exponent() const noexcept {
        return exponent_;
    }

    int type_priority() const override { return 0; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, base_->hash());
        hash_combine(seed, exponent_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const PowerNode&>(other);

        int cmp_exp = exponent_->compare(*o.exponent_);
        if (cmp_exp != 0) return cmp_exp > 0 ? -1 : 1;

        return base_->compare(*o.base_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<PowerNode>(base_->clone(), exponent_->clone());
    }
};
