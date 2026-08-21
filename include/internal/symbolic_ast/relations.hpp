/** @file internal/symbolic_ast/relations.hpp */
#pragma once
#include "functions.hpp"

/**
 * @brief 关系运算节点，表示两个表达式之间的比较关系。
 */
class RelationalNode : public SymbolicNode {
public:
    using Op = lamina::RelationOp;

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> left_;
    const std::shared_ptr<const SymbolicNode> right_;
    const Op op_;

    RelationalNode(std::shared_ptr<const SymbolicNode> l, std::shared_ptr<const SymbolicNode> r, Op o)
        : left_(std::move(l)), right_(std::move(r)), op_(o) {
        if (!left_ || !right_) {
            throw std::invalid_argument("RelationalNode operands cannot be null");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& left() const noexcept { return left_; }
    const std::shared_ptr<const SymbolicNode>& right() const noexcept { return right_; }
    Op op() const noexcept { return op_; }

    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<RelationalNode>(left_->clone(), right_->clone(), op_);
    }

    int type_priority() const override { return 100; }

    std::size_t compute_hash() const override {
        std::size_t h = std::hash<int>{}((int)op_);
        hash_combine(h, left_->hash());
        hash_combine(h, right_->hash());
        return h;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const RelationalNode&>(other);
        if (op_ != o.op_) return (int)op_ < (int)o.op_ ? -1 : 1;
        int cmp = left_->compare(*o.left_);
        if (cmp != 0) return cmp;
        return right_->compare(*o.right_);
    }

    /**
     * @brief 将关系运算符转换为字符串表示。
     * @param op 关系运算符
     * @return 对应的字符串（如 "=", "<", ">=" 等）
     */
    static std::string op_to_string(Op op) {
        switch(op) {
            case Op::EQ: return "=";
            case Op::NEQ: return "!=";
            case Op::LT: return "<";
            case Op::GT: return ">";
            case Op::LEQ: return "<=";
            case Op::GEQ: return ">=";
            default: return "?";
        }
    }
};

/**
 * @brief 逻辑运算节点，表示 And/Or 逻辑组合。
 */
class LogicalNode : public SymbolicNode {
public:
    /** @brief 逻辑运算符类型 */
    enum class Op {
        And,     ///< 逻辑与
        Or,      ///< 逻辑或
        Not,     ///< 逻辑非（一元运算，right 为 nullptr）
        Implies  ///< 逻辑蕴含 A ⇒ B
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> left_;
    const std::shared_ptr<const SymbolicNode> right_;
    const Op op_;

    LogicalNode(std::shared_ptr<const SymbolicNode> l, std::shared_ptr<const SymbolicNode> r, Op o)
        : left_(std::move(l)), right_(std::move(r)), op_(o) {
        if (!left_) {
            throw std::invalid_argument("LogicalNode left operand cannot be null");
        }
        if (op_ != Op::Not && !right_) {
            throw std::invalid_argument("LogicalNode binary right operand cannot be null");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& left() const noexcept { return left_; }
    const std::shared_ptr<const SymbolicNode>& right() const noexcept { return right_; }
    Op op() const noexcept { return op_; }

    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<LogicalNode>(
            left_ ? left_->clone() : nullptr,
            right_ ? right_->clone() : nullptr,
            op_);
    }

    int type_priority() const override { return 101; }

    std::size_t compute_hash() const override {
        std::size_t h = std::hash<int>{}((int)op_);
        if (left_) hash_combine(h, left_->hash());
        if (right_) hash_combine(h, right_->hash());
        return h;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const LogicalNode&>(other);
        if (op_ != o.op_) return (int)op_ < (int)o.op_ ? -1 : 1;
        /// 处理一元 Not 运算（right 为 nullptr）
        bool l_null = !left_;
        bool ol_null = !o.left_;
        if (l_null != ol_null) return l_null ? -1 : 1;
        if (!l_null) {
            int cmp = left_->compare(*o.left_);
            if (cmp != 0) return cmp;
        }
        bool r_null = !right_;
        bool or_null = !o.right_;
        if (r_null != or_null) return r_null ? -1 : 1;
        if (!r_null) {
            return right_->compare(*o.right_);
        }
        return 0;
    }

    /**
     * @brief 将逻辑运算符转换为字符串表示。
     * @param op 逻辑运算符
     * @return "And"、"Or"、"Not" 或 "Implies"
     */
    static std::string op_to_string(Op op) {
        switch(op) {
            case Op::And: return "And";
            case Op::Or: return "Or";
            case Op::Not: return "Not";
            case Op::Implies: return "Implies";
            default: return "?";
        }
    }
};
