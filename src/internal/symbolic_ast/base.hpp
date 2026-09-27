/**
 * @file internal/symbolic_ast/base.hpp
 * @brief 内部 AST 的所有权、节点协议、访问者契约与工厂声明。
 */
#pragma once
#define LMCAS_INTERNAL_AST_INCLUDED 1
#include <memory>
#include <vector>
#include <string>
#include <variant>
#include <algorithm>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include "lmcas_export.hpp"
#include "symbolic.hpp"
#include "rational.hpp"
#include "bigint.hpp"
#include "unit.hpp"
#include "internal/rewrite_budget.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <stdexcept>
#include <atomic>
#include <type_traits>
#include <utility>

namespace LMCAS {

class SymbolicNode;
class NumberNode;
class VariableNode;
class AddNode;
class MultiplyNode;
class PowerNode;
class FunctionNode;
class UninterpretedFunctionNode;
class MatrixNode;
class RelationalNode;
class LogicalNode;
class TransformNode;
class QuantifierNode;
class SetBuilderNode;
class FiniteSetNode;
class IntervalNode;
class MembershipNode;
class QuantityNode;
class PiecewiseNode;
class SummationNode;
class ProductNode;
class ComplexNode;
class IntegralNode;
class LimitNode;
class RootOfNode;

struct SymbolicExpr::Impl {
    explicit Impl(std::shared_ptr<const SymbolicNode> root_node)
        : root(std::move(root_node)) {
        if (!root) {
            throw std::invalid_argument("SymbolicExpr requires a non-null AST root");
        }
    }

    const std::shared_ptr<const SymbolicNode> root;
};

}


namespace LMCAS::detail {

using SymbolicNodePtr = std::shared_ptr<const SymbolicNode>;

struct SymbolicExprAccess {
    static const SymbolicNodePtr& node(const ::LMCAS::SymbolicExpr& expression) noexcept {
        return expression.impl_->root;
    }

    static ::LMCAS::SymbolicExpr expression_from_node(SymbolicNodePtr root) {
        return ::LMCAS::SymbolicExpr(
            std::make_shared<const ::LMCAS::SymbolicExpr::Impl>(std::move(root)));
    }

    static std::shared_ptr<::LMCAS::SymbolicExpr> make_expression_ptr(SymbolicNodePtr root) {
        return std::shared_ptr<::LMCAS::SymbolicExpr>(new ::LMCAS::SymbolicExpr(
            std::make_shared<const ::LMCAS::SymbolicExpr::Impl>(std::move(root))));
    }
};

inline const SymbolicNodePtr& node(const ::LMCAS::SymbolicExpr& expression) noexcept {
    return SymbolicExprAccess::node(expression);
}

inline const SymbolicNodePtr& node(
    const std::shared_ptr<::LMCAS::SymbolicExpr>& expression) noexcept {
    static const SymbolicNodePtr empty;
    return expression ? node(*expression) : empty;
}

inline const SymbolicNodePtr& node(
    const std::shared_ptr<const ::LMCAS::SymbolicExpr>& expression) noexcept {
    static const SymbolicNodePtr empty;
    return expression ? node(*expression) : empty;
}

inline ::LMCAS::SymbolicExpr expression_from_node(SymbolicNodePtr root) {
    return SymbolicExprAccess::expression_from_node(std::move(root));
}

inline std::shared_ptr<::LMCAS::SymbolicExpr> make_expression_ptr(SymbolicNodePtr root) {
    return SymbolicExprAccess::make_expression_ptr(std::move(root));
}

inline std::shared_ptr<::LMCAS::SymbolicExpr> make_expression_ptr(
    const ::LMCAS::SymbolicExpr& expression) {
    return std::make_shared<::LMCAS::SymbolicExpr>(expression);
}

class SymbolicVisitor;

template <typename Node, typename... Args>
std::shared_ptr<const Node> make_node(Args&&... args) {
    static_assert(std::is_base_of<SymbolicNode, Node>::value,
                  "make_node only constructs SymbolicNode implementations");
    return std::shared_ptr<const Node>(
        new Node(std::forward<Args>(args)...));
}

}

#define LMCAS_AST_NODE_FACTORY_FRIEND                                      \
    template <typename Node, typename... Args>                              \
    friend std::shared_ptr<const Node>                                      \
        LMCAS::detail::make_node(Args&&... args)

namespace LMCAS {

/**
 * @brief 将一个哈希值混合到种子中,用于组合多个字段的哈希.
 * @param seed 当前哈希种子,混合后就地更新
 * @param value 要混合的哈希值
 */
inline void hash_combine(std::size_t& seed, std::size_t value) {
    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
class SymbolicNode {
protected:
    mutable std::atomic<std::size_t> cached_hash{0};
    mutable std::atomic<bool> hash_computed{false};

    SymbolicNode() = default;
    SymbolicNode(const SymbolicNode&) : cached_hash(0), hash_computed(false) {}
    SymbolicNode& operator=(const SymbolicNode&) { return *this; }
    virtual std::size_t compute_hash() const = 0;
    virtual int compare_same_type(const SymbolicNode& other) const = 0;

public:
    virtual ~SymbolicNode() = default;
    virtual void accept(LMCAS::detail::SymbolicVisitor& visitor) const = 0;
    virtual std::shared_ptr<const SymbolicNode> clone() const = 0;
    virtual int type_priority() const = 0;

    /**
     * @brief 获取节点哈希值(带缓存).
     * @return 哈希值
     */
    std::size_t hash() const {
        if (!hash_computed.load(std::memory_order_acquire)) {
            std::size_t h = compute_hash();
            cached_hash.store(h, std::memory_order_relaxed);
            hash_computed.store(true, std::memory_order_release);
            return h;
        }
        return cached_hash.load(std::memory_order_relaxed);
    }

    /**
     * @brief 与另一个节点进行全序比较.
     * @param other 待比较的节点
     * @return 小于返回 -1,等于返回 0,大于返回 1
     */
    int compare(const SymbolicNode& other) const {
        if (type_priority() != other.type_priority()) {
            return type_priority() < other.type_priority() ? -1 : 1;
        }
        return compare_same_type(other);
    }

    /**
     * @brief 判断两个节点是否结构相等.
     * @param other 待比较的节点
     * @return 相等返回 true
     */
    bool equals(const SymbolicNode& other) const {
        if (this == &other) return true;
        if (hash() != other.hash()) return false;
        if (type_priority() != other.type_priority()) return false;
        return compare_same_type(other) == 0;
    }
    virtual bool is_number() const { return false; }
    virtual bool is_one() const { return false; }
    virtual bool is_zero() const { return false; }
    virtual bool is_positive() const { return false; }
};

/**
 * @brief AST 访问者基类，由子类实现各节点的 visit 方法。
 * 遍历深度超限时抛出异常。
 */

}

namespace LMCAS::detail {

class SymbolicVisitor {
protected:
    int current_depth = 0;
    static constexpr int MAX_DEPTH = 200;
    RewriteBudget* rewrite_budget() const noexcept { return rewrite_budget_; }

private:
    RewriteBudget* rewrite_budget_ = nullptr;
public:
    SymbolicVisitor() = default;
    explicit SymbolicVisitor(RewriteBudget* budget) noexcept : rewrite_budget_(budget) {}
    struct DepthGuard {
        SymbolicVisitor& visitor;
        DepthGuard(SymbolicVisitor& v) : visitor(v) {
            if (visitor.current_depth >= MAX_DEPTH) {
                throw std::runtime_error("AST traversal depth limit exceeded");
            }
            if (visitor.rewrite_budget_) { visitor.rewrite_budget_->enter(); }
            ++visitor.current_depth;
        }
        ~DepthGuard() {
            --visitor.current_depth;
            if (visitor.rewrite_budget_) { visitor.rewrite_budget_->leave(); }
        }
        DepthGuard(const DepthGuard&) = delete;
        DepthGuard& operator=(const DepthGuard&) = delete;
    };

    virtual ~SymbolicVisitor() = default;

    virtual void visit(const NumberNode& node) = 0;
    virtual void visit(const VariableNode& node) = 0;
    virtual void visit(const AddNode& node) = 0;
    virtual void visit(const MultiplyNode& node) = 0;
    virtual void visit(const PowerNode& node) = 0;
    virtual void visit(const FunctionNode& node) = 0;
    virtual void visit(const UninterpretedFunctionNode& node) = 0;
    virtual void visit(const MatrixNode& node) = 0;
    virtual void visit(const RelationalNode& node) = 0;
    virtual void visit(const LogicalNode& node) = 0;
    virtual void visit(const PiecewiseNode& node) = 0;
    virtual void visit(const SummationNode& node) = 0;
    virtual void visit(const ProductNode& node) = 0;
    virtual void visit(const TransformNode& node) = 0;
    virtual void visit(const QuantifierNode& node) = 0;
    virtual void visit(const SetBuilderNode& node) = 0;
    virtual void visit(const FiniteSetNode& node) = 0;
    virtual void visit(const IntervalNode& node) = 0;
    virtual void visit(const MembershipNode& node) = 0;
    virtual void visit(const QuantityNode& node) = 0;
    virtual void visit(const ComplexNode& node) = 0;
    virtual void visit(const IntegralNode& node) = 0;
    virtual void visit(const LimitNode& node) = 0;
    virtual void visit(const RootOfNode& node) = 0;
};

} namespace LMCAS {
struct NodeHash {
    std::size_t operator()(const std::shared_ptr<const SymbolicNode>& node) const {
        return node ? node->hash() : 0;
    }
};
struct NodeEqual {
    bool operator()(const std::shared_ptr<const SymbolicNode>& lhs, const std::shared_ptr<const SymbolicNode>& rhs) const {
        if (!lhs || !rhs) return lhs == rhs;
        return lhs->equals(*rhs);
    }
};

template<typename T>
using NodeMap = std::unordered_map<std::shared_ptr<const SymbolicNode>, T, NodeHash, NodeEqual>;

using NodeSet = std::unordered_set<std::shared_ptr<const SymbolicNode>, NodeHash, NodeEqual>;
class SymbolicFactory {
public:
    static std::shared_ptr<const SymbolicNode> create_number(const ::LMCAS::BigInt& v);
    static std::shared_ptr<const SymbolicNode> create_number(const ::LMCAS::Rational& v);
    static std::shared_ptr<const SymbolicNode> create_number(lmmc_real_t v);
    static std::shared_ptr<const SymbolicNode> create_variable(const std::string& name);

    /**
     * @brief 创建加法节点,自动扁平化嵌套加法并消除零项.
     * @param ops 操作数列表
     * @return 简化后的节点
     */
    static std::shared_ptr<const SymbolicNode> create_add(
        std::vector<std::shared_ptr<const SymbolicNode>> ops,
        detail::RewriteBudget* budget = nullptr);

    /**
     * @brief 创建乘法节点,自动扁平化嵌套乘法并消除单位元.
     * @param ops 操作数列表
     * @return 简化后的节点
     */
    static std::shared_ptr<const SymbolicNode> create_multiply(
        std::vector<std::shared_ptr<const SymbolicNode>> ops,
        detail::RewriteBudget* budget = nullptr);

    /**
     * @brief 创建幂运算节点,自动处理指数为 0/1 及底数为 0/1 的情况.
     * @param base 底数节点
     * @param exponent 指数节点
     * @return 简化后的节点
     */
    static std::shared_ptr<const SymbolicNode> create_power(std::shared_ptr<const SymbolicNode> base, std::shared_ptr<const SymbolicNode> exponent);

    /**
     * @brief 创建复数节点,自动简化(若虚部为 0,则返回实部).
     */
    static std::shared_ptr<const SymbolicNode> create_complex(std::shared_ptr<const SymbolicNode> real, std::shared_ptr<const SymbolicNode> imag);
};

}
