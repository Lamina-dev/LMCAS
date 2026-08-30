/**
 * @file internal/symbolic_ast/base.hpp
 * @brief Internal AST ownership, node protocol, visitor contract, and factory declarations.
 */
#pragma once
#define LAMINA_INTERNAL_AST_INCLUDED 1
#include <memory>
#include <vector>
#include <string>
#include <variant>
#include <algorithm>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include "lamina_export.hpp"
#include "symbolic.hpp"
#include "rational.hpp"
#include "bigint.hpp"
#include "unit.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <stdexcept>
#include <atomic>
#include <type_traits>
#include <utility>

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

namespace lamina::detail {

using SymbolicNodePtr = std::shared_ptr<const SymbolicNode>;

struct SymbolicExprAccess {
    static const SymbolicNodePtr& node(const ::SymbolicExpr& expression) noexcept {
        return expression.impl_->root;
    }

    static ::SymbolicExpr expression_from_node(SymbolicNodePtr root) {
        return ::SymbolicExpr(
            std::make_shared<const ::SymbolicExpr::Impl>(std::move(root)));
    }

    static std::shared_ptr<::SymbolicExpr> make_expression_ptr(SymbolicNodePtr root) {
        return std::shared_ptr<::SymbolicExpr>(new ::SymbolicExpr(
            std::make_shared<const ::SymbolicExpr::Impl>(std::move(root))));
    }
};

inline const SymbolicNodePtr& node(const ::SymbolicExpr& expression) noexcept {
    return SymbolicExprAccess::node(expression);
}

inline const SymbolicNodePtr& node(
    const std::shared_ptr<::SymbolicExpr>& expression) noexcept {
    static const SymbolicNodePtr empty;
    return expression ? node(*expression) : empty;
}

inline const SymbolicNodePtr& node(
    const std::shared_ptr<const ::SymbolicExpr>& expression) noexcept {
    static const SymbolicNodePtr empty;
    return expression ? node(*expression) : empty;
}

inline ::SymbolicExpr expression_from_node(SymbolicNodePtr root) {
    return SymbolicExprAccess::expression_from_node(std::move(root));
}

inline std::shared_ptr<::SymbolicExpr> make_expression_ptr(SymbolicNodePtr root) {
    return SymbolicExprAccess::make_expression_ptr(std::move(root));
}

inline std::shared_ptr<::SymbolicExpr> make_expression_ptr(
    const ::SymbolicExpr& expression) {
    return std::make_shared<::SymbolicExpr>(expression);
}

class SymbolicVisitor;

template <typename Node, typename... Args>
std::shared_ptr<const Node> make_node(Args&&... args) {
    static_assert(std::is_base_of<SymbolicNode, Node>::value,
                  "make_node only constructs SymbolicNode implementations");
    return std::shared_ptr<const Node>(
        new Node(std::forward<Args>(args)...));
}

} // namespace lamina::detail

#define LAMINA_AST_NODE_FACTORY_FRIEND                                      \
    template <typename Node, typename... Args>                              \
    friend std::shared_ptr<const Node>                                      \
        lamina::detail::make_node(Args&&... args)

/**
 * @brief 将一个哈希值混合到种子中，用于组合多个字段的哈希。
 * @param seed 当前哈希种子，混合后就地更新
 * @param value 要混合的哈希值
 */
inline void hash_combine(std::size_t& seed, std::size_t value) {
    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * @brief 符号表达式 AST 节点基类。
 *
 * 所有具体节点（数值、变量、运算符、函数等）均继承此类。
 * 提供哈希、比较、克隆等通用接口，支持 Visitor 模式遍历。
 */
class SymbolicNode {
protected:
    mutable std::atomic<std::size_t> cached_hash{0};
    mutable std::atomic<bool> hash_computed{false};

    SymbolicNode() = default;
    SymbolicNode(const SymbolicNode&) : cached_hash(0), hash_computed(false) {}
    SymbolicNode& operator=(const SymbolicNode&) { return *this; }

    /** @brief 计算当前节点的哈希值（由子类实现）。 */
    virtual std::size_t compute_hash() const = 0;

    /** @brief 与同类型节点进行比较（由子类实现）。 */
    virtual int compare_same_type(const SymbolicNode& other) const = 0;

public:
    virtual ~SymbolicNode() = default;

    /** @brief 接受 Visitor 访问。 */
    virtual void accept(lamina::detail::SymbolicVisitor& visitor) const = 0;

    /** @brief 深拷贝当前节点及其子树。 */
    virtual std::shared_ptr<const SymbolicNode> clone() const = 0;

    /** @brief 返回节点类型的排序优先级，用于规范化排序。 */
    virtual int type_priority() const = 0;

    /**
     * @brief 获取节点哈希值（带缓存）。
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
     * @brief 与另一个节点进行全序比较。
     * @param other 待比较的节点
     * @return 小于返回 -1，等于返回 0，大于返回 1
     */
    int compare(const SymbolicNode& other) const {
        if (type_priority() != other.type_priority()) {
            return type_priority() < other.type_priority() ? -1 : 1;
        }
        return compare_same_type(other);
    }

    /**
     * @brief 判断两个节点是否结构相等。
     * @param other 待比较的节点
     * @return 相等返回 true
     */
    bool equals(const SymbolicNode& other) const {
        if (this == &other) return true;
        if (hash() != other.hash()) return false;
        if (type_priority() != other.type_priority()) return false;
        return compare_same_type(other) == 0;
    }

    /** @brief 判断节点是否为数值节点。 */
    virtual bool is_number() const { return false; }
    /** @brief 判断节点数值是否为 1。 */
    virtual bool is_one() const { return false; }
    /** @brief 判断节点数值是否为 0。 */
    virtual bool is_zero() const { return false; }
    /** @brief 判断节点是否可证明为严格正数（非数值节点默认返回 false）。 */
    virtual bool is_positive() const { return false; }
};

/**
 * @brief AST 访问者基类（Visitor 模式）。
 *
 * 子类实现各 visit 方法以处理不同节点类型。
 * 深度计数达到资源上限时返回诊断，使递归保持在配置范围内。
 */
namespace lamina::detail {

class SymbolicVisitor {
protected:
    int current_depth = 0;
    static constexpr int MAX_DEPTH = 200;
public:
    /** @brief 深度守卫，进入时递增深度，超限时抛出异常。 */
    struct DepthGuard {
        SymbolicVisitor& visitor;
        DepthGuard(SymbolicVisitor& v) : visitor(v) {
            if (++visitor.current_depth > MAX_DEPTH) {
                throw std::runtime_error("AST traversal depth limit exceeded");
            }
        }
        ~DepthGuard() {
            --visitor.current_depth;
        }
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

} // namespace lamina::detail

/** @brief 基于节点哈希的哈希函数对象，用于无序容器。 */
struct NodeHash {
    std::size_t operator()(const std::shared_ptr<const SymbolicNode>& node) const {
        return node ? node->hash() : 0;
    }
};

/** @brief 基于节点结构相等性的比较函数对象，用于无序容器。 */
struct NodeEqual {
    bool operator()(const std::shared_ptr<const SymbolicNode>& lhs, const std::shared_ptr<const SymbolicNode>& rhs) const {
        if (!lhs || !rhs) return lhs == rhs;
        return lhs->equals(*rhs);
    }
};

template<typename T>
using NodeMap = std::unordered_map<std::shared_ptr<const SymbolicNode>, T, NodeHash, NodeEqual>;

using NodeSet = std::unordered_set<std::shared_ptr<const SymbolicNode>, NodeHash, NodeEqual>;

/**
 * @brief 符号节点工厂类，提供创建常用节点的静态方法。
 *
 * 创建时自动执行基本简化（如零元素消除、单元素折叠、扁平化）。
 */
class SymbolicFactory {
public:
    /** @brief 创建大整数数值节点。 */
    static std::shared_ptr<const SymbolicNode> create_number(const ::BigInt& v);
    /** @brief 创建有理数数值节点。 */
    static std::shared_ptr<const SymbolicNode> create_number(const ::Rational& v);
    /** @brief 创建浮点数值节点。 */
    static std::shared_ptr<const SymbolicNode> create_number(lmmc_real_t v);
    /** @brief 创建变量节点。 */
    static std::shared_ptr<const SymbolicNode> create_variable(const std::string& name);

    /**
     * @brief 创建加法节点，自动扁平化嵌套加法并消除零项。
     * @param ops 操作数列表
     * @return 简化后的节点
     */
    static std::shared_ptr<const SymbolicNode> create_add(std::vector<std::shared_ptr<const SymbolicNode>> ops);

    /**
     * @brief 创建乘法节点，自动扁平化嵌套乘法并消除单位元。
     * @param ops 操作数列表
     * @return 简化后的节点
     */
    static std::shared_ptr<const SymbolicNode> create_multiply(std::vector<std::shared_ptr<const SymbolicNode>> ops);

    /**
     * @brief 创建幂运算节点，自动处理指数为 0/1 及底数为 0/1 的情况。
     * @param base 底数节点
     * @param exponent 指数节点
     * @return 简化后的节点
     */
    static std::shared_ptr<const SymbolicNode> create_power(std::shared_ptr<const SymbolicNode> base, std::shared_ptr<const SymbolicNode> exponent);

    /**
     * @brief 创建复数节点，自动简化（若虚部为 0，则返回实部）。
     */
    static std::shared_ptr<const SymbolicNode> create_complex(std::shared_ptr<const SymbolicNode> real, std::shared_ptr<const SymbolicNode> imag);
};
