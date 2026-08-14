/**
 * @file symbolic_ast.hpp
 * @brief AST 节点定义：NumberNode, VariableNode, AddNode, MultiplyNode, PowerNode, FunctionNode, MatrixNode, RelationalNode, LogicalNode。
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
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <stdexcept>
#include <atomic>
#include <type_traits>
#include <utility>

class SymbolicNode;
class BooleanNode;
class NumberNode;
class VariableNode;
class AddNode;
class MultiplyNode;
class PowerNode;
class FunctionNode;
class MatrixNode;
class RelationalNode;
class LogicalNode;
class FiniteSetNode;
class IntervalNode;
class MembershipNode;
class UninterpretedFunctionNode;
class TransformNode;
class QuantifierNode;
class SetBuilderNode;
class PiecewiseNode;
class SummationNode;
class ProductNode_Op;
class ComplexNode;

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
 * 子类实现各 visit 方法以对不同节点类型执行操作。
 * 内置深度保护，防止递归过深导致栈溢出。
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

    virtual void visit(const BooleanNode&) {
        throw std::runtime_error("BooleanNode is not supported by this visitor");
    }
    virtual void visit(const NumberNode& node) = 0;
    virtual void visit(const VariableNode& node) = 0;
    virtual void visit(const AddNode& node) = 0;
    virtual void visit(const MultiplyNode& node) = 0;
    virtual void visit(const PowerNode& node) = 0;
    virtual void visit(const FunctionNode& node) = 0;
    virtual void visit(const MatrixNode& node) = 0;
    virtual void visit(const RelationalNode& node) = 0;
    virtual void visit(const LogicalNode& node) = 0;
    virtual void visit(const PiecewiseNode& node) = 0;
    virtual void visit(const SummationNode& node) = 0;
    virtual void visit(const ProductNode_Op& node) = 0;
    virtual void visit(const TransformNode& node) = 0;
    virtual void visit(const QuantifierNode& node) = 0;
    virtual void visit(const SetBuilderNode& node) = 0;
    virtual void visit(const ComplexNode& node) = 0;
    virtual void visit(const FiniteSetNode&) {
        throw std::runtime_error("FiniteSetNode is not supported by this visitor");
    }
    virtual void visit(const IntervalNode&) {
        throw std::runtime_error("IntervalNode is not supported by this visitor");
    }
    virtual void visit(const MembershipNode&) {
        throw std::runtime_error("MembershipNode is not supported by this visitor");
    }
    virtual void visit(const UninterpretedFunctionNode&) {
        throw std::runtime_error("UninterpretedFunctionNode is not supported by this visitor");
    }
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

class BooleanNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const bool value_;

    explicit BooleanNode(bool value) : value_(value) {}

public:
    bool value() const noexcept { return value_; }

    int type_priority() const override { return -20; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<std::size_t>(value_));
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const BooleanNode&>(other);
        if (value_ == o.value_) return 0;
        return value_ ? 1 : -1;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<BooleanNode>(value_);
    }
};

/**
 * @brief 数值节点，存储 BigInt、Rational 或浮点数。
 */
class NumberNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::variant<BigInt, Rational, lmmc_real_t> value_;

    explicit NumberNode(const BigInt& v) : value_(v) {}
    explicit NumberNode(const Rational& v) : value_(v) {}
    explicit NumberNode(lmmc_real_t v) : value_(v) {}
    explicit NumberNode(std::variant<BigInt, Rational, lmmc_real_t> v) : value_(std::move(v)) {}

public:
    const std::variant<BigInt, Rational, lmmc_real_t>& value() const noexcept {
        return value_;
    }

    int type_priority() const override { return -10; }

protected:
    std::size_t compute_hash() const override {
        if (std::holds_alternative<lmmc_real_t>(value_)) {
            lmmc_real_t real = std::get<lmmc_real_t>(value_);
            if (real == 0.0) real = 0.0; // canonicalize negative zero
            std::size_t seed = std::hash<int>{}(2); // approximate-number domain
            hash_combine(seed, std::hash<lmmc_real_t>{}(real));
            return seed;
        }

        // BigInt n and Rational n/1 are the same exact number. Approximate
        // reals intentionally occupy a different structural domain.
        if (std::holds_alternative<Rational>(value_)) {
            return std::get<Rational>(value_).hash();
        }
        return Rational(std::get<BigInt>(value_)).hash();
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const NumberNode&>(other);

        bool is_l_real = std::holds_alternative<lmmc_real_t>(value_);
        bool is_r_real = std::holds_alternative<lmmc_real_t>(o.value_);

        if (is_l_real) {
            if (!is_r_real) {
                const Rational lhs = Rational::from_double(std::get<lmmc_real_t>(value_));
                const Rational rhs = std::holds_alternative<Rational>(o.value_)
                    ? std::get<Rational>(o.value_)
                    : Rational(std::get<BigInt>(o.value_));
                if (lhs < rhs) return -1;
                if (rhs < lhs) return 1;
                return 1; // same numeric value, but approximate numbers remain structurally distinct
            }
            const lmmc_real_t lhs = std::get<lmmc_real_t>(value_);
            const lmmc_real_t rhs = std::get<lmmc_real_t>(o.value_);
            const bool lhs_nan = std::isnan(lhs);
            const bool rhs_nan = std::isnan(rhs);
            if (lhs_nan || rhs_nan) {
                if (lhs_nan && rhs_nan) return 0;
                return lhs_nan ? 1 : -1;
            }
            if (lhs < rhs) return -1;
            if (lhs > rhs) return 1;
            return 0;
        }
        if (is_r_real) {
            const Rational lhs = std::holds_alternative<Rational>(value_)
                ? std::get<Rational>(value_)
                : Rational(std::get<BigInt>(value_));
            const Rational rhs = Rational::from_double(std::get<lmmc_real_t>(o.value_));
            if (lhs < rhs) return -1;
            if (rhs < lhs) return 1;
            return -1; // same numeric value, but exact numbers remain structurally distinct
        }

        bool is_l_rat = std::holds_alternative<Rational>(value_);
        bool is_r_rat = std::holds_alternative<Rational>(o.value_);

        if (is_l_rat || is_r_rat) {
            Rational r1 = is_l_rat ? std::get<Rational>(value_) : Rational(std::get<BigInt>(value_));
            Rational r2 = is_r_rat ? std::get<Rational>(o.value_) : Rational(std::get<BigInt>(o.value_));
            if (r1 < r2) return -1;
            if (r1 == r2) return 0;
            return 1;
        }

        const auto& b1 = std::get<BigInt>(value_);
        const auto& b2 = std::get<BigInt>(o.value_);
        if (b1 < b2) return -1;
        if (b1 == b2) return 0;
        return 1;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {

        if (std::holds_alternative<BigInt>(value_)) return lamina::detail::make_node<NumberNode>(std::get<BigInt>(value_));
        if (std::holds_alternative<Rational>(value_)) return lamina::detail::make_node<NumberNode>(std::get<Rational>(value_));
        return lamina::detail::make_node<NumberNode>(std::get<lmmc_real_t>(value_));
    }

    bool is_number() const override { return true; }
    bool is_zero() const override {

        if (std::holds_alternative<BigInt>(value_)) return std::get<BigInt>(value_) == BigInt(0);
        if (std::holds_alternative<Rational>(value_)) return std::get<Rational>(value_) == Rational(0);
        return std::get<lmmc_real_t>(value_) == 0.0;
    }
    bool is_one() const override {
        if (std::holds_alternative<BigInt>(value_)) return std::get<BigInt>(value_) == BigInt(1);
        if (std::holds_alternative<Rational>(value_)) return std::get<Rational>(value_) == Rational(1);
        return std::get<lmmc_real_t>(value_) == 1.0;
    }
    bool is_positive() const override {
        if (std::holds_alternative<BigInt>(value_)) {
            const auto& b = std::get<BigInt>(value_);
            return !b.IsNegative() && !(b == BigInt(0));
        }
        if (std::holds_alternative<Rational>(value_)) {
            const auto& r = std::get<Rational>(value_);
            return r > Rational(0);
        }
        lmmc_real_t v = std::get<lmmc_real_t>(value_);
        return std::isfinite(v) && v > 0.0;
    }
};

/**
 * @brief 复数节点，表示 real + imag * i。
 */
class ComplexNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

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

    int type_priority() const override { return -5; } // 在 NumberNode 和 VariableNode 之间

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
        if (cmp_r != 0) return cmp_r;
        return imag_->compare(*o.imag_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<ComplexNode>(real_->clone(), imag_->clone());
    }
    
    bool is_zero() const override {
        return real_->is_zero() && imag_->is_zero();
    }
};

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

/**
 * @brief 函数节点，表示数学函数调用（三角函数、对数、特殊函数等）。
 */
class FunctionNode : public SymbolicNode {
public:
    /** @brief 函数类型枚举 */
    enum class FuncType {
        Sin, Cos, Tan, Cot, Sec, Csc,       ///< 三角函数
        ArcSin, ArcCos, ArcTan,              ///< 反三角函数
        Sinh, Cosh, Tanh,                    ///< 双曲函数
        Ln, Log, Abs, Sqrt,                  ///< 对数、绝对值、平方根
        Exp,                                 ///< 指数函数
        LambertW,                            ///< Lambert W 函数
        RootOf,                              ///< 多项式根表示
        Atan2,                               ///< 双参数反正切
        Calculus_Integral,                   ///< 积分
        Infinity,                            ///< 无穷大
        Limit,                               ///< 极限
        Erf,                                 ///< 误差函数 erf(x) = (2/√π)∫₀ˣ e^(-t²) dt
        Ei,                                  ///< 指数积分 Ei(x)
        Si,                                  ///< 正弦积分 Si(x) = ∫₀ˣ sin(t)/t dt
        Ci,                                  ///< 余弦积分 Ci(x)
        Li,                                  ///< 对数积分 Li(x) = ∫₀ˣ 1/ln(t) dt
        Max,                                 ///< 最大值 max(a, b)
        Min,                                 ///< 最小值 min(a, b)
        Sgn,                                 ///< 符号函数 sgn(x) ∈ {-1, 0, 1}
        Floor,                               ///< 下取整 ⌊x⌋
        Ceil,                                ///< 上取整 ⌈x⌉
        Round,                               ///< 四舍五入 round(x)
        RealPart,                            ///< 实部 Re(z)
        ImagPart,                            ///< 虚部 Im(z)
        Conjugate,                           ///< 共轭 conj(z)
        ComplexAbs,                          ///< 复数模 |z|
        ComplexArg                           ///< 复数辐角 arg(z)
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const FuncType type_;
    const std::vector<std::shared_ptr<const SymbolicNode>> arguments_;

    FunctionNode(FuncType t, std::vector<std::shared_ptr<const SymbolicNode>> args)
        : type_(t), arguments_(std::move(args)) {
        for (const auto& arg : arguments_) {
            if (!arg) {
                throw std::invalid_argument("FunctionNode argument cannot be null");
            }
        }
    }

public:
    FuncType type() const noexcept { return type_; }
    const std::vector<std::shared_ptr<const SymbolicNode>>& arguments() const noexcept {
        return arguments_;
    }

    int type_priority() const override { return 2; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<size_t>(type_));
        for (const auto& arg : arguments_) {
            hash_combine(seed, arg->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const FunctionNode&>(other);
        if (type_ != o.type_) {
            return static_cast<int>(type_) < static_cast<int>(o.type_) ? -1 : 1;
        }
        if (arguments_.size() != o.arguments_.size()) {
            return arguments_.size() < o.arguments_.size() ? -1 : 1;
        }
        for (size_t i = 0; i < arguments_.size(); ++i) {
            int cmp = arguments_[i]->compare(*o.arguments_[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (const auto& arg : arguments_) new_args.push_back(arg->clone());
        return lamina::detail::make_node<FunctionNode>(type_, std::move(new_args));
    }
};

#include <map>

/**
 * @brief 矩阵节点，支持稠密和稀疏两种存储方式。
 *
 * 当非零元素占比低于 20% 时自动使用稀疏存储，否则使用稠密存储。
 */
class MatrixNode : public SymbolicNode {
public:
    using DenseStorage = std::vector<std::shared_ptr<const SymbolicNode>>;   ///< 稠密存储（按行优先展开）
    using SparseStorage = std::map<size_t, std::shared_ptr<const SymbolicNode>>; ///< 稀疏存储（索引 -> 节点）

    /**
     * @brief 验证并返回网格的列数（要求各行列数一致）。
     * @param grid 二维节点网格
     * @return 列数
     * @throw std::invalid_argument 各行列数不一致时抛出
     */
    static size_t validate_grid_columns(const std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>& grid) {
        if (grid.empty()) {
            throw std::invalid_argument("MatrixNode: matrix must have at least one row");
        }
        size_t ncols = grid[0].size();
        if (ncols == 0) {
            throw std::invalid_argument("MatrixNode: matrix must have at least one column");
        }
        for (const auto& item : grid[0]) {
            if (!item) {
                throw std::invalid_argument("MatrixNode: matrix elements cannot be null");
            }
        }
        for (size_t i = 1; i < grid.size(); ++i) {
            if (grid[i].size() != ncols) {
                throw std::invalid_argument("MatrixNode: all rows must have the same number of columns");
            }
            for (const auto& item : grid[i]) {
                if (!item) {
                    throw std::invalid_argument("MatrixNode: matrix elements cannot be null");
                }
            }
        }
        return ncols;
    }

    static DenseStorage validate_dense_storage(size_t r, size_t c, DenseStorage dense) {
        if (r == 0 || c == 0) {
            throw std::invalid_argument("MatrixNode: matrix dimensions must be non-zero");
        }
        if (r > std::numeric_limits<size_t>::max() / c) {
            throw std::length_error("MatrixNode: matrix dimensions overflow");
        }
        if (dense.size() != r * c) {
            throw std::invalid_argument("MatrixNode: dense storage size does not match dimensions");
        }
        for (const auto& item : dense) {
            if (!item) {
                throw std::invalid_argument("MatrixNode: dense storage elements cannot be null");
            }
        }
        return dense;
    }

    static SparseStorage validate_sparse_storage(size_t r, size_t c, SparseStorage sparse) {
        if (r == 0 || c == 0) {
            throw std::invalid_argument("MatrixNode: matrix dimensions must be non-zero");
        }
        if (r > std::numeric_limits<size_t>::max() / c) {
            throw std::length_error("MatrixNode: matrix dimensions overflow");
        }
        const size_t total = r * c;
        for (const auto& [idx, item] : sparse) {
            if (idx >= total) {
                throw std::invalid_argument("MatrixNode: sparse storage index is out of bounds");
            }
            if (!item) {
                throw std::invalid_argument("MatrixNode: sparse storage elements cannot be null");
            }
        }
        return sparse;
    }

    /**
     * @brief 根据稀疏度选择存储方式，从二维网格创建存储。
     * @param grid 二维节点网格
     * @param total_elements 总元素数
     * @param ncols 列数
     * @return 稠密或稀疏存储
     */
    static std::variant<DenseStorage, SparseStorage> create_storage_from_grid(
        const std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>& grid, size_t total_elements, size_t ncols);

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const size_t rows_;
    const size_t cols_;
    const std::variant<DenseStorage, SparseStorage> storage_;

    MatrixNode(const std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>& grid)
        : rows_(grid.size()),
          cols_(validate_grid_columns(grid)),
          storage_(create_storage_from_grid(grid, rows_ * cols_, cols_)) {}

    MatrixNode(size_t r, size_t c, DenseStorage dense)
        : rows_(r), cols_(c), storage_(validate_dense_storage(r, c, std::move(dense))) {}

    MatrixNode(size_t r, size_t c, SparseStorage sparse)
        : rows_(r), cols_(c), storage_(validate_sparse_storage(r, c, std::move(sparse))) {}

public:
    size_t rows() const noexcept { return rows_; }
    size_t cols() const noexcept { return cols_; }
    const std::variant<DenseStorage, SparseStorage>& storage() const noexcept {
        return storage_;
    }

    int type_priority() const override { return 6; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, rows_);
        hash_combine(seed, cols_);

        if (std::holds_alternative<DenseStorage>(storage_)) {
            const auto& dense = std::get<DenseStorage>(storage_);
            for (size_t i = 0; i < dense.size(); ++i) {
                if (dense[i] && !dense[i]->is_zero()) {
                    hash_combine(seed, i);
                    hash_combine(seed, dense[i]->hash());
                }
            }
        } else {
            const auto& sparse = std::get<SparseStorage>(storage_);
            for (const auto& [idx, val] : sparse) {

                if (val && !val->is_zero()) {
                    hash_combine(seed, idx);
                    hash_combine(seed, val->hash());
                }
            }
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const MatrixNode&>(other);
        if (rows_ != o.rows_) return rows_ < o.rows_ ? -1 : 1;
        if (cols_ != o.cols_) return cols_ < o.cols_ ? -1 : 1;

        for (size_t r = 0; r < rows_; ++r) {
            for (size_t c = 0; c < cols_; ++c) {
                auto v1 = this->get(r, c);
                auto v2 = o.get(r, c);

                bool z1 = !v1 || v1->is_zero();
                bool z2 = !v2 || v2->is_zero();

                if (z1 && z2) continue;
                if (z1) return -1;
                if (z2) return 1;

                int cmp = v1->compare(*v2);
                if (cmp != 0) return cmp;
            }
        }
        return 0;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        if (std::holds_alternative<DenseStorage>(storage_)) {
            const auto& dense = std::get<DenseStorage>(storage_);
            DenseStorage new_dense;
            new_dense.reserve(dense.size());
            for (const auto& e : dense) {
                new_dense.push_back(e->clone());
            }
            return lamina::detail::make_node<MatrixNode>(rows_, cols_, std::move(new_dense));
        } else {
            const auto& sparse = std::get<SparseStorage>(storage_);
            SparseStorage new_sparse;
            for(const auto& [idx, node] : sparse) {
                new_sparse[idx] = node->clone();
            }
            return lamina::detail::make_node<MatrixNode>(rows_, cols_, std::move(new_sparse));
        }
    }

    /**
     * @brief 获取指定位置的矩阵元素。
     * @param r 行索引（从 0 开始）
     * @param c 列索引（从 0 开始）
     * @return 对应节点，越界返回 nullptr，稀疏缺失返回零节点
     */
    std::shared_ptr<const SymbolicNode> get(size_t r, size_t c) const {
        if (r >= rows_ || c >= cols_) return nullptr;

        size_t idx = r * cols_ + c;
        if (std::holds_alternative<DenseStorage>(storage_)) {
            return std::get<DenseStorage>(storage_)[idx];
        } else {
            const auto& sparse = std::get<SparseStorage>(storage_);
            auto it = sparse.find(idx);
            if (it != sparse.end()) return it->second;
            return lamina::detail::make_node<NumberNode>(0.0);
        }
    }

    /** @brief 判断当前矩阵是否使用稀疏存储。 */
    bool is_sparse() const { return std::holds_alternative<SparseStorage>(storage_); }

private:

};

inline std::variant<MatrixNode::DenseStorage, MatrixNode::SparseStorage> MatrixNode::create_storage_from_grid(
    const std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>& grid, size_t total_elements, size_t ncols) {

    size_t non_zeros = 0;
    for (const auto& row : grid) {
        for (const auto& item : row) {
            if (item && !item->is_zero()) non_zeros++;
        }
    }

    bool use_sparse = total_elements > 0 && ((double)non_zeros / total_elements < 0.2);

    if (use_sparse) {
        SparseStorage s;
        size_t r = 0;
        for (const auto& row : grid) {
            size_t c = 0;
            for (const auto& item : row) {
                if (item && !item->is_zero()) {
                    s[r * ncols + c] = item;
                }
                c++;
            }
            r++;
        }
        return s;
    } else {
        DenseStorage d;
        d.reserve(total_elements);
        for (const auto& row : grid) {
            d.insert(d.end(), row.begin(), row.end());
        }
        return d;
    }
}

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

/**
 * @brief 分段函数节点，表示条件分支表达式。
 *
 * 存储有序的 (表达式, 条件) 对列表和可选的默认表达式。
 * 条件应为 RelationalNode 或 LogicalNode 表达式。
 */
class PiecewiseNode : public SymbolicNode {
public:
    /** @brief 分支结构，包含表达式和对应条件。 */
    struct Branch {
        std::shared_ptr<const SymbolicNode> expression; ///< 分支值
        std::shared_ptr<const SymbolicNode> condition;  ///< 条件（RelationalNode 或 LogicalNode）
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::vector<Branch> branches_;
    const std::shared_ptr<const SymbolicNode> default_expr_;

    /**
     * @brief 构造分段函数节点。
     * @param br 分支列表
     * @param def 默认表达式（可为 nullptr）
     */
    PiecewiseNode(std::vector<Branch> br, std::shared_ptr<const SymbolicNode> def = nullptr)
        : branches_(std::move(br)), default_expr_(std::move(def)) {
        if (branches_.empty()) {
            throw std::invalid_argument("PiecewiseNode requires at least one branch");
        }
        for (const auto& branch : branches_) {
            if (!branch.expression || !branch.condition) {
                throw std::invalid_argument("PiecewiseNode branch fields cannot be null");
            }
        }
    }

public:
    const std::vector<Branch>& branches() const noexcept { return branches_; }
    const std::shared_ptr<const SymbolicNode>& default_expr() const noexcept {
        return default_expr_;
    }

    int type_priority() const override { return 7; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        for (const auto& b : branches_) {
            hash_combine(seed, b.expression->hash());
            hash_combine(seed, b.condition->hash());
        }
        if (default_expr_) {
            hash_combine(seed, default_expr_->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const PiecewiseNode&>(other);
        if (branches_.size() != o.branches_.size()) {
            return branches_.size() < o.branches_.size() ? -1 : 1;
        }
        for (size_t i = 0; i < branches_.size(); ++i) {
            int cmp = branches_[i].expression->compare(*o.branches_[i].expression);
            if (cmp != 0) return cmp;
            cmp = branches_[i].condition->compare(*o.branches_[i].condition);
            if (cmp != 0) return cmp;
        }
        bool has_def = (default_expr_ != nullptr);
        bool o_has_def = (o.default_expr_ != nullptr);
        if (has_def != o_has_def) return has_def ? 1 : -1;
        if (has_def && o_has_def) {
            return default_expr_->compare(*o.default_expr_);
        }
        return 0;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<Branch> new_branches;
        new_branches.reserve(branches_.size());
        for (const auto& b : branches_) {
            new_branches.push_back({b.expression->clone(), b.condition->clone()});
        }
        auto new_def = default_expr_ ? default_expr_->clone() : nullptr;
        return lamina::detail::make_node<PiecewiseNode>(std::move(new_branches), std::move(new_def));
    }
};

/**
 * @brief 求和节点，表示符号有限/无限求和 ∑_{k=a}^{b} f(k)。
 */
class SummationNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> body_;
    const std::string index_var_;
    const std::shared_ptr<const SymbolicNode> lower_bound_;
    const std::shared_ptr<const SymbolicNode> upper_bound_;

    /**
     * @brief 构造求和节点。
     * @param b 通项表达式
     * @param idx 指标变量名
     * @param lo 下界
     * @param hi 上界
     */
    SummationNode(std::shared_ptr<const SymbolicNode> b, std::string idx,
                  std::shared_ptr<const SymbolicNode> lo, std::shared_ptr<const SymbolicNode> hi)
        : body_(std::move(b)), index_var_(std::move(idx)),
          lower_bound_(std::move(lo)), upper_bound_(std::move(hi)) {
        if (!body_ || !lower_bound_ || !upper_bound_) {
            throw std::invalid_argument("SummationNode children cannot be null");
        }
        if (index_var_.empty()) {
            throw std::invalid_argument("SummationNode index variable cannot be empty");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& body() const noexcept { return body_; }
    const std::string& index_var() const noexcept { return index_var_; }
    const std::shared_ptr<const SymbolicNode>& lower_bound() const noexcept {
        return lower_bound_;
    }
    const std::shared_ptr<const SymbolicNode>& upper_bound() const noexcept {
        return upper_bound_;
    }

    int type_priority() const override { return 8; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, body_->hash());
        hash_combine(seed, std::hash<std::string>{}(index_var_));
        hash_combine(seed, lower_bound_->hash());
        hash_combine(seed, upper_bound_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const SummationNode&>(other);
        int cmp = index_var_.compare(o.index_var_);
        if (cmp != 0) return cmp;
        cmp = lower_bound_->compare(*o.lower_bound_);
        if (cmp != 0) return cmp;
        cmp = upper_bound_->compare(*o.upper_bound_);
        if (cmp != 0) return cmp;
        return body_->compare(*o.body_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<SummationNode>(
            body_->clone(), index_var_, lower_bound_->clone(), upper_bound_->clone());
    }
};

/**
 * @brief 连乘节点，表示符号有限/无限连乘 ∏_{k=a}^{b} f(k)。
 */
class ProductNode_Op : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> body_;
    const std::string index_var_;
    const std::shared_ptr<const SymbolicNode> lower_bound_;
    const std::shared_ptr<const SymbolicNode> upper_bound_;

    /**
     * @brief 构造连乘节点。
     * @param b 通项表达式
     * @param idx 指标变量名
     * @param lo 下界
     * @param hi 上界
     */
    ProductNode_Op(std::shared_ptr<const SymbolicNode> b, std::string idx,
                   std::shared_ptr<const SymbolicNode> lo, std::shared_ptr<const SymbolicNode> hi)
        : body_(std::move(b)), index_var_(std::move(idx)),
          lower_bound_(std::move(lo)), upper_bound_(std::move(hi)) {
        if (!body_ || !lower_bound_ || !upper_bound_) {
            throw std::invalid_argument("ProductNode children cannot be null");
        }
        if (index_var_.empty()) {
            throw std::invalid_argument("ProductNode index variable cannot be empty");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& body() const noexcept { return body_; }
    const std::string& index_var() const noexcept { return index_var_; }
    const std::shared_ptr<const SymbolicNode>& lower_bound() const noexcept {
        return lower_bound_;
    }
    const std::shared_ptr<const SymbolicNode>& upper_bound() const noexcept {
        return upper_bound_;
    }

    int type_priority() const override { return 9; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, body_->hash());
        hash_combine(seed, std::hash<std::string>{}(index_var_));
        hash_combine(seed, lower_bound_->hash());
        hash_combine(seed, upper_bound_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const ProductNode_Op&>(other);
        int cmp = index_var_.compare(o.index_var_);
        if (cmp != 0) return cmp;
        cmp = lower_bound_->compare(*o.lower_bound_);
        if (cmp != 0) return cmp;
        cmp = upper_bound_->compare(*o.upper_bound_);
        if (cmp != 0) return cmp;
        return body_->compare(*o.body_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<ProductNode_Op>(
            body_->clone(), index_var_, lower_bound_->clone(), upper_bound_->clone());
    }
};

/**
 * @brief 积分变换节点，表示 Laplace、Fourier、Z 变换及其逆变换。
 *
 * 统一表示各类积分变换：L{f(t)}(s)、F{f(t)}(ω)、Z{f[n]}(z) 等。
 */
class TransformNode : public SymbolicNode {
public:
    /** @brief 变换类型枚举 */
    enum class TransformType {
        Laplace,         ///< Laplace 变换 L{f(t)}(s)
        InverseLaplace,  ///< 逆 Laplace 变换 L⁻¹{F(s)}(t)
        Fourier,         ///< Fourier 变换 F{f(t)}(ω)
        InverseFourier,  ///< 逆 Fourier 变换 F⁻¹{F(ω)}(t)
        ZTransform       ///< Z 变换 Z{f[n]}(z)
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const TransformType transform_type_;
    const std::shared_ptr<const SymbolicNode> body_;
    const std::string source_var_;
    const std::string target_var_;

    /**
     * @brief 构造积分变换节点。
     * @param tt 变换类型
     * @param b 被变换的表达式
     * @param src 源变量名
     * @param tgt 目标变量名
     */
    TransformNode(TransformType tt, std::shared_ptr<const SymbolicNode> b,
                  std::string src, std::string tgt)
        : transform_type_(tt), body_(std::move(b)),
          source_var_(std::move(src)), target_var_(std::move(tgt)) {
        if (!body_) {
            throw std::invalid_argument("TransformNode body cannot be null");
        }
        if (source_var_.empty() || target_var_.empty()) {
            throw std::invalid_argument("TransformNode variables cannot be empty");
        }
    }

public:
    TransformType transform_type() const noexcept { return transform_type_; }
    const std::shared_ptr<const SymbolicNode>& body() const noexcept { return body_; }
    const std::string& source_var() const noexcept { return source_var_; }
    const std::string& target_var() const noexcept { return target_var_; }

    int type_priority() const override { return 11; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<std::size_t>(transform_type_));
        hash_combine(seed, body_->hash());
        hash_combine(seed, std::hash<std::string>{}(source_var_));
        hash_combine(seed, std::hash<std::string>{}(target_var_));
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const TransformNode&>(other);
        if (transform_type_ != o.transform_type_) {
            return static_cast<int>(transform_type_) < static_cast<int>(o.transform_type_) ? -1 : 1;
        }
        int cmp = source_var_.compare(o.source_var_);
        if (cmp != 0) return cmp;
        cmp = target_var_.compare(o.target_var_);
        if (cmp != 0) return cmp;
        return body_->compare(*o.body_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<TransformNode>(
            transform_type_, body_->clone(), source_var_, target_var_);
    }
};

/**
 * @brief 量词节点，表示全称量词 (∀) 或存在量词 (∃)。
 *
 * 用于表示逻辑公式中的量化表达式，如 ∀x∈D: P(x) 或 ∃x∈D: P(x)。
 */
class QuantifierNode : public SymbolicNode {
public:
    /** @brief 量词类型枚举 */
    enum class Type {
        ForAll, ///< 全称量词 ∀
        Exists  ///< 存在量词 ∃
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const Type quantifier_type_;
    const std::string bound_var_;
    const std::shared_ptr<const SymbolicNode> domain_;
    const std::shared_ptr<const SymbolicNode> predicate_;

    /**
     * @brief 构造量词节点。
     * @param qt 量词类型
     * @param var 约束变量名
     * @param dom 定义域表达式
     * @param pred 谓词表达式
     */
    QuantifierNode(Type qt, std::string var,
                   std::shared_ptr<const SymbolicNode> dom, std::shared_ptr<const SymbolicNode> pred)
        : quantifier_type_(qt), bound_var_(std::move(var)),
          domain_(std::move(dom)), predicate_(std::move(pred)) {
        if (!domain_ || !predicate_) {
            throw std::invalid_argument("QuantifierNode domain and predicate cannot be null");
        }
        if (bound_var_.empty()) {
            throw std::invalid_argument("QuantifierNode bound variable cannot be empty");
        }
    }

public:
    Type quantifier_type() const noexcept { return quantifier_type_; }
    const std::string& bound_var() const noexcept { return bound_var_; }
    const std::shared_ptr<const SymbolicNode>& domain() const noexcept { return domain_; }
    const std::shared_ptr<const SymbolicNode>& predicate() const noexcept {
        return predicate_;
    }

    int type_priority() const override { return 102; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<std::size_t>(quantifier_type_));
        hash_combine(seed, std::hash<std::string>{}(bound_var_));
        hash_combine(seed, domain_->hash());
        hash_combine(seed, predicate_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const QuantifierNode&>(other);
        if (quantifier_type_ != o.quantifier_type_) {
            return static_cast<int>(quantifier_type_) < static_cast<int>(o.quantifier_type_) ? -1 : 1;
        }
        int cmp = bound_var_.compare(o.bound_var_);
        if (cmp != 0) return cmp;
        cmp = domain_->compare(*o.domain_);
        if (cmp != 0) return cmp;
        return predicate_->compare(*o.predicate_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<QuantifierNode>(
            quantifier_type_, bound_var_, domain_->clone(), predicate_->clone());
    }
};

/**
 * @brief 集合构造器节点，表示集合构造式 {x ∈ D | P(x)}。
 *
 * 用于表示满足特定条件的元素集合。
 */
class SetBuilderNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::string element_var_;
    const std::shared_ptr<const SymbolicNode> domain_;
    const std::shared_ptr<const SymbolicNode> predicate_;

    /**
     * @brief 构造集合构造器节点。
     * @param var 元素变量名
     * @param dom 定义域表达式
     * @param pred 成员条件表达式
     */
    SetBuilderNode(std::string var, std::shared_ptr<const SymbolicNode> dom,
                   std::shared_ptr<const SymbolicNode> pred)
        : element_var_(std::move(var)), domain_(std::move(dom)),
          predicate_(std::move(pred)) {
        if (!domain_ || !predicate_) {
            throw std::invalid_argument("SetBuilderNode domain and predicate cannot be null");
        }
        if (element_var_.empty()) {
            throw std::invalid_argument("SetBuilderNode element variable cannot be empty");
        }
    }

public:
    const std::string& element_var() const noexcept { return element_var_; }
    const std::shared_ptr<const SymbolicNode>& domain() const noexcept { return domain_; }
    const std::shared_ptr<const SymbolicNode>& predicate() const noexcept {
        return predicate_;
    }

    int type_priority() const override { return 103; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, std::hash<std::string>{}(element_var_));
        hash_combine(seed, domain_->hash());
        hash_combine(seed, predicate_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const SetBuilderNode&>(other);
        int cmp = element_var_.compare(o.element_var_);
        if (cmp != 0) return cmp;
        cmp = domain_->compare(*o.domain_);
        if (cmp != 0) return cmp;
        return predicate_->compare(*o.predicate_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<SetBuilderNode>(
            element_var_, domain_->clone(), predicate_->clone());
    }
};

class FiniteSetNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::vector<std::shared_ptr<const SymbolicNode>> elements_;

    explicit FiniteSetNode(std::vector<std::shared_ptr<const SymbolicNode>> elements)
        : elements_(std::move(elements)) {
        for (const auto& element : elements_) {
            if (!element) {
                throw std::invalid_argument("FiniteSetNode element cannot be null");
            }
        }
    }

public:
    const std::vector<std::shared_ptr<const SymbolicNode>>& elements() const noexcept {
        return elements_;
    }

    int type_priority() const override { return 104; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, elements_.size());
        for (const auto& element : elements_) {
            hash_combine(seed, element->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const FiniteSetNode&>(other);
        if (elements_.size() != o.elements_.size()) {
            return elements_.size() < o.elements_.size() ? -1 : 1;
        }
        for (size_t i = 0; i < elements_.size(); ++i) {
            int cmp = elements_[i]->compare(*o.elements_[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_elements;
        new_elements.reserve(elements_.size());
        for (const auto& element : elements_) {
            new_elements.push_back(element->clone());
        }
        return lamina::detail::make_node<FiniteSetNode>(std::move(new_elements));
    }
};

class IntervalNode : public SymbolicNode {
public:
    enum class Bound {
        Open,
        Closed
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> lower_;
    const std::shared_ptr<const SymbolicNode> upper_;
    const Bound lower_bound_;
    const Bound upper_bound_;

    IntervalNode(std::shared_ptr<const SymbolicNode> lower,
                 std::shared_ptr<const SymbolicNode> upper,
                 Bound lower_bound,
                 Bound upper_bound)
        : lower_(std::move(lower)),
          upper_(std::move(upper)),
          lower_bound_(lower_bound),
          upper_bound_(upper_bound) {
        if (!lower_ || !upper_) {
            throw std::invalid_argument("IntervalNode bounds cannot be null");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& lower() const noexcept { return lower_; }
    const std::shared_ptr<const SymbolicNode>& upper() const noexcept { return upper_; }
    Bound lower_bound() const noexcept { return lower_bound_; }
    Bound upper_bound() const noexcept { return upper_bound_; }

    int type_priority() const override { return 105; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, lower_->hash());
        hash_combine(seed, upper_->hash());
        hash_combine(seed, static_cast<std::size_t>(lower_bound_));
        hash_combine(seed, static_cast<std::size_t>(upper_bound_));
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const IntervalNode&>(other);
        if (lower_bound_ != o.lower_bound_) {
            return static_cast<int>(lower_bound_) < static_cast<int>(o.lower_bound_) ? -1 : 1;
        }
        if (upper_bound_ != o.upper_bound_) {
            return static_cast<int>(upper_bound_) < static_cast<int>(o.upper_bound_) ? -1 : 1;
        }
        int cmp = lower_->compare(*o.lower_);
        if (cmp != 0) return cmp;
        return upper_->compare(*o.upper_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<IntervalNode>(
            lower_->clone(), upper_->clone(), lower_bound_, upper_bound_);
    }
};

class MembershipNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> element_;
    const std::shared_ptr<const SymbolicNode> set_;
    const bool negated_;

    MembershipNode(std::shared_ptr<const SymbolicNode> element,
                   std::shared_ptr<const SymbolicNode> set,
                   bool negated)
        : element_(std::move(element)), set_(std::move(set)), negated_(negated) {
        if (!element_ || !set_) {
            throw std::invalid_argument("MembershipNode operands cannot be null");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& element() const noexcept { return element_; }
    const std::shared_ptr<const SymbolicNode>& set() const noexcept { return set_; }
    bool negated() const noexcept { return negated_; }

    int type_priority() const override { return 106; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, element_->hash());
        hash_combine(seed, set_->hash());
        hash_combine(seed, static_cast<std::size_t>(negated_));
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const MembershipNode&>(other);
        if (negated_ != o.negated_) return negated_ ? 1 : -1;
        int cmp = element_->compare(*o.element_);
        if (cmp != 0) return cmp;
        return set_->compare(*o.set_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<MembershipNode>(
            element_->clone(), set_->clone(), negated_);
    }
};

class UninterpretedFunctionNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::string name_;
    const std::vector<std::shared_ptr<const SymbolicNode>> arguments_;

    UninterpretedFunctionNode(std::string name,
                              std::vector<std::shared_ptr<const SymbolicNode>> arguments)
        : name_(std::move(name)), arguments_(std::move(arguments)) {
        if (name_.empty()) {
            throw std::invalid_argument("UninterpretedFunctionNode name cannot be empty");
        }
        for (const auto& argument : arguments_) {
            if (!argument) {
                throw std::invalid_argument("UninterpretedFunctionNode argument cannot be null");
            }
        }
    }

public:
    const std::string& name() const noexcept { return name_; }
    const std::vector<std::shared_ptr<const SymbolicNode>>& arguments() const noexcept {
        return arguments_;
    }

    int type_priority() const override { return 107; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, std::hash<std::string>{}(name_));
        hash_combine(seed, arguments_.size());
        for (const auto& argument : arguments_) {
            hash_combine(seed, argument->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const UninterpretedFunctionNode&>(other);
        int cmp = name_.compare(o.name_);
        if (cmp != 0) return cmp;
        if (arguments_.size() != o.arguments_.size()) {
            return arguments_.size() < o.arguments_.size() ? -1 : 1;
        }
        for (size_t i = 0; i < arguments_.size(); ++i) {
            cmp = arguments_[i]->compare(*o.arguments_[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_arguments;
        new_arguments.reserve(arguments_.size());
        for (const auto& argument : arguments_) {
            new_arguments.push_back(argument->clone());
        }
        return lamina::detail::make_node<UninterpretedFunctionNode>(
            name_, std::move(new_arguments));
    }
};

#undef LAMINA_AST_NODE_FACTORY_FRIEND

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(const ::BigInt& v) { return lamina::detail::make_node<NumberNode>(v); }
inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(const ::Rational& v) { return lamina::detail::make_node<NumberNode>(v); }
inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(lmmc_real_t v) { return lamina::detail::make_node<NumberNode>(v); }
inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_variable(const std::string& name) { return lamina::detail::make_node<VariableNode>(name); }

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_add(std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    if (ops.empty()) return create_number(::BigInt(0));
    for (const auto& op : ops) {
        if (!op) {
            throw std::invalid_argument("create_add operand cannot be null");
        }
    }
    if (ops.size() == 1) return ops[0];

    std::vector<std::shared_ptr<const SymbolicNode>> flat_ops;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_zero()) continue;
        if (auto add = std::dynamic_pointer_cast<const AddNode>(op)) {
            flat_ops.insert(flat_ops.end(), add->operands().begin(), add->operands().end());
        } else {
            flat_ops.push_back(op);
        }
    }

    if (flat_ops.empty()) return create_number(::BigInt(0));
    if (flat_ops.size() == 1) return flat_ops[0];
    return lamina::detail::make_node<AddNode>(std::move(flat_ops));
}

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_multiply(std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    if (ops.empty()) return create_number(::BigInt(1));

    for (const auto& op : ops) {
        if (!op) {
            throw std::invalid_argument("create_multiply operand cannot be null");
        }
        if (op->is_zero()) return op;
    }

    if (ops.size() == 1) return ops[0];

    std::vector<std::shared_ptr<const SymbolicNode>> flat_ops;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_one()) continue;
        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(op)) {
            flat_ops.insert(flat_ops.end(), mul->operands().begin(), mul->operands().end());
        } else {
            flat_ops.push_back(op);
        }
    }

    if (flat_ops.empty()) return create_number(::BigInt(1));
    if (flat_ops.size() == 1) return flat_ops[0];
    return lamina::detail::make_node<MultiplyNode>(std::move(flat_ops));
}

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_power(std::shared_ptr<const SymbolicNode> base, std::shared_ptr<const SymbolicNode> exponent) {
    if (!base || !exponent) {
        throw std::invalid_argument("create_power operands cannot be null");
    }
    if (exponent->is_zero()) {
        // 不能盲目把 x^0 折叠成 1：base 可能为 0 时原表达式有定义域条件。
        if (base->is_zero()) {
            return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
        }
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(base)) {
            if (!num->is_zero()) return create_number(::BigInt(1));
        }
        if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(base)) {
            auto real_num = std::dynamic_pointer_cast<const NumberNode>(complex->real());
            auto imag_num = std::dynamic_pointer_cast<const NumberNode>(complex->imag());
            if ((real_num && !real_num->is_zero()) || (imag_num && !imag_num->is_zero())) {
                return create_number(::BigInt(1));
            }
        }
        return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
    }
    if (exponent->is_one()) return base;
    if (base->is_zero()) {
        // 0^x = 0 仅当指数严格为正时才安全；x<=0 或符号未知时保留 PowerNode。
        if (exponent->is_positive()) return create_number(::BigInt(0));
        return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
    }
    if (base->is_one()) return create_number(::BigInt(1));
    return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
}

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_complex(std::shared_ptr<const SymbolicNode> real, std::shared_ptr<const SymbolicNode> imag) {
    if (!real || !imag) {
        throw std::invalid_argument("create_complex operands cannot be null");
    }
    if (imag->is_zero()) {
        return real;
    }
    return lamina::detail::make_node<ComplexNode>(std::move(real), std::move(imag));
}
