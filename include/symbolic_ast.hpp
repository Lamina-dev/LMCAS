/**
 * @file symbolic_ast.hpp
 * @brief AST 节点定义：NumberNode, VariableNode, AddNode, MultiplyNode, PowerNode, FunctionNode, MatrixNode, RelationalNode, LogicalNode。
 */
#pragma once
#include <memory>
#include <vector>
#include <string>
#include <variant>
#include <algorithm>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include "rational.hpp"
#include "bigint.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <stdexcept>
#include <atomic>

class SymbolicVisitor;
class NumberNode;
class VariableNode;
class AddNode;
class MultiplyNode;
class PowerNode;
class FunctionNode;
class MatrixNode;
class RelationalNode;
class LogicalNode;

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
    virtual void accept(SymbolicVisitor& visitor) = 0;

    /** @brief 深拷贝当前节点及其子树。 */
    virtual std::shared_ptr<SymbolicNode> clone() const = 0;

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
};

/**
 * @brief AST 访问者基类（Visitor 模式）。
 *
 * 子类实现各 visit 方法以对不同节点类型执行操作。
 * 内置深度保护，防止递归过深导致栈溢出。
 */
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

    virtual void visit(NumberNode& node) = 0;
    virtual void visit(VariableNode& node) = 0;
    virtual void visit(AddNode& node) = 0;
    virtual void visit(MultiplyNode& node) = 0;
    virtual void visit(PowerNode& node) = 0;
    virtual void visit(FunctionNode& node) = 0;
    virtual void visit(MatrixNode& node) = 0;
    virtual void visit(RelationalNode& node) {}
    virtual void visit(LogicalNode& node) {}
};

/** @brief 基于节点哈希的哈希函数对象，用于无序容器。 */
struct NodeHash {
    std::size_t operator()(const std::shared_ptr<SymbolicNode>& node) const {
        return node ? node->hash() : 0;
    }
};

/** @brief 基于节点结构相等性的比较函数对象，用于无序容器。 */
struct NodeEqual {
    bool operator()(const std::shared_ptr<SymbolicNode>& lhs, const std::shared_ptr<SymbolicNode>& rhs) const {
        if (!lhs || !rhs) return lhs == rhs;
        return lhs->equals(*rhs);
    }
};

template<typename T>
using NodeMap = std::unordered_map<std::shared_ptr<SymbolicNode>, T, NodeHash, NodeEqual>;

using NodeSet = std::unordered_set<std::shared_ptr<SymbolicNode>, NodeHash, NodeEqual>;

/**
 * @brief 符号节点工厂类，提供创建常用节点的静态方法。
 *
 * 创建时自动执行基本简化（如零元素消除、单元素折叠、扁平化）。
 */
class SymbolicFactory {
public:
    /** @brief 创建大整数数值节点。 */
    static std::shared_ptr<SymbolicNode> create_number(const ::BigInt& v);
    /** @brief 创建有理数数值节点。 */
    static std::shared_ptr<SymbolicNode> create_number(const ::Rational& v);
    /** @brief 创建浮点数值节点。 */
    static std::shared_ptr<SymbolicNode> create_number(lmmc_real_t v);
    /** @brief 创建变量节点。 */
    static std::shared_ptr<SymbolicNode> create_variable(const std::string& name);

    /**
     * @brief 创建加法节点，自动扁平化嵌套加法并消除零项。
     * @param ops 操作数列表
     * @return 简化后的节点
     */
    static std::shared_ptr<SymbolicNode> create_add(std::vector<std::shared_ptr<SymbolicNode>> ops);

    /**
     * @brief 创建乘法节点，自动扁平化嵌套乘法并消除单位元。
     * @param ops 操作数列表
     * @return 简化后的节点
     */
    static std::shared_ptr<SymbolicNode> create_multiply(std::vector<std::shared_ptr<SymbolicNode>> ops);

    /**
     * @brief 创建幂运算节点，自动处理指数为 0/1 及底数为 0/1 的情况。
     * @param base 底数节点
     * @param exponent 指数节点
     * @return 简化后的节点
     */
    static std::shared_ptr<SymbolicNode> create_power(std::shared_ptr<SymbolicNode> base, std::shared_ptr<SymbolicNode> exponent);
};

/**
 * @brief 数值节点，存储 BigInt、Rational 或浮点数。
 */
class NumberNode : public SymbolicNode {
public:
    std::variant<BigInt, Rational, lmmc_real_t> value; ///< 数值，支持大整数、有理数、浮点数

    explicit NumberNode(const BigInt& v) : value(v) {}
    explicit NumberNode(const Rational& v) : value(v) {}
    explicit NumberNode(lmmc_real_t v) : value(v) {}
    explicit NumberNode(std::variant<BigInt, Rational, lmmc_real_t> v) : value(std::move(v)) {}

    int type_priority() const override { return -10; }

protected:
    std::size_t compute_hash() const override {
        if (std::holds_alternative<Rational>(value)) {
            const auto& r = std::get<Rational>(value);
            return r.hash();
        } else if (std::holds_alternative<BigInt>(value)) {
            const auto& b = std::get<BigInt>(value);
            return b.hash();
        } else {
             auto d = std::get<lmmc_real_t>(value);
             return std::hash<lmmc_real_t>{}(d);
        }
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const NumberNode&>(other);

        bool is_l_real = std::holds_alternative<lmmc_real_t>(value);
        bool is_r_real = std::holds_alternative<lmmc_real_t>(o.value);

        if (is_l_real || is_r_real) {
            auto to_lmmc_real = [](const std::variant<BigInt, Rational, lmmc_real_t>& v) {
                if (std::holds_alternative<lmmc_real_t>(v)) return std::get<lmmc_real_t>(v);
                if (std::holds_alternative<Rational>(v)) return (lmmc_real_t)std::get<Rational>(v).to_double();
                return (lmmc_real_t)std::get<BigInt>(v).to_double();
            };

            lmmc_real_t v1 = to_lmmc_real(value);
            lmmc_real_t v2 = to_lmmc_real(o.value);
            return LMMC_REAL_CMP(&v1, &v2);
        }

        bool is_l_rat = std::holds_alternative<Rational>(value);
        bool is_r_rat = std::holds_alternative<Rational>(o.value);

        if (is_l_rat || is_r_rat) {
            Rational r1 = is_l_rat ? std::get<Rational>(value) : Rational(std::get<BigInt>(value));
            Rational r2 = is_r_rat ? std::get<Rational>(o.value) : Rational(std::get<BigInt>(o.value));
            if (r1 < r2) return -1;
            if (r1 == r2) return 0;
            return 1;
        }

        const auto& b1 = std::get<BigInt>(value);
        const auto& b2 = std::get<BigInt>(o.value);
        if (b1 < b2) return -1;
        if (b1 == b2) return 0;
        return 1;
    }

public:
    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {

        if (std::holds_alternative<BigInt>(value)) return std::make_shared<NumberNode>(std::get<BigInt>(value));
        if (std::holds_alternative<Rational>(value)) return std::make_shared<NumberNode>(std::get<Rational>(value));
        return std::make_shared<NumberNode>(std::get<lmmc_real_t>(value));
    }

    bool is_number() const override { return true; }
    bool is_zero() const override {

        if (std::holds_alternative<BigInt>(value)) return std::get<BigInt>(value) == BigInt(0);
        if (std::holds_alternative<Rational>(value)) return std::get<Rational>(value) == Rational(0);
        lmmc_real_t v = std::get<lmmc_real_t>(value);
        int eq;
        lmmc_double_nearly_equal(v, 0.0, &eq);
        return eq != 0;
    }
    bool is_one() const override {
        if (std::holds_alternative<BigInt>(value)) return std::get<BigInt>(value) == BigInt(1);
        if (std::holds_alternative<Rational>(value)) return std::get<Rational>(value) == Rational(1);
        lmmc_real_t v = std::get<lmmc_real_t>(value);
        int eq;
        lmmc_double_nearly_equal(v, 1.0, &eq);
        return eq != 0;
    }
};

/**
 * @brief 变量节点，表示一个符号变量。
 */
class VariableNode : public SymbolicNode {
public:
    std::string name; ///< 变量名

    explicit VariableNode(std::string n) : name(std::move(n)) {}

    int type_priority() const override { return 10; }

protected:
    std::size_t compute_hash() const override {
        return std::hash<std::string>{}(name);
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const VariableNode&>(other);
        return name.compare(o.name);
    }

public:
    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override { return std::make_shared<VariableNode>(name); }
};

/**
 * @brief 加法节点，表示多个操作数的求和。
 *
 * 构造时自动扁平化嵌套的 AddNode，并按规范顺序排序操作数。
 */
class AddNode : public SymbolicNode {
public:
    std::vector<std::shared_ptr<SymbolicNode>> operands; ///< 加法操作数列表

    explicit AddNode(std::vector<std::shared_ptr<SymbolicNode>> ops) {

        for (const auto& op : ops) {
            if (auto add = std::dynamic_pointer_cast<AddNode>(op)) {
                operands.insert(operands.end(), add->operands.begin(), add->operands.end());
            } else {
                operands.push_back(op);
            }
        }

        std::sort(operands.begin(), operands.end(), [](const auto& a, const auto& b) {
            bool a_num = a->is_number();
            bool b_num = b->is_number();
            if (a_num != b_num) return !a_num;
            return a->compare(*b) < 0;
        });
    }

    int type_priority() const override { return 5; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        for (const auto& op : operands) {
            hash_combine(seed, op->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const AddNode&>(other);
        if (operands.size() != o.operands.size()) {
            return operands.size() < o.operands.size() ? -1 : 1;
        }
        for (size_t i = 0; i < operands.size(); ++i) {
            int cmp = operands[i]->compare(*o.operands[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        new_ops.reserve(operands.size());
        for (const auto& op : operands) new_ops.push_back(op->clone());
        return std::make_shared<AddNode>(std::move(new_ops));
    }
};

/**
 * @brief 乘法节点，表示多个操作数的乘积。
 *
 * 构造时自动扁平化嵌套的 MultiplyNode，并按规范顺序排序操作数。
 */
class MultiplyNode : public SymbolicNode {
public:
    std::vector<std::shared_ptr<SymbolicNode>> operands; ///< 乘法操作数列表

    explicit MultiplyNode(std::vector<std::shared_ptr<SymbolicNode>> ops) {

        for (const auto& op : ops) {
            if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {
                operands.insert(operands.end(), mul->operands.begin(), mul->operands.end());
            } else {
                operands.push_back(op);
            }
        }

        std::sort(operands.begin(), operands.end(), [](const auto& a, const auto& b) {
            return a->compare(*b) < 0;
        });
    }

    int type_priority() const override { return 4; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        for (const auto& op : operands) {
            hash_combine(seed, op->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const MultiplyNode&>(other);
        if (operands.size() != o.operands.size()) {
            return operands.size() < o.operands.size() ? -1 : 1;
        }
        for (size_t i = 0; i < operands.size(); ++i) {
            int cmp = operands[i]->compare(*o.operands[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        new_ops.reserve(operands.size());
        for (const auto& op : operands) new_ops.push_back(op->clone());
        return std::make_shared<MultiplyNode>(std::move(new_ops));
    }
};

/**
 * @brief 幂运算节点，表示 base^exponent。
 */
class PowerNode : public SymbolicNode {
public:
    std::shared_ptr<SymbolicNode> base;     ///< 底数
    std::shared_ptr<SymbolicNode> exponent; ///< 指数

    PowerNode(std::shared_ptr<SymbolicNode> b, std::shared_ptr<SymbolicNode> e)
        : base(std::move(b)), exponent(std::move(e)) {}

    int type_priority() const override { return 0; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, base->hash());
        hash_combine(seed, exponent->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const PowerNode&>(other);

        int cmp_exp = exponent->compare(*o.exponent);
        if (cmp_exp != 0) return cmp_exp > 0 ? -1 : 1;

        return base->compare(*o.base);
    }

public:
    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        return std::make_shared<PowerNode>(base->clone(), exponent->clone());
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
        Limit                                ///< 极限
    };

    FuncType type;                                          ///< 函数类型
    std::vector<std::shared_ptr<SymbolicNode>> arguments;   ///< 函数参数列表

    FunctionNode(FuncType t, std::vector<std::shared_ptr<SymbolicNode>> args)
        : type(t), arguments(std::move(args)) {}

    int type_priority() const override { return 2; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<size_t>(type));
        for (const auto& arg : arguments) {
            hash_combine(seed, arg->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const FunctionNode&>(other);
        if (type != o.type) {
            return static_cast<int>(type) < static_cast<int>(o.type) ? -1 : 1;
        }
        if (arguments.size() != o.arguments.size()) {
            return arguments.size() < o.arguments.size() ? -1 : 1;
        }
        for (size_t i = 0; i < arguments.size(); ++i) {
            int cmp = arguments[i]->compare(*o.arguments[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<SymbolicNode>> new_args;
        for (const auto& arg : arguments) new_args.push_back(arg->clone());
        return std::make_shared<FunctionNode>(type, std::move(new_args));
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
    using DenseStorage = std::vector<std::shared_ptr<SymbolicNode>>;   ///< 稠密存储（按行优先展开）
    using SparseStorage = std::map<size_t, std::shared_ptr<SymbolicNode>>; ///< 稀疏存储（索引 -> 节点）

    /**
     * @brief 验证并返回网格的列数（取各行最大长度）。
     * @param grid 二维节点网格
     * @return 列数
     */
    static size_t validate_grid_columns(const std::vector<std::vector<std::shared_ptr<SymbolicNode>>>& grid) {
        if (grid.empty()) return 0;
        size_t ncols = grid[0].size();
        for (size_t i = 1; i < grid.size(); ++i) {
            if (grid[i].size() != ncols) {

                if (grid[i].size() > ncols) ncols = grid[i].size();
            }
        }
        return ncols;
    }

    /**
     * @brief 根据稀疏度选择存储方式，从二维网格创建存储。
     * @param grid 二维节点网格
     * @param total_elements 总元素数
     * @param ncols 列数
     * @return 稠密或稀疏存储
     */
    static std::variant<DenseStorage, SparseStorage> create_storage_from_grid(
        const std::vector<std::vector<std::shared_ptr<SymbolicNode>>>& grid, size_t total_elements, size_t ncols);

    const size_t rows; ///< 行数
    const size_t cols; ///< 列数
    const std::variant<DenseStorage, SparseStorage> storage; ///< 实际存储

    MatrixNode(const std::vector<std::vector<std::shared_ptr<SymbolicNode>>>& grid)
        : rows(grid.size()),
          cols(grid.empty() ? 0 : validate_grid_columns(grid)),
          storage(create_storage_from_grid(grid, rows * cols, cols)) {}

    MatrixNode(size_t r, size_t c, DenseStorage dense)
        : storage(std::move(dense)), rows(r), cols(c) {}

    MatrixNode(size_t r, size_t c, SparseStorage sparse)
        : storage(std::move(sparse)), rows(r), cols(c) {}

    int type_priority() const override { return 6; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, rows);
        hash_combine(seed, cols);

        if (std::holds_alternative<DenseStorage>(storage)) {
            const auto& dense = std::get<DenseStorage>(storage);
            for (size_t i = 0; i < dense.size(); ++i) {
                if (dense[i] && !dense[i]->is_zero()) {
                    hash_combine(seed, i);
                    hash_combine(seed, dense[i]->hash());
                }
            }
        } else {
            const auto& sparse = std::get<SparseStorage>(storage);
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
        if (rows != o.rows) return rows < o.rows ? -1 : 1;
        if (cols != o.cols) return cols < o.cols ? -1 : 1;

        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
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
    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<SymbolicNode> clone() const override {
        if (std::holds_alternative<DenseStorage>(storage)) {
            const auto& dense = std::get<DenseStorage>(storage);
            DenseStorage new_dense;
            new_dense.reserve(dense.size());
            for(const auto& e : dense) new_dense.push_back(e->clone());
            return std::make_shared<MatrixNode>(rows, cols, std::move(new_dense));
        } else {
            const auto& sparse = std::get<SparseStorage>(storage);
            SparseStorage new_sparse;
            for(const auto& [idx, node] : sparse) {
                new_sparse[idx] = node->clone();
            }
            return std::make_shared<MatrixNode>(rows, cols, std::move(new_sparse));
        }
    }

    /**
     * @brief 获取指定位置的矩阵元素。
     * @param r 行索引（从 0 开始）
     * @param c 列索引（从 0 开始）
     * @return 对应节点，越界返回 nullptr，稀疏缺失返回零节点
     */
    std::shared_ptr<SymbolicNode> get(size_t r, size_t c) const {
        if (r >= rows || c >= cols) return nullptr;

        size_t idx = r * cols + c;
        if (std::holds_alternative<DenseStorage>(storage)) {
            return std::get<DenseStorage>(storage)[idx];
        } else {
            const auto& sparse = std::get<SparseStorage>(storage);
            auto it = sparse.find(idx);
            if (it != sparse.end()) return it->second;
            return std::make_shared<NumberNode>(0.0);
        }
    }

    /** @brief 判断当前矩阵是否使用稀疏存储。 */
    bool is_sparse() const { return std::holds_alternative<SparseStorage>(storage); }

private:

};

inline std::variant<MatrixNode::DenseStorage, MatrixNode::SparseStorage> MatrixNode::create_storage_from_grid(
    const std::vector<std::vector<std::shared_ptr<SymbolicNode>>>& grid, size_t total_elements, size_t ncols) {

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
    /** @brief 关系运算符类型 */
    enum class Op {
        EQ,   ///< 等于
        NEQ,  ///< 不等于
        LT,   ///< 小于
        GT,   ///< 大于
        LEQ,  ///< 小于等于
        GEQ   ///< 大于等于
    };

    std::shared_ptr<SymbolicNode> left;  ///< 左操作数
    std::shared_ptr<SymbolicNode> right; ///< 右操作数
    Op op;                               ///< 关系运算符

    RelationalNode(std::shared_ptr<SymbolicNode> l, std::shared_ptr<SymbolicNode> r, Op o)
        : left(std::move(l)), right(std::move(r)), op(o) {}

    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<SymbolicNode> clone() const override {
        return std::make_shared<RelationalNode>(left->clone(), right->clone(), op);
    }

    int type_priority() const override { return 100; }

    std::size_t compute_hash() const override {
        std::size_t h = std::hash<int>{}((int)op);
        hash_combine(h, left->hash());
        hash_combine(h, right->hash());
        return h;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const RelationalNode&>(other);
        if (op != o.op) return (int)op < (int)o.op ? -1 : 1;
        int cmp = left->compare(*o.left);
        if (cmp != 0) return cmp;
        return right->compare(*o.right);
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
        And, ///< 逻辑与
        Or   ///< 逻辑或
    };

    std::shared_ptr<SymbolicNode> left;  ///< 左操作数
    std::shared_ptr<SymbolicNode> right; ///< 右操作数
    Op op;                               ///< 逻辑运算符

    LogicalNode(std::shared_ptr<SymbolicNode> l, std::shared_ptr<SymbolicNode> r, Op o)
        : left(std::move(l)), right(std::move(r)), op(o) {}

    void accept(SymbolicVisitor& visitor) override { SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<SymbolicNode> clone() const override {
        return std::make_shared<LogicalNode>(left->clone(), right->clone(), op);
    }

    int type_priority() const override { return 101; }

    std::size_t compute_hash() const override {
        std::size_t h = std::hash<int>{}((int)op);
        hash_combine(h, left->hash());
        hash_combine(h, right->hash());
        return h;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const LogicalNode&>(other);
        if (op != o.op) return (int)op < (int)o.op ? -1 : 1;
        int cmp = left->compare(*o.left);
        if (cmp != 0) return cmp;
        return right->compare(*o.right);
    }

    /**
     * @brief 将逻辑运算符转换为字符串表示。
     * @param op 逻辑运算符
     * @return "And" 或 "Or"
     */
    static std::string op_to_string(Op op) {
        switch(op) {
            case Op::And: return "And";
            case Op::Or: return "Or";
            default: return "?";
        }
    }
};

inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_number(const ::BigInt& v) { return std::make_shared<NumberNode>(v); }
inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_number(const ::Rational& v) { return std::make_shared<NumberNode>(v); }
inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_number(lmmc_real_t v) { return std::make_shared<NumberNode>(v); }
inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_variable(const std::string& name) { return std::make_shared<VariableNode>(name); }

inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_add(std::vector<std::shared_ptr<SymbolicNode>> ops) {
    if (ops.empty()) return create_number(0.0);
    if (ops.size() == 1) return ops[0];

    std::vector<std::shared_ptr<SymbolicNode>> flat_ops;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_zero()) continue;
        if (auto add = std::dynamic_pointer_cast<AddNode>(op)) {
            flat_ops.insert(flat_ops.end(), add->operands.begin(), add->operands.end());
        } else {
            flat_ops.push_back(op);
        }
    }

    if (flat_ops.empty()) return create_number(0.0);
    if (flat_ops.size() == 1) return flat_ops[0];
    return std::make_shared<AddNode>(std::move(flat_ops));
}

inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_multiply(std::vector<std::shared_ptr<SymbolicNode>> ops) {
    if (ops.empty()) return create_number(1.0);

    for (const auto& op : ops) {
        if (op->is_zero()) return op;
    }

    if (ops.size() == 1) return ops[0];

    std::vector<std::shared_ptr<SymbolicNode>> flat_ops;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_one()) continue;
        if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {
            flat_ops.insert(flat_ops.end(), mul->operands.begin(), mul->operands.end());
        } else {
            flat_ops.push_back(op);
        }
    }

    if (flat_ops.empty()) return create_number(1.0);
    if (flat_ops.size() == 1) return flat_ops[0];
    return std::make_shared<MultiplyNode>(std::move(flat_ops));
}

inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_power(std::shared_ptr<SymbolicNode> base, std::shared_ptr<SymbolicNode> exponent) {
    if (!base || !exponent) return nullptr;
    if (exponent->is_zero()) return create_number(1.0);
    if (exponent->is_one()) return base;
    if (base->is_zero()) return create_number(0.0);
    if (base->is_one()) return create_number(1.0);
    return std::make_shared<PowerNode>(std::move(base), std::move(exponent));
}
