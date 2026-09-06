/** @file internal/symbolic_ast/functions.hpp */
#pragma once
#include "arithmetic.hpp"

namespace LMCAS {

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
        Atan2,                               ///< 双参数反正切
        Infinity,                            ///< 无穷大
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
    LMCAS_AST_NODE_FACTORY_FRIEND;

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
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (const auto& arg : arguments_) new_args.push_back(arg->clone());
        return LMCAS::detail::make_node<FunctionNode>(type_, std::move(new_args));
    }
};

/** A symbolic call whose function name has no built-in semantics. */
class UninterpretedFunctionNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;
    const std::string name_;
    const std::vector<std::shared_ptr<const SymbolicNode>> arguments_;

    UninterpretedFunctionNode(
        std::string name,
        std::vector<std::shared_ptr<const SymbolicNode>> arguments)
        : name_(std::move(name)), arguments_(std::move(arguments)) {
        if (name_.empty()) {
            throw std::invalid_argument("uninterpreted function name cannot be empty");
        }
        for (const auto& argument : arguments_) {
            if (!argument) {
                throw std::invalid_argument("uninterpreted function arguments cannot be null");
            }
        }
    }

public:
    const std::string& name() const noexcept { return name_; }
    const auto& arguments() const noexcept { return arguments_; }
    int type_priority() const override { return 111; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        hash_combine(seed, std::hash<std::string>{}(name_));
        for (const auto& argument : arguments_) hash_combine(seed, argument->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& function = static_cast<const UninterpretedFunctionNode&>(other);
        int comparison = name_.compare(function.name_);
        if (comparison != 0) return comparison;
        const auto count = std::min(arguments_.size(), function.arguments_.size());
        for (std::size_t index = 0; index < count; ++index) {
            comparison = arguments_[index]->compare(*function.arguments_[index]);
            if (comparison != 0) return comparison;
        }
        if (arguments_.size() == function.arguments_.size()) return 0;
        return arguments_.size() < function.arguments_.size() ? -1 : 1;
    }

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override {
        LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }

    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> arguments;
        arguments.reserve(arguments_.size());
        for (const auto& argument : arguments_) arguments.push_back(argument->clone());
        return LMCAS::detail::make_node<UninterpretedFunctionNode>(
            name_, std::move(arguments));
    }
};

} // namespace LMCAS


#include <map>

namespace LMCAS {

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
    LMCAS_AST_NODE_FACTORY_FRIEND;

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
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        if (std::holds_alternative<DenseStorage>(storage_)) {
            const auto& dense = std::get<DenseStorage>(storage_);
            DenseStorage new_dense;
            new_dense.reserve(dense.size());
            for (const auto& e : dense) {
                new_dense.push_back(e->clone());
            }
            return LMCAS::detail::make_node<MatrixNode>(rows_, cols_, std::move(new_dense));
        } else {
            const auto& sparse = std::get<SparseStorage>(storage_);
            SparseStorage new_sparse;
            for(const auto& [idx, node] : sparse) {
                new_sparse[idx] = node->clone();
            }
            return LMCAS::detail::make_node<MatrixNode>(rows_, cols_, std::move(new_sparse));
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
            static const std::shared_ptr<const SymbolicNode> zero =
                LMCAS::detail::make_node<NumberNode>(BigInt(0));
            return zero;
        }
    }

    /** @brief 判断当前矩阵是否使用稀疏存储。 */
    bool is_sparse() const { return std::holds_alternative<SparseStorage>(storage_); }

private:

};

inline bool MultiplyNode::contains_matrix(const SymbolicNode& node) {
    if (dynamic_cast<const MatrixNode*>(&node)) return true;
    if (const auto* power = dynamic_cast<const PowerNode*>(&node)) {
        return contains_matrix(*power->base());
    }
    const std::vector<std::shared_ptr<const SymbolicNode>>* children = nullptr;
    if (const auto* sum = dynamic_cast<const AddNode*>(&node)) {
        children = &sum->operands();
    } else if (const auto* product = dynamic_cast<const MultiplyNode*>(&node)) {
        children = &product->operands();
    } else if (const auto* function = dynamic_cast<const FunctionNode*>(&node)) {
        children = &function->arguments();
    } else if (const auto* function = dynamic_cast<const UninterpretedFunctionNode*>(&node)) {
        children = &function->arguments();
    }
    return children && std::any_of(
        children->begin(), children->end(),
        [](const auto& child) { return contains_matrix(*child); });
}

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

} // namespace LMCAS
