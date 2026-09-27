#pragma once
#include "internal/symbolic_ast/functions.hpp"
#include <map>

namespace LMCAS {

/**
 * @brief 矩阵节点，支持稠密和稀疏两种存储方式。
 *
 * 当非零元素占比低于 20% 时自动使用稀疏存储，否则使用稠密存储。
 */
class MatrixNode : public SymbolicNode {
public:
    using DenseStorage = std::vector<std::shared_ptr<const SymbolicNode>>;   /**< 稠密存储，按行优先展开。 */
    using SparseStorage = std::map<size_t, std::shared_ptr<const SymbolicNode>>; /**< 稀疏存储，索引映射到节点。 */

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

    int compare_same_type(const SymbolicNode& other) const override;

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

}
