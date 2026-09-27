#include "internal/symbolic_ast.hpp"

namespace LMCAS {
namespace {

int compare_matrix_elements(const std::shared_ptr<const SymbolicNode>& left,
                            const std::shared_ptr<const SymbolicNode>& right) {
    const bool left_zero = !left || left->is_zero();
    const bool right_zero = !right || right->is_zero();
    if (left_zero && right_zero) {
        return 0;
    }
    if (left_zero) {
        return -1;
    }
    if (right_zero) {
        return 1;
    }
    return left->compare(*right);
}

}

int MatrixNode::compare_same_type(const SymbolicNode& other) const {
    const auto& matrix = static_cast<const MatrixNode&>(other);
    if (rows_ != matrix.rows_) {
        return rows_ < matrix.rows_ ? -1 : 1;
    }
    if (cols_ != matrix.cols_) {
        return cols_ < matrix.cols_ ? -1 : 1;
    }
    for (size_t row = 0; row < rows_; ++row) {
        for (size_t column = 0; column < cols_; ++column) {
            const int comparison = compare_matrix_elements(get(row, column), matrix.get(row, column));
            if (comparison != 0) {
                return comparison;
            }
        }
    }
    return 0;
}

int IntervalNode::compare_same_type(const SymbolicNode& other) const {
    const auto& interval = static_cast<const IntervalNode&>(other);
    int comparison = lower_->compare(*interval.lower_);
    if (comparison != 0) {
        return comparison;
    }
    comparison = upper_->compare(*interval.upper_);
    if (comparison != 0) {
        return comparison;
    }
    if (lower_closed_ != interval.lower_closed_) {
        return lower_closed_ ? 1 : -1;
    }
    if (upper_closed_ != interval.upper_closed_) {
        return upper_closed_ ? 1 : -1;
    }
    return 0;
}

}
