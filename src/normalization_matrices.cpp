#include "internal/visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {
namespace {
std::shared_ptr<const SymbolicNode> matrix_element(const MatrixNode& matrix, std::size_t index) {
    std::shared_ptr<const SymbolicNode> value;
    if (const auto* dense = std::get_if<MatrixNode::DenseStorage>(&matrix.storage())) {
        if (index < dense->size()) {
            value = (*dense)[index];
        }
    } else {
        const auto& sparse = std::get<MatrixNode::SparseStorage>(matrix.storage());
        const auto found = sparse.find(index);
        if (found != sparse.end()) {
            value = found->second;
        }
    }
    return value ? value : detail::make_node<NumberNode>(BigInt(0));
}

bool matching_matrix_shapes(const std::vector<std::shared_ptr<const SymbolicNode>>& operands,
                            std::size_t rows, std::size_t cols) {
    for (const auto& operand : operands) {
        const auto matrix = std::dynamic_pointer_cast<const MatrixNode>(operand);
        if (!matrix || matrix->rows() != rows || matrix->cols() != cols) {
            return false;
        }
    }
    return true;
}
}

bool NormalizationVisitor::normalize_matrix_sum(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands) {
    if (operands.empty()) {
        return false;
    }
    const auto first = std::dynamic_pointer_cast<const MatrixNode>(operands.front());
    if (!first) {
        return false;
    }
    const auto rows = first->rows();
    const auto cols = first->cols();
    if (!matching_matrix_shapes(operands, rows, cols)) {
        return false;
    }
    MatrixNode::DenseStorage elements;
    const auto count = normalization_check_count(rewrite_budget(), rows, cols);
    std::size_t nodes = 1;
    elements.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        std::vector<std::shared_ptr<const SymbolicNode>> terms;
        std::size_t term_nodes = 0;
        normalization_check_count(rewrite_budget(), operands.size(), 1, false);
        for (const auto& operand : operands) {
            const auto matrix = std::dynamic_pointer_cast<const MatrixNode>(operand);
            normalization_append(rewrite_budget(), term_nodes, terms, matrix_element(*matrix, i));
        }
        normalization_check_arithmetic<AddNode>(rewrite_budget(), terms, term_nodes);
        detail::make_node<AddNode>(terms)->accept(*this);
        normalization_append(rewrite_budget(), nodes, elements, result);
    }
    set_result(detail::make_node<MatrixNode>(rows, cols, elements));
    return true;
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::matrix_product_element(
    const MatrixNode& left, const MatrixNode& right, std::size_t row, std::size_t col) {
    const auto* a = std::get_if<MatrixNode::DenseStorage>(&left.storage());
    const auto* b = std::get_if<MatrixNode::DenseStorage>(&right.storage());
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    std::size_t nodes = 0;
    normalization_check_count(rewrite_budget(), left.cols(), 1, false);
    terms.reserve(left.cols());
    for (std::size_t k = 0; k < left.cols(); ++k) {
        auto first = a ? (*a)[row * left.cols() + k] : matrix_element(left, row * left.cols() + k);
        auto second = b ? (*b)[k * right.cols() + col] : matrix_element(right, k * right.cols() + col);
        normalization_append(rewrite_budget(), nodes, terms,
                             make_normalized_multiply_node({first, second}, rewrite_budget()));
    }
    NormalizationVisitor visitor(context_, facts_, domain_, rewrite_budget());
    normalization_check_arithmetic<AddNode>(rewrite_budget(), terms, nodes);
    detail::make_node<AddNode>(terms)->accept(visitor);
    return visitor.get_result();
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::fuse_matrix_product(
    const std::shared_ptr<const SymbolicNode>& left,
    const std::shared_ptr<const SymbolicNode>& right) {
    const auto a = std::dynamic_pointer_cast<const MatrixNode>(left);
    const auto b = std::dynamic_pointer_cast<const MatrixNode>(right);
    if (!a || !b) {
        return nullptr;
    }
    if (a->cols() != b->rows()) {
        return nullptr;
    }
    MatrixNode::DenseStorage elements;
    const auto count = normalization_check_count(rewrite_budget(), a->rows(), b->cols());
    std::size_t nodes = 1;
    elements.reserve(count);
    for (std::size_t row = 0; row < a->rows(); ++row) {
        for (std::size_t col = 0; col < b->cols(); ++col) {
            normalization_append(rewrite_budget(), nodes, elements,
                                 matrix_product_element(*a, *b, row, col));
        }
    }
    return detail::make_node<MatrixNode>(a->rows(), b->cols(), elements);
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::scale_matrix(
    const MatrixNode& matrix, const std::shared_ptr<const NumberNode>& scalar) {
    MatrixNode::DenseStorage elements;
    const auto count = normalization_check_count(rewrite_budget(), matrix.rows(), matrix.cols());
    std::size_t nodes = 1;
    elements.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        auto element = matrix_element(matrix, i);
        normalization_check_arithmetic<MultiplyNode>(rewrite_budget(), 0, scalar, element);
        detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{scalar, element})->accept(*this);
        normalization_append(rewrite_budget(), nodes, elements, result);
    }
    return detail::make_node<MatrixNode>(matrix.rows(), matrix.cols(), std::move(elements));
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::normalize_matrix_product(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands) {
    std::vector<std::shared_ptr<const SymbolicNode>> nonscalar;
    std::size_t nonscalar_nodes = 0;
    auto scalar = detail::make_node<NumberNode>(BigInt(1));
    for (const auto& operand : operands) {
        if (const auto number = std::dynamic_pointer_cast<const NumberNode>(operand)) {
            scalar = multiply_numbers(scalar, number);
        } else {
            normalization_append(rewrite_budget(), nonscalar_nodes, nonscalar, operand);
        }
    }
    std::vector<std::shared_ptr<const SymbolicNode>> fused;
    std::size_t fused_nodes = 0;
    std::size_t last_nodes = 0;
    if (!nonscalar.empty()) {
        normalization_append(rewrite_budget(), fused_nodes, fused, nonscalar.front());
        last_nodes = fused_nodes;
    }
    for (std::size_t i = 1; i < nonscalar.size(); ++i) {
        if (auto product = fuse_matrix_product(fused.back(), nonscalar[i])) {
            if (rewrite_budget()) {
                const auto replacement_nodes = rewrite_budget()->measure(product);
                fused_nodes = rewrite_budget()->append_size(fused_nodes - last_nodes, replacement_nodes);
                last_nodes = replacement_nodes;
            }
            fused.back() = std::move(product);
        } else {
            const auto previous_nodes = fused_nodes;
            normalization_append(rewrite_budget(), fused_nodes, fused, nonscalar[i]);
            last_nodes = fused_nodes - previous_nodes;
        }
    }
    if (fused.size() == 1 && !scalar->is_one()) {
        if (const auto matrix = std::dynamic_pointer_cast<const MatrixNode>(fused.front())) {
            return scale_matrix(*matrix, scalar);
        }
    }
    if (!scalar->is_one()) {
        if (rewrite_budget()) { fused_nodes = rewrite_budget()->append_size(fused_nodes, 1); }
        fused.insert(fused.begin(), scalar);
    }
    return make_normalized_multiply_node(fused, rewrite_budget());
}
}
