#include "poly_utils.hpp"


namespace LMCAS {
namespace {

using ResultantMatrix = std::vector<std::vector<Rational>>;

ResultantMatrix sylvester_matrix(const Polynomial<Rational>& left,
                                 const Polynomial<Rational>& right) {
    const int left_degree = left.degree();
    const int right_degree = right.degree();
    const int size = left_degree + right_degree;
    ResultantMatrix matrix(size, std::vector<Rational>(size, Rational(0)));
    for (int row = 0; row < right_degree; ++row) {
        for (int column = 0; column <= left_degree; ++column) {
            matrix[row][row + column] = left.coeffs[left_degree - column];
        }
    }
    for (int row = 0; row < left_degree; ++row) {
        for (int column = 0; column <= right_degree; ++column) {
            matrix[right_degree + row][row + column] =
                right.coeffs[right_degree - column];
        }
    }
    return matrix;
}

void eliminate_resultant_column(ResultantMatrix& matrix, int pivot,
                                const Rational& previous) {
    const int size = static_cast<int>(matrix.size());
    for (int row = pivot + 1; row < size; ++row) {
        for (int column = pivot + 1; column < size; ++column) {
            matrix[row][column] =
                (matrix[row][column] * matrix[pivot][pivot] -
                 matrix[row][pivot] * matrix[pivot][column]) / previous;
        }
        matrix[row][pivot] = Rational(0);
    }
}

Rational resultant_determinant(ResultantMatrix& matrix) {
    Rational previous(1);
    int sign = 1;
    const int size = static_cast<int>(matrix.size());
    for (int pivot = 0; pivot < size; ++pivot) {
        if (matrix[pivot][pivot] == Rational(0)) {
            int replacement = pivot + 1;
            while (replacement < size &&
                   matrix[replacement][pivot] == Rational(0)) {
                ++replacement;
            }
            if (replacement == size) { return Rational(0); }
            std::swap(matrix[pivot], matrix[replacement]);
            sign = -sign;
        }
        eliminate_resultant_column(matrix, pivot, previous);
        previous = matrix[pivot][pivot];
    }
    auto determinant = matrix.back().back();
    if (sign < 0) { determinant = Rational(0) - determinant; }
    return determinant;
}

}

std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_resultant(
    const std::shared_ptr<SymbolicExpr>& left,
    const std::shared_ptr<SymbolicExpr>& right,
    const std::string& variable) {
    if (!left || !right || variable.empty()) { return nullptr; }

    try {
        auto left_conversion = LMCAS::symbolic_to_poly<Rational>(left, variable);
        auto right_conversion = LMCAS::symbolic_to_poly<Rational>(right, variable);
        if (!left_conversion || !right_conversion) return nullptr;
        const auto& left_poly = left_conversion.value();
        const auto& right_poly = right_conversion.value();
        const int left_degree = left_poly.degree();
        const int right_degree = right_poly.degree();
        if (left_degree < 0 || right_degree < 0) { return SymbolicExpr::number(0); }
        if (left_degree == 0 && right_degree == 0) { return SymbolicExpr::number(1); }

        auto matrix = sylvester_matrix(left_poly, right_poly);
        return SymbolicExpr::number(resultant_determinant(matrix));
    } catch (const std::exception&) {
        return nullptr;
    }
}

} // namespace LMCAS
