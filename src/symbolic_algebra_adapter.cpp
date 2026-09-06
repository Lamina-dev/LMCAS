#include "poly_utils.hpp"


namespace LMCAS {

std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_resultant(
    const std::shared_ptr<SymbolicExpr>& left,
    const std::shared_ptr<SymbolicExpr>& right,
    const std::string& variable) {
    if (!left || !right || variable.empty()) return nullptr;

    try {
        auto left_poly = LMCAS::symbolic_to_poly<Rational>(left, variable);
        auto right_poly = LMCAS::symbolic_to_poly<Rational>(right, variable);
        const int left_degree = left_poly.degree();
        const int right_degree = right_poly.degree();
        if (left_degree < 0 || right_degree < 0) return SymbolicExpr::number(0);
        if (left_degree == 0 && right_degree == 0) return SymbolicExpr::number(1);

        const int size = left_degree + right_degree;
        std::vector<std::vector<Rational>> matrix(
            size, std::vector<Rational>(size, Rational(0)));
        for (int row = 0; row < right_degree; ++row) {
            for (int column = 0; column <= left_degree; ++column) {
                matrix[row][row + column] =
                    left_poly.coeffs[left_degree - column];
            }
        }
        for (int row = 0; row < left_degree; ++row) {
            for (int column = 0; column <= right_degree; ++column) {
                matrix[right_degree + row][row + column] =
                    right_poly.coeffs[right_degree - column];
            }
        }

        Rational previous(1);
        int sign = 1;
        for (int pivot = 0; pivot < size; ++pivot) {
            if (matrix[pivot][pivot] == Rational(0)) {
                int replacement = pivot + 1;
                while (replacement < size &&
                       matrix[replacement][pivot] == Rational(0)) {
                    ++replacement;
                }
                if (replacement == size) return SymbolicExpr::number(0);
                std::swap(matrix[pivot], matrix[replacement]);
                sign = -sign;
            }
            for (int row = pivot + 1; row < size; ++row) {
                for (int column = pivot + 1; column < size; ++column) {
                    matrix[row][column] =
                        (matrix[row][column] * matrix[pivot][pivot] -
                         matrix[row][pivot] * matrix[pivot][column]) /
                        previous;
                }
                matrix[row][pivot] = Rational(0);
            }
            previous = matrix[pivot][pivot];
        }

        auto determinant = matrix.back().back();
        if (sign < 0) determinant = Rational(0) - determinant;
        return SymbolicExpr::number(determinant);
    } catch (const std::exception&) {
        return nullptr;
    }
}

} // namespace LMCAS
