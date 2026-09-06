#include "internal/exact_matrix.hpp"

#include "assumption_context.hpp"
#include "symbolic_ast.hpp"

#include <algorithm>
#include <limits>
#include <new>
#include <optional>
#include <set>

namespace LMCAS::detail {
namespace {

Result<void> validate_matrix(
    const ExactMatrixData& matrix,
    const std::string& operation) {
    if (matrix.rows == 0 || matrix.cols == 0) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "exact matrix dimensions must be non-zero", operation);
    }
    if (matrix.rows > std::numeric_limits<std::size_t>::max() / matrix.cols ||
        matrix.entries.size() != matrix.rows * matrix.cols) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "exact matrix storage does not match its dimensions", operation);
    }
    for (const auto& entry : matrix.entries) {
        if (!entry || !node(entry)) {
            return Result<void>::failure(
                CasErrc::InvalidArgument,
                "exact matrix entries cannot be null", operation);
        }
    }
    return Result<void>::success();
}

class ApproximateNumberDetector final : public RecursiveSymbolicVisitor {
public:
    bool found = false;

    void visit(const NumberNode& number) override {
        found = found || std::holds_alternative<lmmc_real_t>(number.value());
    }
};

bool contains_approximate_number(const ExprPtr& expression) {
    if (!expression || !node(expression)) return false;
    ApproximateNumberDetector detector;
    node(expression)->accept(detector);
    return detector.found;
}

std::optional<Rational> exact_rational(const ExprPtr& expression) {
    if (!expression || !node(expression)) return std::nullopt;
    auto number = std::dynamic_pointer_cast<const NumberNode>(node(expression));
    if (!number) return std::nullopt;
    if (std::holds_alternative<BigInt>(number->value())) {
        return Rational(std::get<BigInt>(number->value()));
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value());
    }
    return std::nullopt;
}

ExprPtr simplify_expr(ExprPtr expression) {
    if (!expression) return expression;
    auto simplified = expression->simplify();
    return simplified ? simplified : expression;
}

ExprPtr exact_add(const ExprPtr& left, const ExprPtr& right) {
    return simplify_expr(SymbolicExpr::add(left, right));
}

ExprPtr exact_negate(const ExprPtr& expression) {
    return simplify_expr(SymbolicExpr::multiply(
        SymbolicExpr::number(-1), expression));
}

ExprPtr exact_subtract(const ExprPtr& left, const ExprPtr& right) {
    return exact_add(left, exact_negate(right));
}

ExprPtr exact_multiply(const ExprPtr& left, const ExprPtr& right) {
    return simplify_expr(SymbolicExpr::multiply(left, right));
}

ExprPtr exact_divide(const ExprPtr& numerator, const ExprPtr& denominator) {
    return simplify_expr(SymbolicExpr::divide(numerator, denominator));
}

void swap_rows(ExactMatrixData& matrix, std::size_t first, std::size_t second) {
    if (first == second) return;
    for (std::size_t column = 0; column < matrix.cols; ++column) {
        std::swap(matrix.at(first, column), matrix.at(second, column));
    }
}

Result<std::optional<std::size_t>> choose_pivot(
    const ExactMatrixData& matrix,
    std::size_t first_row,
    std::size_t column,
    ComputationContext& context,
    const std::string& operation) {
    bool saw_unknown = false;
    for (std::size_t row = first_row; row < matrix.rows; ++row) {
        auto proof = classify_exact_zero(matrix.at(row, column), context, operation);
        if (!proof) {
            return Result<std::optional<std::size_t>>::failure(proof.error());
        }
        if (proof.value() == ZeroProof::NonZero) {
            return Result<std::optional<std::size_t>>::success(row);
        }
        saw_unknown = saw_unknown || proof.value() == ZeroProof::Unknown;
    }
    if (saw_unknown) {
        return Result<std::optional<std::size_t>>::failure(
            CasErrc::Inconclusive,
            "exact pivot is not provably zero or non-zero", operation);
    }
    return Result<std::optional<std::size_t>>::success(std::nullopt);
}

Result<Rational> rational_determinant_value_impl(
    std::size_t dimension,
    std::vector<Rational> values,
    ComputationContext& context,
    const std::string& operation) {
    if (dimension == 0 ||
        dimension > std::numeric_limits<std::size_t>::max() / dimension ||
        values.size() != dimension * dimension) {
        return Result<Rational>::failure(
            CasErrc::InvalidArgument,
            "rational determinant storage does not match its dimension",
            operation);
    }
    Rational sign(1);
    Rational previous(1);
    if (dimension == 1) return Result<Rational>::success(values[0]);
    for (std::size_t pivot_column = 0;
         pivot_column + 1 < dimension; ++pivot_column) {
        auto budget = context.consume_steps(1, operation);
        if (!budget) return Result<Rational>::failure(budget.error());
        std::size_t pivot_row = pivot_column;
        while (pivot_row < dimension &&
               values[pivot_row * dimension + pivot_column] == Rational(0)) {
            ++pivot_row;
        }
        if (pivot_row == dimension) {
            return Result<Rational>::success(Rational(0));
        }
        if (pivot_row != pivot_column) {
            for (std::size_t column = 0; column < dimension; ++column) {
                std::swap(values[pivot_row * dimension + column],
                          values[pivot_column * dimension + column]);
            }
            sign = -sign;
        }
        const Rational pivot =
            values[pivot_column * dimension + pivot_column];
        for (std::size_t row = pivot_column + 1; row < dimension; ++row) {
            for (std::size_t column = pivot_column + 1;
                 column < dimension; ++column) {
                values[row * dimension + column] =
                    (values[row * dimension + column] * pivot -
                     values[row * dimension + pivot_column] *
                         values[pivot_column * dimension + column]) /
                    previous;
            }
            values[row * dimension + pivot_column] = Rational(0);
        }
        previous = pivot;
    }
    return Result<Rational>::success(
        sign * values[(dimension - 1) * dimension + (dimension - 1)]);
}

Result<ExprPtr> rational_determinant(
    const ExactMatrixData& input,
    ComputationContext& context,
    const std::string& operation) {
    std::vector<Rational> values;
    values.reserve(input.entries.size());
    for (const auto& entry : input.entries) {
        auto value = exact_rational(entry);
        if (!value) {
            return Result<ExprPtr>::failure(
                CasErrc::Inconclusive,
                "matrix is not entirely rational",
                operation);
        }
        values.push_back(*value);
    }
    auto determinant = rational_determinant_value_impl(
        input.rows, std::move(values), context, operation);
    if (!determinant) return Result<ExprPtr>::failure(determinant.error());
    return Result<ExprPtr>::success(
        SymbolicExpr::number(determinant.value()));
}

Result<ExprPtr> determinant_cofactor(
    const ExactMatrixData& matrix,
    ComputationContext& context,
    const std::string& operation) {
    auto budget = context.consume_steps(1, operation);
    if (!budget) return Result<ExprPtr>::failure(budget.error());
    if (matrix.rows == 1) {
        return Result<ExprPtr>::success(matrix.entries[0]);
    }
    auto sum = SymbolicExpr::number(0);
    for (std::size_t column = 0; column < matrix.cols; ++column) {
        if (auto coefficient = exact_rational(matrix.at(0, column));
            coefficient && *coefficient == Rational(0)) {
            continue;
        }
        ExactMatrixData minor{
            matrix.rows - 1, matrix.cols - 1, {}};
        minor.entries.reserve(minor.rows * minor.cols);
        for (std::size_t row = 1; row < matrix.rows; ++row) {
            for (std::size_t source = 0; source < matrix.cols; ++source) {
                if (source != column) {
                    minor.entries.push_back(matrix.at(row, source));
                }
            }
        }
        auto minor_determinant = determinant_cofactor(
            minor, context, operation);
        if (!minor_determinant) {
            return Result<ExprPtr>::failure(minor_determinant.error());
        }
        auto term = exact_multiply(
            matrix.at(0, column), minor_determinant.value());
        if (column % 2 != 0) term = exact_negate(term);
        sum = exact_add(sum, term);
    }
    return Result<ExprPtr>::success(simplify_expr(sum));
}

Result<ExactLinearSolution> rational_linear_solve(
    const ExactMatrixData& augmented,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation) {
    auto budget = context.consume_steps(
        augmented.rows * augmented.cols, operation);
    if (!budget) {
        return Result<ExactLinearSolution>::failure(budget.error());
    }
    std::vector<Rational> matrix;
    matrix.reserve(augmented.entries.size());
    for (const auto& entry : augmented.entries) {
        auto value = exact_rational(entry);
        if (!value) {
            return Result<ExactLinearSolution>::failure(
                CasErrc::Inconclusive,
                "rational specialization received a symbolic entry",
                operation);
        }
        matrix.push_back(*value);
    }
    const std::size_t columns = augmented.cols;
    std::vector<std::size_t> pivots;
    std::size_t pivot_row = 0;
    for (std::size_t column = 0;
         column < coefficient_columns && pivot_row < augmented.rows;
         ++column) {
        std::size_t selected = pivot_row;
        while (selected < augmented.rows &&
               matrix[selected * columns + column] == Rational(0)) {
            ++selected;
        }
        if (selected == augmented.rows) continue;
        if (selected != pivot_row) {
            for (std::size_t target = 0; target < columns; ++target) {
                std::swap(matrix[selected * columns + target],
                          matrix[pivot_row * columns + target]);
            }
        }
        const Rational pivot = matrix[pivot_row * columns + column];
        for (std::size_t target = 0; target < columns; ++target) {
            matrix[pivot_row * columns + target] =
                matrix[pivot_row * columns + target] / pivot;
        }
        for (std::size_t row = 0; row < augmented.rows; ++row) {
            if (row == pivot_row) continue;
            const Rational factor = matrix[row * columns + column];
            if (factor == Rational(0)) continue;
            for (std::size_t target = 0; target < columns; ++target) {
                matrix[row * columns + target] =
                    matrix[row * columns + target] -
                    factor * matrix[pivot_row * columns + target];
            }
        }
        pivots.push_back(column);
        ++pivot_row;
    }
    for (std::size_t row = 0; row < augmented.rows; ++row) {
        bool zero = true;
        for (std::size_t column = 0;
             column < coefficient_columns; ++column) {
            zero = zero &&
                matrix[row * columns + column] == Rational(0);
        }
        if (zero &&
            matrix[row * columns + coefficient_columns] != Rational(0)) {
            return Result<ExactLinearSolution>::success(
                InconsistentLinearSolution{});
        }
    }
    std::vector<std::optional<std::size_t>> pivot_rows(
        coefficient_columns);
    for (std::size_t row = 0; row < pivots.size(); ++row) {
        pivot_rows[pivots[row]] = row;
    }
    std::vector<std::size_t> free_columns;
    for (std::size_t column = 0; column < coefficient_columns; ++column) {
        if (!pivot_rows[column]) free_columns.push_back(column);
    }
    std::vector<ExprPtr> particular(
        coefficient_columns, SymbolicExpr::number(0));
    for (std::size_t column = 0; column < coefficient_columns; ++column) {
        if (pivot_rows[column]) {
            particular[column] = SymbolicExpr::number(
                matrix[*pivot_rows[column] * columns +
                       coefficient_columns]);
        }
    }
    if (free_columns.empty()) {
        return Result<ExactLinearSolution>::success(
            UniqueLinearSolution{std::move(particular)});
    }
    std::vector<std::vector<ExprPtr>> basis;
    for (const auto free_column : free_columns) {
        std::vector<ExprPtr> vector(
            coefficient_columns, SymbolicExpr::number(0));
        vector[free_column] = SymbolicExpr::number(1);
        for (std::size_t column = 0;
             column < coefficient_columns; ++column) {
            if (pivot_rows[column]) {
                vector[column] = SymbolicExpr::number(
                    -matrix[*pivot_rows[column] * columns + free_column]);
            }
        }
        basis.push_back(std::move(vector));
    }
    return Result<ExactLinearSolution>::success(
        ParametricLinearSolution{
            std::move(particular), std::move(basis),
            std::move(free_columns)});
}

} // namespace
Result<Rational> rational_determinant_exact(
    std::size_t dimension,
    std::vector<Rational> values,
    ComputationContext& context,
    const std::string& operation) {
    try {
        return rational_determinant_value_impl(
            dimension, std::move(values), context, operation);
    } catch (const std::bad_alloc&) {
        return Result<Rational>::failure(
            CasErrc::ResourceLimit,
            "rational determinant allocation failed",
            operation);
    } catch (const std::exception& error) {
        return Result<Rational>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}


Result<ZeroProof> classify_exact_zero(
    const ExprPtr& expression,
    ComputationContext& context,
    const std::string& operation) {
    auto access = context.consume_steps(1, operation);
    if (!access) return Result<ZeroProof>::failure(access.error());
    if (!expression || !node(expression)) {
        return Result<ZeroProof>::failure(
            CasErrc::InvalidArgument,
            "zero classification requires an expression", operation);
    }
    try {
        if (contains_approximate_number(expression)) {
            return Result<ZeroProof>::success(ZeroProof::Unknown);
        }
        auto simplified = simplify_expr(expression);
        if (auto value = exact_rational(simplified)) {
            return Result<ZeroProof>::success(
                *value == Rational(0) ? ZeroProof::Zero : ZeroProof::NonZero);
        }
        if (auto power =
                std::dynamic_pointer_cast<const PowerNode>(node(simplified))) {
            auto base = exact_rational(
                make_expression_ptr(power->base()));
            auto exponent = exact_rational(
                make_expression_ptr(power->exponent()));
            if (base && exponent) {
                if (*base > Rational(0)) {
                    return Result<ZeroProof>::success(ZeroProof::NonZero);
                }
                if (*base == Rational(0) && *exponent > Rational(0)) {
                    return Result<ZeroProof>::success(ZeroProof::Zero);
                }
            }
        }
        if (auto function =
                std::dynamic_pointer_cast<const FunctionNode>(
                    node(simplified));
            function && function->type() == FunctionNode::FuncType::Sqrt &&
            function->arguments().size() == 1) {
            auto argument = exact_rational(
                make_expression_ptr(function->arguments()[0]));
            if (argument && *argument >= Rational(0)) {
                return Result<ZeroProof>::success(
                    *argument == Rational(0)
                        ? ZeroProof::Zero : ZeroProof::NonZero);
            }
        }
        if (auto root =
                std::dynamic_pointer_cast<const RootOfNode>(
                    node(simplified));
            root && !root->exact_id().polynomial.coeffs.empty() &&
            root->exact_id().polynomial.coeffs[0] != Rational(0)) {
            return Result<ZeroProof>::success(ZeroProof::NonZero);
        }
        if (!context.assumptions()) {
            return Result<ZeroProof>::success(ZeroProof::Unknown);
        }
        auto assumption = context.assumptions()->is_nonzero_checked(*simplified);
        if (!assumption) return Result<ZeroProof>::failure(assumption.error());
        if (assumption.value() == Tribool::True) {
            return Result<ZeroProof>::success(ZeroProof::NonZero);
        }
        if (assumption.value() == Tribool::False) {
            return Result<ZeroProof>::success(ZeroProof::Zero);
        }
        return Result<ZeroProof>::success(ZeroProof::Unknown);
    } catch (const std::bad_alloc&) {
        return Result<ZeroProof>::failure(
            CasErrc::ResourceLimit,
            "zero classification allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<ZeroProof>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<EliminationResult> eliminate_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    EliminationForm form,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(input, operation);
    if (!valid) return Result<EliminationResult>::failure(valid.error());
    if (coefficient_columns == 0 || coefficient_columns > input.cols) {
        return Result<EliminationResult>::failure(
            CasErrc::InvalidArgument,
            "coefficient column count is outside matrix shape", operation);
    }
    try {
        EliminationResult result;
        result.matrix = std::move(input);
        std::size_t pivot_row = 0;
        for (std::size_t column = 0;
             column < coefficient_columns && pivot_row < result.matrix.rows;
             ++column) {
            auto pivot = choose_pivot(
                result.matrix, pivot_row, column, context, operation);
            if (!pivot) return Result<EliminationResult>::failure(pivot.error());
            if (!pivot.value()) continue;
            if (*pivot.value() != pivot_row) {
                swap_rows(result.matrix, *pivot.value(), pivot_row);
                ++result.row_swaps;
            }
            auto pivot_value = simplify_expr(result.matrix.at(pivot_row, column));
            result.pivot_columns.push_back(column);
            result.pivot_values.push_back(pivot_value);

            for (std::size_t row = pivot_row + 1; row < result.matrix.rows; ++row) {
                auto entry_proof = classify_exact_zero(
                    result.matrix.at(row, column), context, operation);
                if (!entry_proof) {
                    return Result<EliminationResult>::failure(entry_proof.error());
                }
                if (entry_proof.value() == ZeroProof::Zero) continue;
                auto factor = exact_divide(
                    result.matrix.at(row, column), pivot_value);
                for (std::size_t target = 0; target < result.matrix.cols; ++target) {
                    result.matrix.at(row, target) = exact_subtract(
                        result.matrix.at(row, target),
                        exact_multiply(factor,
                                       result.matrix.at(pivot_row, target)));
                }
            }
            ++pivot_row;
        }
        result.rank = result.pivot_columns.size();

        if (form == EliminationForm::ReducedRowEchelon) {
            for (std::size_t pivot_index = result.rank;
                 pivot_index-- > 0;) {
                const std::size_t row = pivot_index;
                const std::size_t column = result.pivot_columns[pivot_index];
                auto pivot_value = simplify_expr(result.matrix.at(row, column));
                for (std::size_t target = 0; target < result.matrix.cols; ++target) {
                    result.matrix.at(row, target) = exact_divide(
                        result.matrix.at(row, target), pivot_value);
                }
                for (std::size_t upper = 0; upper < row; ++upper) {
                    auto proof = classify_exact_zero(
                        result.matrix.at(upper, column), context, operation);
                    if (!proof) {
                        return Result<EliminationResult>::failure(proof.error());
                    }
                    if (proof.value() == ZeroProof::Zero) continue;
                    auto factor = result.matrix.at(upper, column);
                    for (std::size_t target = 0;
                         target < result.matrix.cols; ++target) {
                        result.matrix.at(upper, target) = exact_subtract(
                            result.matrix.at(upper, target),
                            exact_multiply(factor,
                                           result.matrix.at(row, target)));
                    }
                }
            }
        }
        return Result<EliminationResult>::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return Result<EliminationResult>::failure(
            CasErrc::ResourceLimit,
            "exact elimination allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<EliminationResult>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<ExprPtr> determinant_exact(
    const ExactMatrixData& input,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(input, operation);
    if (!valid) return Result<ExprPtr>::failure(valid.error());
    if (input.rows != input.cols) {
        return Result<ExprPtr>::failure(
            CasErrc::InvalidArgument,
            "determinant requires a square matrix", operation);
    }
    for (const auto& entry : input.entries) {
        if (contains_approximate_number(entry)) {
            return Result<ExprPtr>::failure(
                CasErrc::Inconclusive,
                "exact determinant rejects approximate entries", operation);
        }
    }
    bool all_rational = true;
    for (const auto& entry : input.entries) {
        all_rational = all_rational && exact_rational(entry).has_value();
    }
    if (all_rational) return rational_determinant(input, context, operation);

    try {
        ExactMatrixData matrix = input;
        const std::size_t n = matrix.rows;
        if (n == 1) return Result<ExprPtr>::success(matrix.entries[0]);
        auto previous = SymbolicExpr::number(1);
        std::size_t swaps = 0;
        for (std::size_t column = 0; column + 1 < n; ++column) {
            auto pivot = choose_pivot(matrix, column, column, context, operation);
            if (!pivot) {
                if (pivot.error().code == CasErrc::Inconclusive) {
                    return determinant_cofactor(input, context, operation);
                }
                return Result<ExprPtr>::failure(pivot.error());
            }
            if (!pivot.value()) {
                return Result<ExprPtr>::success(SymbolicExpr::number(0));
            }
            if (*pivot.value() != column) {
                swap_rows(matrix, *pivot.value(), column);
                ++swaps;
            }
            auto pivot_value = matrix.at(column, column);
            for (std::size_t row = column + 1; row < n; ++row) {
                for (std::size_t target = column + 1; target < n; ++target) {
                    auto numerator = exact_subtract(
                        exact_multiply(matrix.at(row, target), pivot_value),
                        exact_multiply(matrix.at(row, column),
                                       matrix.at(column, target)));
                    matrix.at(row, target) = exact_divide(numerator, previous);
                }
                matrix.at(row, column) = SymbolicExpr::number(0);
            }
            previous = pivot_value;
        }
        auto determinant = matrix.at(n - 1, n - 1);
        if (swaps % 2 != 0) determinant = exact_negate(determinant);
        return Result<ExprPtr>::success(simplify_expr(determinant));
    } catch (const std::bad_alloc&) {
        return Result<ExprPtr>::failure(
            CasErrc::ResourceLimit,
            "exact determinant allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<ExprPtr>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<ExactMatrixData> rref_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation) {
    auto eliminated = eliminate_exact(
        std::move(input), coefficient_columns,
        EliminationForm::ReducedRowEchelon, context, operation);
    if (!eliminated) return Result<ExactMatrixData>::failure(eliminated.error());
    return Result<ExactMatrixData>::success(
        std::move(eliminated.value().matrix));
}

Result<std::size_t> rank_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation) {
    auto eliminated = eliminate_exact(
        std::move(input), coefficient_columns,
        EliminationForm::Echelon, context, operation);
    if (!eliminated) return Result<std::size_t>::failure(eliminated.error());
    return Result<std::size_t>::success(eliminated.value().rank);
}

Result<ExactMatrixData> inverse_exact(
    const ExactMatrixData& input,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(input, operation);
    if (!valid) return Result<ExactMatrixData>::failure(valid.error());
    if (input.rows != input.cols) {
        return Result<ExactMatrixData>::failure(
            CasErrc::InvalidArgument,
            "inverse requires a square matrix", operation);
    }
    const std::size_t n = input.rows;
    ExactMatrixData augmented{n, 2 * n, {}};
    augmented.entries.reserve(n * 2 * n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
            augmented.entries.push_back(input.at(row, column));
        }
        for (std::size_t column = 0; column < n; ++column) {
            augmented.entries.push_back(
                SymbolicExpr::number(row == column ? 1 : 0));
        }
    }
    auto eliminated = eliminate_exact(
        std::move(augmented), n,
        EliminationForm::ReducedRowEchelon, context, operation);
    if (!eliminated) return Result<ExactMatrixData>::failure(eliminated.error());
    if (eliminated.value().rank != n) {
        return Result<ExactMatrixData>::failure(
            CasErrc::DomainError,
            "matrix is singular and has no inverse", operation);
    }
    ExactMatrixData inverse{n, n, {}};
    inverse.entries.reserve(n * n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
            inverse.entries.push_back(
                eliminated.value().matrix.at(row, n + column));
        }
    }
    return Result<ExactMatrixData>::success(std::move(inverse));
}

Result<std::vector<Rational>> solve_rational_unique(
    std::size_t rows,
    std::size_t coefficient_columns,
    std::vector<Rational> matrix,
    ComputationContext& context,
    const std::string& operation) {
    const std::size_t columns = coefficient_columns + 1;
    if (rows == 0 || coefficient_columns == 0 ||
        matrix.size() != rows * columns) {
        return Result<std::vector<Rational>>::failure(
            CasErrc::InvalidArgument,
            "rational linear-system storage is invalid", operation);
    }
    auto budget = context.consume_steps(rows * columns, operation);
    if (!budget) return Result<std::vector<Rational>>::failure(budget.error());
    std::vector<std::optional<std::size_t>> pivot_rows(
        coefficient_columns);
    std::size_t pivot_row = 0;
    for (std::size_t column = 0;
         column < coefficient_columns && pivot_row < rows; ++column) {
        std::size_t selected = pivot_row;
        while (selected < rows &&
               matrix[selected * columns + column] == Rational(0)) {
            ++selected;
        }
        if (selected == rows) continue;
        if (selected != pivot_row) {
            for (std::size_t target = 0; target < columns; ++target) {
                std::swap(matrix[selected * columns + target],
                          matrix[pivot_row * columns + target]);
            }
        }
        const Rational pivot = matrix[pivot_row * columns + column];
        for (std::size_t target = column; target < columns; ++target) {
            matrix[pivot_row * columns + target] =
                matrix[pivot_row * columns + target] / pivot;
        }
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivot_row) continue;
            const Rational factor = matrix[row * columns + column];
            if (factor == Rational(0)) continue;
            for (std::size_t target = column; target < columns; ++target) {
                matrix[row * columns + target] =
                    matrix[row * columns + target] -
                    factor * matrix[pivot_row * columns + target];
            }
        }
        pivot_rows[column] = pivot_row++;
    }
    for (std::size_t row = 0; row < rows; ++row) {
        bool zero = true;
        for (std::size_t column = 0;
             column < coefficient_columns; ++column) {
            zero = zero &&
                matrix[row * columns + column] == Rational(0);
        }
        if (zero && matrix[row * columns + coefficient_columns] != Rational(0)) {
            return Result<std::vector<Rational>>::failure(
                CasErrc::DomainError,
                "rational linear system is inconsistent", operation);
        }
    }
    std::vector<Rational> solution(
        coefficient_columns, Rational(0));
    for (std::size_t column = 0; column < coefficient_columns; ++column) {
        if (!pivot_rows[column]) {
            return Result<std::vector<Rational>>::failure(
                CasErrc::Inconclusive,
                "rational linear system does not have a unique solution",
                operation);
        }
        solution[column] =
            matrix[*pivot_rows[column] * columns + coefficient_columns];
    }
    return Result<std::vector<Rational>>::success(std::move(solution));
}
Result<ExactLinearSolution> solve_linear_exact(
    ExactMatrixData augmented,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(augmented, operation);

    if (!valid) return Result<ExactLinearSolution>::failure(valid.error());
    if (augmented.cols != coefficient_columns + 1) {
        return Result<ExactLinearSolution>::failure(
            CasErrc::InvalidArgument,
            "linear solve requires one augmented right-hand side", operation);
    }
    bool all_rational = true;
    for (const auto& entry : augmented.entries) {
        all_rational = all_rational && exact_rational(entry).has_value();
    }
    if (all_rational) {
        return rational_linear_solve(
            augmented, coefficient_columns, context, operation);
    }
    auto eliminated = eliminate_exact(
        std::move(augmented), coefficient_columns,
        EliminationForm::ReducedRowEchelon, context, operation);
    if (!eliminated) return Result<ExactLinearSolution>::failure(eliminated.error());
    auto& rref = eliminated.value().matrix;

    for (std::size_t row = 0; row < rref.rows; ++row) {
        bool coefficients_zero = true;
        for (std::size_t column = 0; column < coefficient_columns; ++column) {
            auto proof = classify_exact_zero(rref.at(row, column), context, operation);
            if (!proof) return Result<ExactLinearSolution>::failure(proof.error());
            if (proof.value() == ZeroProof::Unknown) {
                return Result<ExactLinearSolution>::failure(
                    CasErrc::Inconclusive,
                    "linear-system consistency is not provable", operation);
            }
            coefficients_zero = coefficients_zero && proof.value() == ZeroProof::Zero;
        }
        if (!coefficients_zero) continue;
        auto rhs = classify_exact_zero(
            rref.at(row, coefficient_columns), context, operation);
        if (!rhs) return Result<ExactLinearSolution>::failure(rhs.error());
        if (rhs.value() == ZeroProof::Unknown) {
            return Result<ExactLinearSolution>::failure(
                CasErrc::Inconclusive,
                "linear-system right-hand side is not provably zero", operation);
        }
        if (rhs.value() == ZeroProof::NonZero) {
            return Result<ExactLinearSolution>::success(
                InconsistentLinearSolution{});
        }
    }

    std::vector<std::optional<std::size_t>> pivot_row(coefficient_columns);
    for (std::size_t index = 0;
         index < eliminated.value().pivot_columns.size(); ++index) {
        const auto column = eliminated.value().pivot_columns[index];
        if (column < coefficient_columns) pivot_row[column] = index;
    }
    std::vector<std::size_t> free_columns;
    for (std::size_t column = 0; column < coefficient_columns; ++column) {
        if (!pivot_row[column]) free_columns.push_back(column);
    }
    std::vector<ExprPtr> particular(
        coefficient_columns, SymbolicExpr::number(0));
    for (std::size_t column = 0; column < coefficient_columns; ++column) {
        if (pivot_row[column]) {
            particular[column] = rref.at(*pivot_row[column], coefficient_columns);
        }
    }
    if (free_columns.empty()) {
        return Result<ExactLinearSolution>::success(
            UniqueLinearSolution{std::move(particular)});
    }

    std::vector<std::vector<ExprPtr>> basis;
    basis.reserve(free_columns.size());
    for (const auto free_column : free_columns) {
        std::vector<ExprPtr> vector(
            coefficient_columns, SymbolicExpr::number(0));
        vector[free_column] = SymbolicExpr::number(1);
        for (std::size_t column = 0; column < coefficient_columns; ++column) {
            if (pivot_row[column]) {
                vector[column] = exact_negate(
                    rref.at(*pivot_row[column], free_column));
            }
        }
        basis.push_back(std::move(vector));
    }
    return Result<ExactLinearSolution>::success(ParametricLinearSolution{
        std::move(particular), std::move(basis), std::move(free_columns)});
}

Result<std::vector<std::vector<ExprPtr>>> nullspace_exact(
    ExactMatrixData coefficients,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(coefficients, operation);
    if (!valid) {
        return Result<std::vector<std::vector<ExprPtr>>>::failure(valid.error());
    }
    const std::size_t columns = coefficients.cols;
    auto eliminated = eliminate_exact(
        std::move(coefficients), columns,
        EliminationForm::ReducedRowEchelon, context, operation);
    if (!eliminated) {
        return Result<std::vector<std::vector<ExprPtr>>>::failure(
            eliminated.error());
    }
    std::set<std::size_t> pivots(
        eliminated.value().pivot_columns.begin(),
        eliminated.value().pivot_columns.end());
    std::vector<std::vector<ExprPtr>> basis;
    for (std::size_t free_column = 0; free_column < columns; ++free_column) {
        if (pivots.count(free_column) != 0) continue;
        std::vector<ExprPtr> vector(columns, SymbolicExpr::number(0));
        vector[free_column] = SymbolicExpr::number(1);
        for (std::size_t row = 0;
             row < eliminated.value().pivot_columns.size(); ++row) {
            vector[eliminated.value().pivot_columns[row]] = exact_negate(
                eliminated.value().matrix.at(row, free_column));
        }
        basis.push_back(std::move(vector));
    }
    return Result<std::vector<std::vector<ExprPtr>>>::success(
        std::move(basis));
}

} // namespace LMCAS::detail
