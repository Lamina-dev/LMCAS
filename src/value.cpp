#include "value.hpp"
#include "numeric_evaluation.hpp"
#include "lmmc/numeric.h"
#include <cmath>
#include <stdexcept>

namespace LMCAS {

namespace {
void append_array(std::string& result, const std::vector<Value>& values) {
    result += "[";
    for (size_t index = 0; index < values.size(); ++index) {
        if (index) result += ", ";
        result += values[index].to_string();
    }
    result += "]";
}

std::string format_array(const std::vector<Value>& values) {
    std::string result;
    append_array(result, values);
    return result;
}

std::string format_matrix(const std::vector<std::vector<Value>>& matrix) {
    std::string result = "[";
    for (size_t row = 0; row < matrix.size(); ++row) {
        if (row) result += ", ";
        append_array(result, matrix[row]);
    }
    result += "]";
    return result;
}

std::string format_float(lmmc_real_t value) {
    std::string result = std::to_string(value);
    result.erase(result.find_last_not_of('0') + 1, std::string::npos);
    result.erase(result.find_last_not_of('.') + 1, std::string::npos);
    return result;
}

const char* format_infinity(const int* sign) {
    return sign && *sign > 0 ? "inf" : "-inf";
}

template <typename Exact>
LMCAS::Result<lmmc_real_t> exact_number_as_floating(const Exact& value) {
    const auto range_failure = [] {
        return LMCAS::Result<lmmc_real_t>::failure(
            LMCAS::CasErrc::NumericFailure,
            "exact numeric value is outside the floating-point range",
            "Value::as_number_checked");
    };
    try {
        const lmmc_real_t converted = value.to_double();
        if (!std::isfinite(converted) ||
            (converted == 0.0 && !value.is_zero())) {
            return range_failure();
        }
        return LMCAS::Result<lmmc_real_t>::success(converted);
    } catch (const std::overflow_error&) {
        return range_failure();
    } catch (const std::underflow_error&) {
        return range_failure();
    } catch (const std::out_of_range&) {
        return range_failure();
    } catch (const std::range_error&) {
        return range_failure();
    }
}
}

    Value::Value(lmmc_real_t f) : type(Type::Float), data(f) {
        if (std::isnan(f)) {
            throw std::invalid_argument("Value: floating-point value must not be NaN");
        }
        int res;
        lmmc_isinf(f, &res);
        if (res) {
            if (f < 0) res = -1;
            this->type = Type::Infinity;
            this->data = DataType(std::in_place_index<1>, res);
        }
    }

    Value::Value(const std::vector<Value>& arr) {
        bool is_matrix = !arr.empty() && arr[0].is_array();
        if (is_matrix) {
            std::vector<std::vector<Value>> matrix;
            for (const auto& row : arr) {
                if (row.is_array()) {
                    matrix.push_back(std::get<std::vector<Value>>(row.data));
                } else {
                    type = Type::Array;
                    data = arr;
                    return;
                }
            }
            type = Type::Matrix;
            data = validate_matrix(std::move(matrix));
        } else {
            type = Type::Array;
            data = arr;
        }
    }

    bool Value::is_numeric() const {
        if (type == Type::Int || type == Type::Float ||
            type == Type::BigInt || type == Type::Rational ||
            type == Type::Irrational) {
            return true;
        }
        if (type == Type::Symbolic) {
            const auto& sp = std::get<std::shared_ptr<SymbolicExpr>>(data);
            if (!sp) return false;
            try {
                auto simp = sp->simplify();
                return simp && simp->is_number();
            } catch (const std::invalid_argument&) {
                return false;
            } catch (const std::out_of_range&) {
                return false;
            }
        }
        return false;
    }

    LMCAS::Result<lmmc_real_t> Value::as_number_checked() const {
        if (type == Type::Infinity) {
            lmmc_real_t inf;
            lmmc_inf(&inf);
            const auto sign = std::get_if<int>(&data);
            return LMCAS::Result<lmmc_real_t>::success(
                (sign && *sign > 0) ? inf : -inf);
        }
        if (type == Type::Int) {
            return LMCAS::Result<lmmc_real_t>::success(
                static_cast<lmmc_real_t>(std::get<int>(data)));
        }
        if (type == Type::Float) {
            return LMCAS::Result<lmmc_real_t>::success(
                std::get<lmmc_real_t>(data));
        }
        if (type == Type::BigInt) {
            return exact_number_as_floating(
                std::get<::LMCAS::BigInt>(data));
        }
        if (type == Type::Rational) {
            return exact_number_as_floating(
                std::get<::LMCAS::Rational>(data));
        }
        if (type == Type::Irrational) {
            return LMCAS::Result<lmmc_real_t>::success(
                std::get<::LMCAS::Irrational>(data).to_double());
        }
        if (type == Type::Symbolic) {
            const auto& expression =
                std::get<std::shared_ptr<SymbolicExpr>>(data);
            auto evaluated = LMCAS::evaluate_numeric(*expression);
            if (!evaluated) {
                return LMCAS::Result<lmmc_real_t>::failure(
                    evaluated.error());
            }
            if (!evaluated.value().is_finite() ||
                !std::isfinite(evaluated.value().value)) {
                return LMCAS::Result<lmmc_real_t>::failure(
                    LMCAS::CasErrc::NumericFailure,
                    "symbolic value did not evaluate to a finite number",
                    "Value::as_number_checked");
            }
            return LMCAS::Result<lmmc_real_t>::success(
                static_cast<lmmc_real_t>(evaluated.value().value));
        }
        return LMCAS::Result<lmmc_real_t>::failure(
            LMCAS::CasErrc::InvalidArgument,
            "value is not numeric", "Value::as_number_checked");
    }


    LMCAS::Result<::LMCAS::Rational> Value::as_rational_checked() const {
        if (type == Type::Rational) {
            return LMCAS::Result<::LMCAS::Rational>::success(
                std::get<::LMCAS::Rational>(data));
        }
        if (type == Type::Int) {
            return LMCAS::Result<::LMCAS::Rational>::success(
                ::LMCAS::Rational(std::get<int>(data)));
        }
        if (type == Type::Float) {
            return LMCAS::Result<::LMCAS::Rational>::success(
                ::LMCAS::Rational::from_double(std::get<lmmc_real_t>(data)));
        }
        if (type == Type::BigInt) {
            return LMCAS::Result<::LMCAS::Rational>::success(
                ::LMCAS::Rational(std::get<::LMCAS::BigInt>(data)));
        }
        if (type == Type::Irrational) {
            return std::get<::LMCAS::Irrational>(data).as_rational_checked();
        }
        return LMCAS::Result<::LMCAS::Rational>::failure(
            LMCAS::CasErrc::InvalidArgument,
            "value cannot be represented as a rational",
            "Value::as_rational_checked");
    }


    LMCAS::Result<::LMCAS::Irrational> Value::as_irrational_checked() const {
        if (type == Type::Irrational) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                std::get<::LMCAS::Irrational>(data));
        }
        if (type == Type::Int) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                ::LMCAS::Irrational::constant(::LMCAS::Rational(std::get<int>(data))));
        }
        if (type == Type::Float) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                ::LMCAS::Irrational::constant(std::get<lmmc_real_t>(data)));
        }
        if (type == Type::Rational) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                ::LMCAS::Irrational::constant(
                    std::get<::LMCAS::Rational>(data)));
        }
        if (type == Type::BigInt) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                ::LMCAS::Irrational::constant(
                    ::LMCAS::Rational(std::get<::LMCAS::BigInt>(data))));
        }
        return LMCAS::Result<::LMCAS::Irrational>::failure(
            LMCAS::CasErrc::InvalidArgument,
            "value cannot be represented as an irrational wrapper",
            "Value::as_irrational_checked");
    }


    LMCAS::Result<std::shared_ptr<SymbolicExpr>> Value::as_symbolic_checked() const {
        if (type == Type::Infinity) {
            const auto sign = std::get_if<int>(&data);
            return LMCAS::Result<std::shared_ptr<SymbolicExpr>>::success(
                SymbolicExpr::infinity(sign ? *sign : 1));
        }
        if (type == Type::Symbolic) {
            return LMCAS::Result<std::shared_ptr<SymbolicExpr>>::success(
                std::get<std::shared_ptr<SymbolicExpr>>(data));
        }
        if (type == Type::Float) {
            return LMCAS::Result<std::shared_ptr<SymbolicExpr>>::success(
                SymbolicExpr::number(std::get<lmmc_real_t>(data)));
        }
        if (type == Type::Int || type == Type::Rational ||
            type == Type::BigInt) {
            auto rational = as_rational_checked();
            if (!rational) {
                return LMCAS::Result<std::shared_ptr<SymbolicExpr>>::failure(
                    rational.error());
            }
            return LMCAS::Result<std::shared_ptr<SymbolicExpr>>::success(
                SymbolicExpr::number(rational.value()));
        }
        if (type == Type::Irrational) {
            return LMCAS::Result<std::shared_ptr<SymbolicExpr>>::success(
                std::get<::LMCAS::Irrational>(data).to_symbolic());
        }
        if (type == Type::Matrix) {
            const auto& matrix =
                std::get<std::vector<std::vector<Value>>>(data);
            std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>
                symbolic_matrix;
            symbolic_matrix.reserve(matrix.size());
            for (const auto& row : matrix) {
                std::vector<std::shared_ptr<SymbolicExpr>> symbolic_row;
                symbolic_row.reserve(row.size());
                for (const auto& value : row) {
                    auto symbolic = value.as_symbolic_checked();
                    if (!symbolic) {
                        return LMCAS::Result<
                            std::shared_ptr<SymbolicExpr>>::failure(
                                symbolic.error());
                    }
                    symbolic_row.push_back(std::move(symbolic.value()));
                }
                symbolic_matrix.push_back(std::move(symbolic_row));
            }
            return LMCAS::Result<std::shared_ptr<SymbolicExpr>>::success(
                SymbolicExpr::matrix(symbolic_matrix));
        }
        return LMCAS::Result<std::shared_ptr<SymbolicExpr>>::failure(
            LMCAS::CasErrc::InvalidArgument,
            "value cannot be converted to a symbolic expression",
            "Value::as_symbolic_checked");
    }


    bool Value::as_symbolic_compatible() const {
        switch (type) {
            case Type::Symbolic: {
                return true;
            }
            case Type::Infinity: {
                return true;
            }
            case Type::Int: {
                return true;
            }
            case Type::Float: {
                return true;
            }
            case Type::Rational: {
                return true;
            }
            case Type::BigInt: {
                return true;
            }
            case Type::Irrational: {
                return true;
            }
            case Type::Matrix: {
                return true;
            }
            default:
                return false;
        }
    }

    std::string Value::to_string() const {
        if (type == Type::String) {
            return _str_cache;
        }
        if (!_str_cache.empty()) {
            return _str_cache;
        }
        switch (type) {
            case Type::Infinity: {
                return format_infinity(std::get_if<int>(&data));
            }
            case Type::Null: {
                return "null";
            }
            case Type::Int: {
                return std::to_string(std::get<int>(data));
            }
            case Type::Float: {
                return format_float(std::get<lmmc_real_t>(data));
            }
            case Type::Array: {
                return format_array(std::get<std::vector<Value>>(data));
            }
            case Type::Matrix: {
                return format_matrix(std::get<std::vector<std::vector<Value>>>(data));
            }
            case Type::BigInt: {
                return std::get<::LMCAS::BigInt>(data).to_string();
            }
            case Type::Rational: {
                return std::get<::LMCAS::Rational>(data).to_string();
            }
            case Type::Irrational: {
                return std::get<::LMCAS::Irrational>(data).to_string();
            }
            case Type::Symbolic: {
                const auto& expression =
                    std::get<std::shared_ptr<SymbolicExpr>>(data);
                if (!expression) {
                    return "<invalid-symbolic>";
                }
                return expression->to_string();
            }
            default:
                return "<unknown>";
        }
    }

    bool Value::operator==(const Value& other) const {
        if (type != other.type) return false;
        if (type == Type::String) {
            return _str_cache == other._str_cache;
        }
        if (data.index() != other.data.index()) return false;
        return std::visit([&other](const auto& val1) -> bool {
            using T = std::decay_t<decltype(val1)>;
            if constexpr (std::is_same_v<T, std::nullptr_t>) {
                return true;
            } else if constexpr (std::is_same_v<T, std::shared_ptr<SymbolicExpr>>) {
                if (!val1 && !std::get<T>(other.data)) return true;
                if (!val1 || !std::get<T>(other.data)) return false;
                return val1->compare(std::get<T>(other.data)) == 0;
            } else {
                return val1 == std::get<T>(other.data);
            }
        }, data);
    }

    bool Value::operator<(const Value& other) const {
        if (type != other.type) return static_cast<int>(type) < static_cast<int>(other.type);
        if (type == Type::String) {
            return _str_cache < other._str_cache;
        }
        if (data.index() != other.data.index()) return data.index() < other.data.index();
        return std::visit([&other](const auto& val1) -> bool {
            using T = std::decay_t<decltype(val1)>;
            if constexpr (std::is_same_v<T, std::nullptr_t>) {
                return false;
            } else if constexpr (std::is_same_v<T, std::shared_ptr<SymbolicExpr>>) {
                if (!val1 && !std::get<T>(other.data)) return false;
                if (!val1) return true;
                if (!std::get<T>(other.data)) return false;
                return val1->compare(std::get<T>(other.data)) < 0;
            } else {
                return val1 < std::get<T>(other.data);
            }
        }, data);
    }

    LMCAS::Result<Value> Value::vector_add_checked(const Value& other) const {
        if (!is_array() || !other.is_array()) {
            return LMCAS::Result<Value>::failure(
                LMCAS::CasErrc::InvalidArgument,
                "vector addition requires two arrays",
                "Value::vector_add_checked");
        }
        const auto& lhs = std::get<std::vector<Value>>(data);
        const auto& rhs = std::get<std::vector<Value>>(other.data);
        if (lhs.size() != rhs.size()) {
            return LMCAS::Result<Value>::failure(
                LMCAS::CasErrc::DimensionMismatch,
                "vector lengths do not match",
                "Value::vector_add_checked");
        }
        std::vector<Value> values;
        values.reserve(lhs.size());
        for (size_t index = 0; index < lhs.size(); ++index) {
            auto left = lhs[index].as_number_checked();
            if (!left) return LMCAS::Result<Value>::failure(left.error());
            auto right = rhs[index].as_number_checked();
            if (!right) return LMCAS::Result<Value>::failure(right.error());
            values.emplace_back(left.value() + right.value());
        }
        return LMCAS::Result<Value>::success(Value(values));
    }


    LMCAS::Result<Value> Value::dot_product_checked(const Value& other) const {
        if (!is_array() || !other.is_array()) {
            return LMCAS::Result<Value>::failure(
                LMCAS::CasErrc::InvalidArgument,
                "dot product requires two arrays",
                "Value::dot_product_checked");
        }
        const auto& lhs = std::get<std::vector<Value>>(data);
        const auto& rhs = std::get<std::vector<Value>>(other.data);
        if (lhs.size() != rhs.size()) {
            return LMCAS::Result<Value>::failure(
                LMCAS::CasErrc::DimensionMismatch,
                "vector lengths do not match",
                "Value::dot_product_checked");
        }
        lmmc_real_t sum = 0.0;
        for (size_t index = 0; index < lhs.size(); ++index) {
            auto left = lhs[index].as_number_checked();
            if (!left) return LMCAS::Result<Value>::failure(left.error());
            auto right = rhs[index].as_number_checked();
            if (!right) return LMCAS::Result<Value>::failure(right.error());
            sum += left.value() * right.value();
        }
        return LMCAS::Result<Value>::success(Value(sum));
    }


    LMCAS::Result<Value> Value::matrix_multiply_checked(
        const Value& other) const {
        if (!is_matrix() || !other.is_matrix()) {
            return LMCAS::Result<Value>::failure(
                LMCAS::CasErrc::InvalidArgument,
                "matrix multiplication requires two matrices",
                "Value::matrix_multiply_checked");
        }
        const auto& lhs =
            std::get<std::vector<std::vector<Value>>>(data);
        const auto& rhs =
            std::get<std::vector<std::vector<Value>>>(other.data);
        if (lhs.empty() || rhs.empty() ||
            lhs.front().size() != rhs.size()) {
            return LMCAS::Result<Value>::failure(
                LMCAS::CasErrc::DimensionMismatch,
                "matrix dimensions are not multiplicatively compatible",
                "Value::matrix_multiply_checked");
        }
        const size_t rows = lhs.size();
        const size_t columns = rhs.front().size();
        const size_t inner = lhs.front().size();
        std::vector<std::vector<Value>> values(
            rows, std::vector<Value>(columns, Value(0.0)));
        for (size_t row = 0; row < rows; ++row) {
            for (size_t column = 0; column < columns; ++column) {
                lmmc_real_t sum = 0.0;
                for (size_t index = 0; index < inner; ++index) {
                    auto left = lhs[row][index].as_number_checked();
                    if (!left) {
                        return LMCAS::Result<Value>::failure(left.error());
                    }
                    auto right = rhs[index][column].as_number_checked();
                    if (!right) {
                        return LMCAS::Result<Value>::failure(right.error());
                    }
                    sum += left.value() * right.value();
                }
                values[row][column] = Value(sum);
            }
        }
        return LMCAS::Result<Value>::success(Value(values));
    }


    std::shared_ptr<SymbolicExpr> Value::require_symbolic(
        const std::shared_ptr<SymbolicExpr>& expression) {
        if (!expression) {
            throw std::invalid_argument(
                "Value: symbolic expression must not be null");
        }
        return expression;
    }

    std::vector<std::vector<Value>> Value::validate_matrix(
        std::vector<std::vector<Value>> matrix) {
        if (matrix.empty()) return matrix;
        const size_t columns = matrix.front().size();
        for (const auto& row : matrix) {
            if (row.size() != columns) {
                throw std::invalid_argument(
                    "Value: all matrix rows must have the same number of columns");
            }
        }
        return matrix;
    }

}
