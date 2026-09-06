/**
 * @file value.hpp
 * @brief 运行时值类型 Value,支持数值,符号,容器类型.
 */
#pragma once
#include "lmcas_export.hpp"
#include "bigint.hpp"
#include "irrational.hpp"
#include "rational.hpp"
#include "symbolic.hpp"
#include "result.hpp"
#include "numeric_evaluation.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>
#include <sstream>

namespace LMCAS {

/** @brief 运行时统一值类型,支持数值,符号,容器等类型的动态表示 */
class Value final {
public:
    /** @brief 值的类型枚举 */
    enum class Type {
        Null, Int, Float, BigInt,
        Rational, Irrational, Symbolic,
        Infinity, Array, Matrix,
        String
    };
private:
    Type type;  ///< 当前值的类型

public:
    /** @brief 内部存储的 variant 类型 */
    using DataType = std::variant<
        std::nullptr_t,
        int, lmmc_real_t,
        ::LMCAS::BigInt, ::LMCAS::Rational, ::LMCAS::Irrational,
        std::shared_ptr<SymbolicExpr>,
        std::vector<Value>,
        std::vector<std::vector<Value>>>;

private:
    DataType data;  ///< 实际存储的数据

public:
    /** @brief 返回只读类型标记。 */
    Type kind() const noexcept { return type; }

    /** @brief 返回只读底层存储，调用方不能破坏类型不变量。 */
    const DataType& storage() const noexcept { return data; }

    ~Value() = default;

    /** @brief 默认构造,初始化为 Null */
    Value() : type(Type::Null), data(std::in_place_index<0>, nullptr) {}

    Value(std::nullptr_t) : type(Type::Null), data(std::in_place_index<0>, nullptr) {}
    Value(int i) : type(Type::Int), data(i) {}
    Value(lmmc_real_t f) : type(Type::Float), data(f) {
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
    Value(const ::LMCAS::BigInt& bi) : type(Type::BigInt), data(bi) {}
    Value(const ::LMCAS::Rational& r) : type(Type::Rational), data(r) {}
    Value(const ::LMCAS::Irrational& ir) : type(Type::Irrational), data(ir) {}
    Value(const std::shared_ptr<SymbolicExpr>& sym)
        : type(Type::Symbolic), data(require_symbolic(sym)) {}
    Value(const std::vector<Value>& arr) {
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
    Value(const std::vector<std::vector<Value>>& mat)
        : type(Type::Matrix), data(validate_matrix(mat)) {}

    /// 字符串使用独立的 Type::String 标记,与 Null 保持类型区分.
    Value(const std::string& s) : type(Type::String), data(nullptr), _str_cache(s) {}
    Value(const char* s) : type(Type::String), data(nullptr), _str_cache(s ? s : "") {}

    bool is_null() const { return type == Type::Null; }
    bool is_string() const { return type == Type::String; }
    bool is_infinity() const { return type == Type::Infinity; }
    bool is_int() const { return type == Type::Int; }
    bool is_float() const { return type == Type::Float; }
    bool is_array() const { return type == Type::Array; }
    bool is_matrix() const { return type == Type::Matrix; }
    bool is_bigint() const { return type == Type::BigInt; }
    bool is_rational() const { return type == Type::Rational; }
    bool is_irrational() const { return type == Type::Irrational; }
    bool is_symbolic() const { return type == Type::Symbolic; }
    bool is_numeric() const {
        if (type == Type::Int || type == Type::Float ||
            type == Type::BigInt || type == Type::Rational ||
            type == Type::Irrational) {
            return true;
        }
        /// 可化简为数值节点的 Symbolic 表达式进入 numeric 路径;
        /// 其余符号表达式保持 Symbolic 类型.
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

    /** @brief Checked conversion to a floating-point value. */
    LMCAS::Result<lmmc_real_t> as_number_checked() const {
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
            return LMCAS::Result<lmmc_real_t>::success(
                std::get<::LMCAS::BigInt>(data).to_double());
        }
        if (type == Type::Rational) {
            return LMCAS::Result<lmmc_real_t>::success(
                std::get<::LMCAS::Rational>(data).to_double());
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

    /**
     * @brief Legacy numeric conversion.
     * @return Numeric value, or 0.0 for incompatible non-symbolic values.
     */
    lmmc_real_t as_number() const {
        auto result = as_number_checked();
        if (result) return result.value();
        if (type == Type::Symbolic) {
            throw std::runtime_error(
                "numeric evaluation failed: " + result.error().message);
        }
        return 0.0;
    }

    /** @brief Checked conversion to a rational value. */
    LMCAS::Result<::LMCAS::Rational> as_rational_checked() const {
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
            return LMCAS::Result<::LMCAS::Rational>::success(
                ::LMCAS::Rational::from_double(
                    std::get<::LMCAS::Irrational>(data).to_double()));
        }
        return LMCAS::Result<::LMCAS::Rational>::failure(
            LMCAS::CasErrc::InvalidArgument,
            "value cannot be represented as a rational",
            "Value::as_rational_checked");
    }

    /** @brief Legacy rational conversion; incompatible values map to zero. */
    ::LMCAS::Rational as_rational() const {
        auto result = as_rational_checked();
        return result ? result.value() : ::LMCAS::Rational(0);
    }

    /** @brief Checked conversion to an irrational-value wrapper. */
    LMCAS::Result<::LMCAS::Irrational> as_irrational_checked() const {
        if (type == Type::Irrational) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                std::get<::LMCAS::Irrational>(data));
        }
        if (type == Type::Int) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                ::LMCAS::Irrational::constant(std::get<int>(data)));
        }
        if (type == Type::Float) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                ::LMCAS::Irrational::constant(std::get<lmmc_real_t>(data)));
        }
        if (type == Type::Rational) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                ::LMCAS::Irrational::constant(
                    std::get<::LMCAS::Rational>(data).to_double()));
        }
        if (type == Type::BigInt) {
            return LMCAS::Result<::LMCAS::Irrational>::success(
                ::LMCAS::Irrational::constant(
                    std::get<::LMCAS::BigInt>(data).to_double()));
        }
        return LMCAS::Result<::LMCAS::Irrational>::failure(
            LMCAS::CasErrc::InvalidArgument,
            "value cannot be represented as an irrational wrapper",
            "Value::as_irrational_checked");
    }

    /** @brief Legacy irrational conversion; incompatible values map to zero. */
    ::LMCAS::Irrational as_irrational() const {
        auto result = as_irrational_checked();
        return result ? result.value() : ::LMCAS::Irrational::constant(0);
    }

    /** @brief Checked conversion to a symbolic expression. */
    LMCAS::Result<std::shared_ptr<SymbolicExpr>> as_symbolic_checked() const {
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

    /** @brief Legacy symbolic conversion; incompatible values map to zero. */
    std::shared_ptr<SymbolicExpr> as_symbolic() const {
        auto result = as_symbolic_checked();
        return result ? result.value() : SymbolicExpr::number(0);
    }

    /**
     * @brief 判断值是否可转换为符号表达式
     * @return 若可转换则返回 true
     */
    bool as_symbolic_compatible() const {
        if (type == Type::Symbolic || type == Type::Infinity) return true;
        if (type == Type::Int || type == Type::Float ||
            type == Type::Rational || type == Type::BigInt) return true;
        if (type == Type::Irrational) return true;
        if (type == Type::Matrix) return true;
        return false;
    }

    /**
     * @brief 将值转换为可读字符串
     * @return 格式化的字符串表示
     */
    std::string to_string() const {
        if (type == Type::String) return _str_cache;
        if (!_str_cache.empty()) return _str_cache;
        switch (type) {
            case Type::Infinity: {
                auto sign = std::get_if<int>(&data);
                return (sign && *sign > 0) ? "inf" : "-inf";
            }
            case Type::Null:
                return "null";
            case Type::Int:
                return std::to_string(std::get<int>(data));
            case Type::Float: {
                lmmc_real_t val = std::get<lmmc_real_t>(data);
                std::string str = std::to_string(val);
                str.erase(str.find_last_not_of('0') + 1, std::string::npos);
                str.erase(str.find_last_not_of('.') + 1, std::string::npos);
                return str;
            }
            case Type::Array: {
                std::string res = "[";
                const auto& arr = std::get<std::vector<Value>>(data);
                for (size_t i = 0; i < arr.size(); ++i) {
                    if (i) res += ", ";
                    res += arr[i].to_string();
                }
                res += "]";
                return res;
            }
            case Type::Matrix: {
                std::string res = "[";
                const auto& mat = std::get<std::vector<std::vector<Value>>>(data);
                for (size_t i = 0; i < mat.size(); ++i) {
                    if (i) res += ", ";
                    res += "[";
                    for (size_t j = 0; j < mat[i].size(); ++j) {
                        if (j) res += ", ";
                        res += mat[i][j].to_string();
                    }
                    res += "]";
                }
                res += "]";
                return res;
            }
            case Type::BigInt:
                return std::get<::LMCAS::BigInt>(data).to_string();
            case Type::Rational:
                return std::get<::LMCAS::Rational>(data).to_string();
            case Type::Irrational:
                return std::get<::LMCAS::Irrational>(data).to_string();
            case Type::Symbolic: {
                const auto& expression =
                    std::get<std::shared_ptr<SymbolicExpr>>(data);
                if (!expression) return "<invalid-symbolic>";
                return expression->to_string();
            }
            default:
                return "<unknown>";
        }
    }

    bool operator==(const Value& other) const {
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

    bool operator<(const Value& other) const {
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

    /** @brief Checked element-wise vector addition. */
    LMCAS::Result<Value> vector_add_checked(const Value& other) const {
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

    /** @brief Legacy vector addition; failures map to Null. */
    Value vector_add(const Value& other) const {
        auto result = vector_add_checked(other);
        return result ? std::move(result.value()) : Value();
    }

    /** @brief Checked vector dot product. */
    LMCAS::Result<Value> dot_product_checked(const Value& other) const {
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

    /** @brief Legacy dot product; failures map to Null. */
    Value dot_product(const Value& other) const {
        auto result = dot_product_checked(other);
        return result ? std::move(result.value()) : Value();
    }

    /** @brief Checked matrix multiplication. */
    LMCAS::Result<Value> matrix_multiply_checked(
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

    /** @brief Legacy matrix multiplication; failures map to Null. */
    Value matrix_multiply(const Value& other) const {
        auto result = matrix_multiply_checked(other);
        return result ? std::move(result.value()) : Value();
    }

private:
    static std::shared_ptr<SymbolicExpr> require_symbolic(
        const std::shared_ptr<SymbolicExpr>& expression) {
        if (!expression) {
            throw std::invalid_argument(
                "Value: symbolic expression must not be null");
        }
        return expression;
    }

    static std::vector<std::vector<Value>> validate_matrix(
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

    std::string _str_cache;  ///< 字符串构造时的缓存(兼容旧测试)
};

} // namespace LMCAS
