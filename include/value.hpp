/**
 * @file value.hpp
 * @brief 运行时值类型 Value,支持数值,符号,容器类型.
 */
#pragma once
#include "lamina_export.hpp"
#include "bigint.hpp"
#include "irrational.hpp"
#include "rational.hpp"
#include "symbolic.hpp"
#include "numeric_evaluation.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"

#include <iostream>
#include <string>
#include <variant>
#include <vector>
#include <sstream>

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
    Type type;  ///< 当前值的类型

    /** @brief 内部存储的 variant 类型 */
    using DataType = std::variant<
        std::nullptr_t,
        int, lmmc_real_t,
        ::BigInt, ::Rational, ::Irrational,
        std::shared_ptr<SymbolicExpr>,
        std::vector<Value>,
        std::vector<std::vector<Value>>>;

    DataType data;  ///< 实际存储的数据

    ~Value() = default;

    /** @brief 默认构造,初始化为 Null */
    Value() : type(Type::Null), data(std::in_place_index<0>, nullptr) {}

    Value(std::nullptr_t) : type(Type::Null), data(std::in_place_index<0>, nullptr) {}
    Value(int i) : type(Type::Int), data(i) {}
    Value(lmmc_real_t f) : type(Type::Float), data(f) {
        int res;
        lmmc_isinf(f, &res);
        if (res) {
            if (f < 0) res = -1;
            this->type = Type::Infinity;
            this->data = DataType(std::in_place_index<1>, res);
        }
    }
    Value(const ::BigInt& bi) : type(Type::BigInt), data(bi) {}
    Value(const ::Rational& r) : type(Type::Rational), data(r) {}
    Value(const ::Irrational& ir) : type(Type::Irrational), data(ir) {}
    Value(const std::shared_ptr<SymbolicExpr>& sym) : type(Type::Symbolic), data(sym) {}
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
            data = matrix;
        } else {
            type = Type::Array;
            data = arr;
        }
    }
    Value(const std::vector<std::vector<Value>>& mat) : type(Type::Matrix), data(mat) {}

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
            } catch (...) {
                return false;
            }
        }
        return false;
    }

    /**
     * @brief 将值转换为浮点数
     * @return 浮点数近似值,非数值类型返回 0.0
     */
    lmmc_real_t as_number() const {
        if (type == Type::Infinity) {
            lmmc_real_t inf;
            lmmc_inf(&inf);
            auto sign = std::get_if<int>(&data);
            return (sign && *sign > 0) ? inf : -inf;
        }
        if (type == Type::Int) return static_cast<lmmc_real_t>(std::get<int>(data));
        if (type == Type::Float) return std::get<lmmc_real_t>(data);
        if (type == Type::BigInt) {
            /// BigInt::to_double() 直接保留大整数的双精度近似.
            return std::get<::BigInt>(data).to_double();
        }
        if (type == Type::Rational) {
            return std::get<::Rational>(data).to_double();
        }
        if (type == Type::Irrational) {
            return std::get<::Irrational>(data).to_double();
        }
        if (type == Type::Symbolic) {
            const auto& expr = std::get<std::shared_ptr<SymbolicExpr>>(data);
            if (!expr) {
                throw std::runtime_error(
                    "numeric evaluation failed: symbolic value is null");
            }
            auto evaluated = lamina::evaluate_numeric(*expr);
            if (!evaluated || !evaluated.value().is_finite() ||
                !std::isfinite(evaluated.value().value)) {
                throw std::runtime_error(
                    "numeric evaluation failed: " +
                    (evaluated ? std::string("non-finite result")
                               : evaluated.error().message));
            }
            return static_cast<lmmc_real_t>(evaluated.value().value);
        }
        return 0.0;
    }

    /**
     * @brief 将值转换为有理数
     * @return 有理数表示,非数值类型返回 0
     */
    ::Rational as_rational() const {
        if (type == Type::Rational) return std::get<::Rational>(data);
        if (type == Type::Int) return ::Rational(std::get<int>(data));
        if (type == Type::Float) return ::Rational::from_double(std::get<lmmc_real_t>(data));
        if (type == Type::BigInt) {
            /// 直接从 BigInt 构造 Rational,保留全部整数位.
            return ::Rational(std::get<::BigInt>(data));
        }
        if (type == Type::Irrational) {
            return ::Rational::from_double(std::get<::Irrational>(data).to_double());
        }
        return ::Rational(0);
    }

    /**
     * @brief 将值转换为无理数表示
     * @return 无理数对象
     */
    ::Irrational as_irrational() const {
        if (type == Type::Irrational) return std::get<::Irrational>(data);
        if (type == Type::Int) return ::Irrational::constant(std::get<int>(data));
        if (type == Type::Float) return ::Irrational::constant(std::get<lmmc_real_t>(data));
        if (type == Type::Rational) return ::Irrational::constant(std::get<::Rational>(data).to_double());
        if (type == Type::BigInt) {
            /// to_double() 直接生成大整数的浮点近似,保留超出 int 范围的数量级.
            return ::Irrational::constant(std::get<::BigInt>(data).to_double());
        }
        return ::Irrational::constant(0);
    }

    /**
     * @brief 将值转换为符号表达式
     * @return 符号表达式智能指针,不可转换时返回数值 0
     */
    std::shared_ptr<SymbolicExpr> as_symbolic() const {
        if (type == Type::Infinity) {
            auto sign = std::get_if<int>(&data);
            return SymbolicExpr::infinity(sign ? *sign : 1);
        }
        if (type == Type::Symbolic) return std::get<std::shared_ptr<SymbolicExpr>>(data);
        if (type == Type::Int || type == Type::Float || type == Type::Rational || type == Type::BigInt) {
            return SymbolicExpr::number(as_rational());
        }
        if (type == Type::Irrational) {
            return as_irrational().to_symbolic();
        }
        if (type == Type::Matrix) {
            const auto& mat = std::get<std::vector<std::vector<Value>>>(data);
            std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> sym_mat;
            for (const auto& row : mat) {
                std::vector<std::shared_ptr<SymbolicExpr>> sym_row;
                for (const auto& val : row) {
                    sym_row.push_back(val.as_symbolic());
                }
                sym_mat.push_back(sym_row);
            }
            return SymbolicExpr::matrix(sym_mat);
        }
        return SymbolicExpr::number(0);
    }

    /**
     * @brief 判断值是否可转换为符号表达式
     * @return 若可转换则返回 true
     */
    bool as_symbolic_compatible() const {
        if (type == Type::Symbolic) return true;
        if (type == Type::Int || type == Type::Float || type == Type::Rational || type == Type::BigInt) return true;
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
                return std::get<::BigInt>(data).to_string();
            case Type::Rational:
                return std::get<::Rational>(data).to_string();
            case Type::Irrational:
                return std::get<::Irrational>(data).to_string();
            case Type::Symbolic:
                return std::get<std::shared_ptr<SymbolicExpr>>(data)->to_string();
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

    /**
     * @brief 向量加法
     * @param other 另一个向量
     * @return 逐元素相加的结果向量
     */
    Value vector_add(const Value& other) const {
        if (!is_array() || !other.is_array()) return Value();
        const auto& a = std::get<std::vector<Value>>(data);
        const auto& b = std::get<std::vector<Value>>(other.data);
        if (a.size() != b.size()) return Value();
        std::vector<Value> result;
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i].is_numeric() && b[i].is_numeric()) {
                result.push_back(Value(a[i].as_number() + b[i].as_number()));
            } else {
                return Value();
            }
        }
        return Value(result);
    }

    /**
     * @brief 向量点积
     * @param other 另一个向量
     * @return 点积标量值
     */
    Value dot_product(const Value& other) const {
        if (!is_array() || !other.is_array()) return Value();
        const auto& a = std::get<std::vector<Value>>(data);
        const auto& b = std::get<std::vector<Value>>(other.data);
        if (a.size() != b.size()) return Value();
        lmmc_real_t result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i].is_numeric() && b[i].is_numeric()) {
                result += a[i].as_number() * b[i].as_number();
            } else {
                return Value();
            }
        }
        return Value(result);
    }

    /**
     * @brief 矩阵乘法
     * @param other 右矩阵
     * @return 乘积矩阵
     */
    Value matrix_multiply(const Value& other) const {
        if (!is_matrix() || !other.is_matrix()) return Value();
        const auto& a = std::get<std::vector<std::vector<Value>>>(data);
        const auto& b = std::get<std::vector<std::vector<Value>>>(other.data);
        if (a.empty() || b.empty() || a[0].size() != b.size()) return Value();
        size_t rows = a.size();
        size_t cols = b[0].size();
        size_t inner = a[0].size();
        std::vector<std::vector<Value>> result(rows, std::vector<Value>(cols, Value(0.0)));
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                lmmc_real_t sum = 0.0;
                for (size_t k = 0; k < inner; ++k) {
                    if (!a[i][k].is_numeric() || !b[k][j].is_numeric()) return Value();
                    sum += a[i][k].as_number() * b[k][j].as_number();
                }
                result[i][j] = Value(sum);
            }
        }
        return Value(result);
    }

private:
    std::string _str_cache;  ///< 字符串构造时的缓存(兼容旧测试)
};
