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
#include "lmmc/config.h"

#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace LMCAS {

/** @brief 运行时统一值类型,支持数值,符号,容器等类型的动态表示 */
class LMCAS_API Value final {
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
    Value(lmmc_real_t f);
    Value(const ::LMCAS::BigInt& bi) : type(Type::BigInt), data(bi) {}
    Value(const ::LMCAS::Rational& r) : type(Type::Rational), data(r) {}
    Value(const ::LMCAS::Irrational& ir) : type(Type::Irrational), data(ir) {}
    Value(const std::shared_ptr<SymbolicExpr>& sym)
        : type(Type::Symbolic), data(require_symbolic(sym)) {}
    Value(const std::vector<Value>& arr);
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
    bool is_numeric() const;

    /** @brief Checked conversion to a floating-point value. */
    LMCAS::Result<lmmc_real_t> as_number_checked() const;


    /** @brief Checked conversion to a rational value. */
    LMCAS::Result<::LMCAS::Rational> as_rational_checked() const;


    /** @brief Checked conversion to an irrational-value wrapper. */
    LMCAS::Result<::LMCAS::Irrational> as_irrational_checked() const;


    /** @brief Checked conversion to a symbolic expression. */
    LMCAS::Result<std::shared_ptr<SymbolicExpr>> as_symbolic_checked() const;


    /**
     * @brief 判断值是否可转换为符号表达式
     * @return 若可转换则返回 true
     */
    bool as_symbolic_compatible() const;

    /**
     * @brief 将值转换为可读字符串
     * @return 格式化的字符串表示
     */
    std::string to_string() const;

    bool operator==(const Value& other) const;

    bool operator<(const Value& other) const;

    /** @brief Checked element-wise vector addition. */
    LMCAS::Result<Value> vector_add_checked(const Value& other) const;


    /** @brief Checked vector dot product. */
    LMCAS::Result<Value> dot_product_checked(const Value& other) const;


    /** @brief Checked matrix multiplication. */
    LMCAS::Result<Value> matrix_multiply_checked(
        const Value& other) const;


private:
    static std::shared_ptr<SymbolicExpr> require_symbolic(
        const std::shared_ptr<SymbolicExpr>& expression);

    static std::vector<std::vector<Value>> validate_matrix(
        std::vector<std::vector<Value>> matrix);

    std::string _str_cache;  ///< 字符串构造时的缓存(兼容旧测试)
};

} // namespace LMCAS
