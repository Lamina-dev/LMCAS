/**
 * @file symbolic.hpp
 * @brief 符号表达式主类 SymbolicExpr 及其工厂方法、运算接口。
 */
#pragma once
#define _USE_MATH_DEFINES
#include "lamina_export.hpp"
#include "bigint.hpp"
#include "conditional_result.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <variant>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>
#include <functional>
#include <stdexcept>

class SymbolicExpr;

namespace lamina {
class AssumptionContext;

/** Result of an operation that produces a symbolic expression. */
using ExpressionResult = Result<std::shared_ptr<SymbolicExpr>>;

/** @brief Stable relational operator used by public APIs. */
enum class RelationOp {
    EQ,
    NEQ,
    LT,
    GT,
    LEQ,
    GEQ
};

namespace detail {
struct SymbolicExprAccess;
} // namespace detail
} // namespace lamina

/**
 * @brief 符号表达式主类，封装 AST 根节点并提供运算、化简、求解等接口。
 *
 * 通过静态工厂方法创建表达式，支持算术运算、三角函数、微积分、
 * 矩阵运算、方程求解等符号计算功能。
 */
class LAMINA_API SymbolicExpr : public std::enable_shared_from_this<SymbolicExpr> {
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;

    explicit SymbolicExpr(std::shared_ptr<const Impl> impl) : impl_(std::move(impl)) {
        if (!impl_) {
            throw std::invalid_argument("SymbolicExpr requires a non-null implementation");
        }
    }

    friend struct lamina::detail::SymbolicExprAccess;

public:

    /** @brief 表达式类型枚举（已废弃，保留用于兼容） */
    enum class Type {
        Number,      ///< 数值
        Sqrt,        ///< 平方根
        Root,        ///< n 次根
        Power,       ///< 幂运算
        Multiply,    ///< 乘法
        Add,         ///< 加法
        Subtract,    ///< 减法
        Infinity,    ///< 无穷大
        Variable,    ///< 变量

        Sin, Cos, Tan, Cot, Sec, Csc, ///< 三角函数

        ArcSin, ArcCos, ArcTan, Atan2, ///< 反三角函数

        Sinh, Cosh, Tanh, ///< 双曲函数

        Ln, Log, ///< 对数函数

        Abs, Fac, ///< 绝对值、阶乘

        Diff, Integral, Limit, ///< 微分、积分、极限

        Matrix, ///< 矩阵
        Vector, ///< 向量

    };

    SymbolicExpr() = delete;
    SymbolicExpr(Type) = delete;

    /**
     * @brief 与另一个表达式进行全序比较。
     * @param other 待比较的表达式
     * @return 小于返回负数，等于返回 0，大于返回正数
     */
    int compare(const std::shared_ptr<SymbolicExpr>& other) const;

    /**
     * @brief 将表达式中指定变量替换为给定值。
     * @param var_name 变量名
     * @param value 替换值
     * @return 替换后的新表达式
     */
    std::shared_ptr<SymbolicExpr> substitute(const std::string& var_name, const std::shared_ptr<SymbolicExpr>& value) const;

    /**
     * @brief 展开表达式（分配律展开乘法等）。
     * @return 展开后的表达式
     */
    std::shared_ptr<SymbolicExpr> expand() const;

    /** @brief 判断表达式是否为零。 */
    bool is_zero() const;
    /** @brief 判断表达式是否为一。 */
    bool is_one() const;

    /**
     * @brief 计算两个多项式表达式的最大公因式。
     * @param a 第一个多项式
     * @param b 第二个多项式
     * @return 首一 GCD 表达式；非精确有理多项式返回 nullptr
     * @deprecated Use lamina::symbolic_polynomial_gcd to preserve errors and
     *             resource-limit information.
     */
    static std::shared_ptr<SymbolicExpr> poly_gcd(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

    /**
     * @brief 计算两个多项式关于指定变量的结式。
     * @param a 第一个多项式
     * @param b 第二个多项式
     * @param var 消去的变量名
     * @return 结式表达式
     */
    static std::shared_ptr<SymbolicExpr> poly_resultant(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b, const std::string& var);

    /**
     * @brief 计算矩阵行列式。
     * @param mat 矩阵表达式
     * @return 行列式值
     */
    static std::shared_ptr<SymbolicExpr> determinant(const std::shared_ptr<SymbolicExpr>& mat);

    /**
     * @brief 计算矩阵转置。
     * @param mat 矩阵表达式
     * @return 转置矩阵
     */
    static std::shared_ptr<SymbolicExpr> transpose(const std::shared_ptr<SymbolicExpr>& mat);

    /**
     * @brief 计算矩阵的逆。
     * @param mat 矩阵表达式
     * @return 逆矩阵，不可逆时返回 nullptr
     */
    static std::shared_ptr<SymbolicExpr> inverse(const std::shared_ptr<SymbolicExpr>& mat);

    /**
     * @brief 计算矩阵的行最简阶梯形（RREF）。
     * @param mat 矩阵表达式
     * @return RREF 矩阵
     */
    static std::shared_ptr<SymbolicExpr> rref(const std::shared_ptr<SymbolicExpr>& mat);

    /**
     * @brief 计算矩阵的特征多项式。
     * @param mat 矩阵表达式
     * @param lambda 特征值变量名
     * @return 特征多项式表达式
     */
    static std::shared_ptr<SymbolicExpr> charpoly(const std::shared_ptr<SymbolicExpr>& mat, const std::string& lambda);

    /**
     * @brief 计算矩阵的特征值。
     * @param mat 矩阵表达式
     * @return 特征值列表表达式
     */
    static std::shared_ptr<SymbolicExpr> eigenvalues(const std::shared_ptr<SymbolicExpr>& mat);

    /**
     * @brief 计算矩阵的特征向量。
     * @param mat 矩阵表达式
     * @return 特征值与对应特征向量的列表
     */
    static std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> eigenvectors(const std::shared_ptr<SymbolicExpr>& mat);

    /**
     * @brief 求解方程组（无参数变量版本）。
     * @param equations 方程列表（每个表达式等于零）
     * @param vars 未知数变量名列表
     * @return 解的列表，每个解为变量名到表达式的映射
     */
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solve_system(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& vars);

    /**
     * @brief 求解含参数的方程组。
     * @param equations 方程列表（每个表达式等于零）
     * @param unknowns 未知数变量名列表
     * @param parameters 参数变量名列表
     * @return 解的列表，每个解为变量名到表达式的映射
     */
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solve_system(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    /**
     * @brief 计算表达式的极限。
     * @param var 趋近的变量名
     * @param point 趋近的目标值
     * @param direction 方向："+"（右极限）、"-"（左极限）、""（双侧）
     * @return 极限结果表达式
     */
    std::shared_ptr<SymbolicExpr> limit(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, const std::string& direction = "") const;

    /**
     * @brief 计算不定积分。
     * @param var 积分变量名
     * @return 积分结果表达式
     */
    std::shared_ptr<SymbolicExpr> integrate(const std::string& var) const;

    /**
     * @brief 在指定点展开为泰勒级数。
     * @param var 展开变量名
     * @param point 展开中心点
     * @param order 展开阶数
     * @param ctx 可选的假设上下文，用于验证收敛域和选择展开形式
     * @return 级数表达式
     */
    std::shared_ptr<SymbolicExpr> series(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, int order, const lamina::AssumptionContext* ctx = nullptr) const;

    /**
     * @brief 执行符号除法。
     * @param num 被除数
     * @param den 除数
     * @return 商表达式
     */
    static std::shared_ptr<SymbolicExpr> divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den);

    /**
     * @brief 创建积分表达式节点（不求值）。
     * @param expr 被积表达式
     * @param var 积分变量名
     * @return 积分表达式
     */
    static std::shared_ptr<SymbolicExpr> make_integral(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);

    /**
     * @brief 创建极限表达式节点（不求值）。
     * @param expr 表达式
     * @param var 趋近变量名
     * @param point 趋近目标值
     * @return 极限表达式
     */
    static std::shared_ptr<SymbolicExpr> make_limit(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var, const std::shared_ptr<SymbolicExpr>& point);

    /**
     * @brief 对表达式进行因式分解。
     * @return 因式分解后的表达式
     */
    std::shared_ptr<SymbolicExpr> factor() const;

    /**
     * @brief 有理函数约分（消去分子分母公因式）。
     *
     * 将表达式视为有理函数 P(x)/Q(x)，计算 gcd(P, Q) 并约分。
     * 支持多元多项式（逐变量迭代 GCD）。
     * 例如 (x²-1)/(x-1) 约分为 x+1。
     *
     * @return 约分后的表达式
     */
    std::shared_ptr<SymbolicExpr> cancel() const;

    [[deprecated("Inspect expressions through public predicates and operations")]]
    Type get_type() const;

    [[deprecated("Construct and transform expressions through public operations")]]
    std::vector<std::shared_ptr<SymbolicExpr>> get_operands() const;

    [[deprecated("Use numeric evaluation or public numeric predicates")]]
    std::variant<int, ::BigInt, ::Rational> get_number_value() const;

    [[deprecated("Use expression substitution and public symbol predicates")]]
    std::string get_identifier() const;

    /**
     * @brief 创建整数数值表达式。
     * @param n 整数值
     * @return 数值表达式
     */
    static std::shared_ptr<SymbolicExpr> number(int n);

    /**
     * @brief 创建 long long 数值表达式。
     * @param n 整数值
     * @return 数值表达式
     */
    static std::shared_ptr<SymbolicExpr> number(long long n);

    /**
     * @brief 创建浮点数值表达式。
     * @param n 浮点值
     * @return 数值表达式
     */
    static std::shared_ptr<SymbolicExpr> number(double n);

    /**
     * @brief 创建大整数数值表达式。
     * @param bi 大整数
     * @return 数值表达式
     */
    static std::shared_ptr<SymbolicExpr> number(const ::BigInt& bi);

    /**
     * @brief 创建有理数数值表达式。
     * @param r 有理数
     * @return 数值表达式
     */
    static std::shared_ptr<SymbolicExpr> number(const ::Rational& r);

	/**
	 * @brief 创建无穷大表达式。
	 * @param k 正数表示正无穷，负数表示负无穷
	 * @return 无穷大表达式
	 */
	static std::shared_ptr<SymbolicExpr> infinity(int k = 1);

    /**
     * @brief 创建平方根表达式，等价于 operand^(1/2)。
     * @param operand 被开方数
     * @return 平方根表达式
     */
    static std::shared_ptr<SymbolicExpr> sqrt(std::shared_ptr<SymbolicExpr> operand);

    /**
     * @brief 创建乘法表达式。
     * @param left 左操作数
     * @param right 右操作数
     * @return 乘积表达式
     */
    static std::shared_ptr<SymbolicExpr> multiply(std::shared_ptr<SymbolicExpr> left, std::shared_ptr<SymbolicExpr> right);

    /**
     * @brief 创建加法表达式。
     * @param left 左操作数
     * @param right 右操作数
     * @return 和表达式
     */
    static std::shared_ptr<SymbolicExpr> add(std::shared_ptr<SymbolicExpr> left, std::shared_ptr<SymbolicExpr> right);

    /**
     * @brief 创建幂运算表达式。
     * @param base 底数
     * @param exponent 指数
     * @return 幂表达式
     */
    static std::shared_ptr<SymbolicExpr> power(std::shared_ptr<SymbolicExpr> base, std::shared_ptr<SymbolicExpr> exponent);

    /**
     * @brief 创建正弦函数表达式。
     * @param op 参数表达式
     * @return sin(op)
     */
    static std::shared_ptr<SymbolicExpr> sin(std::shared_ptr<SymbolicExpr> op);

    /**
     * @brief 创建余弦函数表达式。
     * @param op 参数表达式
     * @return cos(op)
     */
    static std::shared_ptr<SymbolicExpr> cos(std::shared_ptr<SymbolicExpr> op);

    /**
     * @brief 创建正切函数表达式。
     * @param op 参数表达式
     * @return tan(op)
     */
    static std::shared_ptr<SymbolicExpr> tan(std::shared_ptr<SymbolicExpr> op);

    /**
     * @brief 创建自然对数表达式。
     * @param op 参数表达式
     * @return ln(op)
     */
    static std::shared_ptr<SymbolicExpr> ln(std::shared_ptr<SymbolicExpr> op);

    /**
     * @brief 创建指数函数表达式 e^op。
     * @param op 指数参数
     * @return exp(op)
     */
    static std::shared_ptr<SymbolicExpr> exp(std::shared_ptr<SymbolicExpr> op);

    /**
     * @brief 创建 Lambert W 函数表达式。
     * @param op 参数表达式
     * @return W(op)
     */
    static std::shared_ptr<SymbolicExpr> lambertw(std::shared_ptr<SymbolicExpr> op);

    /**
     * @brief 创建以指定底数的对数表达式。
     * @param val 真数
     * @param base 底数
     * @return log_base(val)
     */
    static std::shared_ptr<SymbolicExpr> log(std::shared_ptr<SymbolicExpr> val, std::shared_ptr<SymbolicExpr> base);

    /**
     * @brief 创建双参数反正切表达式 atan2(y, x)。
     * @param y y 坐标
     * @param x x 坐标
     * @return atan2 表达式
     */
    static std::shared_ptr<SymbolicExpr> atan2(std::shared_ptr<SymbolicExpr> y, std::shared_ptr<SymbolicExpr> x);

    /**
     * @brief 创建 RootOf 表达式，表示多项式的第 index 个根。
     * @param poly 多项式表达式
     * @param var 多项式变量名
     * @param index 根的索引
     * @return RootOf 表达式
     */
    static std::shared_ptr<SymbolicExpr> root_of(std::shared_ptr<SymbolicExpr> poly, const std::string& var, int index);

	/**
	 * @brief 创建等式关系表达式 lhs = rhs。
	 * @param lhs 左侧表达式
	 * @param rhs 右侧表达式
	 * @return 等式表达式
	 */
	static std::shared_ptr<SymbolicExpr> eq(std::shared_ptr<SymbolicExpr> lhs, std::shared_ptr<SymbolicExpr> rhs);

    /**
     * @brief 计算表达式关于指定变量的不定积分。
     * @param op 被积表达式
     * @param var 积分变量名
     * @return 积分结果
     */
    static std::shared_ptr<SymbolicExpr> integral(std::shared_ptr<SymbolicExpr> op, const std::string& var);

    /**
     * @brief 计算表达式的极限。
     * @param op 表达式
     * @param var 趋近变量名
     * @param target 趋近目标值
     * @return 极限结果
     */
    static std::shared_ptr<SymbolicExpr> limit_func(std::shared_ptr<SymbolicExpr> op, const std::string& var, std::shared_ptr<SymbolicExpr> target);

    /**
     * @brief 从二维元素列表创建矩阵表达式。
     * @param elements 二维表达式数组，elements[i][j] 为第 i 行第 j 列元素
     * @return 矩阵表达式
     */
    static std::shared_ptr<SymbolicExpr> matrix(const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& elements);

    /**
     * @brief 创建变量表达式。
     * @param name 变量名
     * @return 变量表达式
     */
    static std::shared_ptr<SymbolicExpr> variable(const std::string& name);

    /**
     * @brief 对表达式进行化简。
     * @return 化简后的表达式
     */
    std::shared_ptr<SymbolicExpr> simplify() const;

    /**
     * @brief 对表达式进行三角恒等式化简。
     * @return 化简后的表达式
     */
    std::shared_ptr<SymbolicExpr> simplify_trig() const;

    /**
     * @brief 对表达式关于指定变量求导。
     * @param var_name 求导变量名
     * @return 导数表达式
     */
    std::shared_ptr<SymbolicExpr> differentiate(const std::string& var_name) const;

    /**
     * @brief 求解方程，返回指定变量的所有解。
     * @param eq 方程表达式（等于零的形式，或 RelationalNode）
     * @param var_name 求解的变量名
     * @return 解的列表
     */
    static std::vector<std::shared_ptr<SymbolicExpr>> solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name);

    /**
     * @brief 将表达式转换为可读字符串。
     * @return 字符串表示
     */
    std::string to_string() const;

    /** @brief 判断表达式是否为数值节点。 */
    bool is_number() const;

    /** @brief 判断表达式数值是否为零。 */
    bool get_number_value_is_zero() const;

    /** @brief 判断表达式是否为大整数类型。 */
    bool is_big_int() const;

    /** @brief 判断表达式是否为有理数类型。 */
    bool is_rational() const;

    /** @brief 判断表达式是否为整数（包括 BigInt、整数 Rational、近似整数浮点）。 */
    bool is_int() const;

    /**
     * @brief 获取数值表达式的值。
     * @return 数值（int、BigInt 或 Rational）
     * @throws std::runtime_error 表达式不是数值时抛出
     */
    std::variant<int, ::BigInt, ::Rational> get_number() const;

    int get_int() const;
    /**
     * @brief 获取大整数值。
     * @return BigInt 值
     * @throws std::runtime_error 表达式不是 BigInt 时抛出
     */
    ::BigInt get_big_int() const;
    /**
     * @brief 获取有理数值。
     * @return Rational 值
     * @throws std::runtime_error 表达式不是 Rational 时抛出
     */
    ::Rational get_rational() const;
    /**
     * @brief 将数值表达式转换为 Rational（BigInt 也会转换）。
     * @return 有理数值，非数值类型返回 Rational(0)
     */
    ::Rational convert_rational() const;

    /**
     * @brief 将表达式求值为浮点数。
     * @return 数值结果
     */
    lmmc_real_t to_numeric() const;
};
