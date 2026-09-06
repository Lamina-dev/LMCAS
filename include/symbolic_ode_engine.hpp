/**
 * @file symbolic_ode_engine.hpp
 * @brief 统一 ODE 求解引擎:分类,分派,多种求解方法.
 *
 * 扩展已有的 symbolic_ode.hpp,提供统一的 ODE 类型检测与分派机制,
 * 支持可分离变量,一阶线性,齐次,Bernoulli,恰当,高阶常系数,Euler 等类型.
 */
#pragma once
#include "computation_context.hpp"
#include "lmcas_export.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include "proof_outcome.hpp"
#include <memory>
#include <string>
#include <vector>
#include <optional>

namespace LMCAS {

/**
 * @brief ODE 类型枚举,标识方程所属的求解类别.
 */
enum class ODEType {
    Separable,              ///< 可分离变量型 dy/dx = f(x)g(y)
    Linear1,                ///< 一阶线性 dy/dx + P(x)y = Q(x)
    Linear2_ConstCoeff,     ///< 二阶常系数线性
    Homogeneous,            ///< 齐次方程 y' = f(y/x)
    Bernoulli,              ///< Bernoulli 方程 y' + P(x)y = Q(x)y^n
    Exact,                  ///< 恰当方程 M(x,y)dx + N(x,y)dy = 0
    HigherOrder_ConstCoeff, ///< 高阶常系数线性
    Euler,                  ///< Euler (Cauchy-Euler) 方程
    System,                 ///< ODE 系统
    LaplaceMethod,          ///< Laplace 变换法
    Frobenius,              ///< Frobenius 级数解法
    Unknown                 ///< 无法识别的类型
};

/**
 * @brief ODE 分类结果,包含类型及提取的结构信息.
 */
struct ODEClassification {
    ODEType type = ODEType::Unknown;  ///< 识别出的 ODE 类型
    int order = 0;                    ///< 方程阶数

    /// 一阶线性: P(x), Q(x) 使得 y' + P(x)y = Q(x)
    std::shared_ptr<SymbolicExpr> P_coeff;
    std::shared_ptr<SymbolicExpr> Q_coeff;

    /// Bernoulli: P(x), Q(x), n 使得 y' + P(x)y = Q(x)y^n
    std::shared_ptr<SymbolicExpr> bernoulli_P;
    std::shared_ptr<SymbolicExpr> bernoulli_Q;
    int bernoulli_n = 0;

    /// 恰当方程: M(x,y), N(x,y)
    std::shared_ptr<SymbolicExpr> exact_M;
    std::shared_ptr<SymbolicExpr> exact_N;

    /// 常系数: 系数列表 [a_n, a_{n-1}, ..., a_1, a_0] 和非齐次项 f(x)
    std::vector<double> const_coeffs;
    std::shared_ptr<SymbolicExpr> forcing_func;

    /// Euler 方程: 系数列表 [a_n, ..., a_0] 对应 x^n*y^(n) + ... + a_0*y = f(x)
    std::vector<double> euler_coeffs;
    std::shared_ptr<SymbolicExpr> euler_forcing;
};

/**
 * @brief ODE 求解结果.
 */
struct ODESolution {
    std::shared_ptr<SymbolicExpr> general_solution;  ///< 通解表达式
    std::vector<std::string> constants;              ///< 任意常数名称列表 (C1, C2, ...)
    ODEType method_used = ODEType::Unknown;          ///< 使用的求解方法
    std::shared_ptr<SymbolicExpr> normalized_residual;
    std::optional<ProofCertificate> certificate;
};

using ODESolutionResult = Result<ODESolution>;

/**
 * @brief 对一阶 ODE 进行类型分类.
 *
 * 给定 dy/dx = rhs(x, y) 形式的方程右端,检测其属于哪种可求解类型.
 * 检测顺序:可分离 -> 一阶线性 -> 齐次 -> Bernoulli -> 恰当.
 * 空 rhs 返回 Unknown。
 *
 * @param[in] rhs 方程右端表达式 (dy/dx = rhs)
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 分类结果
 */
LMCAS_API ODEClassification classify_first_order_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y);

/**
 * @brief 对高阶线性 ODE 进行类型分类.
 *
 * 给定系数列表和非齐次项,检测是否为常系数或 Euler 类型.
 * 系数列表按降阶排列:coeffs[0] 对应最高阶导数.
 * 常系数分类仅在所有系数都能表示为有限 double 时返回；否则返回 Unknown，
 * 避免在 const_coeffs 中伪造数值。
 *
 * @param[in] coeffs 各阶导数的系数表达式列表(降阶排列)
 * @param[in] forcing 非齐次项 f(x)
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 分类结果
 */
LMCAS_API ODEClassification classify_higher_order_ode(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y);

/**
 * @brief 检测表达式是否为可分离变量形式 f(x)*g(y).
 *
 * 空 rhs、空变量名或相同的 x/y 变量名返回 false。
 *
 * @param[in] rhs 方程右端表达式
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 若可分离返回 true
 */
LMCAS_API bool is_separable(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y);

/**
 * @brief 检测一阶 ODE 是否为线性形式 y' + P(x)y = Q(x).
 *
 * 空 rhs、空变量名或相同的 x/y 变量名返回 false；返回 false 时 P、Q
 * 均置空。
 *
 * @param[in] rhs 方程右端表达式 (dy/dx = rhs)
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @param[out] P 若为线性,输出 P(x)
 * @param[out] Q 若为线性,输出 Q(x)
 * @return 若为一阶线性返回 true
 */
LMCAS_API bool is_linear_first_order(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y,
    std::shared_ptr<SymbolicExpr>& P,
    std::shared_ptr<SymbolicExpr>& Q);

/**
 * @brief 检测一阶 ODE 是否为齐次形式 y' = f(y/x).
 *
 * 通过检查 f(tx, ty) = f(x, y) 来判断齐次性.
 * 判定只接受符号等价性证明；有限采样点相等不会被当作恒等证明。
 * 空变量名或相同的 x/y 变量名返回 false。
 *
 * @param[in] rhs 方程右端表达式
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 若为齐次返回 true
 */
LMCAS_API bool is_homogeneous_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y);

/**
 * @brief 检测一阶 ODE 是否为 Bernoulli 形式 y' + P(x)y = Q(x)y^n.
 *
 * 支持可由 AST 精确分解的任意 int 幂次，包括负幂；不通过有限点采样猜测
 * 幂次。
 * 空 rhs、空变量名或相同的 x/y 变量名返回 false；返回 false 时 P、Q
 * 均置空且 n 置零。
 *
 * @param[in] rhs 方程右端表达式 (dy/dx = rhs)
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @param[out] P 若为 Bernoulli,输出 P(x)
 * @param[out] Q 若为 Bernoulli,输出 Q(x)
 * @param[out] n 若为 Bernoulli,输出幂次 n (n != 0, 1)
 * @return 若为 Bernoulli 返回 true
 */
LMCAS_API bool is_bernoulli_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y,
    std::shared_ptr<SymbolicExpr>& P,
    std::shared_ptr<SymbolicExpr>& Q,
    int& n);

/**
 * @brief 检测方程是否为恰当形式 M(x,y) + N(x,y)y' = 0,即 partialM/partialy = partialN/partialx.
 * 判定只接受符号等价性证明；有限采样点相等不会被当作恒等证明。
 * 空表达式、空变量名或相同的 x/y 变量名返回 false。
 *
 * @param[in] M M(x,y) 表达式
 * @param[in] N N(x,y) 表达式
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 若恰当返回 true
 */
LMCAS_API bool is_exact_ode(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y);

/**
 * @brief 检测高阶 ODE 系数是否全为常数.
 *
 * 系数列表和自变量名必须非空，否则返回 false。
 * @param[in] coeffs 各阶导数的系数表达式列表
 * @param[in] x 自变量名
 * @return 若所有系数不依赖 x 返回 true
 */
LMCAS_API bool is_constant_coefficient(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& x);

/**
 * @brief 检测高阶 ODE 是否为 Euler (Cauchy-Euler) 形式.
 *
 * Euler 方程形如 x^n*y^(n) + a_{n-1}*x^{n-1}*y^{n-1} + ... + a_0*y = f(x),
 * 即第 k 阶导数的系数为常数乘以 x^k.
 * 系数列表和自变量名必须非空，最高阶系数必须非零。判定要求每个比值
 * coefficient_k/x^k 可被符号化简为有限常数，不用有限采样代替恒等证明。
 * 返回 false 时 euler_consts 置空。
 *
 * @param[in] coeffs 各阶导数的系数表达式列表(降阶排列)
 * @param[in] x 自变量名
 * @param[out] euler_consts 提取的常数系数列表
 * @return 若为 Euler 方程返回 true
 */
LMCAS_API bool is_euler_equation(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& x,
    std::vector<double>& euler_consts);


/**
 * @brief 求解齐次 ODE y' = f(y/x).
 *
 * 算法:令 v = y/x,则 y = vx,y' = v + xv'.
 * 代入得 v + xv' = f(v),即 xv' = f(v) - v,
 * 分离变量后积分,最后回代 v = y/x.
 *
 * @param[in] rhs 方程右端表达式 f(x, y)(满足 f(tx,ty)=f(x,y))
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 求解结果
 */
LMCAS_API ODESolutionResult solve_homogeneous_ode_checked(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文求解齐次 ODE,并显式报告无效输入和未覆盖域.
 */
LMCAS_API ODESolutionResult solve_homogeneous_ode_checked(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y);

/**
 * @brief 求解 Bernoulli 方程 y' + P(x)y = Q(x)y^n.
 *
 * 算法:令 v = y^(1-n),则 v' = (1-n)y^(-n)y'.
 * 代入得线性方程 v' + (1-n)P(x)v = (1-n)Q(x),
 * 用积分因子法求解后回代 y = v^(1/(1-n)).
 *
 * @param[in] P P(x) 系数
 * @param[in] Q Q(x) 系数
 * @param[in] n 幂次 (n != 0, 1)
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 求解结果
 */
LMCAS_API ODESolutionResult solve_bernoulli_ode_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    int n,
    const std::string& x,
    const std::string& y,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文求解 Bernoulli 方程,并显式报告无效输入和未覆盖域.
 */
LMCAS_API ODESolutionResult solve_bernoulli_ode_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    int n,
    const std::string& x,
    const std::string& y);

/**
 * @brief 求解恰当方程 M(x,y)dx + N(x,y)dy = 0.
 *
 * 算法:当 partialM/partialy = partialN/partialx 时,存在势函数 F(x,y) 使得
 * partialF/partialx = M, partialF/partialy = N.通过对 M 关于 x 积分并从 N 确定
 * y 相关部分来构造 F,解为 F(x,y) = C.
 *
 * 当方程不恰当时,尝试寻找积分因子 mu(x) 或 mu(y).
 *
 * @param[in] M M(x,y) 表达式
 * @param[in] N N(x,y) 表达式
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 求解结果(隐式解 F(x,y) = C)
 */
LMCAS_API ODESolutionResult solve_exact_ode_checked(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文求解恰当方程,并显式报告无效输入和未覆盖域.
 */
LMCAS_API ODESolutionResult solve_exact_ode_checked(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y);

/**
 * @brief 寻找积分因子使非恰当方程变为恰当方程.
 *
 * 尝试以下形式的积分因子:
 * 1. mu = mu(x):当 (partialM/partialy - partialN/partialx)/N 仅依赖 x 时
 * 2. mu = mu(y):当 (partialN/partialx - partialM/partialy)/M 仅依赖 y 时
 *
 * @param[in] M M(x,y) 表达式
 * @param[in] N N(x,y) 表达式
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 积分因子表达式,若找不到返回 nullptr
 */
LMCAS_API std::shared_ptr<SymbolicExpr> find_integrating_factor(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y);


/**
 * @brief 求解高阶常系数线性 ODE.
 *
 * 算法:构造特征多项式；二次实根使用避免相消的稳定公式；三至六次
 * 多项式先用变量尺度平衡系数，再用多项式与其导数的近似 GCD 提取
 * 高阶重根；最后根据根的类型构造通解:
 * - 实根 r(重数 m):e^(rx), x*e^(rx), ..., x^(m-1)*e^(rx)
 * - 复根 alpha+/-betai(重数 m):e^(alphax)(cos(betax), sin(betax)), x*e^(alphax)(cos(betax), sin(betax)), ...
 *
 * 非齐次支持域当前限于常数 forcing 且 y 的系数非零。
 * 支持阶数最高为 6。统一乘以任意有限非零常数不改变输入方程。
 * 特征根迭代必须收敛且通过后向误差检查；实根分类验证实轴投影相对于
 * 候选根的误差，复根必须形成经尺度平衡后多项式验证的共轭对，不使用
 * 会吞掉小非零虚部的固定绝对容差。归一化、尺度恢复、迭代或验证超出
 * 有限 double 域时返回 CasErrc::NumericFailure。
 *
 * @param[in] coeffs 各阶导数的数值系数(降阶排列:coeffs[0] 对应最高阶)
 * @param[in] forcing 非齐次项 f(x),齐次时为 nullptr 或零
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 求解结果
 */
LMCAS_API ODESolutionResult solve_higher_order_ode_checked(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文求解高阶常系数 ODE,并显式报告无效输入和未覆盖域.
 */
LMCAS_API ODESolutionResult solve_higher_order_ode_checked(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y);

/**
 * @brief 求解 Euler (Cauchy-Euler) 方程.
 *
 * 算法:对 Euler 方程 a_n*x^n*y^(n) + ... + a_0*y = f(x),
 * 令 x = e^t(即 t = ln(x)),将方程转化为关于 t 的常系数 ODE,
 * 求解后回代 t = ln(x).
 *
 * 支持二阶和三阶 Euler 方程.
 *
 * @param[in] euler_coeffs Euler 方程的常数系数列表(降阶排列)
 * @param[in] forcing 非齐次项 f(x),齐次时为 nullptr 或零
 * @param[in] x 自变量名
 * @param[in] y 因变量名
 * @return 求解结果
 */
LMCAS_API ODESolutionResult solve_euler_ode_checked(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文求解 Euler 方程,并显式报告无效输入和未覆盖域.
 */
LMCAS_API ODESolutionResult solve_euler_ode_checked(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y);


/**
 * @brief 奇点类型枚举.
 */
enum class ODESingularityType {
    Ordinary,         ///< 常点
    RegularSingular,  ///< 正则奇点
    IrregularSingular ///< 非正则奇点
};
using ODESingularityResult = Result<ODESingularityType>;

/**
 * @brief Frobenius 级数解结果.
 */
struct FrobeniusSolution {
    std::shared_ptr<SymbolicExpr> series_solution;  ///< 截断级数解
    ODESingularityType point_type;                     ///< 展开点的奇点类型
    std::vector<double> indicial_roots;             ///< 指标方程的根(正则奇点时)
    int truncation_order;                           ///< 截断阶数
    std::shared_ptr<SymbolicExpr> normalized_residual;
    std::optional<ProofCertificate> certificate;
};

using FrobeniusSolutionResult = Result<FrobeniusSolution>;

/**
 * @brief 用参数变分法求特解.
 *
 * @see Boyce and DiPrima, Elementary Differential Equations, Section3.6.
 *
 * 给定二阶线性 ODE y'' + p(x)y' + q(x)y = g(x) 的两个齐次解 y_1, y_2,
 * 计算 Wronskian W = y_1y_2' - y_2y_1',然后:
 *   u_1' = -y_2*g(x)/W,  u_2' = y_1*g(x)/W
 * 积分得 u_1, u_2,特解为 y_p = u_1*y_1 + u_2*y_2.
 *
 * @param[in] y1 第一个齐次解
 * @param[in] y2 第二个齐次解
 * @param[in] g 非齐次项 g(x)(方程已归一化为首项系数 1)
 * @param[in] x 自变量名
 * @return 求解结果,general_solution 为特解 y_p
 */
LMCAS_API ODESolutionResult solve_variation_of_parameters_checked(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文执行参数变分法,并显式报告无效输入和未覆盖域.
 */
LMCAS_API ODESolutionResult solve_variation_of_parameters_checked(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x);

/**
 * @brief 对 ODE 的指定点进行奇点分类.
 * 若代入点后仍含无法数值化的符号参数，则返回 CasErrc::Inconclusive，
 * 而不是把未知值当作系数极点。
 */
LMCAS_API ODESingularityResult classify_singular_point_checked(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    ComputationContext& context);

LMCAS_API ODESingularityResult classify_singular_point_checked(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x);

/**
 * @brief 用 Frobenius 方法求 ODE 的级数解.
 *
 * @see Bender and Orszag, Advanced Mathematical Methods for Scientists and Engineers, Chapter 3.
 *
 * 给定 y'' + p(x)y' + q(x)y = 0:
 * - 常点:假设 y = suma_n(x-x_0)ⁿ,代入匹配系数
 * - 正则奇点:假设 y = (x-x_0)^r*suma_n(x-x_0)ⁿ,
 *   先求指标方程确定 r,再递推确定 a_n
 *
 * 当前 checked 实现要求在截断阶数内使用的系数及其导数均可转换为有限
 * double；否则返回 CasErrc::Inconclusive，而不会把未知系数替换为零。
 *
 * @param[in] p 系数 p(x)(y' 的系数,方程已归一化)
 * @param[in] q 系数 q(x)(y 的系数,方程已归一化)
 * @param[in] x0 展开点
 * @param[in] x 自变量名
 * @param[in] order 截断阶数(默认 6)
 * @return Frobenius 解结果
 */
LMCAS_API FrobeniusSolutionResult solve_frobenius_checked(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    int order,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文求 Frobenius 级数解,并显式报告无效输入和未覆盖域.
 */
LMCAS_API FrobeniusSolutionResult solve_frobenius_checked(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    int order = 6);

} // namespace LMCAS
