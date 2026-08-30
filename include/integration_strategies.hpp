/** @file integration_strategies.hpp */ #pragma once
#include "integration_table.hpp"
#include "computation_context.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <variant>

namespace lamina {

class Integrator;

struct IntegrationNotApplicable {};

struct IntegrationCandidate {
    std::shared_ptr<SymbolicExpr> expression;
    std::string strategy_name;
};

using IntegrationStrategyOutcome =
    std::variant<IntegrationNotApplicable, IntegrationCandidate>;
using IntegrationStrategyResult = Result<IntegrationStrategyOutcome>;


/** @brief 积分策略基类，定义策略接口 */
class LAMINA_API IntegrationStrategy {
public:
    virtual ~IntegrationStrategy() = default;

    IntegrationStrategyResult try_integrate(
        const SymbolicExpr& expr,
        const std::string& var,
        Integrator& ctx,
        ComputationContext& computation,
        int depth = 0);

protected:
    /**
     * @brief 尝试对表达式进行积分
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param ctx 积分器上下文
     * @param depth 当前递归深度
     * @return 积分结果，失败返回 nullptr
     */
    virtual std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr,
        const std::string& var,
        Integrator& ctx,
        ComputationContext& computation,
        int depth = 0) = 0;

public:
    /**
     * @brief 获取策略名称
     * @return 策略名称字符串
     */
    virtual std::string name() const = 0;
    virtual bool requires_residual_verification() const noexcept {
        return true;
    }
};

class BuiltInIntegrationStrategy : public IntegrationStrategy {
public:
    bool requires_residual_verification() const noexcept override {
        return false;
    }
};

/** @brief 查表积分策略 */
class LAMINA_API TableLookupStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "TableLookup"; }
};

/** @brief 幂律积分策略 */
class LAMINA_API PowerRuleStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "PowerRule"; }
};

/** @brief 换元积分策略 */
class LAMINA_API SubstitutionStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "Substitution"; }
};

/**
 * @brief 线性换元积分策略：识别 f(a*x + b) 形式并归约到表项。
 *
 * 当被积函数形如 FunctionNode(arg) 或 PowerNode(base, exp)，且其内部参数可拆解为
 * a*var + b（其中 a、b 关于 var 为常数且 a ≠ 0）时，将参数替换为 var 进行查表，
 * 再代回 a*var + b 并乘以 1/a 得到结果。该策略直接调用积分表
 * （TableLookupStrategy），使递归边界保持在线性换元层。
 */
class LAMINA_API LinearSubstitutionStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "LinearSubstitution"; }

    /**
     * @brief 将表达式 arg 解析为 a*var + b 的线性形式。
     * @param arg 被解析的表达式
     * @param var 积分变量名
     * @param a_out 输出参数：线性系数（关于 var 为常数且非零）
     * @param b_out 输出参数：常数项（关于 var 为常数）
     * @return 解析成功返回 true；若 arg 不依赖 var、依赖关系为非线性、或 a 为零则返回 false
     */
    bool extract_linear_arg(const SymbolicExpr& arg,
                            const std::string& var,
                            std::shared_ptr<SymbolicExpr>& a_out,
                            std::shared_ptr<SymbolicExpr>& b_out);
};

/**
 * @brief 三角组合积分策略：处理 sin^m(x)·cos^n(x)、tan^n(x)、sec^n(x) 形式。
 *
 * 当被积函数形如：
 *   - sin^m(var)·cos^n(var)（m,n 为非负整数，m+n ≤ 8），或
 *   - tan^n(var)（n 为整数，2 ≤ n ≤ 8），或
 *   - sec^n(var)（n 为偶数，2 ≤ n ≤ 8）
 * 时，使用 Pythagorean 恒等式与降幂公式得到闭式原函数：
 *   - 奇数次幂情形（m 或 n 为奇数）通过 u = cos(x) 或 u = sin(x) 换元求解；
 *   - 双偶数情形（m,n 均为偶数）通过半角恒等式 sin² = (1-cos2x)/2 / cos² = (1+cos2x)/2 递归降阶；
 *   - tan^n 通过 tan^n = tan^(n-2)·(sec²-1) 递归化简；
 *   - 偶数 sec^n 通过分部积分得到的递推公式化简。
 *
 * sin/cos/tan/sec 的内部参数必须恰为积分变量本身（如 sin(2x) 由 LinearSubstitutionStrategy 处理）。
 * 若指数为非整数、负数或超出 m+n > 8 的总度数限制，则返回 nullptr 让链上后续策略尝试。
 */
class LAMINA_API TrigCombinationStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "TrigCombination"; }

    /**
     * @brief 检测被积表达式是否为 sin^m(var)·cos^n(var) 形式。
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param m_out sin 的幂次（非负整数）
     * @param n_out cos 的幂次（非负整数）
     * @return 匹配成功返回 true；若包含其他因子、参数不为 var、或幂次为非整数/负数则返回 false
     */
    bool extract_sin_cos_powers(const SymbolicExpr& expr,
                                const std::string& var,
                                int& m_out, int& n_out);

    /**
     * @brief 检测被积表达式是否为 tan^n(var) 形式。
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param n_out tan 的幂次（正整数）
     * @return 匹配成功返回 true
     */
    bool extract_tan_power(const SymbolicExpr& expr,
                           const std::string& var, int& n_out);

    /**
     * @brief 检测被积表达式是否为 sec^n(var) 形式。
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param n_out sec 的幂次（正整数）
     * @return 匹配成功返回 true
     */
    bool extract_sec_power(const SymbolicExpr& expr,
                           const std::string& var, int& n_out);

private:
    /**
     * @brief 根据 m,n 奇偶性分派到 odd / even 情形的内部分发器（带尺度参数）。
     * @param m sin 幂次
     * @param n cos 幂次
     * @param scale sin/cos 内部的尺度系数 c，使被积表达式为 sin^m(c·var)·cos^n(c·var)。
     *              顶层调用使用 scale=1，半角恒等式递归时变为 2、4、8。
     */
    std::shared_ptr<SymbolicExpr> integrate_sin_m_cos_n(
        int m, int n, long long scale,
        const std::string& var, Integrator& ctx, int depth);

    /**
     * @brief 处理 m 或 n 为奇数的情形：剥离一个 sin/cos，使用毕达哥拉斯恒等式
     *        与 u=cos(c·x) 或 u=sin(c·x) 换元，展开成 u 的多项式后逐项积分。
     */
    std::shared_ptr<SymbolicExpr> integrate_odd_case(
        int m, int n, long long scale,
        const std::string& var, Integrator& ctx, int depth);

    /**
     * @brief 处理 m,n 均为偶数的情形：使用半角恒等式
     *        sin²(c·x) = (1 - cos(2c·x))/2，cos²(c·x) = (1 + cos(2c·x))/2，
     *        将被积函数展开为 cos^k(2c·x) 的多项式，并对每个 k 递归调用
     *        integrate_sin_m_cos_n(0, k, 2·scale, ...)。每次递归总度数减半。
     */
    std::shared_ptr<SymbolicExpr> integrate_even_case(
        int m, int n, long long scale,
        const std::string& var, Integrator& ctx, int depth);

    /**
     * @brief 处理 tan^n(var) 形式：递归地应用 tan^n = tan^(n-2)·(sec²-1)。
     *        基础情形：n=0 → x，n=1 → -ln(cos(x))。
     */
    std::shared_ptr<SymbolicExpr> integrate_tan_power(
        int n, const std::string& var, Integrator& ctx, int depth);

    /**
     * @brief 处理 sec^n(var)（偶数 n）形式：使用归约公式
     *        ∫sec^n = sec^(n-2)·tan/(n-1) + (n-2)/(n-1)·∫sec^(n-2)。
     *        基础情形：n=2 → tan(x)。奇数 n 不在本策略范围内，返回 nullptr。
     */
    std::shared_ptr<SymbolicExpr> integrate_sec_power(
        int n, const std::string& var, Integrator& ctx, int depth);
};

/**
 * @brief 高次有理函数部分分式积分策略：处理分母次数 ≥ 3 的有理函数 P(x)/Q(x)。
 *
 * 工作流程：
 *   1. extract_rational：将被积表达式分解为 P(x)/Q(x)，要求 P、Q 均为关于积分变量
 *      的多项式，系数为有理数；非有理形式时返回 false。
 *   2. poly_divide：当 deg(P) ≥ deg(Q) 时进行多项式长除法，得到商和余式；商部分用幂律
 *      逐项积分，剩余的真分式进入部分分式分解。
 *   3. factor_denominator：用 square_free_factorization 加上有理根定理，把 Q 分解为
 *      线性因子和有理域上的二次整体因子；高次整体因子使分解阶段返回 false。
 *   4. solve_coefficients：对每个因子按重数引入待定系数（线性因子贡献 A_l/(x-r)^l，
 *      二次因子贡献 (B_l·x + C_l)/q(x)^l），通过比较多项式系数列出线性方程组，使用
 *      gaussian_eliminate 求解。
 *   5. integrate_term：逐项积分：
 *      - A/(x-r)^l：l=1 生成 ln|x-r|，l>1 生成幂函数；
 *      - (Bx+C)/q(x)：分裂为 ln 和 arctan 两部分（完全平方）；
 *      - q^l（l>=2）映射为保留原语义的未求值积分节点。
 *
 * 支持域映射：
 *   - 零分母、因式分解未决与线性系统未决映射为未求值积分节点；
 *   - 其他表达式返回 nullptr，由策略链继续匹配。
 */
class LAMINA_API RationalDecompositionStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "RationalDecomposition"; }

    /**
     * @brief 把被积表达式拆解为 P(x)/Q(x)（P、Q 均为有理系数多项式）。
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param P_out 输出参数：分子多项式
     * @param Q_out 输出参数：分母多项式
     * @return 成功识别为有理函数返回 true；表达式含 sin/exp/log/非整数指数等
     *         非多项式部分时返回 false
     */
    bool extract_rational(const SymbolicExpr& expr, const std::string& var,
                          Polynomial<Rational>& P_out,
                          Polynomial<Rational>& Q_out);

    /**
     * @brief 多项式长除法：P = Q·quotient + remainder。
     * @param P 被除多项式
     * @param Q 除多项式（前置条件：Q 为非零多项式）
     * @param quotient_out 输出：商多项式
     * @param remainder_out 输出：余式多项式（deg(remainder) < deg(Q)）
     */
    void poly_divide(const Polynomial<Rational>& P, const Polynomial<Rational>& Q,
                     Polynomial<Rational>& quotient_out,
                     Polynomial<Rational>& remainder_out);

    /**
     * @brief 把 Q 分解为线性 (x-r) 与不可约二次 (x²+px+q) 因子的集合。
     * @param Q 分母多项式
     * @param factors_out 输出：(因子, 重数) 列表；每个因子是首一线性或首一不可约二次
     * @return 线性/二次因子完整覆盖时返回 true；高次整体因子或剩余因子使结果为 false
     */
    bool factor_denominator(const Polynomial<Rational>& Q,
                            std::vector<std::pair<Polynomial<Rational>, int>>& factors_out);

    /**
     * @brief 求解部分分式各项的待定系数。
     *
     * 对线性因子 (x-r) 重数 m，引入 m 个常数系数 A_1..A_m；对不可约二次因子 q(x)
     * 重数 m，引入 m 对线性系数 (B_l, C_l)。通过比较多项式系数构造方阵并用
     * gaussian_eliminate 求解。
     *
     * @param P 部分分式分解的目标多项式（真分式的分子，deg(P) < deg(Q)）
     * @param Q 部分分式分解的目标分母
     * @param factors 因子分解（来自 factor_denominator）
     * @param numerators_out 输出：与 (factor_index, l) 一一对应的分子多项式列表，
     *        每个分子要么是常数 (A)，要么是一次多项式 (B·x + C)
     * @return 系数解唯一存在时返回 true；系统无解或退化时返回 false
     */
    bool solve_coefficients(
        const Polynomial<Rational>& P,
        const Polynomial<Rational>& Q,
        const std::vector<std::pair<Polynomial<Rational>, int>>& factors,
        std::vector<Polynomial<Rational>>& numerators_out,
        ComputationContext& context);

    /**
     * @brief 对单个部分分式项 numerator(x)/factor(x)^power 求积分。
     *
     * 支持的形式：
     *   - factor 是 x-r：power=1 → A·ln(x-r)；power≥2 → -A/((power-1)·(x-r)^(power-1))
     *   - factor 是不可约二次 x²+px+q：power=1 → 拆为 (B/2)·ln(q(x)) +
     *     (D-Bp/2)·(2/√(4q-p²))·arctan((2x+p)/√(4q-p²))；power≥2 暂以未求值积分节点返回
     *
     * @param numerator 分子多项式（线性因子时为常数；二次因子时为 B·x+C）
     * @param factor 单个因子多项式（首一线性或首一不可约二次）
     * @param power 因子的幂指数（≥1）
     * @param var 积分变量名
     * @return 该项的原函数表达式
     */
    std::shared_ptr<SymbolicExpr> integrate_term(
        const Polynomial<Rational>& numerator,
        const Polynomial<Rational>& factor,
        int power, const std::string& var);
};

/**
 * @brief 特殊函数积分策略：识别原函数为非初等函数的常见形式。
 *
 * 当被积表达式匹配下列模式之一时，返回以特殊函数（erf, Ei, Si, Ci, Li）表示的结果：
 *   - exp(-x²)        → (√π/2)·erf(x)
 *   - exp(-a·x²)      → (√π/(2·√a))·erf(√a·x)，其中 a 是不依赖积分变量的正常数
 *   - exp(x)/x        → Ei(x)
 *   - sin(x)/x        → Si(x)
 *   - cos(x)/x        → Ci(x)
 *   - 1/ln(x)         → Li(x)
 *
 * 内部参数严格要求为积分变量本身（exp(-x²) 中的 x 必须就是 var；含尺度的二次幂
 * 通过提取 var 的多项式系数识别），其他形式由换元/线性代换等策略处理。
 *
 * 当前规则集之外的输入返回 nullptr，由策略链继续匹配。
 */
class LAMINA_API SpecialFunctionStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "SpecialFunction"; }
};

/** @brief 部分分式积分策略 */
class LAMINA_API PartialFractionStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "PartialFraction"; }
};

/** @brief 分部积分策略 */
class LAMINA_API IBPStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "IntegrationByParts"; }
};

/**
 * @brief 三角换元积分策略（任务 15.1）。
 *
 * 识别含二次根式 √(a²-x²)、√(a²+x²)、√(x²-a²) 的被积函数，应用相应的
 * 三角/双曲换元：
 *   - √(a²-x²): x = a·sin(θ)
 *   - √(a²+x²): x = a·tan(θ)（结果以 arcsinh / ln 形式表达）
 *   - √(x²-a²): x = a·sec(θ)
 * 换元后积分，再用反函数与直角三角关系回代为原变量。
 * 在策略链中位于 SubstitutionStrategy 之后、IBPStrategy 之前。
 */
class LAMINA_API TrigSubstitutionStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "TrigSubstitution"; }
};

/**
 * @brief 万能代换（Weierstrass）积分策略（任务 15.2）。
 *
 * 识别 sin(x)、cos(x) 的有理函数，应用半角代换 t = tan(x/2)：
 *   sin(x) = 2t/(1+t²), cos(x) = (1-t²)/(1+t²), dx = 2/(1+t²) dt
 * 将被积函数转化为 t 的有理函数后递归积分，再回代 t = tan(x/2)。
 * 在策略链中位于 TrigCombinationStrategy 之后。
 */
class LAMINA_API WeierstrassStrategy : public BuiltInIntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate_raw(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
        ComputationContext& computation, int depth = 0) override;
    std::string name() const override { return "Weierstrass"; }
};

} // namespace lamina
