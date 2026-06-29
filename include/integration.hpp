/**
 * @file integration.hpp
 * @brief 符号积分引擎：策略模式，支持查表、幂律、换元、部分分式、分部积分。
 */
#pragma once
#include "symbolic.hpp"
#include "matcher.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <unordered_map>

namespace lamina {

// Forward declaration for optional assumption context integration
class AssumptionContext;

/** @brief 积分表条目，描述一条积分规则（模式 → 结果） */
struct LAMINA_API IntegrationEntry {

    std::string name;                    ///< 规则名称

    SymbolicExpr pattern;                ///< 被积表达式的匹配模式

    SymbolicExpr result;                 ///< 积分结果模板

    std::unordered_set<std::string> wildcards; ///< 模式中的通配符集合

    std::function<bool(const MatchMap&, const std::string& var)> condition; ///< 附加匹配条件

    int priority = 100;                  ///< 优先级，数值越小越优先

    IntegrationEntry() = default;
    IntegrationEntry(std::string name, SymbolicExpr pat, SymbolicExpr res,
                     std::unordered_set<std::string> wc,
                     std::function<bool(const MatchMap&, const std::string& var)> cond = nullptr,
                     int prio = 100)
        : name(std::move(name)), pattern(std::move(pat)), result(std::move(res)),
          wildcards(std::move(wc)), condition(std::move(cond)), priority(prio) {}
};

/** @brief 积分查找表，按类别管理积分规则 */
class LAMINA_API IntegrationTable {
public:
    /** @brief 积分规则类别 */
    enum class Category {
        Polynomial,    ///< 多项式
        Exponential,   ///< 指数函数
        Logarithmic,   ///< 对数函数
        Trigonometric, ///< 三角函数
        InverseTrig,   ///< 反三角函数
        Hyperbolic,    ///< 双曲函数
        Algebraic,     ///< 代数函数
        Special,       ///< 特殊函数
        UserDefined    ///< 用户自定义
    };

    IntegrationTable();

    /**
     * @brief 向指定类别添加积分规则
     * @param cat 规则类别
     * @param entry 积分条目
     */
    void add_entry(Category cat, const IntegrationEntry& entry);

    /**
     * @brief 清空指定类别的所有规则
     * @param cat 规则类别
     */
    void clear_category(Category cat);

    /**
     * @brief 获取指定类别的所有规则
     * @param cat 规则类别
     * @return 该类别下的积分条目列表
     */
    const std::vector<IntegrationEntry>& get_entries(Category cat) const;

    /**
     * @brief 获取所有规则，按优先级排序
     * @return 排序后的积分条目指针列表
     */
    std::vector<const IntegrationEntry*> get_all_sorted() const;

    /** @brief 加载默认积分规则表 */
    void load_defaults();

private:
    std::unordered_map<int, std::vector<IntegrationEntry>> entries_;
    static const std::vector<IntegrationEntry> empty_entries_;
};

class Integrator;

/** @brief 积分策略基类，定义策略接口 */
class LAMINA_API IntegrationStrategy {
public:
    virtual ~IntegrationStrategy() = default;

    /**
     * @brief 尝试对表达式进行积分
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param ctx 积分器上下文
     * @param depth 当前递归深度
     * @return 积分结果，失败返回 nullptr
     */
    virtual std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr,
        const std::string& var,
        Integrator& ctx,
        int depth = 0) = 0;

    /**
     * @brief 获取策略名称
     * @return 策略名称字符串
     */
    virtual std::string name() const = 0;
};

/** @brief 查表积分策略 */
class LAMINA_API TableLookupStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "TableLookup"; }
};

/** @brief 幂律积分策略 */
class LAMINA_API PowerRuleStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "PowerRule"; }
};

/** @brief 换元积分策略 */
class LAMINA_API SubstitutionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "Substitution"; }
};

/**
 * @brief 线性换元积分策略：识别 f(a*x + b) 形式并归约到表项。
 *
 * 当被积函数形如 FunctionNode(arg) 或 PowerNode(base, exp)，且其内部参数可拆解为
 * a*var + b（其中 a、b 关于 var 为常数且 a ≠ 0）时，将参数替换为 var 进行查表，
 * 再代回 a*var + b 并乘以 1/a 得到结果。该策略仅依赖积分表（TableLookupStrategy），
 * 不再次进入完整的策略链，避免与一般换元策略产生递归。
 */
class LAMINA_API LinearSubstitutionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
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
class LAMINA_API TrigCombinationStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
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
 *      线性因子和不可约二次因子（在有理数域上）。无法在有理数上完成分解时返回 false。
 *   4. solve_coefficients：对每个因子按重数引入待定系数（线性因子贡献 A_l/(x-r)^l，
 *      二次因子贡献 (B_l·x + C_l)/q(x)^l），通过比较多项式系数列出线性方程组，使用
 *      gaussian_eliminate 求解。
 *   5. integrate_term：逐项积分：
 *        - A/(x-r)：A·ln|x-r|；A/(x-r)^l (l≥2)：-A/((l-1)(x-r)^(l-1))；
 *        - (Bx+C)/q(x)：分裂为 ln 和 arctan 两部分（完全平方）；
 *        - 不可约二次因子的高次幂 q^l (l≥2) 暂以未求值积分节点保留。
 *
 * 失败保护：
 *   - 分母 Q 为零多项式时直接返回未求值积分节点；
 *   - 表达式不是有理函数时返回 nullptr，让链上后续策略尝试；
 *   - 因式分解或线性方程组无解时返回未求值积分节点保留原积分形式。
 */
class LAMINA_API RationalDecompositionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
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
     * @param Q 除多项式（不能为零）
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
     * @return 全部分解为线性/不可约二次因子返回 true；含有更高次不可约因子或
     *         有理根定理无法穷尽线性因子时返回 false
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
    bool solve_coefficients(const Polynomial<Rational>& P,
                            const Polynomial<Rational>& Q,
                            const std::vector<std::pair<Polynomial<Rational>, int>>& factors,
                            std::vector<Polynomial<Rational>>& numerators_out);

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
 * 任何无法匹配的输入返回 nullptr，让策略链上的后续策略尝试。
 */
class LAMINA_API SpecialFunctionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "SpecialFunction"; }
};

/** @brief 部分分式积分策略 */
class LAMINA_API PartialFractionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "PartialFraction"; }
};

/** @brief 分部积分策略 */
class LAMINA_API IBPStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
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
class LAMINA_API TrigSubstitutionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
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
class LAMINA_API WeierstrassStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "Weierstrass"; }
};

/** @brief 符号积分器，协调各策略完成积分运算 */
class LAMINA_API Integrator {
public:
    Integrator();

    /**
     * @brief 计算不定积分
     * @param expr 被积表达式
     * @param var_name 积分变量名
     * @return 积分结果表达式
     */
    SymbolicExpr integrate(const SymbolicExpr& expr, const std::string& var_name);

    /**
     * @brief 计算定积分
     * @param expr 被积表达式
     * @param var_name 积分变量名
     * @param lower 积分下限
     * @param upper 积分上限
     * @return 定积分结果表达式
     */
    SymbolicExpr integrate_def(const SymbolicExpr& expr, const std::string& var_name,
                               const SymbolicExpr& lower, const SymbolicExpr& upper);

    /**
     * @brief 添加积分策略
     * @param strategy 策略对象
     * @param position 插入位置，-1 表示追加到末尾
     */
    void add_strategy(std::unique_ptr<IntegrationStrategy> strategy, int position = -1);

    /** @brief 获取积分表（可修改） */
    IntegrationTable& table() { return table_; }
    /** @brief 获取积分表（只读） */
    const IntegrationTable& table() const { return table_; }

    /**
     * @brief 递归积分入口
     * @param expr 被积表达式
     * @param var 积分变量名
     * @param depth 当前递归深度
     * @return 积分结果，失败返回 nullptr
     */
    std::shared_ptr<SymbolicExpr> integrate_recursive(
        const SymbolicExpr& expr, const std::string& var, int depth = 0);

    /**
     * @brief 判断表达式是否依赖指定变量
     * @param expr 表达式
     * @param var 变量名
     * @return 依赖返回 true
     */
    static bool depends_on(const SymbolicExpr& expr, const std::string& var);

    /** @brief 获取最大递归深度 */
    int max_depth() const { return max_depth_; }
    /**
     * @brief 设置最大递归深度
     * @param d 深度值
     */
    void set_max_depth(int d) { max_depth_ = d; }

    /**
     * @brief 设置可选的假设上下文，用于指导积分策略选择和被积函数化简。
     *
     * 当提供非空上下文时，积分器可利用变量的符号属性（如正性、非零性）
     * 来简化被积函数（例如将 |x| 简化为 x）或跳过不必要的分情况讨论。
     * 传入 nullptr 恢复默认行为（与未设置上下文时完全一致）。
     *
     * @param ctx 指向 AssumptionContext 的常量指针，nullptr 表示不使用假设
     */
    void set_assumption_context(const AssumptionContext* ctx) { assumption_ctx_ = ctx; }

    /** @brief 获取当前假设上下文（可能为 nullptr） */
    const AssumptionContext* assumption_context() const { return assumption_ctx_; }

private:
    IntegrationTable table_;
    std::vector<std::unique_ptr<IntegrationStrategy>> strategies_;

    struct CycleState {
        std::vector<SymbolicExpr> history;
    };
    CycleState cycle_state_;

    int max_depth_ = 8;

    const AssumptionContext* assumption_ctx_ = nullptr;

    std::shared_ptr<SymbolicExpr> apply_linearity(
        const SymbolicExpr& expr, const std::string& var);

    std::shared_ptr<SymbolicExpr> dispatch_strategies(
        const SymbolicExpr& expr, const std::string& var, int depth);

    static std::shared_ptr<SymbolicExpr> make_integral_node(
        const SymbolicExpr& expr, const std::string& var);

    std::shared_ptr<SymbolicExpr> check_cycle(
        const SymbolicExpr& expr, const std::string& var);
    void resolve_cycle(std::shared_ptr<SymbolicExpr>& result, size_t cycle_idx);
};

/**
 * @brief 多重（迭代）积分引擎：按由内到外的顺序依次对若干积分变量求积分。
 *
 * 给定一个被积表达式和一组按内层（index 0）到外层（最后一个 index）排序的
 * IntegrationStep（每一步包含变量名、可选的下界与上界），引擎会依次对每个
 * 变量调用 Integrator::integrate（不定积分，下/上界均为 nullptr）或
 * Integrator::integrate_def（定积分，下/上界均给出）。如果某一步产生未求值
 * 的积分节点（Calculus_Integral），引擎会停止后续步骤，返回当前的部分结果。
 *
 * 当被积表达式不依赖于某个具有定积分边界 [a,b] 的变量时，引擎会直接将当前
 * 表达式乘以 (b - a)，避免不必要的形式积分。
 *
 * 输入约束：步骤数为 1~3；变量名不能重复；每一步的 lower/upper 必须同时为
 * nullptr（不定积分）或同时非空（定积分）。任何违反约束的输入返回 nullptr。
 */
class LAMINA_API MultipleIntegralEngine {
public:
    /** @brief 单步迭代积分参数：变量与可选的定积分上下界。 */
    struct IntegrationStep {
        std::string variable;                ///< 积分变量名
        std::shared_ptr<SymbolicExpr> lower; ///< 下界，nullptr 表示不定积分
        std::shared_ptr<SymbolicExpr> upper; ///< 上界，nullptr 表示不定积分
    };

    /**
     * @brief 计算迭代积分。
     * @param integrand 被积表达式
     * @param steps 由内层到外层的积分步骤列表（1~3 个，无重复变量）
     * @param integrator 用于实际执行单变量积分的 Integrator 实例
     * @return 迭代积分结果；若输入非法或任一步骤未能求值则可能返回 nullptr 或
     *         保留未求值积分节点的部分结果
     */
    std::shared_ptr<SymbolicExpr> evaluate(
        const SymbolicExpr& integrand,
        const std::vector<IntegrationStep>& steps,
        Integrator& integrator);

private:
    /**
     * @brief 校验步骤列表是否合法。
     *
     * 当满足以下全部条件时返回 true：
     *   - 步骤数为 1~3；
     *   - 变量名互不相同；
     *   - 每一步的 lower 和 upper 同时为 nullptr 或同时非空。
     */
    bool validate(const std::vector<IntegrationStep>& steps) const;
};

}
