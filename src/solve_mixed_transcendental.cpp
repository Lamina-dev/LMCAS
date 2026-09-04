/**
 * @file solve_mixed_transcendental.cpp
 * @brief 混合超越方程求解模块实现.
 */
#include "solve_mixed_transcendental.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "solve_polynomial.hpp"
#include "solve_transcendental.hpp"
#include "symbolic_ast.hpp"
#include "transcendental_factor.hpp"
#include "visitors/differentiation_visitor.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>
#include <new>
#include <stdexcept>

namespace lamina {


/**
 * @internal
 * @brief 判断 FunctionNode::FuncType 是否为超越函数(sin/cos/tan/exp/ln).
 */
static bool is_transcendental_func(FunctionNode::FuncType t) {
    switch (t) {
    case FunctionNode::FuncType::Sin:
    case FunctionNode::FuncType::Cos:
    case FunctionNode::FuncType::Tan:
    case FunctionNode::FuncType::Exp:
    case FunctionNode::FuncType::Ln:
        return true;
    default:
        return false;
    }
}

bool contains_transcendental_of_var(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var)
{
    if (!expr || !lamina::detail::node(expr)) return false;

    struct TranscendentalDetector : public lamina::detail::RecursiveSymbolicVisitor {
        const std::string& target_var;
        bool found = false;

        explicit TranscendentalDetector(const std::string& v) : target_var(v) {}

        void visit(const NumberNode&) override {}
        void visit(const VariableNode&) override {}
        void visit(const MatrixNode&) override {}
        void visit(const RelationalNode& n) override {
            if (found) return;
            if (n.left()) n.left()->accept(*this);
            if (!found && n.right()) n.right()->accept(*this);
        }
        void visit(const LogicalNode& n) override {
            if (found) return;
            if (n.left()) n.left()->accept(*this);
            if (!found && n.right()) n.right()->accept(*this);
        }
        void visit(const PiecewiseNode& n) override {
            for (const auto& branch : n.branches()) {
                if (found) return;
                branch.expression->accept(*this);
                if (!found) branch.condition->accept(*this);
            }
            if (!found && n.default_expr()) n.default_expr()->accept(*this);
        }
        void visit(const SummationNode& n) override {
            n.body()->accept(*this);
            if (!found) n.lower_bound()->accept(*this);
            if (!found) n.upper_bound()->accept(*this);
        }
        void visit(const ProductNode& n) override {
            n.body()->accept(*this);
            if (!found) n.lower_bound()->accept(*this);
            if (!found) n.upper_bound()->accept(*this);
        }
        void visit(const TransformNode& n) override {
            n.body()->accept(*this);
        }
        void visit(const QuantifierNode& n) override {
            n.domain()->accept(*this);
            if (!found) n.predicate()->accept(*this);
        }
        void visit(const SetBuilderNode& n) override {
            n.domain()->accept(*this);
            if (!found) n.predicate()->accept(*this);
        }
        void visit(const ComplexNode& n) override {
            n.real()->accept(*this);
            if (!found) n.imag()->accept(*this);
        }
        void visit(const FiniteSetNode& n) override { for (const auto& e : n.elements()) { if (found) return; e->accept(*this); } }
        void visit(const IntervalNode& n) override { n.lower()->accept(*this); if (!found) n.upper()->accept(*this); }
        void visit(const MembershipNode& n) override { n.element()->accept(*this); if (!found) n.set()->accept(*this); }
        void visit(const QuantityNode& n) override { n.value()->accept(*this); }

        void visit(const AddNode& n) override {
            for (auto& op : n.operands()) {
                if (found) return;
                op->accept(*this);
            }
        }

        void visit(const MultiplyNode& n) override {
            for (auto& op : n.operands()) {
                if (found) return;
                op->accept(*this);
            }
        }

        void visit(const PowerNode& n) override {
            if (found) return;
            n.base()->accept(*this);
            if (!found) n.exponent()->accept(*this);
        }

        void visit(const FunctionNode& n) override {
            if (found) return;

            if (is_transcendental_func(n.type())) {
                /// 检查参数是否依赖目标变量
                for (auto& arg : n.arguments()) {
                    if (expression_depends_on_variable(arg, target_var)) {
                        found = true;
                        return;
                    }
                }
            }

            /// 递归检查参数子树(可能嵌套超越函数)
            for (auto& arg : n.arguments()) {
                if (found) return;
                arg->accept(*this);
            }
        }
    } detector(var);

    lamina::detail::node(expr)->accept(detector);
    return detector.found;
}


/**
 * @internal
 * @brief 计算节点关于指定变量的多项式次数.
 *
 * 与 var 无关的子表达式次数为 0;-1 表示当前结构位于多项式次数支持域之外.
 *
 * @param[in] node 符号节点
 * @param[in] var  变量名
 * @return 关于 var 的次数;-1 表示非多项式结构
 */
static int degree_in_var(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) return 0;

    if (!expression_depends_on_variable(node, var)) return 0;

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return (v->name() == var) ? 1 : 0;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        int max_deg = 0;
        for (const auto& op : add->operands()) {
            int d = degree_in_var(op, var);
            if (d < 0) return -1;
            max_deg = std::max(max_deg, d);
        }
        return max_deg;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int total_deg = 0;
        for (const auto& op : mul->operands()) {
            int d = degree_in_var(op, var);
            if (d < 0) return -1;
            total_deg += d;
        }
        return total_deg;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
        if (!exp_num) return -1;

        int e_val = -1;
        if (std::holds_alternative<BigInt>(exp_num->value())) {
            const auto& bi = std::get<BigInt>(exp_num->value());
            if (!bi.IsNegative()) e_val = bi.to_int();
        } else if (std::holds_alternative<Rational>(exp_num->value())) {
            const auto& r = std::get<Rational>(exp_num->value());
            if (r.is_integer()) {
                BigInt bi = r.to_BigInt();
                if (!bi.IsNegative()) e_val = bi.to_int();
            }
        } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
            lmmc_real_t d = std::get<lmmc_real_t>(exp_num->value());
            if (std::isfinite(d) && d >= 0 && d == std::floor(d) && d < 1000.0) {
                e_val = static_cast<int>(d);
            }
        }

        if (e_val < 0) return -1;

        int base_deg = degree_in_var(pow->base(), var);
        if (base_deg < 0) return -1;
        return base_deg * e_val;
    }

    /// FunctionNode 依赖 var 意味着非多项式结构
    if (std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return -1;
    }

    return 0;
}

bool is_polynomial_after_substitution(
    const TransSubstitutionResult& sub_result)
{
    if (sub_result.mappings.empty()) return false;
    if (!sub_result.poly_expr || !lamina::detail::node(sub_result.poly_expr)) return false;

    const auto& root = lamina::detail::node(sub_result.poly_expr);

    for (const auto& m : sub_result.mappings) {
        if (expression_depends_on_variable(root, m.indeterminate)) {
            int deg = degree_in_var(root, m.indeterminate);
            if (deg < 0) return false;
        }
    }

    return true;
}



/**
 * @internal
 * @brief 从 NumberNode 中提取浮点数值.
 * @param[in] node 数值节点
 * @return 浮点值;非数值节点返回 NaN
 */
static lmmc_real_t extract_real_value(const std::shared_ptr<const SymbolicNode>& node) {
    auto num = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!num) return std::numeric_limits<lmmc_real_t>::quiet_NaN();

    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        return std::get<lmmc_real_t>(num->value());
    }
    if (std::holds_alternative<Rational>(num->value())) {
        return static_cast<lmmc_real_t>(std::get<Rational>(num->value()).to_double());
    }
    if (std::holds_alternative<BigInt>(num->value())) {
        return static_cast<lmmc_real_t>(std::get<BigInt>(num->value()).to_double());
    }
    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

/**
 * @internal
 * @brief 尝试从表达式中提取关于变量的线性系数 k(即 k*x + c 形式中的 k).
 *
 * 支持的模式:
 * - 变量本身 x -> k = 1
 * - 乘法 k*x(数值 * 变量)-> 提取数值
 * - 加法 k*x + c -> 在加法项中查找含变量的项并提取系数
 * - 非线性参数(x^2, sin(x) 等)-> 返回 NaN
 *
 * @param[in] node 参数表达式节点
 * @param[in] var  变量名
 * @return 线性系数 k;非线性时返回 NaN
 */
static lmmc_real_t extract_linear_coefficient(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var)
{
    if (!node || !expression_depends_on_variable(node, var)) return 0.0;

    /// 模式 1: 变量本身 -> k = 1
    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (v->name() == var) return 1.0;
        return std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }

    /// 模式 2: 乘法节点 k*x
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        /// 检查是否恰好是 数值 * 变量 的形式
        if (mul->operands().size() == 2) {
            std::shared_ptr<const SymbolicNode> num_part = nullptr;
            std::shared_ptr<const SymbolicNode> var_part = nullptr;

            for (const auto& op : mul->operands()) {
                if (op->is_number()) {
                    num_part = op;
                } else if (auto vn = std::dynamic_pointer_cast<const VariableNode>(op)) {
                    if (vn->name() == var) var_part = op;
                }
            }

            if (num_part && var_part) {
                return extract_real_value(num_part);
            }
        }
        /// 非简单 k*x 形式 -> 非线性
        return std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }

    /// 模式 3: 加法节点 k*x + c
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        lmmc_real_t coeff = std::numeric_limits<lmmc_real_t>::quiet_NaN();
        bool found_var_term = false;

        for (const auto& op : add->operands()) {
            if (!expression_depends_on_variable(op, var)) {
                /// 常数项,跳过
                continue;
            }
            /// 含变量的项 - 只允许一个
            if (found_var_term) {
                /// 多个含变量的项 -> 非线性
                return std::numeric_limits<lmmc_real_t>::quiet_NaN();
            }
            found_var_term = true;
            coeff = extract_linear_coefficient(op, var);
        }

        return found_var_term ? coeff : std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }

    /// 其他节点类型(PowerNode, FunctionNode 等)-> 非线性
    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

/**
 * @internal
 * @brief 判断 FunctionNode::FuncType 是否为周期函数(sin/cos/tan).
 */
static bool is_periodic_func(FunctionNode::FuncType t) {
    switch (t) {
    case FunctionNode::FuncType::Sin:
    case FunctionNode::FuncType::Cos:
    case FunctionNode::FuncType::Tan:
        return true;
    default:
        return false;
    }
}

/**
 * @internal
 * @brief 遍历 AST 收集所有周期函数的周期信息.
 *
 * 对每个 sin/cos/tan 节点,若其参数为 k*x + c 形式,
 * 计算周期(sin/cos: 2pi/|k|, tan: pi/|k|).
 */
struct PeriodicCollector : public lamina::detail::RecursiveSymbolicVisitor {
    const std::string& target_var;
    lmmc_real_t max_period = 0.0;
    bool found_periodic = false;
    bool has_nonlinear_periodic = false;

    explicit PeriodicCollector(const std::string& v) : target_var(v) {}

    void visit(const NumberNode&) override {}
    void visit(const VariableNode&) override {}
    void visit(const MatrixNode&) override {}
    void visit(const RelationalNode& n) override {
        if (n.left()) n.left()->accept(*this);
        if (n.right()) n.right()->accept(*this);
    }
    void visit(const LogicalNode& n) override {
        if (n.left()) n.left()->accept(*this);
        if (n.right()) n.right()->accept(*this);
    }
    void visit(const PiecewiseNode& n) override {
        for (const auto& branch : n.branches()) {
            branch.expression->accept(*this);
            branch.condition->accept(*this);
        }
        if (n.default_expr()) n.default_expr()->accept(*this);
    }
    void visit(const SummationNode& n) override {
        n.body()->accept(*this);
        n.lower_bound()->accept(*this);
        n.upper_bound()->accept(*this);
    }
    void visit(const ProductNode& n) override {
        n.body()->accept(*this);
        n.lower_bound()->accept(*this);
        n.upper_bound()->accept(*this);
    }
    void visit(const TransformNode& n) override {
        n.body()->accept(*this);
    }
    void visit(const QuantifierNode& n) override {
        n.domain()->accept(*this);
        n.predicate()->accept(*this);
    }
    void visit(const SetBuilderNode& n) override {
        n.domain()->accept(*this);
        n.predicate()->accept(*this);
    }
    void visit(const ComplexNode& n) override {
        n.real()->accept(*this);
        n.imag()->accept(*this);
    }
    void visit(const FiniteSetNode& n) override { for (const auto& e : n.elements()) e->accept(*this); }
    void visit(const IntervalNode& n) override { n.lower()->accept(*this); n.upper()->accept(*this); }
    void visit(const MembershipNode& n) override { n.element()->accept(*this); n.set()->accept(*this); }
    void visit(const QuantityNode& n) override { n.value()->accept(*this); }

    void visit(const AddNode& n) override {
        for (auto& op : n.operands()) op->accept(*this);
    }

    void visit(const MultiplyNode& n) override {
        for (auto& op : n.operands()) op->accept(*this);
    }

    void visit(const PowerNode& n) override {
        n.base()->accept(*this);
        n.exponent()->accept(*this);
    }

    void visit(const FunctionNode& n) override {
        if (is_periodic_func(n.type())) {
            /// 检查参数是否依赖目标变量
            if (!n.arguments().empty() && expression_depends_on_variable(n.arguments()[0], target_var)) {
                found_periodic = true;

                lmmc_real_t k = extract_linear_coefficient(n.arguments()[0], target_var);
                if (std::isnan(k) || k == 0.0) {
                    /// 非线性参数
                    has_nonlinear_periodic = true;
                } else {
                    lmmc_real_t period = 0.0;
                    if (n.type() == FunctionNode::FuncType::Tan) {
                        period = LMMC_CONST_PI / std::fabs(k);
                    } else {
                        /// sin/cos
                        period = 2.0 * LMMC_CONST_PI / std::fabs(k);
                    }
                    if (period > max_period) {
                        max_period = period;
                    }
                }
            }
        }

        /// 递归检查参数子树
        for (auto& arg : n.arguments()) {
            arg->accept(*this);
        }
    }
};

std::optional<SearchInterval> determine_search_interval(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts)
{
    /// 优先级 1: 用户指定区间
    if (opts.has_search_interval) {
        if (opts.search_lo >= opts.search_hi) return std::nullopt;
        if ((opts.search_hi - opts.search_lo) <= opts.tolerance) return std::nullopt;
        return SearchInterval{opts.search_lo, opts.search_hi};
    }

    /// 优先级 2: 周期扩展
    lmmc_real_t lo = -10.0;
    lmmc_real_t hi = 10.0;

    if (expr && lamina::detail::node(expr)) {
        PeriodicCollector collector(var);
        lamina::detail::node(expr)->accept(collector);

        if (collector.found_periodic && !collector.has_nonlinear_periodic && collector.max_period > 0.0) {
            /// 扩展区间覆盖 2 个完整周期(对称于 0)
            lmmc_real_t half_span = collector.max_period;  // 2 periods / 2 = period (symmetric)
            if (half_span > lo * -1.0 || half_span > hi) {
                lo = -half_span;
                hi = half_span;
            }
            /// 钳制到 [-100, 100]
            if (lo < -100.0) lo = -100.0;
            if (hi > 100.0) hi = 100.0;
        }
    }

    /// 最终验证
    if (lo >= hi) return std::nullopt;
    if ((hi - lo) <= opts.tolerance) return std::nullopt;

    return SearchInterval{lo, hi};
}


/**
 * @internal
 * @brief 在指定点对表达式进行数值求值.
 *
 * 将变量替换为数值后调用 to_numeric() 获取浮点结果.
 * 若求值产生 NaN 或无穷大,返回 NaN.
 *
 * @param[in] expr 待求值表达式
 * @param[in] var  变量名
 * @param[in] x    求值点
 * @return 数值结果;异常时返回 NaN
 */
/**
 * @internal
 * @brief 递归数值求值:遍历 AST 计算表达式的浮点值.
 *
 * 处理 NumberNode,FunctionNode,AddNode,MultiplyNode,PowerNode;
 * NaN 表示节点位于当前数值求值支持域之外.
 */
static lmmc_real_t recursive_eval(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return 0.0;

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (std::holds_alternative<lmmc_real_t>(num->value()))
            return std::get<lmmc_real_t>(num->value());
        if (std::holds_alternative<BigInt>(num->value()))
            return static_cast<lmmc_real_t>(std::get<BigInt>(num->value()).to_double());
        if (std::holds_alternative<Rational>(num->value()))
            return static_cast<lmmc_real_t>(std::get<Rational>(num->value()).to_double());
        return 0.0;
    }

    if (std::dynamic_pointer_cast<const VariableNode>(node)) {
        return std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        lmmc_real_t sum = 0.0;
        for (const auto& op : add->operands()) {
            lmmc_real_t v = recursive_eval(op);
            if (std::isnan(v)) return v;
            sum += v;
        }
        return sum;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        lmmc_real_t product = 1.0;
        for (const auto& op : mul->operands()) {
            lmmc_real_t v = recursive_eval(op);
            if (std::isnan(v)) return v;
            product *= v;
        }
        return product;
    }

    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        lmmc_real_t base = recursive_eval(pw->base());
        lmmc_real_t exp = recursive_eval(pw->exponent());
        if (std::isnan(base) || std::isnan(exp))
            return std::numeric_limits<lmmc_real_t>::quiet_NaN();
        if (base == 0.0 && exp < 0.0)
            return std::numeric_limits<lmmc_real_t>::quiet_NaN();
        return std::pow(base, exp);
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->arguments().size() == 1) {
            lmmc_real_t arg = recursive_eval(func->arguments()[0]);
            if (std::isnan(arg)) return arg;
            switch (func->type()) {
                case FunctionNode::FuncType::Sin: return std::sin(arg);
                case FunctionNode::FuncType::Cos: return std::cos(arg);
                case FunctionNode::FuncType::Tan: return std::tan(arg);
                case FunctionNode::FuncType::Exp: return std::exp(arg);
                case FunctionNode::FuncType::Ln:  return std::log(arg);
                case FunctionNode::FuncType::Sqrt: return std::sqrt(arg);
                case FunctionNode::FuncType::Abs: return std::abs(arg);
                case FunctionNode::FuncType::ArcSin:
                    if (arg < -1.0 || arg > 1.0) return std::numeric_limits<lmmc_real_t>::quiet_NaN();
                    return std::asin(arg);
                case FunctionNode::FuncType::ArcCos:
                    if (arg < -1.0 || arg > 1.0) return std::numeric_limits<lmmc_real_t>::quiet_NaN();
                    return std::acos(arg);
                case FunctionNode::FuncType::ArcTan: return std::atan(arg);
                case FunctionNode::FuncType::Sinh: return std::sinh(arg);
                case FunctionNode::FuncType::Cosh: return std::cosh(arg);
                case FunctionNode::FuncType::Tanh: return std::tanh(arg);
                case FunctionNode::FuncType::Sec: {
                    lmmc_real_t c = std::cos(arg);
                    if (std::fabs(c) < 1e-15) return std::numeric_limits<lmmc_real_t>::quiet_NaN();
                    return 1.0 / c;
                }
                case FunctionNode::FuncType::Csc: {
                    lmmc_real_t s = std::sin(arg);
                    if (std::fabs(s) < 1e-15) return std::numeric_limits<lmmc_real_t>::quiet_NaN();
                    return 1.0 / s;
                }
                case FunctionNode::FuncType::Cot: {
                    lmmc_real_t s = std::sin(arg);
                    if (std::fabs(s) < 1e-15) return std::numeric_limits<lmmc_real_t>::quiet_NaN();
                    return std::cos(arg) / s;
                }
                default: break;
            }
        }
        if (func->arguments().size() == 2 && func->type() == FunctionNode::FuncType::Atan2) {
            lmmc_real_t y = recursive_eval(func->arguments()[0]);
            lmmc_real_t x_val = recursive_eval(func->arguments()[1]);
            if (std::isnan(y) || std::isnan(x_val)) return std::numeric_limits<lmmc_real_t>::quiet_NaN();
            return std::atan2(y, x_val);
        }
    }

    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

static lmmc_real_t evaluate_at(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    lmmc_real_t x)
{
    try {
        auto substituted = expr->substitute(var, SymbolicExpr::number(static_cast<double>(x)));
        if (!substituted || !lamina::detail::node(substituted))
            return std::numeric_limits<lmmc_real_t>::quiet_NaN();
        lmmc_real_t val = recursive_eval(lamina::detail::node(substituted));
        if (!std::isfinite(val)) return std::numeric_limits<lmmc_real_t>::quiet_NaN();
        return val;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::exception&) {
        return std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }
}

/**
 * @internal
 * @brief 在指定点求值,若结果为 NaN 则尝试在附近重新采样.
 *
 * 最多重试 max_retries 次,每次将偏移量减半并向两侧尝试.
 *
 * @param[in] expr        表达式
 * @param[in] var         变量名
 * @param[in] x           原始采样点
 * @param[in] half_width  初始偏移量(子区间宽度的一半)
 * @param[in] max_retries 最大重试次数
 * @return 有效数值结果;所有重试失败时返回 NaN
 */
static lmmc_real_t evaluate_with_retry(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    lmmc_real_t x,
    lmmc_real_t half_width,
    int max_retries = 3)
{
    lmmc_real_t val = evaluate_at(expr, var, x);
    if (!std::isnan(val)) return val;

    lmmc_real_t offset = half_width;
    for (int i = 0; i < max_retries; ++i) {
        offset *= 0.5;
        /// 尝试右侧
        val = evaluate_at(expr, var, x + offset);
        if (!std::isnan(val)) return val;
        /// 尝试左侧
        val = evaluate_at(expr, var, x - offset);
        if (!std::isnan(val)) return val;
    }
    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

static bool consume_mixed_step(
    ComputationContext* context,
    std::optional<CasError>* failure)
{
    if (!context) return true;
    auto step = context->consume_steps(
        1, "solve_mixed_transcendental_checked");
    if (step) return true;
    if (failure) *failure = step.error();
    return false;
}


static std::vector<IsolatedInterval> isolate_roots_with_context(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var,
    const SearchInterval& interval,
    const SolveOptions& opts,
    ComputationContext* context,
    std::optional<CasError>* failure,
    bool* complete)
{
    std::vector<IsolatedInterval> result;

    /// 无效区间检查
    if (interval.lo >= interval.hi) return result;

    int max_roots_limit = opts.max_roots;
    constexpr lmmc_real_t MIN_WIDTH = 1e-6;

    /// 使用显式栈执行自适应细分,使递归深度保持恒定.
    struct SubInterval {
        lmmc_real_t lo;
        lmmc_real_t hi;
        lmmc_real_t f_lo;
        lmmc_real_t f_hi;
    };

    /// 初始均匀划分为 N 个子区间
    constexpr int INITIAL_DIVISIONS = 64;
    lmmc_real_t total_width = interval.hi - interval.lo;
    lmmc_real_t step = total_width / INITIAL_DIVISIONS;

    /// 预计算初始采样点的函数值
    std::vector<lmmc_real_t> sample_vals(INITIAL_DIVISIONS + 1);
    for (int i = 0; i <= INITIAL_DIVISIONS; ++i) {
        if (!consume_mixed_step(context, failure)) {
            if (complete) *complete = false;
            return result;
        }
        lmmc_real_t x = interval.lo + i * step;
        sample_vals[i] = evaluate_with_retry(expr, var, x, step * 0.5);
    }

    /// 收集需要进一步处理的子区间(含符号变化或需要细分的)
    std::vector<SubInterval> work_stack;
    for (int i = 0; i < INITIAL_DIVISIONS; ++i) {
        lmmc_real_t a = interval.lo + i * step;
        lmmc_real_t b = interval.lo + (i + 1) * step;
        lmmc_real_t fa = sample_vals[i];
        lmmc_real_t fb = sample_vals[i + 1];

        /// 跳过两端都是 NaN 的子区间
        if (std::isnan(fa) && std::isnan(fb)) continue;

        /// 若一端为 NaN,尝试缩小区间
        if (std::isnan(fa) || std::isnan(fb)) {
            /// 将此子区间加入工作栈进行进一步细分
            work_stack.push_back({a, b, fa, fb});
            continue;
        }

        /// 检测符号变化
        if (fa * fb < 0.0) {
            work_stack.push_back({a, b, fa, fb});
        }
    }

    /// 自适应细分处理工作栈
    while (!work_stack.empty()) {
        if (!consume_mixed_step(context, failure)) {
            if (complete) *complete = false;
            break;
        }
        /// 检查 max_roots 限制
        if (max_roots_limit > 0 &&
            static_cast<int>(result.size()) >= max_roots_limit) {
            if (complete) *complete = false;
            break;
        }

        SubInterval current = work_stack.back();
        work_stack.pop_back();

        lmmc_real_t width = current.hi - current.lo;
        lmmc_real_t fa = current.f_lo;
        lmmc_real_t fb = current.f_hi;

        /// 处理 NaN 端点:尝试重新采样
        if (std::isnan(fa)) {
            fa = evaluate_with_retry(expr, var, current.lo, width * 0.25);
            if (std::isnan(fa)) {
                /// 若宽度足够,细分后重试
                if (width > MIN_WIDTH * 2.0) {
                    lmmc_real_t mid = (current.lo + current.hi) * 0.5;
                    lmmc_real_t fm = evaluate_with_retry(expr, var, mid, width * 0.125);
                    work_stack.push_back({mid, current.hi, fm, fb});
                }
                continue;
            }
        }
        if (std::isnan(fb)) {
            fb = evaluate_with_retry(expr, var, current.hi, width * 0.25);
            if (std::isnan(fb)) {
                if (width > MIN_WIDTH * 2.0) {
                    lmmc_real_t mid = (current.lo + current.hi) * 0.5;
                    lmmc_real_t fm = evaluate_with_retry(expr, var, mid, width * 0.125);
                    work_stack.push_back({current.lo, mid, fa, fm});
                }
                continue;
            }
        }

        /// 无符号变化 -> 跳过
        if (fa * fb >= 0.0) continue;

        /// 符号变化确认 - 检查是否需要进一步细分
        /// 若宽度已达最小限制,直接接受
        if (width <= MIN_WIDTH) {
            bool confirmed = false;
            if (derivative) {
                lmmc_real_t da = evaluate_at(derivative, var, current.lo);
                lmmc_real_t db = evaluate_at(derivative, var, current.hi);
                if (!std::isnan(da) && !std::isnan(db) && da * db > 0.0) {
                    confirmed = true;
                }
            }
            result.push_back({current.lo, current.hi, confirmed});
            continue;
        }

        /// 单调性确认:若导数在区间内不变号,则恰含一个根
        if (derivative) {
            lmmc_real_t da = evaluate_at(derivative, var, current.lo);
            lmmc_real_t db = evaluate_at(derivative, var, current.hi);

            if (!std::isnan(da) && !std::isnan(db) && da * db > 0.0) {
                /// 导数不变号 -> 单调 -> 恰含一个根
                result.push_back({current.lo, current.hi, true});
                continue;
            }

            /// 导数变号 -> 可能含多个根,需要细分
            if (!std::isnan(da) && !std::isnan(db) && da * db < 0.0) {
                /// 仅在宽度允许时细分
                if (width > MIN_WIDTH * 2.0) {
                    lmmc_real_t mid = (current.lo + current.hi) * 0.5;
                    lmmc_real_t fm = evaluate_with_retry(expr, var, mid, width * 0.125);

                    if (std::isnan(fm)) {
                        /// 中点求值失败,接受当前区间
                        result.push_back({current.lo, current.hi, false});
                    } else {
                        /// 将两半加入工作栈
                        if (fa * fm < 0.0) {
                            work_stack.push_back({current.lo, mid, fa, fm});
                        }
                        if (fm * fb < 0.0) {
                            work_stack.push_back({mid, current.hi, fm, fb});
                        }
                        /// 若两半都无符号变化但原区间有,说明根恰在中点附近
                        if (fa * fm >= 0.0 && fm * fb >= 0.0) {
                            result.push_back({current.lo, current.hi, false});
                        }
                    }
                    continue;
                }
            }

            /// 导数包含 NaN 时保留区间,并将单调性标记为未确认.
            result.push_back({current.lo, current.hi, false});
        } else {
            /// 导数信息未决时继续细分,以缩小候选区间.
            if (width > MIN_WIDTH * 4.0) {
                lmmc_real_t mid = (current.lo + current.hi) * 0.5;
                lmmc_real_t fm = evaluate_with_retry(expr, var, mid, width * 0.125);

                if (std::isnan(fm)) {
                    /// 中点求值失败,接受当前区间
                    result.push_back({current.lo, current.hi, false});
                } else {
                    if (fa * fm < 0.0) {
                        work_stack.push_back({current.lo, mid, fa, fm});
                    }
                    if (fm * fb < 0.0) {
                        work_stack.push_back({mid, current.hi, fm, fb});
                    }
                    if (fa * fm >= 0.0 && fm * fb >= 0.0) {
                        result.push_back({current.lo, current.hi, false});
                    }
                }
            } else {
                /// 宽度已足够小,接受
                result.push_back({current.lo, current.hi, false});
            }
        }
    }

    /// 按区间下界排序
    std::sort(result.begin(), result.end(),
        [](const IsolatedInterval& a, const IsolatedInterval& b) {
            return a.lo < b.lo;
        });

    if (max_roots_limit > 0 &&
        static_cast<int>(result.size()) > max_roots_limit) {
        if (complete) *complete = false;
        result.resize(static_cast<size_t>(max_roots_limit));
    }

    return result;
}

std::vector<IsolatedInterval> isolate_roots(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var,
    const SearchInterval& interval,
    const SolveOptions& opts)
{
    return isolate_roots_with_context(
        expr, derivative, var, interval, opts, nullptr, nullptr, nullptr);
}


std::vector<lmmc_real_t> deduplicate_roots(
    std::vector<NumericRoot>& roots,
    lmmc_real_t tolerance,
    int max_roots)
{
    if (roots.empty()) return {};

    /// 按根值升序排序
    std::sort(roots.begin(), roots.end(),
        [](const NumericRoot& a, const NumericRoot& b) {
            return a.value < b.value;
        });

    /// 去重:相邻比较,差值 < 10 * tolerance 视为重复,保留残差较小者
    lmmc_real_t threshold = 10.0 * tolerance;
    std::vector<NumericRoot> survivors;
    survivors.reserve(roots.size());
    survivors.push_back(roots[0]);

    for (size_t i = 1; i < roots.size(); ++i) {
        NumericRoot& last = survivors.back();
        if (std::fabs(roots[i].value - last.value) < threshold) {
            /// 重复:保留残差更小者(相等时保留先出现的,即 last)
            if (roots[i].residual < last.residual) {
                last = roots[i];
            }
        } else {
            survivors.push_back(roots[i]);
        }
    }

    /// 收集存活根的值(已按升序排列)
    std::vector<lmmc_real_t> result;
    result.reserve(survivors.size());
    for (const auto& r : survivors) {
        result.push_back(r.value);
    }

    /// 应用 max_roots 限制
    if (max_roots > 0 && static_cast<int>(result.size()) > max_roots) {
        result.resize(static_cast<size_t>(max_roots));
    }

    return result;
}


/**
 * @internal
 * @brief 对隔离区间执行根精化,返回高精度数值根.
 *
 * 算法流程:
 * 1. 以区间中点为初始猜测;
 * 2. 若导数可用,执行带区间约束的 Newton-Raphson 迭代:
 *    - 每步计算 x_new = x - f(x)/f'(x);
 *    - 若 x_new 超出 [lo, hi] 或残差增大,丢弃该迭代并执行一步二分法;
 *    - 若 f'(x) == 0,执行一步二分法;
 * 3. 若导数不可用(为 nullptr),执行纯二分法;
 * 4. 收敛条件:|f(x)| < tolerance;
 * 5. 达到最大迭代次数后,返回残差最小的迭代结果(仅当残差 < tolerance).
 *
 * @param[in] expr       待求解表达式
 * @param[in] derivative 表达式的导数(可为 nullptr,此时使用纯二分法)
 * @param[in] var        求解变量名
 * @param[in] interval   隔离子区间
 * @param[in] opts       求解选项(tolerance, max_newton_iterations)
 * @return 精化后的数值根;未收敛或残差超限时返回 nullopt
 */
static std::optional<NumericRoot> refine_root_with_context(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var,
    const IsolatedInterval& interval,
    const SolveOptions& opts,
    ComputationContext* context,
    std::optional<CasError>* failure,
    bool* complete)
{
    lmmc_real_t lo = interval.lo;
    lmmc_real_t hi = interval.hi;
    lmmc_real_t x = (lo + hi) * 0.5;

    lmmc_real_t f_lo = evaluate_at(expr, var, lo);
    lmmc_real_t f_hi = evaluate_at(expr, var, hi);

    /// 端点已是根的情形
    if (!std::isnan(f_lo) && std::fabs(f_lo) < opts.tolerance) {
        return NumericRoot{lo, std::fabs(f_lo), 0};
    }
    if (!std::isnan(f_hi) && std::fabs(f_hi) < opts.tolerance) {
        return NumericRoot{hi, std::fabs(f_hi), 0};
    }

    /// 跟踪最佳迭代结果(残差最小)
    lmmc_real_t best_x = x;
    lmmc_real_t best_residual = std::numeric_limits<lmmc_real_t>::max();
    int best_iter = 0;

    int max_iter = opts.max_newton_iterations;

    /// 纯二分法路径:导数不可用
    if (!derivative) {
        /// 确保端点有效且异号
        if (std::isnan(f_lo) || std::isnan(f_hi) || f_lo * f_hi > 0.0) {
            /// 尝试中点求值
            lmmc_real_t fx = evaluate_at(expr, var, x);
            if (!std::isnan(fx) && std::fabs(fx) < opts.tolerance) {
                return NumericRoot{x, std::fabs(fx), 1};
            }
            return std::nullopt;
        }

        for (int i = 1; i <= max_iter; ++i) {
            if (!consume_mixed_step(context, failure)) {
                if (complete) *complete = false;
                return std::nullopt;
            }
            lmmc_real_t mid = (lo + hi) * 0.5;
            lmmc_real_t f_mid = evaluate_at(expr, var, mid);

            if (std::isnan(f_mid)) {
                /// NaN 时缩小区间继续
                hi = mid;
                continue;
            }

            lmmc_real_t residual = std::fabs(f_mid);
            if (residual < best_residual) {
                best_x = mid;
                best_residual = residual;
                best_iter = i;
            }

            if (residual < opts.tolerance) {
                return NumericRoot{mid, residual, i};
            }

            /// 保持符号变化的半区间
            if (f_lo * f_mid < 0.0) {
                hi = mid;
            } else {
                lo = mid;
                f_lo = f_mid;
            }

            /// 区间宽度收敛
            if ((hi - lo) < opts.tolerance) {
                lmmc_real_t final_x = (lo + hi) * 0.5;
                lmmc_real_t final_res = std::fabs(evaluate_at(expr, var, final_x));
                if (!std::isnan(final_res) && final_res < best_residual) {
                    best_x = final_x;
                    best_residual = final_res;
                    best_iter = i;
                }
                break;
            }
        }

        /// 返回最佳结果(仅当残差 < tolerance)
        if (best_residual < opts.tolerance) {
            return NumericRoot{best_x, best_residual, best_iter};
        }
        return std::nullopt;
    }

    /// 带区间约束的 Newton-Raphson 路径
    /// 确保端点函数值有效
    if (std::isnan(f_lo) || std::isnan(f_hi)) {
        /// 尝试重新采样端点
        if (std::isnan(f_lo)) {
            f_lo = evaluate_at(expr, var, lo + (hi - lo) * 0.01);
            lo = lo + (hi - lo) * 0.01;
        }
        if (std::isnan(f_hi)) {
            f_hi = evaluate_at(expr, var, hi - (hi - lo) * 0.01);
            hi = hi - (hi - lo) * 0.01;
        }
        if (std::isnan(f_lo) || std::isnan(f_hi)) {
            return std::nullopt;
        }
    }

    for (int i = 1; i <= max_iter; ++i) {
        if (!consume_mixed_step(context, failure)) {
            if (complete) *complete = false;
            return std::nullopt;
        }
        lmmc_real_t fx = evaluate_at(expr, var, x);

        if (std::isnan(fx)) {
            /// NaN 时执行二分步
            x = (lo + hi) * 0.5;
            continue;
        }

        lmmc_real_t residual = std::fabs(fx);

        /// 更新最佳结果
        if (residual < best_residual) {
            best_x = x;
            best_residual = residual;
            best_iter = i;
        }

        /// 收敛判定
        if (residual < opts.tolerance) {
            return NumericRoot{x, residual, i};
        }

        /// 计算导数值
        lmmc_real_t dfx = evaluate_at(derivative, var, x);

        /// 导数为零或 NaN -> 执行二分步
        if (std::isnan(dfx) || std::fabs(dfx) < 1e-15) {
            /// 二分步:保持符号变化
            lmmc_real_t mid = (lo + hi) * 0.5;
            lmmc_real_t f_mid = evaluate_at(expr, var, mid);

            if (!std::isnan(f_mid)) {
                if (!std::isnan(f_lo) && f_lo * f_mid < 0.0) {
                    hi = mid;
                    f_hi = f_mid;
                } else {
                    lo = mid;
                    f_lo = f_mid;
                }
            }
            x = (lo + hi) * 0.5;
            continue;
        }

        /// Newton 更新
        lmmc_real_t x_new = x - fx / dfx;

        /// 检查 Newton 迭代是否有效:
        /// 1. x_new 位于 [lo, hi] 内
        /// 2. 新残差小于等于当前残差
        bool do_bisection = false;

        if (x_new < lo || x_new > hi) {
            do_bisection = true;
        } else {
            /// 检查新点的残差是否增大
            lmmc_real_t f_new = evaluate_at(expr, var, x_new);
            if (std::isnan(f_new) || std::fabs(f_new) > residual) {
                do_bisection = true;
            } else {
                /// Newton 步有效,更新状态
                /// 缩小区间:根据 fx 的符号更新端点
                if (!std::isnan(f_lo) && fx * f_lo < 0.0) {
                    hi = x;
                    f_hi = fx;
                } else if (!std::isnan(f_hi) && fx * f_hi < 0.0) {
                    lo = x;
                    f_lo = fx;
                }
                x = x_new;
                continue;
            }
        }

        if (do_bisection) {
            /// 丢弃 Newton 迭代,执行一步二分法
            lmmc_real_t mid = (lo + hi) * 0.5;
            lmmc_real_t f_mid = evaluate_at(expr, var, mid);

            if (!std::isnan(f_mid)) {
                if (!std::isnan(f_lo) && f_lo * f_mid < 0.0) {
                    hi = mid;
                    f_hi = f_mid;
                } else {
                    lo = mid;
                    f_lo = f_mid;
                }
            }
            x = (lo + hi) * 0.5;
        }

        /// 区间宽度收敛
        if ((hi - lo) < opts.tolerance) {
            lmmc_real_t final_x = (lo + hi) * 0.5;
            lmmc_real_t final_fx = evaluate_at(expr, var, final_x);
            if (!std::isnan(final_fx) && std::fabs(final_fx) < opts.tolerance) {
                return NumericRoot{final_x, std::fabs(final_fx), i};
            }
            break;
        }
    }

    /// 达到最大迭代次数:返回最佳结果(仅当残差 < tolerance)
    if (best_residual < opts.tolerance) {
        return NumericRoot{best_x, best_residual, best_iter};
    }
    return std::nullopt;
}

std::optional<NumericRoot> refine_root(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var,
    const IsolatedInterval& interval,
    const SolveOptions& opts)
{
    return refine_root_with_context(
        expr, derivative, var, interval, opts, nullptr, nullptr, nullptr);
}


/**
 * @internal
 * @brief 将去重排序后的数值根转换为 NumberNode 符号表达式向量.
 *
 * 对每个根值调用 SymbolicExpr::number(static_cast<double>(value)) 构造
 * NumberNode 表达式.输入为空时返回空向量.
 *
 * @param[in] root_values 去重,升序排列的根值列表
 * @return 对应的 NumberNode 表达式向量;输入为空时返回空向量
 */
static std::vector<std::shared_ptr<SymbolicExpr>> assemble_results(
    const std::vector<lmmc_real_t>& root_values)
{
    if (root_values.empty()) return {};

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    results.reserve(root_values.size());

    for (const auto& value : root_values) {
        results.push_back(SymbolicExpr::number(static_cast<double>(value)));
    }

    return results;
}


/**
 * @internal
 * @brief 对单个因子执行数值路径:根隔离 + 精化.
 *
 * 计算因子的导数,调用 isolate_roots 隔离根,再逐一调用 refine_root 精化.
 *
 * @param[in] factor   待求解因子表达式
 * @param[in] var      变量名
 * @param[in] interval 搜索区间
 * @param[in] opts     求解选项
 * @return 精化后的数值根列表
 */
static Result<std::vector<NumericRoot>> numerical_path(
    const std::shared_ptr<SymbolicExpr>& factor,
    const std::string& var,
    const SearchInterval& interval,
    const SolveOptions& opts,
    ComputationContext& context,
    bool& complete)
{
    std::vector<NumericRoot> roots;
    std::shared_ptr<SymbolicExpr> derivative;
    try {
        if (factor && lamina::detail::node(factor)) {
            DifferentiationVisitor visitor(var);
            lamina::detail::node(factor)->accept(visitor);
            auto derivative_node = visitor.get_result();
            if (derivative_node) {
                derivative =
                    lamina::detail::make_expression_ptr(derivative_node);
            }
        }
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::logic_error&) {
        // 缺少受支持的导数时仍可使用二分路径,但不能证明搜索完备.
        complete = false;
    }

    // 当前符号变化隔离无法证明排除偶重根,因此数值路径始终未决.
    complete = false;
    std::optional<CasError> failure;
    auto intervals = isolate_roots_with_context(
        factor, derivative, var, interval, opts,
        &context, &failure, &complete);
    if (failure) {
        return Result<std::vector<NumericRoot>>::failure(*failure);
    }

    for (const auto& isolated : intervals) {
        auto refined = refine_root_with_context(
            factor, derivative, var, isolated, opts,
            &context, &failure, &complete);
        if (failure) {
            return Result<std::vector<NumericRoot>>::failure(*failure);
        }
        if (refined) roots.push_back(*refined);
    }
    return Result<std::vector<NumericRoot>>::success(std::move(roots));
}

MixedTranscendentalResult solve_mixed_transcendental_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts,
    ComputationContext& context)
{
    constexpr const char* operation = "solve_mixed_transcendental_checked";
    if (!expr || !lamina::detail::node(expr)) {
        return MixedTranscendentalResult::failure(
            CasErrc::InvalidArgument, "待求解表达式不能为空", operation);
    }
    if (var.empty()) {
        return MixedTranscendentalResult::failure(
            CasErrc::InvalidArgument, "求解变量不能为空", operation);
    }
    if (!std::isfinite(opts.tolerance) || opts.tolerance <= 0.0 ||
        opts.max_newton_iterations <= 0 || opts.max_roots == 0 ||
        opts.max_roots < -1) {
        return MixedTranscendentalResult::failure(
            CasErrc::InvalidArgument, "求解选项包含无效界限", operation);
    }
    if (opts.has_search_interval &&
        (!std::isfinite(opts.search_lo) ||
         !std::isfinite(opts.search_hi) ||
         opts.search_lo >= opts.search_hi)) {
        return MixedTranscendentalResult::failure(
            CasErrc::InvalidArgument, "搜索区间必须为有限递增区间", operation);
    }

    try {
        auto initial_step = context.consume_steps(1, operation);
        if (!initial_step) {
            return MixedTranscendentalResult::failure(initial_step.error());
        }
        const auto variables =
            free_variables(lamina::detail::node(expr));
        if (variables.find(var) == variables.end()) {
            if (!variables.empty()) {
                return MixedTranscendentalResult::success(
                    MathResult<std::vector<std::shared_ptr<SymbolicExpr>>>{
                        {}, Completeness::Inconclusive,
                        "equation depends on unresolved parameters but not the solve variable"});
            }

            const auto constant_value =
                recursive_eval(lamina::detail::node(expr));
            if (!std::isfinite(constant_value)) {
                return MixedTranscendentalResult::success(
                    MathResult<std::vector<std::shared_ptr<SymbolicExpr>>>{
                        {}, Completeness::Inconclusive,
                        "constant equation could not be evaluated in the real domain"});
            }
            if (std::abs(constant_value) <= opts.tolerance) {
                return MixedTranscendentalResult::success(
                    MathResult<std::vector<std::shared_ptr<SymbolicExpr>>>{
                        {}, Completeness::Inconclusive,
                        "equation is identically zero; every real value satisfies it"});
            }
            return MixedTranscendentalResult::success(
                MathResult<std::vector<std::shared_ptr<SymbolicExpr>>>{
                    {}, Completeness::Complete, {}});
        }

        auto interval_opt = determine_search_interval(expr, var, opts);
        if (!interval_opt) {
            return MixedTranscendentalResult::failure(
                CasErrc::InvalidArgument, "无法确定有效搜索区间", operation);
        }
        const SearchInterval interval = *interval_opt;
        if (!std::isfinite(interval.lo) || !std::isfinite(interval.hi) ||
            interval.lo >= interval.hi) {
            return MixedTranscendentalResult::failure(
                CasErrc::InvalidArgument, "搜索区间必须为有限递增区间",
                operation);
        }

        std::vector<NumericRoot> all_roots;
        bool complete = true;
        auto accept_value = [&](lmmc_real_t value, int iterations) {
            if (!std::isfinite(value) ||
                value < interval.lo || value > interval.hi) {
                complete = false;
                return;
            }
            const lmmc_real_t residual =
                std::fabs(evaluate_at(expr, var, value));
            if (!std::isfinite(residual) || residual > opts.tolerance) {
                complete = false;
                return;
            }
            all_roots.push_back({value, residual, iterations});
        };
        auto accept_symbolic = [&](const auto& roots) {
            for (const auto& root : roots) {
                if (!root) {
                    complete = false;
                    continue;
                }
                auto simplified = root->simplify();
                if (!simplified || !lamina::detail::node(simplified)) {
                    complete = false;
                    continue;
                }
                accept_value(
                    recursive_eval(lamina::detail::node(simplified)), 0);
            }
        };

        auto factors = factor_transcendental(expr, var);
        if (factors.empty()) {
            factors.push_back(expr);
            complete = false;
        }
        for (const auto& factor : factors) {
            auto factor_step = context.consume_steps(1, operation);
            if (!factor_step) {
                return MixedTranscendentalResult::failure(
                    factor_step.error());
            }
            if (!factor || !lamina::detail::node(factor)) {
                complete = false;
                continue;
            }
            if (!expression_depends_on_variable(
                    lamina::detail::node(factor), var)) {
                continue;
            }

            bool solved_as_polynomial = false;
            if (!contains_transcendental_of_var(factor, var)) {
                try {
                    auto polynomial =
                        symbolic_to_poly<SymbolicPolyCoeff>(factor, var);
                    const int degree = polynomial.degree();
                    if (degree >= 1 && degree <= 4) {
                        accept_symbolic(solve_by_factoring(polynomial, var));
                        solved_as_polynomial = true;
                    }
                } catch (const std::invalid_argument&) {
                    complete = false;
                } catch (const std::domain_error&) {
                    complete = false;
                }
            }
            if (solved_as_polynomial) continue;

            accept_symbolic(solve_transcendental(factor, var));

            auto numeric =
                numerical_path(factor, var, interval, opts, context, complete);
            if (!numeric) {
                return MixedTranscendentalResult::failure(numeric.error());
            }
            for (const auto& root : numeric.value()) {
                accept_value(root.value, root.iterations);
            }
        }

        auto values =
            deduplicate_roots(all_roots, opts.tolerance, -1);
        if (opts.max_roots > 0 &&
            static_cast<int>(values.size()) > opts.max_roots) {
            values.resize(static_cast<std::size_t>(opts.max_roots));
            complete = false;
        }
        auto results = assemble_results(values);
        return MixedTranscendentalResult::success(
            MathResult<std::vector<std::shared_ptr<SymbolicExpr>>>{
                std::move(results),
                complete ? Completeness::Complete
                         : Completeness::Inconclusive,
                complete ? std::string{}
                         : "有界数值隔离未能证明排除全部未检测根"});
    } catch (const std::bad_alloc&) {
        return MixedTranscendentalResult::failure(
            CasErrc::ResourceLimit, "混合求解分配失败", operation);
    } catch (const std::exception& error) {
        return MixedTranscendentalResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

MixedTranscendentalResult solve_mixed_transcendental_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts)
{
    ComputationContext context;
    return solve_mixed_transcendental_checked(expr, var, opts, context);
}

} // namespace lamina
