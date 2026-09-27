/** @file integrator.hpp */ #pragma once
#include "integration_strategies.hpp"
#include "computation_context.hpp"
#include <memory>
#include <vector>

namespace LMCAS {

/** @brief 符号积分器，协调各策略完成积分运算 */
class LMCAS_API Integrator {
public:
    Integrator();

    /**
     * @brief 计算不定积分。
     *
     * @param expr 被积表达式
     * @param var_name 积分变量名
     * @return 积分结果表达式；原函数尚未证明时可成功返回未求值的 IntegralNode，闭式是否存在仍待定。
     * @note 已求值原函数遵循原始公共定义域约束，可去间断点的逐点定义另行处理。
     * 传入上下文提供假设、资源及取消限制；计算错误返回失败。
     */
    Result<SymbolicExpr> integrate(
        const SymbolicExpr& expr, const std::string& var_name);
    Result<SymbolicExpr> integrate_checked(
        const SymbolicExpr& expr,
        const std::string& var_name,
        ComputationContext& context);
    /**
     * @brief 计算普通实广义定积分。
     *
     * @param expr 被积表达式
     * @param var_name 积分变量名
     * @param lower 积分下限
     * @param upper 积分上限
     * @return 经证明的有限表达式；各已证明发散部分同号时返回带符号 Infinity。
     * 异号无穷或已证明积分不存在时返回 DomainError（integrate.definite.convergence）。
     * @note 采用普通广义积分语义，Cauchy 主值另论。约消前划分原始定义域；
     * 可去间断点仅以单侧极限衔接，保留函数在该点未定义的状态。
     * 定义域、原函数或收敛性尚未证明时，以单个 IntegralNode 保留完整原始被积式及上下限。
     * checked 调用与边界极限证明共享传入上下文，并传播取消、资源耗尽及不变量失败。
     */
    Result<SymbolicExpr> integrate_def(
        const SymbolicExpr& expr, const std::string& var_name,
        const SymbolicExpr& lower, const SymbolicExpr& upper);
    Result<SymbolicExpr> integrate_def_checked(
        const SymbolicExpr& expr,
        const std::string& var_name,
        const SymbolicExpr& lower,
        const SymbolicExpr& upper,
        ComputationContext& context);

    /**
     * @brief 添加积分策略
     * @param strategy 策略对象
     * @param position 插入位置，-1 表示追加到末尾
     */
    Result<void> add_strategy(
        std::unique_ptr<IntegrationStrategy> strategy, int position = -1);

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
    Result<std::shared_ptr<SymbolicExpr>> integrate_recursive(
        const SymbolicExpr& expr, const std::string& var,
        ComputationContext& context, int depth = 0);

    /**
     * @brief 判断表达式是否依赖指定变量
     * @param expr 表达式
     * @param var 变量名
     * @return 依赖返回 true
     */
    static bool depends_on(const SymbolicExpr& expr, const std::string& var);

private:
    IntegrationTable table_;
    std::vector<std::unique_ptr<IntegrationStrategy>> strategies_;

    struct CycleState {
        std::vector<SymbolicExpr> history;
    };
    CycleState cycle_state_;
    std::size_t query_depth_ = 0;

    Result<std::shared_ptr<SymbolicExpr>> apply_linearity(
        const SymbolicExpr& expr, const std::string& var,
        ComputationContext& context, int depth);

    Result<std::shared_ptr<SymbolicExpr>> dispatch_strategies(
        const SymbolicExpr& expr, const std::string& var,
        ComputationContext& context, int depth);

    static std::shared_ptr<SymbolicExpr> make_integral_node(
        const SymbolicExpr& expr, const std::string& var);

    std::shared_ptr<SymbolicExpr> check_cycle(
        const SymbolicExpr& expr, const std::string& var);
    void resolve_cycle(std::shared_ptr<SymbolicExpr>& result, size_t cycle_idx);
};

} // namespace LMCAS
