#include "internal/mixed_transcendental_support.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include <new>
#include <stdexcept>

namespace LMCAS::detail {

/**
 * @internal
 * @brief 隔离并精化单个因子的数值根。
 *
 * 计算导数后用 isolate_roots 隔离根，再用 refine_root 逐根精化。
 *
 * @param[in] factor   待求解因子表达式
 * @param[in] var      变量名
 * @param[in] interval 搜索区间
 * @param[in] opts     求解选项
 * @return 精化后的数值根列表
 */
Result<std::vector<NumericRoot>> mixed_numerical_path(
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
        if (factor && LMCAS::detail::node(factor)) {
            DifferentiationVisitor visitor(var);
            LMCAS::detail::node(factor)->accept(visitor);
            auto derivative_node = visitor.get_result();
            if (derivative_node) {
                derivative =
                    LMCAS::detail::make_expression_ptr(derivative_node);
            }
        }
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::logic_error&) {
        complete = false;
    }
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
        if (refined) {
            roots.push_back(*refined);
        }
    }
    return Result<std::vector<NumericRoot>>::success(std::move(roots));
}

}
