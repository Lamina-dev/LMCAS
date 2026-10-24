#pragma once

#include "lmcas_export.hpp"
#include "result.hpp"

#include <cstddef>
#include <string>

namespace LMCAS {

enum class EqvProfile {
    Core,
    TrigBasic,
    ExpLogBasic
};

/** @brief 各限额均为非零值；深度与规模在重写运行中检查，超限返回 ResourceLimit。 */
struct EqvBudget {
    /** @brief 一次证明及 profile 切换共享的重写遍历、候选追加与阶段步数。
     * 形状查询和缓存清空的实际工作计入 ResourceLimits::max_steps。
     */
    std::size_t max_rewrite_steps = 256;
    /** @brief 同时活动的重写及嵌套 visitor 层数，根为第 1 层。 */
    std::size_t max_rewrite_depth = 64;
    /**
     * @brief 候选树及未封装孩子列表的规模上限取 context.max_ast_nodes
     * 与原始两输入节点总数乘此因子的较小值。
     * 每次子树引用分别计数；整个证明及 profile 切换共用原始基数。
     * 乘法饱和为 SIZE_MAX；真实输入节点和溢出返回 ResourceLimit。
     */
    std::size_t max_node_growth_factor = 4;
};

struct EqvOptions {
    EqvProfile profile = EqvProfile::Core;
    EqvBudget budget = {};
};

LMCAS_API Result<EqvProfile> eqv_profile_from_name(const std::string& name);
LMCAS_API Result<void> set_eqv_profile(EqvOptions& options,
                                       const std::string& name);
LMCAS_API Result<void> set_eqv_budget(EqvOptions& options,
                                      std::size_t steps,
                                      std::size_t depth,
                                      std::size_t growth);

}
