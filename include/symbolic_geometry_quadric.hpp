#pragma once
#include "symbolic_geometry_vector.hpp"

namespace LMCAS {

/** @brief 隐式曲面 F(x,y,z) = 0。 */
struct SurfaceSymbolic {
    std::shared_ptr<SymbolicExpr> F;                 /**< 曲面方程左端，等于零。 */
    std::vector<std::string> vars;                   /**< 坐标变量名 {x,y,z}。 */
};

/**
 * @brief 对二次曲面进行分类。
 * 对称特征分解将含混合项的二次型旋转到主轴，再完成平方，
 * 依据中心常数和惯性判定实数轨迹。零空间线性分量采用逐主轴点积误差界，
 * 使抛物方向的判定独立于其他轴的大系数；秩一非中心二次型归为抛物柱面。
 * 非零特征值落入特征分解反向误差界时返回 `Inconclusive`，保持原秩；
 * 空集、单点、直线和平面对按退化情形处理。
 * @param surf 隐式曲面。
 * @return 分类字符串："sphere"、"ellipsoid"、"paraboloid"、"hyperboloid"、
 *         "cone"、"cylinder"、"unknown"。
 */

LMCAS_API VectorStringResult classify_quadric_checked(
    const SurfaceSymbolic& surf,
    ComputationContext& context
);

LMCAS_API VectorStringResult classify_quadric_checked(
    const SurfaceSymbolic& surf
);

}
