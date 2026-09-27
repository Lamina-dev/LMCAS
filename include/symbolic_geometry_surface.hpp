#pragma once
#include "symbolic_geometry_quadric.hpp"
#include "symbolic_geometry_line_plane.hpp"

namespace LMCAS {

/**
 * @brief 计算曲面在某点的单位法向量 ∇F/|∇F|。
 *
 * 数值梯度采用按最大分量缩放的二范数归一化，因此有限梯度的平方
 * 即使超出 `double` 范围，只要单位方向可表示仍能返回结果。
 * @param surf 隐式曲面
 * @param point 曲面上的点（变量名到值的映射）
 * @return 单位法向量分量
 */

LMCAS_API VectorExprListResult surface_normal_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    ComputationContext& context
);

LMCAS_API VectorExprListResult surface_normal_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point
);

/**
 * @brief 计算曲面在某点的切平面。
 *
 * 奇异点判断直接检查有限梯度分量，不构造可能溢出的平方和。有限数值
 * 梯度与切点的点积存在中间溢出风险时，等比例缩放平面系数后计算常数。
 * @param surf 隐式曲面
 * @param point 切点
 * @return 切平面（法向量与 ∇F(point) 同向）
 */

LMCAS_API PlaneSymbolicResult tangent_plane_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    ComputationContext& context
);

LMCAS_API PlaneSymbolicResult tangent_plane_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point
);

}
