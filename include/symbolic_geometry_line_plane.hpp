#pragma once
#include "symbolic_geometry_vector.hpp"

namespace LMCAS {

/** @brief 由一点和方向向量定义的符号直线。 */
struct LineSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> point;      /**< 直线上一点。 */
    std::vector<std::shared_ptr<SymbolicExpr>> direction;  /**< 方向向量。 */
};

using LineSymbolicResult = Result<LineSymbolic>;

/** @brief 由法向量和常数 d 定义的符号平面 ax + by + cz = d。 */
struct PlaneSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> normal;  /**< 法向量。 */
    std::shared_ptr<SymbolicExpr> d;                    /**< 常数项。 */
};

using PlaneSymbolicResult = Result<PlaneSymbolic>;

/**
 * @brief 由一般方程系数构造平面 ax + by + cz = d。
 * @param a x 方向法向量分量
 * @param b y 方向法向量分量
 * @param c z 方向法向量分量
 * @param d 常数项
 * @return 构造的 PlaneSymbolic 对象
 */
inline PlaneSymbolic plane_general(
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c,
    std::shared_ptr<SymbolicExpr> d
) {
    return PlaneSymbolic{{a, b, c}, d};
}

/**
 * @brief 计算直线与平面的交点。
 *
 * 有限向量的点积或交点坐标更新存在溢出、下溢风险时，按最大分量
 * 缩放相关向量和坐标，并同步缩放平面常数及直线参数化。
 * 点积按维数预留累加余量；非零尺度变换保持交点不变。
 * @param line 直线
 * @param plane 平面
 * @return 交点坐标向量
 */

LMCAS_API VectorExprListResult line_plane_intersection_checked(
    const LineSymbolic& line,
    const PlaneSymbolic& plane,
    ComputationContext& context
);

LMCAS_API VectorExprListResult line_plane_intersection_checked(
    const LineSymbolic& line,
    const PlaneSymbolic& plane
);

/**
 * @brief 计算点到平面的距离。
 *
 * 有限法向量的平方或点积部分和存在溢出、下溢风险时，同尺度缩放
 * 法向量与常数项，保持平面和距离不变。点积按维数预留累加余量，
 * 使最终可表示的抵消结果免受部分和溢出影响。
 * @param point 空间点坐标
 * @param plane 平面
 * @return 距离的符号表达式
 */

LMCAS_API ExpressionResult point_plane_distance_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane,
    ComputationContext& context
);

LMCAS_API ExpressionResult point_plane_distance_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane
);

/**
 * @brief 计算两条异面直线之间的距离。
 *
 * 有限方向向量分别按自身最大分量归一化后计算叉积。
 * 有限点坐标差可能溢出时，先按公共最大分量缩放两点，
 * 求出距离比值后恢复点尺度；方向及临时点尺度保持距离不变。
 * @param l1 第一条直线
 * @param l2 第二条直线
 * @return 距离的符号表达式
 */

LMCAS_API ExpressionResult skew_lines_distance_checked(
    const LineSymbolic& l1,
    const LineSymbolic& l2,
    ComputationContext& context
);

LMCAS_API ExpressionResult skew_lines_distance_checked(
    const LineSymbolic& l1,
    const LineSymbolic& l2
);

/**
 * @brief 由两点构造直线。
 *
 * 有限点坐标差可能溢出时，先按两点的公共最大分量缩放再求差；
 * 公共方向尺度保持直线不变。
 * @param p1 第一个点
 * @param p2 第二个点
 * @return 直线（point = p1, direction = p2 - p1）
 */

LMCAS_API LineSymbolicResult line_from_two_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    ComputationContext& context
);

LMCAS_API LineSymbolicResult line_from_two_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2
);

/**
 * @brief 由三点构造平面（法向量 = (p2-p1)×(p3-p1)）。
 *
 * 有限点坐标差可能溢出时，按三点的公共尺度构造边；
 * 叉积或其范数平方仍有溢出、下溢风险时，分别缩放两条边。
 * 这些尺度变换保持平面不变。
 * @param p1 第一个点
 * @param p2 第二个点
 * @param p3 第三个点
 * @return 平面（法向量 + 常数 d = n·p1）
 */
LMCAS_API PlaneSymbolicResult plane_from_three_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3,
    ComputationContext& context
);

LMCAS_API PlaneSymbolicResult plane_from_three_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3
);

/**
 * @brief 计算两平面之间的二面角 arccos(|n₁·n₂| / (|n₁||n₂|))。
 *
 * 有限法向量的非零验证和角度计算均避开直接分量平方和，
 * 支持方向可表示但平方会溢出的法向量。
 * @param p1 第一个平面
 * @param p2 第二个平面
 * @return 二面角表达式（弧度）
 */

LMCAS_API ExpressionResult dihedral_angle_checked(
    const PlaneSymbolic& p1,
    const PlaneSymbolic& p2,
    ComputationContext& context
);

LMCAS_API ExpressionResult dihedral_angle_checked(
    const PlaneSymbolic& p1,
    const PlaneSymbolic& p2
);

}
