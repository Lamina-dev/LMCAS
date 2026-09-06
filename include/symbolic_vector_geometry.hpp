/**
 * @file symbolic_vector_geometry.hpp
 * @brief 向量几何：点积、叉积、夹角、直线与平面交点、异面直线距离。
 */
#pragma once
#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <memory>
#include <string>
#include <vector>

namespace LMCAS {

using VectorExprListResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;
using VectorAngleResult = Result<double>;
using VectorStringResult = Result<std::string>;

/**
 * @brief 计算两个符号向量的点积
 * @param a 向量 a 的各分量
 * @param b 向量 b 的各分量
 * @return 点积结果的符号表达式
 */
LMCAS_API ExpressionResult vector_dot_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
);

LMCAS_API ExpressionResult vector_dot_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

LMCAS_API std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/**
 * @brief 计算两个三维符号向量的叉积
 * @param a 向量 a 的各分量
 * @param b 向量 b 的各分量
 * @return 叉积结果向量的各分量
 */
LMCAS_API VectorExprListResult vector_cross_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
);

LMCAS_API VectorExprListResult vector_cross_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> vector_cross(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/**
 * @brief 采用独立向量尺度归一化计算两个有限数值向量的夹角（弧度）.
 *
 * 缩放范数和补偿点积避免分量平方、乘积的可避免溢出及相消误差。
 * 零向量返回 `CasErrc::DomainError`;非数值或非有限分量返回数值错误。
 * @param a 向量 a 的各分量
 * @param b 向量 b 的各分量
 * @return `[0,pi]` 内的有限夹角
 */
LMCAS_API VectorAngleResult vector_angle_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
);

LMCAS_API VectorAngleResult vector_angle_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);


/** @brief 符号直线，由一点和方向向量定义 */
struct LineSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> point;      ///< 直线上一点
    std::vector<std::shared_ptr<SymbolicExpr>> direction;  ///< 方向向量
};

using LineSymbolicResult = Result<LineSymbolic>;

/** @brief 符号平面，由法向量和常数 d 定义（ax + by + cz = d） */
struct PlaneSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> normal;  ///< 法向量
    std::shared_ptr<SymbolicExpr> d;                    ///< 常数项
};

using PlaneSymbolicResult = Result<PlaneSymbolic>;

/**
 * @brief 由一般方程系数构造平面 ax + by + cz = d
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
 * 当有限法向量与方向向量、直线点的点积或最终坐标更新可能溢出、
 * 下溢时，相关向量和空间坐标按最大分量缩放，并同步缩放平面常数
 * 及直线参数化。点积路径按维数保留累加余量；这些非零尺度变换不
 * 改变几何交点。
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
 * @brief 计算点到平面的距离
 *
 * 当有限平面法向量的直接平方或点积部分和可能溢出、下溢时，法向量
 * 和常数项按同一尺度缩放；该缩放不改变平面或距离。点积风险路径还
 * 按维数留出累加余量，避免最终可表示的抵消结果在部分和处溢出。
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
 * @brief 计算两条异面直线之间的距离
 *
 * 两条有限数值方向向量分别按自身最大分量归一化后再计算叉积。
 * 若有限数值点的直接坐标差会溢出，则两点先按公共最大分量缩放，
 * 距离比值求出后再恢复点尺度；方向尺度和临时点尺度均不改变距离。
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
 * @brief 由两点构造直线（点 + 方向向量）。
 *
 * 若有限点坐标的直接差会溢出，则先按两点的公共最大分量缩放后求差；
 * 方向向量的公共尺度不改变所表示的直线。
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
 * 对有限数值点，先在必要时用三点的公共尺度形成边，避免坐标差溢出；
 * 若后续叉积或其范数平方仍有溢出/下溢风险，再分别缩放两条边。
 * 这些公共尺度均不改变所表示的平面。
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
 * 有限数值法向量的非零验证与角度计算均避免直接构造分量平方和，
 * 因而接受方向可表示但平方会溢出的法向量。
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

/** @brief 隐式曲面 F(x,y,z) = 0 */
struct SurfaceSymbolic {
    std::shared_ptr<SymbolicExpr> F;                 ///< 曲面方程左端（= 0）
    std::vector<std::string> vars;                   ///< 坐标变量名 {x,y,z}
};

/**
 * @brief 对二次曲面进行分类。
 *
 * 二次型通过对称特征分解旋转到主轴坐标，再完成平方并依据中心常数和
 * 惯性确认实数轨迹。该过程支持混合二次项；零空间中的线性分量使用
 * 逐主轴点积误差界判定，避免其他轴的大系数掩盖抛物方向。秩一非中心
 * 二次型识别为抛物柱面。落在特征分解反向误差界内的非零特征值返回
 * `Inconclusive`，不会静默降秩；空集、单点、直线及平面对同样不会
 * 冒充非退化曲面。
 * @param surf 隐式曲面
 * @return 分类字符串："sphere"、"ellipsoid"、"paraboloid"、"hyperboloid"、
 *         "cone"、"cylinder"、"unknown"
 */

LMCAS_API VectorStringResult classify_quadric_checked(
    const SurfaceSymbolic& surf,
    ComputationContext& context
);

LMCAS_API VectorStringResult classify_quadric_checked(
    const SurfaceSymbolic& surf
);

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
