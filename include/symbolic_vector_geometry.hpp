/**
 * @file symbolic_vector_geometry.hpp
 * @brief 向量几何：点积、叉积、夹角、直线与平面交点、异面直线距离。
 */
#pragma once
#include "symbolic.hpp"
#include <memory>
#include <string>
#include <vector>

namespace lamina {

/**
 * @brief 计算两个符号向量的点积
 * @param a 向量 a 的各分量
 * @param b 向量 b 的各分量
 * @return 点积结果的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/**
 * @brief 计算两个三维符号向量的叉积
 * @param a 向量 a 的各分量
 * @param b 向量 b 的各分量
 * @return 叉积结果向量的各分量
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> vector_cross(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/**
 * @brief 计算两个向量的夹角（弧度）
 * @param a 向量 a 的各分量
 * @param b 向量 b 的各分量
 * @return 夹角的数值（弧度）
 */
LAMINA_API double vector_angle(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/** @brief 符号直线，由一点和方向向量定义 */
struct LineSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> point;      ///< 直线上一点
    std::vector<std::shared_ptr<SymbolicExpr>> direction;  ///< 方向向量
};

/** @brief 符号平面，由法向量和常数 d 定义（ax + by + cz = d） */
struct PlaneSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> normal;  ///< 法向量
    std::shared_ptr<SymbolicExpr> d;                    ///< 常数项
};

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
 * @brief 计算直线与平面的交点
 * @param line 直线
 * @param plane 平面
 * @return 交点坐标向量
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> line_plane_intersection(
    const LineSymbolic& line,
    const PlaneSymbolic& plane
);

/**
 * @brief 计算点到平面的距离
 * @param point 空间点坐标
 * @param plane 平面
 * @return 距离的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> point_plane_distance(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane
);

/**
 * @brief 计算两条异面直线之间的距离
 * @param l1 第一条直线
 * @param l2 第二条直线
 * @return 距离的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> skew_lines_distance(
    const LineSymbolic& l1,
    const LineSymbolic& l2
);

/**
 * @brief 由两点构造直线（点 + 方向向量）。
 * @param p1 第一个点
 * @param p2 第二个点
 * @return 直线（point = p1, direction = p2 - p1）
 */
LAMINA_API LineSymbolic line_from_two_points(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2
);

/**
 * @brief 由三点构造平面（法向量 = (p2-p1)×(p3-p1)）。
 * @param p1 第一个点
 * @param p2 第二个点
 * @param p3 第三个点
 * @return 平面（法向量 + 常数 d = n·p1）
 */
LAMINA_API PlaneSymbolic plane_from_three_points(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3
);

/**
 * @brief 计算两平面之间的二面角 arccos(|n₁·n₂| / (|n₁||n₂|))。
 * @param p1 第一个平面
 * @param p2 第二个平面
 * @return 二面角表达式（弧度）
 */
LAMINA_API std::shared_ptr<SymbolicExpr> dihedral_angle(
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
 * @param surf 隐式曲面
 * @return 分类字符串："sphere"、"ellipsoid"、"paraboloid"、"hyperboloid"、
 *         "cone"、"cylinder"、"unknown"
 */
LAMINA_API std::string classify_quadric(const SurfaceSymbolic& surf);

/**
 * @brief 计算曲面在某点的单位法向量 ∇F/|∇F|。
 * @param surf 隐式曲面
 * @param point 曲面上的点（变量名到值的映射）
 * @return 单位法向量分量
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> surface_normal(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point
);

/**
 * @brief 计算曲面在某点的切平面。
 * @param surf 隐式曲面
 * @param point 切点
 * @return 切平面（法向量 = ∇F(point)）
 */
LAMINA_API PlaneSymbolic tangent_plane(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point
);

}
