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
std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/**
 * @brief 计算两个三维符号向量的叉积
 * @param a 向量 a 的各分量
 * @param b 向量 b 的各分量
 * @return 叉积结果向量的各分量
 */
std::vector<std::shared_ptr<SymbolicExpr>> vector_cross(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/**
 * @brief 计算两个向量的夹角（弧度）
 * @param a 向量 a 的各分量
 * @param b 向量 b 的各分量
 * @return 夹角的数值（弧度）
 */
double vector_angle(
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
std::vector<std::shared_ptr<SymbolicExpr>> line_plane_intersection(
    const LineSymbolic& line,
    const PlaneSymbolic& plane
);

/**
 * @brief 计算点到平面的距离
 * @param point 空间点坐标
 * @param plane 平面
 * @return 距离的符号表达式
 */
std::shared_ptr<SymbolicExpr> point_plane_distance(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane
);

/**
 * @brief 计算两条异面直线之间的距离
 * @param l1 第一条直线
 * @param l2 第二条直线
 * @return 距离的符号表达式
 */
std::shared_ptr<SymbolicExpr> skew_lines_distance(
    const LineSymbolic& l1,
    const LineSymbolic& l2
);

}
