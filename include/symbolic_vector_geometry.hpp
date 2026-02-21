#pragma once
#include "symbolic.hpp"
#include <memory>
#include <string>
#include <vector>

namespace lamina {
// 向量点积
std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

// 向量叉积（仅支持三维）
std::vector<std::shared_ptr<SymbolicExpr>> vector_cross(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

// 求夹角（弧度）
double vector_angle(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

// 直线符号类 r = a + t b
struct LineSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> point;
    std::vector<std::shared_ptr<SymbolicExpr>> direction;
};

// 平面符号类 r ⋅ n = d
struct PlaneSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> normal;
    std::shared_ptr<SymbolicExpr> d;
};

// 平面一般式 ax+by+cz=d
inline PlaneSymbolic plane_general(
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c,
    std::shared_ptr<SymbolicExpr> d
) {
    return PlaneSymbolic{{a, b, c}, d};
}

// 求直线与平面交点
std::vector<std::shared_ptr<SymbolicExpr>> line_plane_intersection(
    const LineSymbolic& line,
    const PlaneSymbolic& plane
);

// 点到平面距离
std::shared_ptr<SymbolicExpr> point_plane_distance(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane
);

// 异面直线距离
std::shared_ptr<SymbolicExpr> skew_lines_distance(
    const LineSymbolic& l1,
    const LineSymbolic& l2
);

}
