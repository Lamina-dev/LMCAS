#pragma once
#include "symbolic.hpp"
#include <memory>
#include <string>
#include <vector>

namespace lamina {

std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

std::vector<std::shared_ptr<SymbolicExpr>> vector_cross(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

double vector_angle(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

struct LineSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> point;
    std::vector<std::shared_ptr<SymbolicExpr>> direction;
};

struct PlaneSymbolic {
    std::vector<std::shared_ptr<SymbolicExpr>> normal;
    std::shared_ptr<SymbolicExpr> d;
};

inline PlaneSymbolic plane_general(
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c,
    std::shared_ptr<SymbolicExpr> d
) {
    return PlaneSymbolic{{a, b, c}, d};
}

std::vector<std::shared_ptr<SymbolicExpr>> line_plane_intersection(
    const LineSymbolic& line,
    const PlaneSymbolic& plane
);

std::shared_ptr<SymbolicExpr> point_plane_distance(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane
);

std::shared_ptr<SymbolicExpr> skew_lines_distance(
    const LineSymbolic& l1,
    const LineSymbolic& l2
);

}
