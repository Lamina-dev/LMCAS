#include "../include/symbolic_vector_geometry.hpp"
#include "../include/symbolic.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lamina {

std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("vector_dot: operand dimensions do not match");
    }
    std::vector<std::shared_ptr<SymbolicNode>> sum_terms;
    for (size_t i = 0; i < a.size(); ++i) {
        sum_terms.push_back(SymbolicExpr::multiply(a[i], b[i])->root);
    }
    return std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(sum_terms));
}

std::vector<std::shared_ptr<SymbolicExpr>> vector_cross(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {
    if (a.size() != 3 || b.size() != 3) {
        throw std::invalid_argument("vector_cross: requires 3-dimensional vectors");
    }
    auto x = SymbolicExpr::add(SymbolicExpr::multiply(a[1], b[2]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[2], b[1])));
    auto y = SymbolicExpr::add(SymbolicExpr::multiply(a[2], b[0]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[0], b[2])));
    auto z = SymbolicExpr::add(SymbolicExpr::multiply(a[0], b[1]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[1], b[0])));
    return {x, y, z};
}

double vector_angle(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {

    if (a.size() != b.size() || a.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double norm_a = 0, norm_b = 0, dot = 0;
    bool numeric = true;
    for (size_t i = 0; i < a.size(); ++i) {
        if (!a[i]->is_number() || !b[i]->is_number()) {
            numeric = false;
            break;
        }
        auto na_var = a[i]->get_number_value();
        auto nb_var = b[i]->get_number_value();

        auto to_double = [](const auto& v) -> double {
            using T = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<T, int>) return static_cast<double>(v);
            else if constexpr (std::is_same_v<T, BigInt>) return v.to_double();
            else if constexpr (std::is_same_v<T, Rational>) return v.to_double();
            else return 0.0;
        };
        double na = std::visit(to_double, na_var);
        double nb = std::visit(to_double, nb_var);
        norm_a += na * na;
        norm_b += nb * nb;
        dot += na * nb;
    }
    if (!numeric) return std::numeric_limits<double>::quiet_NaN();
    if (norm_a == 0.0 || norm_b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double cosv = dot / (std::sqrt(norm_a) * std::sqrt(norm_b));
    // Clamp into [-1, 1] to avoid std::acos domain errors caused by rounding.
    if (cosv > 1.0) cosv = 1.0;
    else if (cosv < -1.0) cosv = -1.0;
    return std::acos(cosv);
}

std::vector<std::shared_ptr<SymbolicExpr>> line_plane_intersection(
    const LineSymbolic& line,
    const PlaneSymbolic& plane
) {

    auto n_dot_a = vector_dot(plane.normal, line.point);
    auto n_dot_b = vector_dot(plane.normal, line.direction);
    // If n . direction == 0 the line is parallel to the plane: either lies in
    // the plane (infinite intersections) or has no intersection. In both cases
    // we cannot return a single point, so signal "no unique intersection" via
    // an empty vector instead of dividing by zero.
    auto n_dot_b_simpl = n_dot_b ? n_dot_b->simplify() : nullptr;
    if (!n_dot_b_simpl || n_dot_b_simpl->is_zero()) {
        return {};
    }
    auto t = SymbolicExpr::divide(SymbolicExpr::add(plane.d, SymbolicExpr::multiply(SymbolicExpr::number(-1), n_dot_a)), n_dot_b);
    std::vector<std::shared_ptr<SymbolicExpr>> intersection;
    for (size_t i = 0; i < line.point.size(); ++i) {
        intersection.push_back(SymbolicExpr::add(line.point[i], SymbolicExpr::multiply(t, line.direction[i])));
    }
    return intersection;
}

std::shared_ptr<SymbolicExpr> point_plane_distance(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane
) {

    auto n_dot_r = vector_dot(plane.normal, point);
    auto diff = SymbolicExpr::add(n_dot_r, SymbolicExpr::multiply(SymbolicExpr::number(-1), plane.d));
    auto norm_n = SymbolicExpr::sqrt(vector_dot(plane.normal, plane.normal));

    auto abs_diff = std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Abs, std::vector<std::shared_ptr<SymbolicNode>>{diff->root}));
    return SymbolicExpr::divide(abs_diff, norm_n);
}

std::shared_ptr<SymbolicExpr> skew_lines_distance(
    const LineSymbolic& l1,
    const LineSymbolic& l2
) {

    std::vector<std::shared_ptr<SymbolicExpr>> a2_minus_a1;
    for (size_t i = 0; i < l1.point.size(); ++i) {
        a2_minus_a1.push_back(SymbolicExpr::add(l2.point[i], SymbolicExpr::multiply(SymbolicExpr::number(-1), l1.point[i])));
    }
    auto cross = vector_cross(l1.direction, l2.direction);
    auto cross_norm = SymbolicExpr::sqrt(vector_dot(cross, cross));
    auto numerator = vector_dot(a2_minus_a1, cross);

    auto abs_num = std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Abs, std::vector<std::shared_ptr<SymbolicNode>>{numerator->root}));
    return SymbolicExpr::divide(abs_num, cross_norm);
}

LineSymbolic line_from_two_points(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2) {
    LineSymbolic line;
    line.point = p1;
    for (size_t i = 0; i < p1.size() && i < p2.size(); ++i) {
        line.direction.push_back(
            SymbolicExpr::add(p2[i], SymbolicExpr::multiply(SymbolicExpr::number(-1), p1[i]))->simplify());
    }
    return line;
}

PlaneSymbolic plane_from_three_points(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3) {
    std::vector<std::shared_ptr<SymbolicExpr>> v1, v2;
    for (size_t i = 0; i < 3; ++i) {
        v1.push_back(SymbolicExpr::add(p2[i], SymbolicExpr::multiply(SymbolicExpr::number(-1), p1[i])));
        v2.push_back(SymbolicExpr::add(p3[i], SymbolicExpr::multiply(SymbolicExpr::number(-1), p1[i])));
    }
    auto n = vector_cross(v1, v2);
    PlaneSymbolic plane;
    plane.normal = {n[0]->simplify(), n[1]->simplify(), n[2]->simplify()};
    // d = n · p1
    plane.d = vector_dot(plane.normal, p1)->simplify();
    return plane;
}

std::shared_ptr<SymbolicExpr> dihedral_angle(
    const PlaneSymbolic& p1, const PlaneSymbolic& p2) {
    auto dot = vector_dot(p1.normal, p2.normal);
    auto n1 = SymbolicExpr::sqrt(vector_dot(p1.normal, p1.normal));
    auto n2 = SymbolicExpr::sqrt(vector_dot(p2.normal, p2.normal));
    auto abs_dot = std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(
        FunctionNode::FuncType::Abs, std::vector<std::shared_ptr<SymbolicNode>>{dot->root}));
    auto cos_theta = SymbolicExpr::divide(abs_dot, SymbolicExpr::multiply(n1, n2));
    auto arccos = std::make_shared<FunctionNode>(FunctionNode::FuncType::ArcCos,
        std::vector<std::shared_ptr<SymbolicNode>>{cos_theta->root});
    return std::make_shared<SymbolicExpr>(arccos)->simplify();
}

std::string classify_quadric(const SurfaceSymbolic& surf) {
    if (!surf.F || surf.vars.size() < 3) return "unknown";
    auto& vars = surf.vars;
    auto F = surf.F->expand();
    if (!F) F = surf.F;

    // 提取二次项系数 a*x^2、线性项、常数。
    auto coeff_of_sq = [&](const std::string& v) -> double {
        auto d2 = F->differentiate(v)->differentiate(v)->simplify();
        if (d2->is_number()) return d2->to_numeric() / 2.0;
        return 0.0;
    };
    auto coeff_of_lin = [&](const std::string& v) -> double {
        // 在所有变量=0 处的一阶偏导
        auto d = F->differentiate(v);
        for (auto& w : vars) d = d->substitute(w, SymbolicExpr::number(0));
        d = d->simplify();
        if (d->is_number()) return d->to_numeric();
        return 0.0;
    };
    double a = coeff_of_sq(vars[0]);
    double b = coeff_of_sq(vars[1]);
    double c = coeff_of_sq(vars[2]);
    double lz = coeff_of_lin(vars[2]);

    int n_sq = (a != 0) + (b != 0) + (c != 0);
    int n_pos = (a > 0) + (b > 0) + (c > 0);
    int n_neg = (a < 0) + (b < 0) + (c < 0);

    if (n_sq == 2 && lz != 0 && c == 0) return "paraboloid";
    if (n_sq == 2) return "cylinder";
    if (n_sq == 3) {
        if (n_neg == 0) {
            // 全正：球或椭球。系数相等→球
            if (a == b && b == c) return "sphere";
            return "ellipsoid";
        }
        if (n_pos == 0) return "ellipsoid"; // 全负，移项后等价
        // 混合符号：锥面（常数项0）或双曲面
        for (auto& w : vars) F = F->substitute(w, SymbolicExpr::number(0));
        auto cst = F->simplify();
        double cval = cst->is_number() ? cst->to_numeric() : 1.0;
        if (std::abs(cval) < 1e-9) return "cone";
        return "hyperboloid";
    }
    return "unknown";
}

std::vector<std::shared_ptr<SymbolicExpr>> surface_normal(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point) {
    std::vector<std::shared_ptr<SymbolicExpr>> grad;
    for (size_t i = 0; i < surf.vars.size(); ++i) {
        auto d = surf.F->differentiate(surf.vars[i]);
        for (size_t j = 0; j < surf.vars.size() && j < point.size(); ++j) {
            d = d->substitute(surf.vars[j], point[j]);
        }
        grad.push_back(d->simplify());
    }
    // 归一化
    auto norm_sq = SymbolicExpr::number(0);
    for (auto& g : grad) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(g, g));
    auto norm = SymbolicExpr::sqrt(norm_sq);
    std::vector<std::shared_ptr<SymbolicExpr>> result;
    for (auto& g : grad) result.push_back(SymbolicExpr::divide(g, norm)->simplify());
    return result;
}

PlaneSymbolic tangent_plane(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point) {
    PlaneSymbolic plane;
    for (size_t i = 0; i < surf.vars.size(); ++i) {
        auto d = surf.F->differentiate(surf.vars[i]);
        for (size_t j = 0; j < surf.vars.size() && j < point.size(); ++j) {
            d = d->substitute(surf.vars[j], point[j]);
        }
        plane.normal.push_back(d->simplify());
    }
    // d = n · point
    plane.d = vector_dot(plane.normal, point)->simplify();
    return plane;
}

}
