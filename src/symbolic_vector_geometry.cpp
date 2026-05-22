#include "../include/symbolic_vector_geometry.hpp"
#include "../include/symbolic.hpp"
#include <cmath>
#include <vector>

namespace lamina {

std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {
    if (a.size() != b.size()) return nullptr;
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
    if (a.size() != 3 || b.size() != 3) return {};
    auto x = SymbolicExpr::add(SymbolicExpr::multiply(a[1], b[2]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[2], b[1])));
    auto y = SymbolicExpr::add(SymbolicExpr::multiply(a[2], b[0]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[0], b[2])));
    auto z = SymbolicExpr::add(SymbolicExpr::multiply(a[0], b[1]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[1], b[0])));
    return {x, y, z};
}

double vector_angle(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {

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
    if (!numeric) return NAN;
    double angle = std::acos(dot / (std::sqrt(norm_a) * std::sqrt(norm_b)));
    return angle;
}

std::vector<std::shared_ptr<SymbolicExpr>> line_plane_intersection(
    const LineSymbolic& line,
    const PlaneSymbolic& plane
) {

    auto n_dot_a = vector_dot(plane.normal, line.point);
    auto n_dot_b = vector_dot(plane.normal, line.direction);
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

}
