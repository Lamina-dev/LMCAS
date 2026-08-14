#include "internal/vector_calculus_support.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "symbolic_ast.hpp"

#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace lamina {

using namespace vector_calculus_detail;

std::shared_ptr<SymbolicExpr> dot_product(const VectorField& a, const VectorField& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("dot_product: vectors must have the same dimension");
    }
    std::shared_ptr<SymbolicExpr> sum = SymbolicExpr::number(0);
    for (size_t i = 0; i < a.size(); ++i) {
        if (!a[i] || !b[i]) continue;
        sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(a[i], b[i]));
    }
    return sum->simplify();
}

VectorField cross_product(const VectorField& a, const VectorField& b)
{
    if (a.size() != 3 || b.size() != 3) {
        throw std::invalid_argument("cross_product: vectors must be 3-dimensional");
    }
    auto sub = [](const std::shared_ptr<SymbolicExpr>& p,
                  const std::shared_ptr<SymbolicExpr>& q) {
        return SymbolicExpr::add(p, SymbolicExpr::multiply(SymbolicExpr::number(-1), q));
    };
    VectorField result(3);
    result[0] = sub(SymbolicExpr::multiply(a[1], b[2]), SymbolicExpr::multiply(a[2], b[1]))->simplify();
    result[1] = sub(SymbolicExpr::multiply(a[2], b[0]), SymbolicExpr::multiply(a[0], b[2]))->simplify();
    result[2] = sub(SymbolicExpr::multiply(a[0], b[1]), SymbolicExpr::multiply(a[1], b[0]))->simplify();
    return result;
}

VectorField vector_project(const VectorField& a, const VectorField& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("vector_project: vectors must have the same dimension");
    }
    auto bb = dot_product(b, b);
    if (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero()) {
        return VectorField(a.size(), SymbolicExpr::number(0));
    }
    auto ab = dot_product(a, b);
    auto coeff = SymbolicExpr::divide(ab, bb);
    VectorField result;
    result.reserve(b.size());
    for (const auto& comp : b) {
        result.push_back(SymbolicExpr::multiply(coeff, comp)->simplify());
    }
    return result;
}

std::shared_ptr<SymbolicExpr> scalar_project(const VectorField& a, const VectorField& b)
{
    auto bb = dot_product(b, b);
    if (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero()) {
        return nullptr;
    }
    auto ab = dot_product(a, b);
    return SymbolicExpr::divide(ab, SymbolicExpr::sqrt(bb))->simplify();
}

std::shared_ptr<SymbolicExpr> vector_angle_symbolic(const VectorField& a, const VectorField& b)
{
    auto aa = dot_product(a, a);
    auto bb = dot_product(b, b);
    if ((lamina::detail::node(aa) && lamina::detail::node(aa)->is_zero()) || (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero())) {
        return nullptr;
    }
    auto ab = dot_product(a, b);
    auto denom = SymbolicExpr::multiply(SymbolicExpr::sqrt(aa), SymbolicExpr::sqrt(bb));
    auto cos_theta = SymbolicExpr::divide(ab, denom);
    auto arccos_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::ArcCos,
        std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(cos_theta)});
    return lamina::detail::make_expression_ptr(arccos_node)->simplify();
}

std::shared_ptr<SymbolicExpr> mixed_product(const VectorField& a, const VectorField& b,
    const VectorField& c)
{
    return dot_product(a, cross_product(b, c));
}


} // namespace lamina
