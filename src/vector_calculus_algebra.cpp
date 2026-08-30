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

VectorCalculusExprResult dot_product(const VectorField& a, const VectorField& b)
{
    if (a.size() != b.size()) {
        return VectorCalculusExprResult::failure(
            CasErrc::DimensionMismatch,
            "vectors must have the same dimension",
            "dot_product");
    }
    std::shared_ptr<SymbolicExpr> sum = SymbolicExpr::number(0);
    for (size_t i = 0; i < a.size(); ++i) {
        if (!a[i] || !b[i]) continue;
        sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(a[i], b[i]));
    }
    return VectorCalculusExprResult::success(sum->simplify());
}

VectorCalculusFieldResult cross_product(const VectorField& a, const VectorField& b)
{
    if (a.size() != 3 || b.size() != 3) {
        return VectorCalculusFieldResult::failure(
            CasErrc::DimensionMismatch,
            "vectors must be three-dimensional",
            "cross_product");
    }
    auto sub = [](const std::shared_ptr<SymbolicExpr>& p,
                  const std::shared_ptr<SymbolicExpr>& q) {
        return SymbolicExpr::add(p, SymbolicExpr::multiply(SymbolicExpr::number(-1), q));
    };
    VectorField result(3);
    result[0] = sub(SymbolicExpr::multiply(a[1], b[2]), SymbolicExpr::multiply(a[2], b[1]))->simplify();
    result[1] = sub(SymbolicExpr::multiply(a[2], b[0]), SymbolicExpr::multiply(a[0], b[2]))->simplify();
    result[2] = sub(SymbolicExpr::multiply(a[0], b[1]), SymbolicExpr::multiply(a[1], b[0]))->simplify();
    return VectorCalculusFieldResult::success(std::move(result));
}

VectorCalculusFieldResult vector_project(const VectorField& a, const VectorField& b)
{
    if (a.size() != b.size()) {
        return VectorCalculusFieldResult::failure(
            CasErrc::DimensionMismatch,
            "vectors must have the same dimension",
            "vector_project");
    }
    auto bb_result = dot_product(b, b);
    if (!bb_result) return VectorCalculusFieldResult::failure(bb_result.error());
    auto bb = bb_result.value();
    if (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero()) {
        return VectorCalculusFieldResult::success(
            VectorField(a.size(), SymbolicExpr::number(0)));
    }
    auto ab_result = dot_product(a, b);
    if (!ab_result) return VectorCalculusFieldResult::failure(ab_result.error());
    auto ab = ab_result.value();
    auto coeff = SymbolicExpr::divide(ab, bb);
    VectorField result;
    result.reserve(b.size());
    for (const auto& comp : b) {
        result.push_back(SymbolicExpr::multiply(coeff, comp)->simplify());
    }
    return VectorCalculusFieldResult::success(std::move(result));
}

VectorCalculusExprResult scalar_project(const VectorField& a, const VectorField& b)
{
    auto bb_result = dot_product(b, b);
    if (!bb_result) return bb_result;
    auto bb = bb_result.value();
    if (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero()) {
        return VectorCalculusExprResult::success(nullptr);
    }
    auto ab_result = dot_product(a, b);
    if (!ab_result) return ab_result;
    return VectorCalculusExprResult::success(
        SymbolicExpr::divide(ab_result.value(), SymbolicExpr::sqrt(bb))->simplify());
}

VectorCalculusExprResult vector_angle_symbolic(const VectorField& a, const VectorField& b)
{
    auto aa_result = dot_product(a, a);
    if (!aa_result) return aa_result;
    auto bb_result = dot_product(b, b);
    if (!bb_result) return bb_result;
    auto aa = aa_result.value();
    auto bb = bb_result.value();
    if ((lamina::detail::node(aa) && lamina::detail::node(aa)->is_zero()) ||
        (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero())) {
        return VectorCalculusExprResult::success(nullptr);
    }
    auto ab_result = dot_product(a, b);
    if (!ab_result) return ab_result;
    auto denom = SymbolicExpr::multiply(SymbolicExpr::sqrt(aa), SymbolicExpr::sqrt(bb));
    auto cos_theta = SymbolicExpr::divide(ab_result.value(), denom);
    auto arccos_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::ArcCos,
        std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(cos_theta)});
    return VectorCalculusExprResult::success(
        lamina::detail::make_expression_ptr(arccos_node)->simplify());
}

VectorCalculusExprResult mixed_product(
    const VectorField& a, const VectorField& b, const VectorField& c)
{
    auto cross = cross_product(b, c);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());
    return dot_product(a, cross.value());
}


} // namespace lamina
