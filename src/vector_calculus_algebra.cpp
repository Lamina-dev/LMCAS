#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_vectors.hpp"
#include "internal/symbolic_ast.hpp"


namespace LMCAS {

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

std::shared_ptr<SymbolicExpr>
vector_calculus_detail::vector_calculus_cross_component(
    const VectorField& left, const VectorField& right, std::size_t index) {
    static constexpr std::size_t first[] = {1, 2, 0};
    static constexpr std::size_t second[] = {2, 0, 1};
    return SymbolicExpr::add(
        SymbolicExpr::multiply(left[first[index]], right[second[index]]),
        SymbolicExpr::multiply(
            SymbolicExpr::number(-1),
            SymbolicExpr::multiply(
                left[second[index]], right[first[index]])));
}

VectorCalculusFieldResult cross_product(
    const VectorField& a, const VectorField& b) {
    if (a.size() != 3 || b.size() != 3) {
        return VectorCalculusFieldResult::failure(
            CasErrc::DimensionMismatch,
            "vectors must be three-dimensional",
            "cross_product");
    }
    VectorField result;
    result.reserve(3);
    for (std::size_t index = 0; index < 3; ++index) {
        result.push_back(
            vector_calculus_cross_component(a, b, index)->simplify());
    }
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
    if (LMCAS::detail::node(bb) && LMCAS::detail::node(bb)->is_zero()) {
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
    if (LMCAS::detail::node(bb) && LMCAS::detail::node(bb)->is_zero()) {
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
    if (!aa_result) {
        return aa_result;
    }
    auto bb_result = dot_product(b, b);
    if (!bb_result) {
        return bb_result;
    }
    auto aa = aa_result.value();
    auto bb = bb_result.value();
    if ((LMCAS::detail::node(aa) && LMCAS::detail::node(aa)->is_zero()) ||
        (LMCAS::detail::node(bb) && LMCAS::detail::node(bb)->is_zero())) {
        return VectorCalculusExprResult::success(nullptr);
    }
    auto ab_result = dot_product(a, b);
    if (!ab_result) {
        return ab_result;
    }
    auto denom = SymbolicExpr::multiply(SymbolicExpr::sqrt(aa), SymbolicExpr::sqrt(bb));
    auto cos_theta = SymbolicExpr::divide(ab_result.value(), denom);
    auto arccos_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::ArcCos,
        std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(cos_theta)});
    return VectorCalculusExprResult::success(
        LMCAS::detail::make_expression_ptr(arccos_node)->simplify());
}

VectorCalculusExprResult mixed_product(
    const VectorField& a, const VectorField& b, const VectorField& c)
{
    auto cross = cross_product(b, c);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());
    return dot_product(a, cross.value());
}

}
