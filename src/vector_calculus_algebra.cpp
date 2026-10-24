#include "internal/expression_transform.hpp"
#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_vectors.hpp"
#include "internal/symbolic_ast.hpp"


namespace LMCAS {

using namespace vector_calculus_detail;

static VectorCalculusExprResult dot_product_impl(const VectorField& a, const VectorField& b, detail::RewriteBudget& transform_budget)
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
    return VectorCalculusExprResult::success(detail::make_expression_ptr(detail::normalize_node(detail::node(sum), transform_budget)));
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

static VectorCalculusFieldResult cross_product_impl(
    const VectorField& a, const VectorField& b, detail::RewriteBudget& transform_budget) {
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
            detail::make_expression_ptr(detail::normalize_node(detail::node(vector_calculus_cross_component(a, b, index)), transform_budget)));
    }
    return VectorCalculusFieldResult::success(std::move(result));
}

static VectorCalculusFieldResult vector_project_impl(const VectorField& a, const VectorField& b, detail::RewriteBudget& transform_budget)
{
    if (a.size() != b.size()) {
        return VectorCalculusFieldResult::failure(
            CasErrc::DimensionMismatch,
            "vectors must have the same dimension",
            "vector_project");
    }
    auto bb_result = dot_product_impl(b, b, transform_budget);
    if (!bb_result) return VectorCalculusFieldResult::failure(bb_result.error());
    auto bb = bb_result.value();
    if (LMCAS::detail::node(bb) && LMCAS::detail::node(bb)->is_zero()) {
        return VectorCalculusFieldResult::success(
            VectorField(a.size(), SymbolicExpr::number(0)));
    }
    auto ab_result = dot_product_impl(a, b, transform_budget);
    if (!ab_result) return VectorCalculusFieldResult::failure(ab_result.error());
    auto ab = ab_result.value();
    auto coeff = SymbolicExpr::divide(ab, bb);
    VectorField result;
    result.reserve(b.size());
    for (const auto& comp : b) {
        result.push_back(detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(coeff, comp)), transform_budget)));
    }
    return VectorCalculusFieldResult::success(std::move(result));
}

static VectorCalculusExprResult scalar_project_impl(const VectorField& a, const VectorField& b, detail::RewriteBudget& transform_budget)
{
    auto bb_result = dot_product_impl(b, b, transform_budget);
    if (!bb_result) return bb_result;
    auto bb = bb_result.value();
    if (LMCAS::detail::node(bb) && LMCAS::detail::node(bb)->is_zero()) {
        return VectorCalculusExprResult::success(nullptr);
    }
    auto ab_result = dot_product_impl(a, b, transform_budget);
    if (!ab_result) return ab_result;
    return VectorCalculusExprResult::success(
        detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(ab_result.value(), SymbolicExpr::sqrt(bb))), transform_budget)));
}

static VectorCalculusExprResult vector_angle_symbolic_impl(const VectorField& a, const VectorField& b, detail::RewriteBudget& transform_budget)
{
    auto aa_result = dot_product_impl(a, a, transform_budget);
    if (!aa_result) {
        return aa_result;
    }
    auto bb_result = dot_product_impl(b, b, transform_budget);
    if (!bb_result) {
        return bb_result;
    }
    auto aa = aa_result.value();
    auto bb = bb_result.value();
    if ((LMCAS::detail::node(aa) && LMCAS::detail::node(aa)->is_zero()) ||
        (LMCAS::detail::node(bb) && LMCAS::detail::node(bb)->is_zero())) {
        return VectorCalculusExprResult::success(nullptr);
    }
    auto ab_result = dot_product_impl(a, b, transform_budget);
    if (!ab_result) {
        return ab_result;
    }
    auto denom = SymbolicExpr::multiply(SymbolicExpr::sqrt(aa), SymbolicExpr::sqrt(bb));
    auto cos_theta = SymbolicExpr::divide(ab_result.value(), denom);
    auto arccos_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::ArcCos,
        std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(cos_theta)});
    return VectorCalculusExprResult::success(
        detail::make_expression_ptr(detail::normalize_node(detail::node(LMCAS::detail::make_expression_ptr(arccos_node)), transform_budget)));
}

static VectorCalculusExprResult mixed_product_impl(
    const VectorField& a, const VectorField& b, const VectorField& c, detail::RewriteBudget& transform_budget)
{
    auto cross = cross_product_impl(b, c, transform_budget);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());
    return dot_product_impl(a, cross.value(), transform_budget);
}


VectorCalculusExprResult dot_product(const VectorField& a, const VectorField& b) try {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "dot_product");
    return dot_product_impl(a, b, transform_budget);
} catch (const CasError& error) {
    return VectorCalculusExprResult::failure(error);
} catch (const std::bad_alloc&) {
    return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
        "vector calculus allocation failed", "dot_product");
} catch (const std::length_error& error) {
    return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), "dot_product");
} catch (const std::exception& error) {
    return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, error.what(), "dot_product");
}

VectorCalculusFieldResult cross_product(
    const VectorField& a, const VectorField& b) try {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "cross_product");
    return cross_product_impl(a, b, transform_budget);
} catch (const CasError& error) {
    return VectorCalculusFieldResult::failure(error);
} catch (const std::bad_alloc&) {
    return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit,
        "vector calculus allocation failed", "cross_product");
} catch (const std::length_error& error) {
    return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit, error.what(), "cross_product");
} catch (const std::exception& error) {
    return VectorCalculusFieldResult::failure(CasErrc::InternalInvariant, error.what(), "cross_product");
}

VectorCalculusFieldResult vector_project(const VectorField& a, const VectorField& b) try {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "vector_project");
    return vector_project_impl(a, b, transform_budget);
} catch (const CasError& error) {
    return VectorCalculusFieldResult::failure(error);
} catch (const std::bad_alloc&) {
    return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit,
        "vector calculus allocation failed", "vector_project");
} catch (const std::length_error& error) {
    return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit, error.what(), "vector_project");
} catch (const std::exception& error) {
    return VectorCalculusFieldResult::failure(CasErrc::InternalInvariant, error.what(), "vector_project");
}

VectorCalculusExprResult scalar_project(const VectorField& a, const VectorField& b) try {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "scalar_project");
    return scalar_project_impl(a, b, transform_budget);
} catch (const CasError& error) {
    return VectorCalculusExprResult::failure(error);
} catch (const std::bad_alloc&) {
    return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
        "vector calculus allocation failed", "scalar_project");
} catch (const std::length_error& error) {
    return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), "scalar_project");
} catch (const std::exception& error) {
    return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, error.what(), "scalar_project");
}

VectorCalculusExprResult vector_angle_symbolic(const VectorField& a, const VectorField& b) try {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "vector_angle_symbolic");
    return vector_angle_symbolic_impl(a, b, transform_budget);
} catch (const CasError& error) {
    return VectorCalculusExprResult::failure(error);
} catch (const std::bad_alloc&) {
    return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
        "vector calculus allocation failed", "vector_angle_symbolic");
} catch (const std::length_error& error) {
    return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), "vector_angle_symbolic");
} catch (const std::exception& error) {
    return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, error.what(), "vector_angle_symbolic");
}

VectorCalculusExprResult mixed_product(
    const VectorField& a, const VectorField& b, const VectorField& c) try {
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "mixed_product");
    return mixed_product_impl(a, b, c, transform_budget);
} catch (const CasError& error) {
    return VectorCalculusExprResult::failure(error);
} catch (const std::bad_alloc&) {
    return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
        "vector calculus allocation failed", "mixed_product");
} catch (const std::length_error& error) {
    return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), "mixed_product");
} catch (const std::exception& error) {
    return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, error.what(), "mixed_product");
}

}
