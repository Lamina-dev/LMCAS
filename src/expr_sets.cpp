#include "expr.hpp"

#include <exception>
#include <utility>
#include <vector>

#include "internal/assumption_simplification.hpp"
#include "internal/expr_internal.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS {

ExprSet::ExprSet()
    : ExprSet(std::make_shared<const SymbolicExpr>(detail::expression_from_node(
          detail::make_node<FiniteSetNode>(
              std::vector<std::shared_ptr<const SymbolicNode>>{})))) {}

ExprSet::ExprSet(ConstExprPtr expression)
    : expression_(std::move(expression)) {
    const auto& node = static_cast<const FiniteSetNode&>(*detail::node(*expression_));
    elements_.reserve(node.elements().size());
    for (const auto& element : node.elements()) {
        elements_.push_back(std::make_shared<const SymbolicExpr>(
            detail::expression_from_node(element)));
    }
}

Result<ExprSet> ExprSet::make(std::vector<ExprPtr> elements) {
    ComputationContext context;
    for (const auto& element : elements) {
        if (!element) {
            return expr_set_failure(CasErrc::InvalidArgument,
                                    "set<Expr> elements cannot be null",
                                    kExprSetOperation);
        }
    }
    auto expression = make_finite_set(std::move(elements), context);
    if (!expression) { return Result<ExprSet>::failure(expression.error()); }
    return Result<ExprSet>::success(
        ExprSet(std::make_shared<const SymbolicExpr>(*expression.value())));
}

bool ExprSet::contains(const SymbolicExpr& expression) const {
    const auto& node = static_cast<const FiniteSetNode&>(*detail::node(*expression_));
    return node.contains(*detail::node(expression));
}

bool ExprSet::subset_of(const ExprSet& other) const {
    for (const auto& element : elements_) {
        if (!other.contains(*element)) {
            return false;
        }
    }
    return true;
}

ExprSet ExprSet::set_union(const ExprSet& other) const {
    ComputationContext context;
    auto result = finite_set_union(*expression_, *other.expression_, context);
    return ExprSet(std::make_shared<const SymbolicExpr>(*result.value()));
}

ExprSet ExprSet::intersection(const ExprSet& other) const {
    ComputationContext context;
    auto result = finite_set_intersection(*expression_, *other.expression_, context);
    return ExprSet(std::make_shared<const SymbolicExpr>(*result.value()));
}

ExprSet ExprSet::difference(const ExprSet& other) const {
    ComputationContext context;
    auto result = finite_set_difference(*expression_, *other.expression_, context);
    return ExprSet(std::make_shared<const SymbolicExpr>(*result.value()));
}

ExprSet ExprSet::symmetric_difference(const ExprSet& other) const {
    ComputationContext context;
    auto result = finite_set_symmetric_difference(*expression_, *other.expression_, context);
    return ExprSet(std::make_shared<const SymbolicExpr>(*result.value()));
}


ExprSetResult expr_set(std::vector<ExprPtr> elements) {
    try {
        return ExprSet::make(std::move(elements));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(CasErrc::ResourceLimit,
                                "set<Expr> allocation failed",
                                kExprSetOperation);
    } catch (const std::exception& error) {
        return expr_set_failure(CasErrc::InternalInvariant, error.what(),
                                kExprSetOperation);
    }
}


Result<bool> expr_set_contains(const ExprSet& set,
                               const ExprPtr& element) {
    if (!element) {
        return Result<bool>::failure(CasErrc::InvalidArgument,
                                     "set<Expr> membership element cannot be null",
                                     kExprSetOperation);
    }
    return Result<bool>::success(set.contains(*element));
}

Result<bool> expr_set_not_contains(const ExprSet& set,
                                   const ExprPtr& element) {
    auto result = expr_set_contains(set, element);
    if (!result) { return result; }
    return Result<bool>::success(!result.value());
}

Result<bool> expr_set_subset(const ExprSet& lhs,
                             const ExprSet& rhs) {
    return Result<bool>::success(lhs.subset_of(rhs));
}

Result<bool> expr_set_subset_domain(const ExprSet& set,
                                    const NumberDomainSet& domain) {
    for (const auto& element : set.elements()) {
        auto contains = domain.contains(element);
        if (!contains) {
            return contains;
        }
        if (!contains.value()) {
            return Result<bool>::success(false);
        }
    }
    return Result<bool>::success(true);
}

ExprSetResult expr_set_union(
    const ExprSet& lhs, const ExprSet& rhs)
{
    try {
        ComputationContext context;
        auto expression =
            finite_set_union(*lhs.expression_, *rhs.expression_, context);
        if (!expression) {
            return ExprSetResult::failure(expression.error());
        }
        return ExprSetResult::success(ExprSet(
            std::make_shared<const SymbolicExpr>(*expression.value())));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(
            CasErrc::ResourceLimit, "set<Expr> union allocation failed",
            kExprSetOperation);
    }
}

ExprSetResult expr_set_intersection(
    const ExprSet& lhs, const ExprSet& rhs)
{
    try {
        ComputationContext context;
        auto expression = finite_set_intersection(
            *lhs.expression_, *rhs.expression_, context);
        if (!expression) {
            return ExprSetResult::failure(expression.error());
        }
        return ExprSetResult::success(ExprSet(
            std::make_shared<const SymbolicExpr>(*expression.value())));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(
            CasErrc::ResourceLimit,
            "set<Expr> intersection allocation failed",
            kExprSetOperation);
    }
}

ExprSetResult expr_set_difference(
    const ExprSet& lhs, const ExprSet& rhs)
{
    try {
        ComputationContext context;
        auto expression =
            finite_set_difference(*lhs.expression_, *rhs.expression_, context);
        if (!expression) {
            return ExprSetResult::failure(expression.error());
        }
        return ExprSetResult::success(ExprSet(
            std::make_shared<const SymbolicExpr>(*expression.value())));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(
            CasErrc::ResourceLimit,
            "set<Expr> difference allocation failed",
            kExprSetOperation);
    }
}

ExprSetResult expr_set_symmetric_difference(
    const ExprSet& lhs, const ExprSet& rhs)
{
    try {
        ComputationContext context;
        auto expression = finite_set_symmetric_difference(
            *lhs.expression_, *rhs.expression_, context);
        if (!expression) {
            return ExprSetResult::failure(expression.error());
        }
        return ExprSetResult::success(ExprSet(
            std::make_shared<const SymbolicExpr>(*expression.value())));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(
            CasErrc::ResourceLimit,
            "set<Expr> symmetric difference allocation failed",
            kExprSetOperation);
    }
}

}
