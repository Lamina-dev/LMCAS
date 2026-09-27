#pragma once

#include "computation_context.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include <exception>
#include <memory>
#include <string>
#include <utility>

namespace LMCAS::detail {

template <typename Transform>
ExpressionResult checked_transform_expr(const ExprPtr& expression,
                                        ComputationContext& context,
                                        const char* operation,
                                        const char* transform_name,
                                        Transform&& transform) {
    auto step = context.consume_steps(1, operation);
    if (!step) return ExpressionResult::failure(step.error());
    if (!expression || !LMCAS::detail::node(expression)) {
        return ExpressionResult::failure(CasErrc::InvalidArgument,
                                         std::string(transform_name) +
                                             " argument cannot be null",
                                         operation);
    }
    try {
        auto result = std::forward<Transform>(transform)(*expression);
        if (!result || !LMCAS::detail::node(result)) {
            return ExpressionResult::failure(CasErrc::InternalInvariant,
                                             std::string(transform_name) +
                                                 " returned null",
                                             operation);
        }
        return ExpressionResult::success(std::move(result));
    } catch (const CasError& error) {
        return ExpressionResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                         std::string(transform_name) +
                                             " allocation failed",
                                         operation);
    } catch (const detail::UnsupportedDifferentiation& error) {
        return ExpressionResult::failure(CasErrc::UnsupportedExpression, error.what(),
                                         operation);
    } catch (const std::exception& error) {
        return ExpressionResult::failure(CasErrc::InvalidArgument, error.what(),
                                         operation);
    }
}

}
