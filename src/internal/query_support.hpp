#pragma once

#include <utility>

#include "query_interface.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS::query_detail {


template <typename T, typename Query>
Result<T> checked_expression_result(
    const SymbolicExpr& expr,
    const std::string& operation,
    Query&& query) {
    if (!LMCAS::detail::node(expr)) {
        return Result<T>::failure(
            CasErrc::InvalidArgument, "query expression must not be null", operation);
    }
    try {
        return Result<T>::success(std::forward<Query>(query)());
    } catch (const std::bad_alloc&) {
        return Result<T>::failure(
            CasErrc::ResourceLimit, "query allocation failed", operation);
    } catch (const std::exception& ex) {
        return Result<T>::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
}

template <typename Query>
auto checked_expression_result(
    const SymbolicExpr& expr,
    const std::string& operation,
    Query&& query) -> decltype(std::forward<Query>(query)()) {
    using ResultType = decltype(std::forward<Query>(query)());
    if (!LMCAS::detail::node(expr)) {
        return ResultType::failure(
            CasErrc::InvalidArgument, "query expression must not be null", operation);
    }
    try {
        return std::forward<Query>(query)();
    } catch (const std::bad_alloc&) {
        return ResultType::failure(
            CasErrc::ResourceLimit, "query allocation failed", operation);
    } catch (const std::exception& ex) {
        return ResultType::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

}
