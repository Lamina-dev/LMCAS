#pragma once

#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "interval.hpp"
#include "numeric_evaluation.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "symbolic_ast.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <unordered_set>

namespace lamina {

bool is_integer_number(const NumberNode& number);
bool is_even_integer_number(const NumberNode& number);
bool is_positive_integer_number(const NumberNode& number);
bool is_zero_number(const NumberNode& number);
bool is_exponent_neg_one(const std::shared_ptr<const SymbolicNode>& node);
std::optional<double> endpoint_to_double(const Endpoint& endpoint);

inline int infinity_sign(const SymbolicExpr& expression) {
    const auto& node = detail::node(expression);
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return function->type() == FunctionNode::FuncType::Infinity ? 1 : 0;
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!product) return 0;

    bool has_infinity = false;
    int sign = 1;
    for (const auto& operand : product->operands()) {
        if (auto function = std::dynamic_pointer_cast<const FunctionNode>(operand)) {
            has_infinity = has_infinity ||
                function->type() == FunctionNode::FuncType::Infinity;
        }
        if (auto number = std::dynamic_pointer_cast<const NumberNode>(operand)) {
            if (std::holds_alternative<BigInt>(number->value()) &&
                std::get<BigInt>(number->value()) < BigInt(0)) {
                sign = -sign;
            } else if (std::holds_alternative<Rational>(number->value()) &&
                       std::get<Rational>(number->value()) < Rational(0)) {
                sign = -sign;
            } else if (std::holds_alternative<lmmc_real_t>(number->value()) &&
                       std::get<lmmc_real_t>(number->value()) < 0) {
                sign = -sign;
            }
        }
    }
    return has_infinity ? sign : 0;
}

template <typename T>
Result<T> checked_inference_result(
    const SymbolicExpr& expr,
    const std::string& operation,
    const std::function<T()>& query) {
    if (!detail::node(expr)) {
        return Result<T>::failure(
            CasErrc::InvalidArgument, "inference expression must not be null", operation);
    }
    try {
        return Result<T>::success(query());
    } catch (const detail::ResultPropagation& error) {
        return Result<T>::failure(error.error());
    } catch (const std::bad_alloc&) {
        return Result<T>::failure(
            CasErrc::ResourceLimit, "inference query allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<T>::failure(CasErrc::InternalInvariant, error.what(), operation);
    }
}

inline Tribool inference_query_or_unknown(InferenceTriboolResult result) {
    return result ? result.value() : Tribool::Unknown;
}

inline std::optional<SymbolicExpr> inference_period_or_empty(InferencePeriodResult result) {
    return result ? result.value() : std::nullopt;
}

struct InferenceEngine::Impl {
    explicit Impl(const AssumptionContext& context)
        : ctx(context), max_depth(context.get_max_query_depth()) {}

    const AssumptionContext& ctx;
    int max_depth;
    mutable std::unordered_set<const void*> visited;
    mutable int current_depth = 0;
};

class InferenceEngine::DepthGuard {
public:
    DepthGuard(const InferenceEngine& engine, const void* node);
    ~DepthGuard();

    bool should_abort() const { return abort_; }

    DepthGuard(const DepthGuard&) = delete;
    DepthGuard& operator=(const DepthGuard&) = delete;

private:
    const InferenceEngine& engine_;
    const void* node_;
    bool abort_ = false;
    bool inserted_ = false;
};

} // namespace lamina
