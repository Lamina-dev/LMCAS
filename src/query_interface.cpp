/**
 * @file query_interface.cpp
 * @brief Implementation of the QueryInterface class.
 *
 * The QueryInterface is the single public entry point for property queries.
 * It handles special cases (null, NaN, Infinity, Matrix, Relational, Logical)
 * and delegates to the InferenceEngine for all other node types.
 *
 * Query results are cached by (expression_hash, property_type) to avoid
 * redundant inference. The cache is invalidated on any AssumptionContext mutation.
 */

#include "query_interface.hpp"
#include "symbolic_ast.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include <cmath>

namespace LMCAS {

namespace {

Tribool query_or_unknown(QueryTriboolResult result) {
    if (!result) {
        return Tribool::Unknown;
    }
    return result.value();
}


QueryTriboolResult checked_query_result(
    const SymbolicExpr& expr,
    const std::string& operation,
    const std::function<QueryTriboolResult()>& query) {
    if (!LMCAS::detail::node(expr)) {
        return QueryTriboolResult::failure(
            CasErrc::InvalidArgument,
            "query expression must not be null", operation);
    }
    try {
        return query();
    } catch (const std::bad_alloc&) {
        return QueryTriboolResult::failure(
            CasErrc::ResourceLimit, "query allocation failed", operation);
    } catch (const std::exception& ex) {
        return QueryTriboolResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

template <typename T>
Result<T> checked_expression_result(
    const SymbolicExpr& expr,
    const std::string& operation,
    const std::function<T()>& query) {
    if (!LMCAS::detail::node(expr)) {
        return Result<T>::failure(
            CasErrc::InvalidArgument, "query expression must not be null", operation);
    }
    try {
        return Result<T>::success(query());
    } catch (const std::bad_alloc&) {
        return Result<T>::failure(
            CasErrc::ResourceLimit, "query allocation failed", operation);
    } catch (const std::exception& ex) {
        return Result<T>::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
}

} // anonymous namespace

// Construction

QueryInterface::QueryInterface(const AssumptionContext& ctx)
    : ctx_(ctx), observed_generation_(ctx.cache_generation()) {}

// Cache management

void QueryInterface::invalidate_cache() const {
    cache_.clear();
    observed_generation_ = ctx_.cache_generation();
}

Tribool QueryInterface::cached_query(const SymbolicExpr& expr, PropType prop,
                                     const std::function<Tribool()>& compute) const {
    return query_or_unknown(cached_query_checked(expr, prop, "cached_query", compute));
}

QueryTriboolResult QueryInterface::cached_query_checked(
    const SymbolicExpr& expr,
    PropType prop,
    const std::string& operation,
    const std::function<QueryTriboolResult()>& compute) const {
    if (!LMCAS::detail::node(expr)) {
        return QueryTriboolResult::failure(
            CasErrc::InvalidArgument, "query expression must not be null", operation);
    }

    // Check if the context has been mutated since we last validated the cache.
    // If so, invalidate and update our observed generation.
    uint64_t current_gen = ctx_.cache_generation();
    if (current_gen != observed_generation_) {
        cache_.clear();
        observed_generation_ = current_gen;
    }

    CacheKey key{LMCAS::detail::node(expr)->hash(), prop};
    auto it = cache_.find(key);
    if (it != cache_.end()) {
        for (const auto& entry : it->second) {
            if (LMCAS::detail::node(entry.expression)->equals(*LMCAS::detail::node(expr))) {
                return QueryTriboolResult::success(entry.result);
            }
        }
    }

    auto result = compute();
    if (!result) return result;
    cache_[key].push_back(CacheEntry{expr, result.value()});
    return result;
}

// Private helpers

bool QueryInterface::is_unhandled_type(const SymbolicExpr& expression) const {
    const auto& node = LMCAS::detail::node(expression);
    if (!node) return true;
    if (std::dynamic_pointer_cast<const MatrixNode>(node)) return true;
    if (std::dynamic_pointer_cast<const RelationalNode>(node)) return true;
    if (std::dynamic_pointer_cast<const LogicalNode>(node)) return true;
    return false;
}

bool QueryInterface::is_nan_number(const SymbolicExpr& expression) const {
    const auto node = std::dynamic_pointer_cast<const NumberNode>(
        LMCAS::detail::node(expression));
    if (node && std::holds_alternative<lmmc_real_t>(node->value())) {
        lmmc_real_t v = std::get<lmmc_real_t>(node->value());
        return std::isnan(v);
    }
    return false;
}

bool QueryInterface::is_infinity_node(const SymbolicExpr& expression) const {
    const auto& node = LMCAS::detail::node(expression);
    // Direct infinity: FunctionNode with FuncType::Infinity
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return func->type() == FunctionNode::FuncType::Infinity;
    }
    // Negative infinity: MultiplyNode(-1, Infinity)
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        if (mul->operands().size() == 2) {
            for (const auto& op : mul->operands()) {
                if (auto func = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                    if (func->type() == FunctionNode::FuncType::Infinity) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

int QueryInterface::get_infinity_sign(const SymbolicExpr& expression) const {
    const auto& node = LMCAS::detail::node(expression);
    // Direct infinity node -> positive infinity
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->type() == FunctionNode::FuncType::Infinity) {
            return +1;
        }
    }
    // Negative infinity: MultiplyNode containing -1 and Infinity
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        bool has_infinity = false;
        bool has_neg_one = false;
        for (const auto& op : mul->operands()) {
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                if (func->type() == FunctionNode::FuncType::Infinity) {
                    has_infinity = true;
                }
            }
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(op)) {
                // Check if it's -1
                if (std::holds_alternative<BigInt>(num->value())) {
                    if (std::get<BigInt>(num->value()) == BigInt(-1)) {
                        has_neg_one = true;
                    }
                } else if (std::holds_alternative<Rational>(num->value())) {
                    if (std::get<Rational>(num->value()) == Rational(-1)) {
                        has_neg_one = true;
                    }
                } else if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    if (v == -1.0) {
                        has_neg_one = true;
                    }
                }
            }
        }
        if (has_infinity && has_neg_one) return -1;
        if (has_infinity) return +1; // Infinity multiplied by something else positive
    }
    return 0; // indeterminate
}

// Public query methods

QueryTriboolResult QueryInterface::query_positive(const SymbolicExpr& expr) const {
    return query_positive_checked(expr);
}

QueryTriboolResult QueryInterface::query_positive_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_positive_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Positive, "query_positive_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        int sign = get_infinity_sign(expr);
                        if (sign > 0) return Tribool::True;
                        if (sign < 0) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_positive_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_negative(const SymbolicExpr& expr) const {
    return query_negative_checked(expr);
}

QueryTriboolResult QueryInterface::query_negative_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_negative_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Negative, "query_negative_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        int sign = get_infinity_sign(expr);
                        if (sign < 0) return Tribool::True;
                        if (sign > 0) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_negative_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_nonnegative(const SymbolicExpr& expr) const {
    return query_nonnegative_checked(expr);
}

QueryTriboolResult QueryInterface::query_nonnegative_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_nonnegative_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::NonNegative, "query_nonnegative_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        int sign = get_infinity_sign(expr);
                        if (sign > 0) return Tribool::True;
                        if (sign < 0) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_nonnegative_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_real(const SymbolicExpr& expr) const {
    return query_real_checked(expr);
}

QueryTriboolResult QueryInterface::query_real_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_real_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Real, "query_real_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::Unknown;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_real_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_integer(const SymbolicExpr& expr) const {
    return query_integer_checked(expr);
}

QueryTriboolResult QueryInterface::query_integer_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_integer_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Integer, "query_integer_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::False;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_integer_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_nonzero(const SymbolicExpr& expr) const {
    return query_nonzero_checked(expr);
}

QueryTriboolResult QueryInterface::query_nonzero_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_nonzero_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::NonZero, "query_nonzero_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::True;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_nonzero_checked(expr);
                });
            return result;
        });
}

// Extended property queries

QueryTriboolResult QueryInterface::query_algebraic(const SymbolicExpr& expr) const {
    return query_algebraic_checked(expr);
}

QueryTriboolResult QueryInterface::query_algebraic_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_algebraic_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Algebraic, "query_algebraic_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_algebraic_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_transcendental(const SymbolicExpr& expr) const {
    return query_transcendental_checked(expr);
}

QueryTriboolResult QueryInterface::query_transcendental_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_transcendental_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Transcendental, "query_transcendental_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_transcendental_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_finite(const SymbolicExpr& expr) const {
    return query_finite_checked(expr);
}

QueryTriboolResult QueryInterface::query_finite_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_finite_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Finite, "query_finite_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_finite_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_divergent(const SymbolicExpr& expr) const {
    return query_divergent_checked(expr);
}

QueryTriboolResult QueryInterface::query_divergent_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_divergent_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Divergent, "query_divergent_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::True;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_divergent_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_periodic(const SymbolicExpr& expr) const {
    return query_periodic_checked(expr);
}

QueryTriboolResult QueryInterface::query_periodic_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_periodic_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Periodic, "query_periodic_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::False;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_periodic_checked(expr);
                });
            return result;
        });
}

QueryPeriodResult QueryInterface::get_period(const SymbolicExpr& expr) const {
    return get_period_checked(expr);
}


QueryPeriodResult QueryInterface::get_period_checked(
    const SymbolicExpr& expr) const {
    if (!LMCAS::detail::node(expr)) {
        return QueryPeriodResult::failure(
            CasErrc::InvalidArgument,
            "query expression must not be null",
            "get_period");
    }
    if (is_unhandled_type(expr) || is_infinity_node(expr)) {
        return QueryPeriodResult::success(std::nullopt);
    }
    InferenceEngine engine(ctx_);
    return engine.infer_period_checked(expr);
}

QueryTriboolResult QueryInterface::query_positive_definite(
    const SymbolicExpr& expr) const {
    return query_positive_definite_checked(expr);
}

QueryTriboolResult QueryInterface::query_positive_definite_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_positive_definite_checked",
        [&]() {
            auto result = cached_query_checked(
                expr,
                PropType::PositiveDefinite,
                "query_positive_definite_checked",
                [&]() -> QueryTriboolResult {
                    if (auto var = std::dynamic_pointer_cast<const VariableNode>(
                            LMCAS::detail::node(expr))) {
                        const auto& props = ctx_.current_properties();
                        Definiteness d = props.get_definiteness(var->name());
                        if (d == Definiteness::PositiveDefinite) return Tribool::True;
                        if (d == Definiteness::NegativeDefinite ||
                            d == Definiteness::NegativeSemiDefinite ||
                            d == Definiteness::Indefinite) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    return Tribool::Unknown;
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_positive_semidefinite(
    const SymbolicExpr& expr) const {
    return query_positive_semidefinite_checked(expr);
}

QueryTriboolResult QueryInterface::query_positive_semidefinite_checked(const SymbolicExpr& expr) const {
    return checked_query_result(expr, "query_positive_semidefinite_checked",
        [&]() {
            auto result = cached_query_checked(
                expr,
                PropType::PositiveSemiDefinite,
                "query_positive_semidefinite_checked",
                [&]() -> QueryTriboolResult {
                    if (auto var = std::dynamic_pointer_cast<const VariableNode>(
                            LMCAS::detail::node(expr))) {
                        const auto& props = ctx_.current_properties();
                        Definiteness d = props.get_definiteness(var->name());
                        if (d == Definiteness::PositiveDefinite ||
                            d == Definiteness::PositiveSemiDefinite) return Tribool::True;
                        if (d == Definiteness::NegativeDefinite ||
                            d == Definiteness::NegativeSemiDefinite ||
                            d == Definiteness::Indefinite) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    return Tribool::Unknown;
                });
            return result;
        });
}

// Query conditions (query mode)

QueryInterface::QueryConditionSetsResult QueryInterface::query_conditions(
    const SymbolicExpr& expr, Sign target) const {
    return query_conditions_checked(expr, target);
}

std::vector<QueryInterface::ConditionSet> QueryInterface::query_conditions_impl(
    const SymbolicExpr& expr, Sign target) const {
    if (!LMCAS::detail::node(expr)) {
        return {};
    }
    if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
        ConditionSet cs;
        cs.sign_conditions.emplace_back(var->name(), target);
        return {cs};
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        if (add->operands().size() == 2) {
            std::shared_ptr<const SymbolicNode> pos_operand = nullptr;
            std::shared_ptr<const SymbolicNode> neg_operand = nullptr; // the term being subtracted

            for (size_t i = 0; i < 2; ++i) {
                auto mul = std::dynamic_pointer_cast<const MultiplyNode>(add->operands()[i]);
                if (mul && mul->operands().size() == 2) {
                    for (size_t j = 0; j < 2; ++j) {
                        auto num = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[j]);
                        if (num) {
                            bool is_neg_one = false;
                            if (std::holds_alternative<BigInt>(num->value())) {
                                is_neg_one = (std::get<BigInt>(num->value()) == BigInt(-1));
                            } else if (std::holds_alternative<Rational>(num->value())) {
                                is_neg_one = (std::get<Rational>(num->value()) == Rational(-1));
                            } else if (std::holds_alternative<lmmc_real_t>(num->value())) {
                                is_neg_one = (std::get<lmmc_real_t>(num->value()) == -1.0);
                            }
                            if (is_neg_one) {
                                neg_operand = mul->operands()[1 - j];
                                pos_operand = add->operands()[1 - i];
                                break;
                            }
                        }
                    }
                }
                if (neg_operand) break;
            }
            if (pos_operand && neg_operand) {
                auto lhs_var = std::dynamic_pointer_cast<const VariableNode>(pos_operand);
                auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(neg_operand);

                if (lhs_var && rhs_var && target == Sign::Positive) {
                    ConditionSet cs1;
                    cs1.sign_conditions.emplace_back(lhs_var->name(), Sign::Positive);
                    cs1.sign_conditions.emplace_back(rhs_var->name(), Sign::Negative);
                    ConditionSet cs2;
                    auto lhs_expr = LMCAS::detail::expression_from_node(pos_operand);
                    auto rhs_expr = LMCAS::detail::expression_from_node(neg_operand);
                    Relation rel{lhs_expr, rhs_expr, RelationalNode::Op::GT};
                    cs2.relational_conditions.push_back(rel);
                    cs2.sign_conditions.emplace_back(rhs_var->name(), Sign::NonNegative);

                    return {cs1, cs2};
                } else if (lhs_var && rhs_var && target == Sign::Negative) {
                    ConditionSet cs1;
                    cs1.sign_conditions.emplace_back(lhs_var->name(), Sign::Negative);
                    cs1.sign_conditions.emplace_back(rhs_var->name(), Sign::Positive);

                    ConditionSet cs2;
                    auto lhs_expr = LMCAS::detail::expression_from_node(pos_operand);
                    auto rhs_expr = LMCAS::detail::expression_from_node(neg_operand);
                    Relation rel{rhs_expr, lhs_expr, RelationalNode::Op::GT};
                    cs2.relational_conditions.push_back(rel);
                    cs2.sign_conditions.emplace_back(lhs_var->name(), Sign::NonNegative);

                    return {cs1, cs2};
                } else if (lhs_var && rhs_var && target == Sign::NonNegative) {
                    ConditionSet cs1;
                    cs1.sign_conditions.emplace_back(lhs_var->name(), Sign::NonNegative);
                    cs1.sign_conditions.emplace_back(rhs_var->name(), Sign::NonPositive);

                    ConditionSet cs2;
                    auto lhs_expr = LMCAS::detail::expression_from_node(pos_operand);
                    auto rhs_expr = LMCAS::detail::expression_from_node(neg_operand);
                    Relation rel{lhs_expr, rhs_expr, RelationalNode::Op::GEQ};
                    cs2.relational_conditions.push_back(rel);

                    return {cs1, cs2};
                }
            }
        }
        if (target == Sign::Positive || target == Sign::NonNegative) {
            bool all_variables = true;
            ConditionSet cs;
            for (const auto& op : add->operands()) {
                if (auto var = std::dynamic_pointer_cast<const VariableNode>(op)) {
                    cs.sign_conditions.emplace_back(var->name(), target);
                } else {
                    all_variables = false;
                    break;
                }
            }
            if (all_variables && !cs.sign_conditions.empty()) {
                return {cs};
            }
        }
        return {};
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        if (mul->operands().size() == 2) {
            std::shared_ptr<const SymbolicNode> numerator = nullptr;
            std::shared_ptr<const SymbolicNode> denominator = nullptr;

            for (size_t i = 0; i < 2; ++i) {
                auto pow = std::dynamic_pointer_cast<const PowerNode>(mul->operands()[i]);
                if (pow) {
                    auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
                    if (exp_num) {
                        bool is_neg_one = false;
                        if (std::holds_alternative<BigInt>(exp_num->value())) {
                            is_neg_one = (std::get<BigInt>(exp_num->value()) == BigInt(-1));
                        } else if (std::holds_alternative<Rational>(exp_num->value())) {
                            is_neg_one = (std::get<Rational>(exp_num->value()) == Rational(-1));
                        } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
                            is_neg_one = (std::get<lmmc_real_t>(exp_num->value()) == -1.0);
                        }
                        if (is_neg_one) {
                            denominator = pow->base();
                            numerator = mul->operands()[1 - i];
                            break;
                        }
                    }
                }
            }
            if (numerator && denominator) {
                auto num_var = std::dynamic_pointer_cast<const VariableNode>(numerator);
                auto den_var = std::dynamic_pointer_cast<const VariableNode>(denominator);

                if (num_var && den_var) {
                    if (target == Sign::Positive) {
                        // x/y > 0: both positive OR both negative
                        ConditionSet cs1;
                        cs1.sign_conditions.emplace_back(num_var->name(), Sign::Positive);
                        cs1.sign_conditions.emplace_back(den_var->name(), Sign::Positive);

                        ConditionSet cs2;
                        cs2.sign_conditions.emplace_back(num_var->name(), Sign::Negative);
                        cs2.sign_conditions.emplace_back(den_var->name(), Sign::Negative);

                        return {cs1, cs2};
                    } else if (target == Sign::Negative) {
                        // x/y < 0: one positive, one negative
                        ConditionSet cs1;
                        cs1.sign_conditions.emplace_back(num_var->name(), Sign::Positive);
                        cs1.sign_conditions.emplace_back(den_var->name(), Sign::Negative);

                        ConditionSet cs2;
                        cs2.sign_conditions.emplace_back(num_var->name(), Sign::Negative);
                        cs2.sign_conditions.emplace_back(den_var->name(), Sign::Positive);

                        return {cs1, cs2};
                    }
                }
                return {};
            }
        }
        if (target == Sign::Positive) {
            bool all_variables = true;
            ConditionSet cs;
            for (const auto& op : mul->operands()) {
                if (auto var = std::dynamic_pointer_cast<const VariableNode>(op)) {
                    cs.sign_conditions.emplace_back(var->name(), Sign::Positive);
                } else {
                    all_variables = false;
                    break;
                }
            }
            if (all_variables && !cs.sign_conditions.empty()) {
                return {cs};
            }
        }

        return {};
    }
    return {};
}

QueryInterface::QueryConditionSetsResult QueryInterface::query_conditions_checked(
    const SymbolicExpr& expr, Sign target) const {
    return checked_expression_result<std::vector<ConditionSet>>(
        expr, "query_conditions", [&]() {
            return query_conditions_impl(expr, target);
        });
}

} // namespace LMCAS
