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
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include <cmath>

namespace lamina {

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
    // Null root cannot be hashed — skip cache
    if (!expr.root) {
        return compute();
    }

    // Check if the context has been mutated since we last validated the cache.
    // If so, invalidate and update our observed generation.
    uint64_t current_gen = ctx_.cache_generation();
    if (current_gen != observed_generation_) {
        cache_.clear();
        observed_generation_ = current_gen;
    }

    CacheKey key{expr.root->hash(), prop};
    auto it = cache_.find(key);
    if (it != cache_.end()) {
        return it->second;
    }

    Tribool result = compute();
    cache_.emplace(key, result);
    return result;
}

// Private helpers

bool QueryInterface::is_unhandled_type(const std::shared_ptr<SymbolicNode>& node) const {
    if (!node) return true;
    if (std::dynamic_pointer_cast<MatrixNode>(node)) return true;
    if (std::dynamic_pointer_cast<RelationalNode>(node)) return true;
    if (std::dynamic_pointer_cast<LogicalNode>(node)) return true;
    return false;
}

bool QueryInterface::is_nan_number(const NumberNode& node) const {
    if (std::holds_alternative<lmmc_real_t>(node.value)) {
        lmmc_real_t v = std::get<lmmc_real_t>(node.value);
        return std::isnan(v);
    }
    return false;
}

bool QueryInterface::is_infinity_node(const std::shared_ptr<SymbolicNode>& node) const {
    // Direct infinity: FunctionNode with FuncType::Infinity
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        return func->type == FunctionNode::FuncType::Infinity;
    }
    // Negative infinity: MultiplyNode(-1, Infinity)
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        if (mul->operands.size() == 2) {
            for (const auto& op : mul->operands) {
                if (auto func = std::dynamic_pointer_cast<FunctionNode>(op)) {
                    if (func->type == FunctionNode::FuncType::Infinity) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

int QueryInterface::get_infinity_sign(const std::shared_ptr<SymbolicNode>& node) const {
    // Direct infinity node → positive infinity
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (func->type == FunctionNode::FuncType::Infinity) {
            return +1;
        }
    }
    // Negative infinity: MultiplyNode containing -1 and Infinity
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        bool has_infinity = false;
        bool has_neg_one = false;
        for (const auto& op : mul->operands) {
            if (auto func = std::dynamic_pointer_cast<FunctionNode>(op)) {
                if (func->type == FunctionNode::FuncType::Infinity) {
                    has_infinity = true;
                }
            }
            if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
                // Check if it's -1
                if (std::holds_alternative<BigInt>(num->value)) {
                    if (std::get<BigInt>(num->value) == BigInt(-1)) {
                        has_neg_one = true;
                    }
                } else if (std::holds_alternative<Rational>(num->value)) {
                    if (std::get<Rational>(num->value) == Rational(-1)) {
                        has_neg_one = true;
                    }
                } else if (std::holds_alternative<lmmc_real_t>(num->value)) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value);
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

Tribool QueryInterface::query_positive(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Positive, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }
        if (is_infinity_node(expr.root)) {
            int sign = get_infinity_sign(expr.root);
            if (sign > 0) return Tribool::True;
            if (sign < 0) return Tribool::False;
            return Tribool::Unknown;
        }
        InferenceEngine engine(ctx_);
        return engine.query_positive(expr);
    });
}

Tribool QueryInterface::query_negative(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Negative, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }
        if (is_infinity_node(expr.root)) {
            int sign = get_infinity_sign(expr.root);
            if (sign < 0) return Tribool::True;
            if (sign > 0) return Tribool::False;
            return Tribool::Unknown;
        }
        InferenceEngine engine(ctx_);
        return engine.query_negative(expr);
    });
}

Tribool QueryInterface::query_nonnegative(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::NonNegative, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }
        if (is_infinity_node(expr.root)) {
            int sign = get_infinity_sign(expr.root);
            if (sign > 0) return Tribool::True;   // +∞ is nonnegative
            if (sign < 0) return Tribool::False;   // -∞ is not nonnegative
            return Tribool::Unknown;
        }
        InferenceEngine engine(ctx_);
        return engine.query_nonnegative(expr);
    });
}

Tribool QueryInterface::query_real(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Real, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }

        // NaN is not a well-defined real value → Unknown
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }

        // Infinity is not a real number (it's extended real)
        if (is_infinity_node(expr.root)) {
            return Tribool::Unknown;
        }
        InferenceEngine engine(ctx_);
        return engine.query_real(expr);
    });
}

Tribool QueryInterface::query_integer(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Integer, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }

        // NaN: False for integer
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::False;
            }
        }

        // Infinity: False for integer
        if (is_infinity_node(expr.root)) {
            return Tribool::False;
        }
        InferenceEngine engine(ctx_);
        return engine.query_integer(expr);
    });
}

Tribool QueryInterface::query_nonzero(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::NonZero, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }

        // NaN: Unknown for sign queries
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }

        // Infinity is nonzero
        if (is_infinity_node(expr.root)) {
            return Tribool::True;
        }
        InferenceEngine engine(ctx_);
        return engine.query_nonzero(expr);
    });
}

// Extended property queries

Tribool QueryInterface::query_algebraic(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Algebraic, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }
        if (is_infinity_node(expr.root)) {
            return Tribool::False;
        }
        InferenceEngine engine(ctx_);
        return engine.query_algebraic(expr);
    });
}

Tribool QueryInterface::query_transcendental(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Transcendental, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }
        if (is_infinity_node(expr.root)) {
            return Tribool::False;
        }
        InferenceEngine engine(ctx_);
        return engine.query_transcendental(expr);
    });
}

Tribool QueryInterface::query_finite(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Finite, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }
        if (is_infinity_node(expr.root)) {
            return Tribool::False;
        }
        InferenceEngine engine(ctx_);
        return engine.query_finite(expr);
    });
}

Tribool QueryInterface::query_divergent(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Divergent, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::Unknown;
            }
        }
        if (is_infinity_node(expr.root)) {
            return Tribool::True;
        }
        InferenceEngine engine(ctx_);
        return engine.query_divergent(expr);
    });
}

Tribool QueryInterface::query_periodic(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::Periodic, [&]() -> Tribool {
        if (is_unhandled_type(expr.root)) {
            return Tribool::Unknown;
        }
        if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
            if (is_nan_number(*num)) {
                return Tribool::False;
            }
        }
        if (is_infinity_node(expr.root)) {
            return Tribool::False;
        }
        InferenceEngine engine(ctx_);
        return engine.query_periodic(expr);
    });
}

std::optional<SymbolicExpr> QueryInterface::get_period(const SymbolicExpr& expr) const {
    if (!expr.root) return std::nullopt;
    if (is_unhandled_type(expr.root)) return std::nullopt;
    if (is_infinity_node(expr.root)) return std::nullopt;

    InferenceEngine engine(ctx_);
    return engine.infer_period(expr);
}

Tribool QueryInterface::query_positive_definite(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::PositiveDefinite, [&]() -> Tribool {
        if (!expr.root) return Tribool::Unknown;
        if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
            const auto& props = ctx_.current_properties();
            Definiteness d = props.get_definiteness(var->name);
            if (d == Definiteness::PositiveDefinite) return Tribool::True;
            if (d == Definiteness::NegativeDefinite ||
                d == Definiteness::NegativeSemiDefinite ||
                d == Definiteness::Indefinite) return Tribool::False;
            return Tribool::Unknown;
        }
        return Tribool::Unknown;
    });
}

Tribool QueryInterface::query_positive_semidefinite(const SymbolicExpr& expr) const {
    return cached_query(expr, PropType::PositiveSemiDefinite, [&]() -> Tribool {
        if (!expr.root) return Tribool::Unknown;
        if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
            const auto& props = ctx_.current_properties();
            Definiteness d = props.get_definiteness(var->name);
            if (d == Definiteness::PositiveDefinite ||
                d == Definiteness::PositiveSemiDefinite) return Tribool::True;
            if (d == Definiteness::NegativeDefinite ||
                d == Definiteness::NegativeSemiDefinite ||
                d == Definiteness::Indefinite) return Tribool::False;
            return Tribool::Unknown;
        }
        return Tribool::Unknown;
    });
}

// Query conditions (query mode)

std::vector<QueryInterface::ConditionSet> QueryInterface::query_conditions(
    const SymbolicExpr& expr, Sign target) const {
    if (!expr.root) {
        return {};
    }
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        ConditionSet cs;
        cs.sign_conditions.emplace_back(var->name, target);
        return {cs};
    }
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        if (add->operands.size() == 2) {
            std::shared_ptr<SymbolicNode> pos_operand = nullptr;
            std::shared_ptr<SymbolicNode> neg_operand = nullptr; // the term being subtracted

            for (size_t i = 0; i < 2; ++i) {
                auto mul = std::dynamic_pointer_cast<MultiplyNode>(add->operands[i]);
                if (mul && mul->operands.size() == 2) {
                    for (size_t j = 0; j < 2; ++j) {
                        auto num = std::dynamic_pointer_cast<NumberNode>(mul->operands[j]);
                        if (num) {
                            bool is_neg_one = false;
                            if (std::holds_alternative<BigInt>(num->value)) {
                                is_neg_one = (std::get<BigInt>(num->value) == BigInt(-1));
                            } else if (std::holds_alternative<Rational>(num->value)) {
                                is_neg_one = (std::get<Rational>(num->value) == Rational(-1));
                            } else if (std::holds_alternative<lmmc_real_t>(num->value)) {
                                is_neg_one = (std::get<lmmc_real_t>(num->value) == -1.0);
                            }
                            if (is_neg_one) {
                                neg_operand = mul->operands[1 - j];
                                pos_operand = add->operands[1 - i];
                                break;
                            }
                        }
                    }
                }
                if (neg_operand) break;
            }
            if (pos_operand && neg_operand) {
                auto lhs_var = std::dynamic_pointer_cast<VariableNode>(pos_operand);
                auto rhs_var = std::dynamic_pointer_cast<VariableNode>(neg_operand);

                if (lhs_var && rhs_var && target == Sign::Positive) {
                    ConditionSet cs1;
                    cs1.sign_conditions.emplace_back(lhs_var->name, Sign::Positive);
                    cs1.sign_conditions.emplace_back(rhs_var->name, Sign::Negative);
                    ConditionSet cs2;
                    SymbolicExpr lhs_expr(pos_operand);
                    SymbolicExpr rhs_expr(neg_operand);
                    Relation rel{lhs_expr, rhs_expr, RelationalNode::Op::GT};
                    cs2.relational_conditions.push_back(rel);
                    cs2.sign_conditions.emplace_back(rhs_var->name, Sign::NonNegative);

                    return {cs1, cs2};
                } else if (lhs_var && rhs_var && target == Sign::Negative) {
                    ConditionSet cs1;
                    cs1.sign_conditions.emplace_back(lhs_var->name, Sign::Negative);
                    cs1.sign_conditions.emplace_back(rhs_var->name, Sign::Positive);

                    ConditionSet cs2;
                    SymbolicExpr lhs_expr(pos_operand);
                    SymbolicExpr rhs_expr(neg_operand);
                    Relation rel{rhs_expr, lhs_expr, RelationalNode::Op::GT};
                    cs2.relational_conditions.push_back(rel);
                    cs2.sign_conditions.emplace_back(lhs_var->name, Sign::NonNegative);

                    return {cs1, cs2};
                } else if (lhs_var && rhs_var && target == Sign::NonNegative) {
                    ConditionSet cs1;
                    cs1.sign_conditions.emplace_back(lhs_var->name, Sign::NonNegative);
                    cs1.sign_conditions.emplace_back(rhs_var->name, Sign::NonPositive);

                    ConditionSet cs2;
                    SymbolicExpr lhs_expr(pos_operand);
                    SymbolicExpr rhs_expr(neg_operand);
                    Relation rel{lhs_expr, rhs_expr, RelationalNode::Op::GEQ};
                    cs2.relational_conditions.push_back(rel);

                    return {cs1, cs2};
                }
            }
        }
        if (target == Sign::Positive || target == Sign::NonNegative) {
            bool all_variables = true;
            ConditionSet cs;
            for (const auto& op : add->operands) {
                if (auto var = std::dynamic_pointer_cast<VariableNode>(op)) {
                    cs.sign_conditions.emplace_back(var->name, target);
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
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        if (mul->operands.size() == 2) {
            std::shared_ptr<SymbolicNode> numerator = nullptr;
            std::shared_ptr<SymbolicNode> denominator = nullptr;

            for (size_t i = 0; i < 2; ++i) {
                auto pow = std::dynamic_pointer_cast<PowerNode>(mul->operands[i]);
                if (pow) {
                    auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent);
                    if (exp_num) {
                        bool is_neg_one = false;
                        if (std::holds_alternative<BigInt>(exp_num->value)) {
                            is_neg_one = (std::get<BigInt>(exp_num->value) == BigInt(-1));
                        } else if (std::holds_alternative<Rational>(exp_num->value)) {
                            is_neg_one = (std::get<Rational>(exp_num->value) == Rational(-1));
                        } else if (std::holds_alternative<lmmc_real_t>(exp_num->value)) {
                            is_neg_one = (std::get<lmmc_real_t>(exp_num->value) == -1.0);
                        }
                        if (is_neg_one) {
                            denominator = pow->base;
                            numerator = mul->operands[1 - i];
                            break;
                        }
                    }
                }
            }
            if (numerator && denominator) {
                auto num_var = std::dynamic_pointer_cast<VariableNode>(numerator);
                auto den_var = std::dynamic_pointer_cast<VariableNode>(denominator);

                if (num_var && den_var) {
                    if (target == Sign::Positive) {
                        // x/y > 0: both positive OR both negative
                        ConditionSet cs1;
                        cs1.sign_conditions.emplace_back(num_var->name, Sign::Positive);
                        cs1.sign_conditions.emplace_back(den_var->name, Sign::Positive);

                        ConditionSet cs2;
                        cs2.sign_conditions.emplace_back(num_var->name, Sign::Negative);
                        cs2.sign_conditions.emplace_back(den_var->name, Sign::Negative);

                        return {cs1, cs2};
                    } else if (target == Sign::Negative) {
                        // x/y < 0: one positive, one negative
                        ConditionSet cs1;
                        cs1.sign_conditions.emplace_back(num_var->name, Sign::Positive);
                        cs1.sign_conditions.emplace_back(den_var->name, Sign::Negative);

                        ConditionSet cs2;
                        cs2.sign_conditions.emplace_back(num_var->name, Sign::Negative);
                        cs2.sign_conditions.emplace_back(den_var->name, Sign::Positive);

                        return {cs1, cs2};
                    }
                }
                return {};
            }
        }
        if (target == Sign::Positive) {
            bool all_variables = true;
            ConditionSet cs;
            for (const auto& op : mul->operands) {
                if (auto var = std::dynamic_pointer_cast<VariableNode>(op)) {
                    cs.sign_conditions.emplace_back(var->name, Sign::Positive);
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

} // namespace lamina
