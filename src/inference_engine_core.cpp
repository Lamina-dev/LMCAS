#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"

namespace lamina {

InferenceEngine::InferenceEngine(const AssumptionContext& ctx)
    : impl_(std::make_unique<Impl>(ctx)) {}

InferenceEngine::~InferenceEngine() = default;
InferenceEngine::InferenceEngine(InferenceEngine&&) noexcept = default;
InferenceEngine& InferenceEngine::operator=(InferenceEngine&&) noexcept = default;

void InferenceEngine::set_max_depth(int depth) {
    if (depth > 0) {
        impl_->max_depth = depth;
    }
}

int InferenceEngine::get_max_depth() const {
    return impl_->max_depth;
}

InferenceEngine::DepthGuard::DepthGuard(const InferenceEngine& engine, const void* node)
    : engine_(engine), node_(node) {
    engine_.impl_->current_depth++;

    if (engine_.impl_->current_depth > engine_.impl_->max_depth) {
        abort_ = true;
        return;
    }

    // Insert and detect cycles atomically. If allocation throws, restore the
    // depth counter because a throwing constructor has no matching destructor.
    if (node_) {
        try {
            inserted_ = engine_.impl_->visited.insert(node_).second;
        } catch (...) {
            engine_.impl_->current_depth--;
            if (engine_.impl_->current_depth == 0) {
                engine_.impl_->visited.clear();
            }
            throw;
        }
        if (!inserted_) {
            abort_ = true;
        }
    }
}

InferenceEngine::DepthGuard::~DepthGuard() {
    if (inserted_ && node_) {
        engine_.impl_->visited.erase(node_);
    }

    engine_.impl_->current_depth--;

    if (engine_.impl_->current_depth == 0) {
        engine_.impl_->visited.clear();
    }
}


InferenceTriboolResult InferenceEngine::query_sign_of_checked(const SymbolicExpr& expr, Sign sign) const {
    switch (sign) {
        case Sign::Positive:    return query_positive_checked(expr);
        case Sign::Negative:    return query_negative_checked(expr);
        case Sign::NonNegative: return query_nonnegative_checked(expr);
        case Sign::NonPositive: return query_nonpositive_checked(expr);
        case Sign::Zero:
            // Zero means both NonNegative and NonPositive
            {
                auto nn = query_nonnegative_checked(expr);
                if (!nn) return nn;
                auto np = query_nonpositive_checked(expr);
                if (!np) return np;
                if (nn.value() == Tribool::True && np.value() == Tribool::True)
                    return InferenceTriboolResult::success(Tribool::True);
                if (nn.value() == Tribool::False || np.value() == Tribool::False)
                    return InferenceTriboolResult::success(Tribool::False);
                return InferenceTriboolResult::success(Tribool::Unknown);
            }
        case Sign::NonZero:     return query_nonzero_checked(expr);
    }
    return InferenceTriboolResult::success(Tribool::Unknown);
}


InferenceTriboolResult InferenceEngine::query_domain_of_checked(const SymbolicExpr& expr, Domain domain) const {
    switch (domain) {
        case Domain::Integer:  return query_integer_checked(expr);
        case Domain::Real:     return query_real_checked(expr);
        case Domain::Rational: {
            if (!lamina::detail::node(expr)) {
                return InferenceTriboolResult::failure(
                    CasErrc::InvalidArgument,
                    "inference expression must not be null",
                    "query_domain_of_checked");
            }
            try {
                // Rational: Integer subset Rational, so Integer implies Rational
                auto integer = query_integer_checked(expr);
                if (!integer) return InferenceTriboolResult::failure(integer.error());
                if (integer.value() == Tribool::True) {
                    return InferenceTriboolResult::success(Tribool::True);
                }
                // Check if it's a Rational number literal
                if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                    if (std::holds_alternative<Rational>(num->value())) {
                        return InferenceTriboolResult::success(Tribool::True);
                    }
                }
                // Check variable domain
                if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                    const auto& props = impl_->ctx.current_properties();
                    if (props.has_domain(var->name(), Domain::Rational)) {
                        return InferenceTriboolResult::success(Tribool::True);
                    }
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            } catch (const std::bad_alloc&) {
                return InferenceTriboolResult::failure(
                    CasErrc::ResourceLimit,
                    "inference query allocation failed",
                    "query_domain_of_checked");
            } catch (const std::exception& ex) {
                return InferenceTriboolResult::failure(
                    CasErrc::InternalInvariant,
                    ex.what(),
                    "query_domain_of_checked");
            }
        }
        case Domain::Natural: {
            return checked_inference_result<Tribool>(expr, "query_domain_of_checked",
                [&]() {
                    // Natural: non-negative integers (0, 1, 2, ...)
                    // Check if it's a non-negative integer literal
                    if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                        if (is_integer_number(*num)) {
                            if (std::holds_alternative<BigInt>(num->value())) {
                                const auto& b = std::get<BigInt>(num->value());
                                if (!b.IsNegative()) return Tribool::True;
                            } else if (std::holds_alternative<Rational>(num->value())) {
                                BigInt n = std::get<Rational>(num->value()).get_numerator();
                                if (!n.IsNegative()) return Tribool::True;
                            } else {
                                double v = std::get<lmmc_real_t>(num->value());
                                if (std::isfinite(v) && v >= 0.0 && v == std::floor(v))
                                    return Tribool::True;
                            }
                        }
                        return Tribool::False;
                    }
                    // Check variable domain
                    if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                        const auto& props = impl_->ctx.current_properties();
                        if (props.has_domain(var->name(), Domain::Natural)) return Tribool::True;
                    }
                    return Tribool::Unknown;
                });
        }
        default:
            return InferenceTriboolResult::success(Tribool::Unknown);
    }
}

// Public query methods — dispatch based on node type

// These public methods are called by the QueryInterface for composite nodes.
// They inspect the root node type and dispatch to the appropriate inference method.


InferenceTriboolResult InferenceEngine::query_positive_checked(const SymbolicExpr& expr) const {
    const int infinity = infinity_sign(expr);
    if (infinity != 0) {
        return InferenceTriboolResult::success(
            infinity > 0 ? Tribool::True : Tribool::False);
    }
    return checked_inference_result<Tribool>(expr, "query_positive_checked",
        [&]() {
            // Depth guard: detect cycles and enforce depth limit
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // Handle NumberNode: determine sign directly from numeric value
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (num->is_positive()) return Tribool::True;
                if (num->is_zero()) return Tribool::False;
                // Check if negative; if not negative and not zero, it must be positive.
                if (std::holds_alternative<BigInt>(num->value())) {
                    return std::get<BigInt>(num->value()).IsNegative() ? Tribool::False : Tribool::True;
                }
                if (std::holds_alternative<Rational>(num->value())) {
                    return std::get<Rational>(num->value()) < Rational(0) ? Tribool::False : Tribool::True;
                }
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    if (!std::isfinite(v)) return Tribool::Unknown;
                    if (v > 0.0) return Tribool::True;
                    if (v < 0.0) return Tribool::False;
                    return Tribool::False; // zero
                }
                return Tribool::Unknown;
            }

            // Handle VariableNode: check PropertyStore in the AssumptionContext
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.has_sign(var->name(), Sign::Positive)) return Tribool::True;
                // If NonPositive or Negative or Zero, then not positive
                if (props.has_sign(var->name(), Sign::Negative) ||
                    props.has_sign(var->name(), Sign::Zero) ||
                    props.has_sign(var->name(), Sign::NonPositive)) return Tribool::False;
                return Tribool::Unknown;
            }

            if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
                auto result = detail::propagate_result(infer_add_sign_checked(add.get(), Sign::Positive));
                if (result != Tribool::Unknown) return result;
                return detail::propagate_result(infer_sign_from_relations_checked(expr, Sign::Positive));
            }
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
                auto result = detail::propagate_result(infer_multiply_sign_checked(mul.get(), Sign::Positive));
                if (result != Tribool::Unknown) return result;
                return detail::propagate_result(infer_sign_from_relations_checked(expr, Sign::Positive));
            }
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_power_property_checked(pow.get(), Sign::Positive));
            }
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_function_property_checked(func.get(), Sign::Positive));
            }
            // For other expression types (e.g., variables handled above), check relations
            return detail::propagate_result(infer_sign_from_relations_checked(expr, Sign::Positive));
        });
}


InferenceTriboolResult InferenceEngine::query_negative_checked(const SymbolicExpr& expr) const {
    const int infinity = infinity_sign(expr);
    if (infinity != 0) {
        return InferenceTriboolResult::success(
            infinity < 0 ? Tribool::True : Tribool::False);
    }
    return checked_inference_result<Tribool>(expr, "query_negative_checked",
        [&]() {
            // Depth guard: detect cycles and enforce depth limit
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // Handle NumberNode
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (num->is_zero()) return Tribool::False;
                if (std::holds_alternative<BigInt>(num->value())) {
                    return std::get<BigInt>(num->value()).IsNegative() ? Tribool::True : Tribool::False;
                }
                if (std::holds_alternative<Rational>(num->value())) {
                    const auto& r = std::get<Rational>(num->value());
                    if (r < Rational(0)) return Tribool::True;
                    return Tribool::False;
                }
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    if (!std::isfinite(v)) return Tribool::Unknown;
                    if (v < 0.0) return Tribool::True;
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }

            // Handle VariableNode
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.has_sign(var->name(), Sign::Negative)) return Tribool::True;
                if (props.has_sign(var->name(), Sign::Positive) ||
                    props.has_sign(var->name(), Sign::Zero) ||
                    props.has_sign(var->name(), Sign::NonNegative)) return Tribool::False;
                return Tribool::Unknown;
            }

            if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_add_sign_checked(add.get(), Sign::Negative));
            }
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_multiply_sign_checked(mul.get(), Sign::Negative));
            }
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_power_property_checked(pow.get(), Sign::Negative));
            }
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_function_property_checked(func.get(), Sign::Negative));
            }
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_nonnegative_checked(const SymbolicExpr& expr) const {
    const int infinity = infinity_sign(expr);
    if (infinity != 0) {
        return InferenceTriboolResult::success(
            infinity > 0 ? Tribool::True : Tribool::False);
    }
    return checked_inference_result<Tribool>(expr, "query_nonnegative_checked",
        [&]() {
            // Depth guard: detect cycles and enforce depth limit
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // Handle NumberNode
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (num->is_zero()) return Tribool::True;
                if (num->is_positive()) return Tribool::True;
                if (std::holds_alternative<BigInt>(num->value())) {
                    return std::get<BigInt>(num->value()).IsNegative() ? Tribool::False : Tribool::True;
                }
                if (std::holds_alternative<Rational>(num->value())) {
                    return std::get<Rational>(num->value()) < Rational(0) ? Tribool::False : Tribool::True;
                }
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    if (!std::isfinite(v)) return Tribool::Unknown;
                    return v >= 0.0 ? Tribool::True : Tribool::False;
                }
                return Tribool::Unknown;
            }

            // Handle VariableNode
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.has_sign(var->name(), Sign::NonNegative)) return Tribool::True;
                if (props.has_sign(var->name(), Sign::Negative)) return Tribool::False;
                return Tribool::Unknown;
            }

            if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
                auto result = detail::propagate_result(infer_add_sign_checked(add.get(), Sign::NonNegative));
                if (result != Tribool::Unknown) return result;
                return detail::propagate_result(infer_sign_from_relations_checked(expr, Sign::NonNegative));
            }
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
                auto result = detail::propagate_result(infer_multiply_sign_checked(mul.get(), Sign::NonNegative));
                if (result != Tribool::Unknown) return result;
                return detail::propagate_result(infer_sign_from_relations_checked(expr, Sign::NonNegative));
            }
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_power_property_checked(pow.get(), Sign::NonNegative));
            }
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_function_property_checked(func.get(), Sign::NonNegative));
            }
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_nonpositive_checked(const SymbolicExpr& expr) const {
    return checked_inference_result<Tribool>(expr, "query_nonpositive_checked",
        [&]() {
            // Depth guard: detect cycles and enforce depth limit
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // Handle NumberNode
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (num->is_zero()) return Tribool::True;
                if (num->is_positive()) return Tribool::False;
                if (std::holds_alternative<BigInt>(num->value())) {
                    return std::get<BigInt>(num->value()).IsNegative() ? Tribool::True : Tribool::False;
                }
                if (std::holds_alternative<Rational>(num->value())) {
                    return std::get<Rational>(num->value()) > Rational(0) ? Tribool::False : Tribool::True;
                }
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    if (!std::isfinite(v)) return Tribool::Unknown;
                    return v <= 0.0 ? Tribool::True : Tribool::False;
                }
                return Tribool::Unknown;
            }

            // Handle VariableNode
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.has_sign(var->name(), Sign::NonPositive)) return Tribool::True;
                if (props.has_sign(var->name(), Sign::Positive)) return Tribool::False;
                return Tribool::Unknown;
            }

            if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_add_sign_checked(add.get(), Sign::NonPositive));
            }
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_multiply_sign_checked(mul.get(), Sign::NonPositive));
            }
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_power_property_checked(pow.get(), Sign::NonPositive));
            }
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_function_property_checked(func.get(), Sign::NonPositive));
            }
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_real_checked(const SymbolicExpr& expr) const {
    return checked_inference_result<Tribool>(expr, "query_real_checked",
        [&]() {
            // Depth guard: detect cycles and enforce depth limit
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // Handle NumberNode: finite numbers are Real
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (std::holds_alternative<BigInt>(num->value())) return Tribool::True;
                if (std::holds_alternative<Rational>(num->value())) return Tribool::True;
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    return std::isfinite(v) ? Tribool::True : Tribool::Unknown;
                }
                return Tribool::Unknown;
            }

            // Handle VariableNode
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.has_domain(var->name(), Domain::Real)) return Tribool::True;
                return Tribool::Unknown;
            }

            if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_add_domain_checked(add.get(), Domain::Real));
            }
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_multiply_domain_checked(mul.get(), Domain::Real));
            }
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_power_domain_checked(pow.get(), Domain::Real));
            }
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_function_domain_checked(func.get(), Domain::Real));
            }
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_integer_checked(const SymbolicExpr& expr) const {
    if (infinity_sign(expr) != 0) {
        return InferenceTriboolResult::success(Tribool::False);
    }
    return checked_inference_result<Tribool>(expr, "query_integer_checked",
        [&]() {
            // Depth guard: detect cycles and enforce depth limit
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // Handle NumberNode: BigInt values are Integer
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (std::holds_alternative<BigInt>(num->value())) return Tribool::True;
                if (std::holds_alternative<Rational>(num->value())) {
                    const auto& r = std::get<Rational>(num->value());
                    // Integer if denominator is 1
                    return r.is_integer() ? Tribool::True : Tribool::False;
                }
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    if (!std::isfinite(v)) return Tribool::False;
                    return (v == std::floor(v)) ? Tribool::True : Tribool::False;
                }
                return Tribool::Unknown;
            }

            // Handle VariableNode
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.has_domain(var->name(), Domain::Integer)) return Tribool::True;
                return Tribool::Unknown;
            }

            if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_add_domain_checked(add.get(), Domain::Integer));
            }
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_multiply_domain_checked(mul.get(), Domain::Integer));
            }
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_power_domain_checked(pow.get(), Domain::Integer));
            }
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_function_domain_checked(func.get(), Domain::Integer));
            }
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_nonzero_checked(const SymbolicExpr& expr) const {
    if (infinity_sign(expr) != 0) {
        return InferenceTriboolResult::success(Tribool::True);
    }
    return checked_inference_result<Tribool>(expr, "query_nonzero_checked",
        [&]() {
            // Depth guard: detect cycles and enforce depth limit
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // Handle NumberNode
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (num->is_zero()) return Tribool::False;
                if (num->is_positive()) return Tribool::True;
                // Check if it's a non-zero number
                if (std::holds_alternative<BigInt>(num->value())) {
                    return (std::get<BigInt>(num->value()) == BigInt(0)) ? Tribool::False : Tribool::True;
                }
                if (std::holds_alternative<Rational>(num->value())) {
                    return (std::get<Rational>(num->value()) == Rational(0)) ? Tribool::False : Tribool::True;
                }
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    if (!std::isfinite(v)) return Tribool::Unknown;
                    return (v != 0.0) ? Tribool::True : Tribool::False;
                }
                return Tribool::Unknown;
            }

            // Handle VariableNode
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.has_sign(var->name(), Sign::NonZero)) return Tribool::True;
                if (props.has_sign(var->name(), Sign::Zero)) return Tribool::False;
                return Tribool::Unknown;
            }

            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_multiply_sign_checked(mul.get(), Sign::NonZero));
            }
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_power_property_checked(pow.get(), Sign::NonZero));
            }
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                return detail::propagate_result(infer_function_property_checked(func.get(), Sign::NonZero));
            }
            // For addition, NonZero is hard to determine in general
            // (positive + positive = nonzero, but that's covered by positive inference)
            if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
                // If the sum is positive or negative, it's nonzero
                auto pos = detail::propagate_result(infer_add_sign_checked(add.get(), Sign::Positive));
                if (pos == Tribool::True) return Tribool::True;
                auto neg = detail::propagate_result(infer_add_sign_checked(add.get(), Sign::Negative));
                if (neg == Tribool::True) return Tribool::True;
                return Tribool::Unknown;
            }
            return Tribool::Unknown;
        });
}

// Subtraction sign inference


InferenceTriboolResult InferenceEngine::infer_subtraction_sign_checked(const void* opaque_node, Sign target) const {
    if (!opaque_node) {
        return InferenceTriboolResult::failure(
            CasErrc::InvalidArgument,
            "subtraction inference node must not be null",
            "infer_subtraction_sign_checked");
    }

    try {
        const auto& node = *static_cast<const AddNode*>(opaque_node);
        if (node.operands().size() != 2) return InferenceTriboolResult::success(Tribool::Unknown);

        std::shared_ptr<const SymbolicNode> minuend_node;
        std::shared_ptr<const SymbolicNode> subtrahend_node;

        for (const auto& operand : node.operands()) {
            auto mul = std::dynamic_pointer_cast<const MultiplyNode>(operand);
            if (mul && mul->operands().size() == 2) {
                bool found_neg_one = false;
                std::shared_ptr<const SymbolicNode> other_operand;
                for (const auto& mul_op : mul->operands()) {
                    auto num = std::dynamic_pointer_cast<const NumberNode>(mul_op);
                    if (num) {
                        bool is_neg_one = false;
                        if (std::holds_alternative<BigInt>(num->value())) {
                            is_neg_one = (std::get<BigInt>(num->value()) == BigInt(-1));
                        } else if (std::holds_alternative<Rational>(num->value())) {
                            is_neg_one = (std::get<Rational>(num->value()) == Rational(-1));
                        } else if (std::holds_alternative<lmmc_real_t>(num->value())) {
                            lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                            is_neg_one = (std::isfinite(v) && v == -1.0);
                        }
                        if (is_neg_one) {
                            found_neg_one = true;
                            continue;
                        }
                    }
                    other_operand = mul_op;
                }
                if (found_neg_one && other_operand) {
                    subtrahend_node = other_operand;
                } else {
                    minuend_node = operand;
                }
            } else {
                minuend_node = operand;
            }
        }

        if (!minuend_node || !subtrahend_node) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }

        auto minuend_expr = lamina::detail::expression_from_node(minuend_node);
        auto subtrahend_expr = lamina::detail::expression_from_node(subtrahend_node);

        Tribool min_pos = detail::propagate_result(query_sign_of_checked(minuend_expr, Sign::Positive));
        Tribool min_neg = detail::propagate_result(query_sign_of_checked(minuend_expr, Sign::Negative));
        Tribool sub_pos = detail::propagate_result(query_sign_of_checked(subtrahend_expr, Sign::Positive));
        Tribool sub_neg = detail::propagate_result(query_sign_of_checked(subtrahend_expr, Sign::Negative));

        if (min_pos == Tribool::True && sub_neg == Tribool::True) {
            switch (target) {
                case Sign::Positive:    return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::True);
                case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::False);
                case Sign::Zero:        return InferenceTriboolResult::success(Tribool::False);
            }
        }

        if (min_neg == Tribool::True && sub_pos == Tribool::True) {
            switch (target) {
                case Sign::Negative:    return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::True);
                case Sign::Positive:    return InferenceTriboolResult::success(Tribool::False);
                case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::False);
                case Sign::Zero:        return InferenceTriboolResult::success(Tribool::False);
            }
        }

        Tribool min_nn = detail::propagate_result(query_sign_of_checked(minuend_expr, Sign::NonNegative));
        Tribool sub_np = detail::propagate_result(query_sign_of_checked(subtrahend_expr, Sign::NonPositive));
        if (min_nn == Tribool::True && sub_np == Tribool::True) {
            switch (target) {
                case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                default: break;
            }
        }

        Tribool min_np = detail::propagate_result(query_sign_of_checked(minuend_expr, Sign::NonPositive));
        Tribool sub_nn = detail::propagate_result(query_sign_of_checked(subtrahend_expr, Sign::NonNegative));
        if (min_np == Tribool::True && sub_nn == Tribool::True) {
            switch (target) {
                case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::True);
                case Sign::Positive:    return InferenceTriboolResult::success(Tribool::False);
                default: break;
            }
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const detail::ResultPropagation& ex) {
        return InferenceTriboolResult::failure(ex.error());
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "subtraction sign inference allocation failed",
            "infer_subtraction_sign_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_subtraction_sign_checked");
    }
}

// Algebraic / Transcendental / Finite / Divergent queries


InferenceTriboolResult InferenceEngine::query_algebraic_checked(const SymbolicExpr& expr) const {
    return checked_inference_result<Tribool>(expr, "query_algebraic_checked",
        [&]() {
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // NumberNode: BigInt and Rational are algebraic; finite doubles are Unknown
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (std::holds_alternative<BigInt>(num->value())) return Tribool::True;
                if (std::holds_alternative<Rational>(num->value())) return Tribool::True;
                return Tribool::Unknown;
            }

            // VariableNode: check if domain is Algebraic or more specific
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                // Domain hierarchy: Algebraic ⊃ Rational ⊃ Integer ⊃ Natural ⊃ PositiveInt
                if (props.has_domain(var->name(), Domain::Algebraic)) return Tribool::True;
                // If transcendental, definitely not algebraic
                if (props.is_transcendental(var->name())) return Tribool::False;
                return Tribool::Unknown;
            }

            // Composite expressions: return Unknown for now
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_transcendental_checked(const SymbolicExpr& expr) const {
    return checked_inference_result<Tribool>(expr, "query_transcendental_checked",
        [&]() {
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // NumberNode: integers and rationals are algebraic, not transcendental
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (std::holds_alternative<BigInt>(num->value())) return Tribool::False;
                if (std::holds_alternative<Rational>(num->value())) return Tribool::False;
                return Tribool::Unknown;
            }

            // VariableNode: check the transcendental flag
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.is_transcendental(var->name())) return Tribool::True;
                // If domain is Algebraic or more specific, not transcendental
                if (props.has_domain(var->name(), Domain::Algebraic)) return Tribool::False;
                return Tribool::Unknown;
            }

            // Composite expressions: return Unknown for now
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_finite_checked(const SymbolicExpr& expr) const {
    return checked_inference_result<Tribool>(expr, "query_finite_checked",
        [&]() {
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // NumberNode: finite numeric values are Finite
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (std::holds_alternative<BigInt>(num->value())) return Tribool::True;
                if (std::holds_alternative<Rational>(num->value())) return Tribool::True;
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    return std::isfinite(v) ? Tribool::True : Tribool::False;
                }
                return Tribool::Unknown;
            }

            // VariableNode: check Finiteness in PropertyStore
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                Finiteness f = props.get_finiteness(var->name());
                if (f == Finiteness::Finite) return Tribool::True;
                if (f == Finiteness::Divergent) return Tribool::False;
                return Tribool::Unknown;
            }

            // Composite expressions: return Unknown for now
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_divergent_checked(const SymbolicExpr& expr) const {
    return checked_inference_result<Tribool>(expr, "query_divergent_checked",
        [&]() {
            DepthGuard guard(*this, lamina::detail::node(expr).get());
            if (guard.should_abort()) return Tribool::Unknown;

            // NumberNode: finite numeric values are not divergent
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
                if (std::holds_alternative<BigInt>(num->value())) return Tribool::False;
                if (std::holds_alternative<Rational>(num->value())) return Tribool::False;
                if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    return std::isfinite(v) ? Tribool::False : Tribool::Unknown;
                }
                return Tribool::Unknown;
            }

            // VariableNode: check Finiteness in PropertyStore
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                Finiteness f = props.get_finiteness(var->name());
                if (f == Finiteness::Divergent) return Tribool::True;
                if (f == Finiteness::Finite) return Tribool::False;
                return Tribool::Unknown;
            }

            // Composite expressions: return Unknown for now
            return Tribool::Unknown;
        });
}

// Addition sign inference


} // namespace lamina
