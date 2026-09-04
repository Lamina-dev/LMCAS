#include "internal/integration_support.hpp"

namespace lamina {

std::shared_ptr<SymbolicExpr> TableLookupStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int) {
    (void)computation;

    auto all_entries = ctx.table().get_all_sorted();

    for (const auto* entry : all_entries) {
        MatchMap bindings;
        if (Matcher::match(entry->pattern, expr, entry->wildcards, bindings)) {

            // The matcher allows partial matches on commutative operations
            // (Add/Multiply): unmatched operands are bound to __Add_REST__ /
            // __Mul_REST__ for use by rewrite rules. For integration table
            // lookup that behavior is wrong - a pattern like (1 + _u^2)^-1
            // would otherwise match (1 + x^2 + x + x^3)^-1 with _u=x and the
            // rest x + x^3 silently dropped, producing arctan(x) for an
            // integrand that has nothing to do with arctan. Require an exact
            // match (no leftover) so the result genuinely equals the
            // integrand under this binding.
            if (bindings.find("__Add_REST__") != bindings.end() ||
                bindings.find("__Mul_REST__") != bindings.end()) {
                continue;
            }
            const bool complete_binding = std::all_of(
                entry->wildcards.begin(), entry->wildcards.end(),
                [&](const std::string& wildcard_name) {
                    return bindings.find(wildcard_name) != bindings.end();
                });
            if (!complete_binding) continue;

            if (entry->condition && !entry->condition(bindings, var)) {
                continue;
            }

            SymbolicExpr result = Matcher::replace(entry->result, bindings, false);
            auto simplified = result.simplify();
            if (simplified && !contains_unevaluated_integral(lamina::detail::node(simplified))) {
                return simplified;
            }
        }
    }
    return nullptr;
}

std::shared_ptr<SymbolicExpr> PowerRuleStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext&, int) {

    if (auto v_node = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
        if (v_node->name() == var) {
            return SymbolicExpr::multiply(
                SymbolicExpr::power(make_expr_ptr(expr), SymbolicExpr::number(2)),
                sym_rational(1, 2));
        }
    }

    if (auto p_node = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
        auto base = lamina::detail::expression_from_node(p_node->base());
        auto exp_expr = lamina::detail::expression_from_node(p_node->exponent());
        if (auto b_var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(base))) {
            if (b_var->name() == var && !depends_on_integration_variable(exp_expr, var)) {
                auto n_plus_1 = SymbolicExpr::add(make_expr_ptr(exp_expr), SymbolicExpr::number(1))->simplify();
                if (n_plus_1->is_zero()) {

                    return SymbolicExpr::ln(make_expr_ptr(base));
                }
                return SymbolicExpr::divide(
                    SymbolicExpr::power(make_expr_ptr(base), n_plus_1),
                    n_plus_1);
            }
        }
    }
    return nullptr;
}

std::shared_ptr<SymbolicExpr> SubstitutionStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext& computation, int) {

    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        ops = mul->operands();
    } else {
        ops.push_back(lamina::detail::node(expr));
    }

    for (size_t i = 0; i < ops.size(); ++i) {
        auto candidate_term = lamina::detail::expression_from_node(ops[i]);
        std::optional<SymbolicExpr> u;

        if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(candidate_term))) {
            u = lamina::detail::expression_from_node(pow->base());
        } else if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(candidate_term))) {
            if (!func->arguments().empty()) {
                u = lamina::detail::expression_from_node(func->arguments()[0]);
            }
        }

        if (u && depends_on_integration_variable(*u, var)) {
            auto d_ptr = u->differentiate(var);
            if (!d_ptr) continue;
            auto du = d_ptr->simplify();
            if (du->is_zero()) continue;

            auto f_u = candidate_term;
            auto term_times_du = SymbolicExpr::multiply(make_expr_ptr(f_u), make_expr_ptr(*du));
            auto ratio = SymbolicExpr::divide(make_expr_ptr(expr), term_times_du)->simplify();
            bool ratio_independent = !depends_on_integration_variable(*ratio, var);
            if (!ratio_independent && computation.assumptions()) {
                Tribool du_nonzero = detail::propagate_result(
                    computation.assumptions()->is_nonzero(*du));
                if (du_nonzero == Tribool::True) {
                    ratio_independent = !depends_on_integration_variable(*ratio, var);
                }
            }

            if (ratio_independent) {
                std::shared_ptr<SymbolicExpr> prim = nullptr;

                if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(candidate_term))) {
                    auto n = lamina::detail::expression_from_node(pow->exponent());
                    auto np1 = SymbolicExpr::add(make_expr_ptr(n), SymbolicExpr::number(1))->simplify();
                    if (np1->is_zero()) {
                        prim = SymbolicExpr::ln(make_expr_ptr(*u));
                    } else {
                        prim = SymbolicExpr::divide(
                            SymbolicExpr::power(make_expr_ptr(*u), np1), np1);
                    }
                } else if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(candidate_term))) {
                    if (func->type() == FunctionNode::FuncType::Cos) {
                        prim = SymbolicExpr::sin(make_expr_ptr(*u));
                    } else if (func->type() == FunctionNode::FuncType::Sin) {
                        prim = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(make_expr_ptr(*u)));
                    } else if (func->type() == FunctionNode::FuncType::Exp) {
                        prim = SymbolicExpr::exp(make_expr_ptr(*u));
                    }
                }

                if (prim) {
                    return SymbolicExpr::multiply(ratio, prim);
                }
            }
        }
    }
    return nullptr;
}


bool LinearSubstitutionStrategy::extract_linear_arg(
    const SymbolicExpr& arg,
    const std::string& var,
    std::shared_ptr<SymbolicExpr>& a_out,
    std::shared_ptr<SymbolicExpr>& b_out) {

    // Argument must actually depend on the integration variable.
    if (!expression_depends_on_variable(lamina::detail::node(arg), var)) return false;

    Polynomial<SymbolicPolyCoeff> poly;
    try {
        poly = symbolic_to_poly<SymbolicPolyCoeff>(make_expr_ptr(arg), var);
    } catch (const std::invalid_argument&) {
        return false;
    } catch (const std::out_of_range&) {
        return false;
    } catch (const std::runtime_error&) {
        return false;
    }

    // Need to be exactly degree 1 in var: arg = a * var + b with a, b constant.
    if (poly.degree() != 1) return false;
    if (poly.coeffs.size() < 2) return false;

    auto a_expr = poly.coeffs[1].val;
    auto b_expr = poly.coeffs[0].val;
    if (!a_expr || !b_expr) return false;

    // Coefficients must not depend on the integration variable.
    if (expression_depends_on_variable(lamina::detail::node(a_expr), var)) return false;
    if (expression_depends_on_variable(lamina::detail::node(b_expr), var)) return false;

    auto a_simp = a_expr->simplify();
    auto b_simp = b_expr->simplify();
    if (!a_simp) a_simp = a_expr;
    if (!b_simp) b_simp = b_expr;

    // a must be demonstrably non-zero (refuse if simplification yields zero).
    if (a_simp->is_zero()) return false;

    a_out = a_simp;
    b_out = b_simp;
    return true;
}

std::shared_ptr<SymbolicExpr> LinearSubstitutionStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int depth) {

    using FT = FunctionNode::FuncType;

    enum class Wrapper { None, Function, PowerOfFunction, PowerOfLinear };
    Wrapper kind = Wrapper::None;

    std::shared_ptr<const FunctionNode> fn;       // outer FunctionNode (Function case)
    std::shared_ptr<const PowerNode>    pn;       // outer PowerNode (Power* cases)
    std::shared_ptr<const FunctionNode> base_fn;  // inner FunctionNode of a PowerNode
    std::optional<SymbolicExpr> arg_expr;   // the linear-candidate sub-expression

    if ((fn = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr)))) {
        // Skip wrappers that are not "real" single-argument functions.
        switch (fn->type()) {
            case FT::Infinity:
            case FT::Atan2:
            case FT::Log:
                return nullptr;
            default:
                break;
        }
        if (fn->arguments().size() != 1) return nullptr;
        arg_expr = lamina::detail::expression_from_node(fn->arguments()[0]);
        kind = Wrapper::Function;
    } else if ((pn = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr)))) {
        // Exponent must be independent of the integration variable so that
        // substituting the linear argument back is sound.
        auto exp_expr = lamina::detail::expression_from_node(pn->exponent());
        if (expression_depends_on_variable(lamina::detail::node(exp_expr), var)) return nullptr;

        if ((base_fn = std::dynamic_pointer_cast<const FunctionNode>(pn->base()))) {
            switch (base_fn->type()) {
                case FT::Infinity:
                case FT::Atan2:
                case FT::Log:
                    return nullptr;
                default:
                    break;
            }
            if (base_fn->arguments().size() != 1) return nullptr;
            arg_expr = lamina::detail::expression_from_node(base_fn->arguments()[0]);
            kind = Wrapper::PowerOfFunction;
        } else {
            // Treat the whole base as the linear-candidate (e.g. (2x+1)^n).
            arg_expr = lamina::detail::expression_from_node(pn->base());
            kind = Wrapper::PowerOfLinear;
        }
    } else {
        return nullptr;
    }

    std::shared_ptr<SymbolicExpr> a_coeff, b_coeff;
    if (!arg_expr || !extract_linear_arg(*arg_expr, var, a_coeff, b_coeff)) return nullptr;

    // Degenerate case: arg == var (a = 1, b = 0). The plain TableLookup /
    // PowerRule strategies would have handled this earlier in the chain;
    // returning nullptr lets the next strategy try.
    bool a_is_one  = a_coeff && a_coeff->is_one();
    bool b_is_zero = b_coeff && b_coeff->is_zero();
    if (a_is_one && b_is_zero) return nullptr;

    const std::string dummy_name = "__lin_sub_u__";
    auto dummy_var = SymbolicExpr::variable(dummy_name);

    std::shared_ptr<SymbolicExpr> test_expr;
    if (kind == Wrapper::Function) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args{lamina::detail::node(dummy_var)};
        test_expr = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(fn->type(), new_args));
    } else if (kind == Wrapper::PowerOfFunction) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args{lamina::detail::node(dummy_var)};
        auto new_fn = lamina::detail::make_node<FunctionNode>(base_fn->type(), new_args);
        test_expr = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<PowerNode>(new_fn, pn->exponent()));
    } else { // Wrapper::PowerOfLinear
        test_expr = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<PowerNode>(lamina::detail::node(dummy_var), pn->exponent()));
    }
    if (!test_expr) return nullptr;

    TableLookupStrategy table_only;
    auto table_attempt = table_only.try_integrate(
        *test_expr, dummy_name, ctx, computation, depth);
    if (!table_attempt) {
        throw detail::ResultPropagation(table_attempt.error());
    }
    auto* table_candidate =
        std::get_if<IntegrationCandidate>(&table_attempt.value());
    if (!table_candidate || !table_candidate->expression) return nullptr;
    auto F_substituted = table_candidate->expression->substitute(
        dummy_name, make_expr_ptr(*arg_expr));
    if (!F_substituted) return nullptr;

    auto inv_a = SymbolicExpr::power(a_coeff, SymbolicExpr::number(-1));
    auto result = SymbolicExpr::multiply(inv_a, F_substituted);

    auto simplified = result->simplify();
    if (simplified && !contains_unevaluated_integral(lamina::detail::node(simplified))) {
        return simplified;
    }
    return result;
}

std::shared_ptr<SymbolicExpr> PartialFractionStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext& computation, int) {

    std::shared_ptr<SymbolicExpr> den = nullptr;

    if (auto p = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
        double exp_val = 0;
        bool is_inv = false;
        if (auto num_node = std::dynamic_pointer_cast<const NumberNode>(p->exponent())) {
            if (std::holds_alternative<lmmc_real_t>(num_node->value()))
                exp_val = std::get<lmmc_real_t>(num_node->value());
            else if (std::holds_alternative<Rational>(num_node->value()))
                exp_val = std::get<Rational>(num_node->value()).to_double();
            else if (std::holds_alternative<BigInt>(num_node->value()))
                exp_val = std::get<BigInt>(num_node->value()).to_double();
            int eq;
            lmmc_double_nearly_equal_tol(exp_val, -1.0, 1e-9, 1e-9, &eq);
            if (eq) is_inv = true;
        }
        if (is_inv) {
            den = make_expr_ptr(lamina::detail::expression_from_node(p->base()));
        }
    }

    if (!den) return nullptr;

    try {
        Polynomial<SymbolicPolyCoeff> Q = symbolic_to_poly<SymbolicPolyCoeff>(den, var);

        if (Q.degree() == 2) {
            SymbolicExpr c_expr = *(Q.coeffs[0].val);
            SymbolicExpr b_expr = *(Q.coeffs[1].val);
            SymbolicExpr a_expr = *(Q.coeffs[2].val);

            if (!a_expr.is_number() || !b_expr.is_number() || !c_expr.is_number()) {
                return nullptr;
            }

            auto a_checked = try_checked_numeric_constant(a_expr, computation);
            auto b_checked = try_checked_numeric_constant(b_expr, computation);
            auto c_checked = try_checked_numeric_constant(c_expr, computation);
            if (!a_checked || !b_checked || !c_checked) {
                return nullptr;
            }

            double a = *a_checked;
            double b = *b_checked;
            double c = *c_checked;

            int eq_a;
            lmmc_double_nearly_equal_tol(a, 0.0, 1e-9, 1e-9, &eq_a);
            if (eq_a) return nullptr;

            double delta = b * b - 4 * a * c;
            int eq_delta;
            lmmc_double_nearly_equal_tol(delta, 0.0, 1e-9, 1e-9, &eq_delta);

            if (!eq_delta && delta > 0) {
                double sqrt_delta = std::sqrt(delta);
                auto scalar = SymbolicExpr::number(1.0 / sqrt_delta);
                auto two_a = SymbolicExpr::number(2.0 * a);
                auto b_num = SymbolicExpr::number(b);
                auto two_a_x = SymbolicExpr::multiply(make_expr_ptr(*two_a), SymbolicExpr::variable(var));
                auto two_a_x_plus_b = SymbolicExpr::add(make_expr_ptr(*two_a_x), make_expr_ptr(*b_num));
                auto term1_arg = sym_sub(*two_a_x_plus_b, *SymbolicExpr::number(sqrt_delta));
                auto term2_arg = SymbolicExpr::add(make_expr_ptr(*two_a_x_plus_b), SymbolicExpr::number(sqrt_delta));
                auto term1 = SymbolicExpr::ln(make_expr_ptr(*term1_arg));
                auto term2 = SymbolicExpr::ln(make_expr_ptr(*term2_arg));
                return SymbolicExpr::multiply(scalar, sym_sub(*term1, *term2));
            } else if (!eq_delta && delta < 0) {
                double sqrt_neg_delta = std::sqrt(-delta);
                auto scalar = SymbolicExpr::number(2.0 / sqrt_neg_delta);
                auto two_a = SymbolicExpr::number(2.0 * a);
                auto b_num = SymbolicExpr::number(b);
                auto num = SymbolicExpr::add(
                    SymbolicExpr::multiply(make_expr_ptr(*two_a), SymbolicExpr::variable(var)),
                    make_expr_ptr(*b_num));
                auto inner = SymbolicExpr::divide(make_expr_ptr(*num), SymbolicExpr::number(sqrt_neg_delta));
                return SymbolicExpr::multiply(scalar, make_arctan(inner));
            }
        }
    } catch (const std::invalid_argument&) {
    } catch (const std::out_of_range&) {
    } catch (const std::runtime_error&) {
    }

    return nullptr;
}

std::shared_ptr<SymbolicExpr> IBPStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int depth) {

    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        ops = mul->operands();
    } else {
        ops.push_back(lamina::detail::node(expr));
    }

    int best_u_idx = -1;
    int best_score = 100;

    auto get_score = [&](const std::shared_ptr<const SymbolicNode>& node) -> int {
        auto e = lamina::detail::expression_from_node(node);
        if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
            if (fn->type() == FunctionNode::FuncType::Ln || fn->type() == FunctionNode::FuncType::Log) return 1;
            if (fn->type() == FunctionNode::FuncType::ArcSin || fn->type() == FunctionNode::FuncType::ArcTan) return 2;
            if (fn->type() == FunctionNode::FuncType::Sin || fn->type() == FunctionNode::FuncType::Cos) return 4;
            if (fn->type() == FunctionNode::FuncType::Exp) return 5;
        }
        if (!depends_on_integration_variable(e, var)) return 10;
        if (std::dynamic_pointer_cast<const VariableNode>(node)) return 3;
        if (std::dynamic_pointer_cast<const PowerNode>(node)) return 3;
        return 10;
    };

    if (ops.size() == 1) {
        int s = get_score(ops[0]);
        if (s <= 2) best_u_idx = 0;
    } else {
        for (size_t i = 0; i < ops.size(); ++i) {
            if (!depends_on_integration_variable(lamina::detail::expression_from_node(ops[i]), var)) continue;
            int s = get_score(ops[i]);
            if (s < best_score) {
                best_score = s;
                best_u_idx = (int)i;
            }
        }
    }

    if (best_u_idx == -1) return nullptr;

    auto u = lamina::detail::expression_from_node(ops[best_u_idx]);
    std::vector<std::shared_ptr<const SymbolicNode>> dv_ops;
    for (size_t i = 0; i < ops.size(); ++i) {
        if ((int)i != best_u_idx) dv_ops.push_back(ops[i]);
    }

    std::shared_ptr<SymbolicExpr> dv;
    if (dv_ops.empty()) dv = SymbolicExpr::number(1);
    else if (dv_ops.size() == 1) {
        dv = make_expr_ptr(lamina::detail::expression_from_node(dv_ops[0]));
    }
    else dv = lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(dv_ops));

    auto v = detail::propagate_result(
        ctx.integrate_recursive(*dv, var, computation, depth + 1));
    if (!v || !lamina::detail::node(v)) return nullptr;
    if (contains_unevaluated_integral(lamina::detail::node(v))) return nullptr;

    auto du_ptr = u.differentiate(var);
    if (!du_ptr) return nullptr;
    auto simplified_du = du_ptr->simplify();
    if (!simplified_du || !lamina::detail::node(simplified_du)) return nullptr;
    auto du = make_expr_ptr(*simplified_du);

    auto uv = SymbolicExpr::multiply(make_expr_ptr(u), v);
    auto vdu = SymbolicExpr::multiply(v, du);
    if (!uv || !vdu) return nullptr;
    vdu = vdu->cancel()->simplify();
    auto int_vdu = detail::propagate_result(
        ctx.integrate_recursive(*vdu, var, computation, depth + 1));
    if (!int_vdu || !lamina::detail::node(int_vdu) ||
        contains_unevaluated_integral(lamina::detail::node(int_vdu))) {
        return nullptr;
    }

    return sym_sub(*uv, *int_vdu);
}

} // namespace lamina
