#pragma once
#include "internal/integration_support.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS::integration_table_detail {
using FT = FunctionNode::FuncType;

inline auto u_is_var(const std::string& wc) {
    return [wc](const MatchMap& m, const std::string& var, ComputationContext&) {
        auto it = m.find(wc);
        if (it == m.end()) { return false; }
        auto value = std::dynamic_pointer_cast<const VariableNode>(detail::node(it->second));
        return value && !value->is_constant() && value->name() == var;
    };
}

inline auto u_is_var_a_indep(const std::string& u_wc, const std::string& a_wc) {
    return [u_matches = u_is_var(u_wc), a_wc](
        const MatchMap& m, const std::string& var, ComputationContext& context) {
        if (!u_matches(m, var, context)) { return false; }
        auto a = m.find(a_wc);
        return a != m.end() && !expression_depends_on_variable(detail::node(a->second), var);
    };
}

inline auto u_is_var_a_nonzero(const std::string& u_wc, const std::string& a_wc,
                               bool real_scale = false) {
    return [independent = u_is_var_a_indep(u_wc, a_wc), u_wc, a_wc, real_scale](
        const MatchMap& m, const std::string& var, ComputationContext& context) {
        if (!independent(m, var, context)) { return false; }
        std::optional<detail::AssumptionFacts> assumptions;
        if (context.assumptions()) { assumptions.emplace(*context.assumptions()); }
        const FactsQuery& facts = assumptions ? static_cast<const FactsQuery&>(*assumptions) : detail::no_facts();
        auto a = detail::node(m.at(a_wc));
        return (!real_scale ||
                (normalization_proved(detail::query_real_value(detail::node(m.at(u_wc)), facts, context)) &&
                 normalization_proved(detail::query_real_value(a, facts, context)))) &&
            normalization_proved(detail::query_nonzero_value(a, facts, Domain::Real, context));
    };
}

inline auto u_is_real_var(const std::string& wc) {
    return [u_matches = u_is_var(wc), wc](
        const MatchMap& m, const std::string& var, ComputationContext& context) {
        if (!u_matches(m, var, context)) { return false; }
        std::optional<detail::AssumptionFacts> assumptions;
        if (context.assumptions()) { assumptions.emplace(*context.assumptions()); }
        return normalization_proved(detail::query_real_value(detail::node(m.at(wc)),
            assumptions ? static_cast<const FactsQuery&>(*assumptions) : detail::no_facts(), context));
    };
}

void load_exponential_rules(std::vector<detail::IntegrationEntry>& table);
void load_trig_derivatives_rules(std::vector<detail::IntegrationEntry>& table);
void load_trig_squares_rules(std::vector<detail::IntegrationEntry>& table);
void load_trig_primitives_rules(std::vector<detail::IntegrationEntry>& table);
void load_inverse_trig_rules(std::vector<detail::IntegrationEntry>& table);
void load_hyperbolic_rules(std::vector<detail::IntegrationEntry>& table);
void load_algebraic_elementary_rules(std::vector<detail::IntegrationEntry>& table);
void load_algebraic_radicals_rules(std::vector<detail::IntegrationEntry>& table);
void load_algebraic_reciprocals_rules(std::vector<detail::IntegrationEntry>& table);
void load_polynomial_rules(std::vector<detail::IntegrationEntry>& table);
void load_logarithmic_rules(std::vector<detail::IntegrationEntry>& table);
void load_exponential_base_rules(std::vector<detail::IntegrationEntry>& table);
void load_scaled_radical_rules(std::vector<detail::IntegrationEntry>& table);
void load_scaled_trig_rules(std::vector<detail::IntegrationEntry>& table);
void load_exponential_trig_rules(std::vector<detail::IntegrationEntry>& table);
}
