#pragma once
#include "internal/integration_support.hpp"

namespace LMCAS::integration_table_detail {
using Category = IntegrationTable::Category;
using FT = FunctionNode::FuncType;
inline auto u_is_var(const std::string& wc) {
        return [wc](const MatchMap& m, const std::string& var) -> bool {
            auto it = m.find(wc);
            if (it == m.end()) { return false; }
            auto v = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(it->second));
            return v && !v->is_constant() && v->name() == var;
        };
}
inline auto u_is_var_a_indep(const std::string& u_wc, const std::string& a_wc) {
        return [u_wc, a_wc](const MatchMap& m, const std::string& var) -> bool {
            auto it_u = m.find(u_wc);
            if (it_u == m.end()) { return false; }
            auto v = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(it_u->second));
            if (!v || v->is_constant() || v->name() != var) { return false; }
            auto it_a = m.find(a_wc);
            if (it_a == m.end()) { return false; }
            return !expression_depends_on_variable(LMCAS::detail::node(it_a->second), var);
        };
}

void load_exponential_rules(IntegrationTable& table);
void load_trig_derivatives_rules(IntegrationTable& table);
void load_trig_squares_rules(IntegrationTable& table);
void load_trig_primitives_rules(IntegrationTable& table);
void load_inverse_trig_rules(IntegrationTable& table);
void load_hyperbolic_rules(IntegrationTable& table);
void load_algebraic_elementary_rules(IntegrationTable& table);
void load_algebraic_radicals_rules(IntegrationTable& table);
void load_algebraic_reciprocals_rules(IntegrationTable& table);
void load_polynomial_rules(IntegrationTable& table);
void load_logarithmic_rules(IntegrationTable& table);
void load_exponential_base_rules(IntegrationTable& table);
void load_scaled_radical_rules(IntegrationTable& table);
void load_scaled_trig_rules(IntegrationTable& table);
void load_exponential_trig_rules(IntegrationTable& table);
}
