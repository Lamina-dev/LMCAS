#include "internal/integration_table_rules.hpp"

namespace LMCAS::detail {

std::vector<IntegrationEntry> make_integration_entries() {
    std::vector<IntegrationEntry> entries;
    integration_table_detail::load_exponential_rules(entries);
    integration_table_detail::load_trig_derivatives_rules(entries);
    integration_table_detail::load_trig_squares_rules(entries);
    integration_table_detail::load_trig_primitives_rules(entries);
    integration_table_detail::load_inverse_trig_rules(entries);
    integration_table_detail::load_hyperbolic_rules(entries);
    integration_table_detail::load_algebraic_elementary_rules(entries);
    integration_table_detail::load_algebraic_radicals_rules(entries);
    integration_table_detail::load_algebraic_reciprocals_rules(entries);
    integration_table_detail::load_polynomial_rules(entries);
    integration_table_detail::load_logarithmic_rules(entries);
    integration_table_detail::load_exponential_base_rules(entries);
    integration_table_detail::load_scaled_radical_rules(entries);
    integration_table_detail::load_scaled_trig_rules(entries);
    integration_table_detail::load_exponential_trig_rules(entries);
    std::stable_sort(entries.begin(), entries.end(),
        [](const IntegrationEntry& a, const IntegrationEntry& b) {
            return a.priority < b.priority;
        });
    return entries;
}

}
