#include "internal/integration_table_rules.hpp"

namespace LMCAS {

const std::vector<IntegrationEntry> IntegrationTable::empty_entries_;

IntegrationTable::IntegrationTable() {
    load_defaults();
}

void IntegrationTable::add_entry(Category cat, const IntegrationEntry& entry) {
    entries_[static_cast<int>(cat)].push_back(entry);

    auto& vec = entries_[static_cast<int>(cat)];
    std::sort(vec.begin(), vec.end(), [](const IntegrationEntry& a, const IntegrationEntry& b) {
        return a.priority < b.priority;
    });
}

void IntegrationTable::clear_category(Category cat) {
    entries_[static_cast<int>(cat)].clear();
}

const std::vector<IntegrationEntry>& IntegrationTable::get_entries(Category cat) const {
    auto it = entries_.find(static_cast<int>(cat));
    if (it == entries_.end()) { return empty_entries_; }
    return it->second;
}

std::vector<const IntegrationEntry*> IntegrationTable::get_all_sorted() const {
    std::vector<const IntegrationEntry*> all;
    for (const auto& [cat, vec] : entries_) {
        for (const auto& entry : vec) {
            all.push_back(&entry);
        }
    }
    std::sort(all.begin(), all.end(), [](const IntegrationEntry* a, const IntegrationEntry* b) {
        return a->priority < b->priority;
    });
    return all;
}

void IntegrationTable::load_defaults() {
    integration_table_detail::load_exponential_rules(*this);
    integration_table_detail::load_trig_derivatives_rules(*this);
    integration_table_detail::load_trig_squares_rules(*this);
    integration_table_detail::load_trig_primitives_rules(*this);
    integration_table_detail::load_inverse_trig_rules(*this);
    integration_table_detail::load_hyperbolic_rules(*this);
    integration_table_detail::load_algebraic_elementary_rules(*this);
    integration_table_detail::load_algebraic_radicals_rules(*this);
    integration_table_detail::load_algebraic_reciprocals_rules(*this);
    integration_table_detail::load_polynomial_rules(*this);
    integration_table_detail::load_logarithmic_rules(*this);
    integration_table_detail::load_exponential_base_rules(*this);
    integration_table_detail::load_scaled_radical_rules(*this);
    integration_table_detail::load_scaled_trig_rules(*this);
    integration_table_detail::load_exponential_trig_rules(*this);
}

}
