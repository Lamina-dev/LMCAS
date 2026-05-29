/**
 * @file property_store.cpp
 * @brief Implementation of the PropertyStore class.
 */

#include "property_store.hpp"
#include <algorithm>

namespace lamina {

// ============================================================
// Domain specificity ordering
// ============================================================

int PropertyStore::domain_specificity(Domain domain) const {
    switch (domain) {
        case Domain::Complex:     return 0;
        case Domain::Real:        return 1;
        case Domain::Rational:    return 2;
        case Domain::Integer:     return 3;
        case Domain::Natural:     return 4;
        case Domain::PositiveInt: return 5;
    }
    return 0;
}

std::vector<Domain> PropertyStore::get_ancestor_domains(Domain domain) const {
    std::vector<Domain> ancestors;
    switch (domain) {
        case Domain::PositiveInt:
            ancestors.push_back(Domain::Natural);
            ancestors.push_back(Domain::Integer);
            ancestors.push_back(Domain::Rational);
            ancestors.push_back(Domain::Real);
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Natural:
            ancestors.push_back(Domain::Integer);
            ancestors.push_back(Domain::Rational);
            ancestors.push_back(Domain::Real);
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Integer:
            ancestors.push_back(Domain::Rational);
            ancestors.push_back(Domain::Real);
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Rational:
            ancestors.push_back(Domain::Real);
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Real:
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Complex:
            break;
    }
    return ancestors;
}

// ============================================================
// Sign implication and contradiction
// ============================================================

std::unordered_set<Sign, SignHash> PropertyStore::get_implied_signs(Sign sign) const {
    std::unordered_set<Sign, SignHash> implied;
    switch (sign) {
        case Sign::Positive:
            implied.insert(Sign::NonNegative);
            implied.insert(Sign::NonZero);
            break;
        case Sign::Negative:
            implied.insert(Sign::NonPositive);
            implied.insert(Sign::NonZero);
            break;
        case Sign::Zero:
            implied.insert(Sign::NonNegative);
            implied.insert(Sign::NonPositive);
            break;
        case Sign::NonNegative:
        case Sign::NonPositive:
        case Sign::NonZero:
            break;
    }
    return implied;
}

void PropertyStore::check_sign_contradiction(const std::string& symbol,
                                             const std::unordered_set<Sign, SignHash>& existing,
                                             Sign new_sign) {
    auto contradicts = [](Sign a, Sign b) -> bool {
        if (a == Sign::Positive && (b == Sign::Negative || b == Sign::Zero || b == Sign::NonPositive))
            return true;
        if (a == Sign::Negative && (b == Sign::Positive || b == Sign::Zero || b == Sign::NonNegative))
            return true;
        if (a == Sign::NonNegative && b == Sign::Negative)
            return true;
        if (a == Sign::NonPositive && b == Sign::Positive)
            return true;
        if (a == Sign::Zero && (b == Sign::Positive || b == Sign::Negative || b == Sign::NonZero))
            return true;
        if (a == Sign::NonZero && b == Sign::Zero)
            return true;
        return false;
    };

    for (Sign existing_sign : existing) {
        if (contradicts(new_sign, existing_sign)) {
            throw std::invalid_argument(
                "Contradiction for symbol '" + symbol + "': sign conflict between declared signs");
        }
    }
}

void PropertyStore::check_domain_sign_consistency(const std::string& symbol,
                                                  const SymbolProperties& props,
                                                  Domain new_domain) {
    // Natural is incompatible with Negative
    if (new_domain == Domain::Natural || new_domain == Domain::PositiveInt) {
        if (props.signs.count(Sign::Negative)) {
            throw std::invalid_argument(
                "Contradiction for symbol '" + symbol +
                "': domain Natural/PositiveInt incompatible with Negative sign");
        }
    }

    // PositiveInt is additionally incompatible with Zero and NonPositive
    if (new_domain == Domain::PositiveInt) {
        if (props.signs.count(Sign::Zero)) {
            throw std::invalid_argument(
                "Contradiction for symbol '" + symbol +
                "': domain PositiveInt incompatible with Zero sign");
        }
        if (props.signs.count(Sign::NonPositive)) {
            throw std::invalid_argument(
                "Contradiction for symbol '" + symbol +
                "': domain PositiveInt incompatible with NonPositive sign");
        }
    }
}

// ============================================================
// declare_domain
// ============================================================

void PropertyStore::declare_domain(const std::string& symbol, Domain domain) {
    auto& props = properties_[symbol];

    // Idempotent: same domain already set -> no-op
    if (props.most_specific_domain == domain) {
        return;
    }

    // Specificity preservation: if the new domain is an ancestor (less specific)
    // of the current most-specific domain, it's a no-op.
    if (domain_specificity(domain) < domain_specificity(props.most_specific_domain)) {
        return;
    }

    // Check domain-sign cross-constraints before applying
    check_domain_sign_consistency(symbol, props, domain);

    // Set the new most-specific domain
    props.most_specific_domain = domain;
}

// ============================================================
// declare_sign
// ============================================================

void PropertyStore::declare_sign(const std::string& symbol, Sign sign) {
    auto& props = properties_[symbol];

    // Idempotent: if sign already present, no-op
    if (props.signs.count(sign)) {
        return;
    }

    // Collect all signs to add (the declared sign + its implications)
    std::unordered_set<Sign, SignHash> to_add;
    to_add.insert(sign);
    auto implied = get_implied_signs(sign);
    to_add.insert(implied.begin(), implied.end());

    // Check each new sign against existing signs for contradictions
    for (Sign s : to_add) {
        if (!props.signs.count(s)) {
            check_sign_contradiction(symbol, props.signs, s);
        }
    }

    // Also check domain-sign cross-constraints for the new signs
    if (props.most_specific_domain == Domain::Natural ||
        props.most_specific_domain == Domain::PositiveInt) {
        for (Sign s : to_add) {
            if (s == Sign::Negative) {
                throw std::invalid_argument(
                    "Contradiction for symbol '" + symbol +
                    "': domain Natural/PositiveInt incompatible with Negative sign");
            }
        }
    }
    if (props.most_specific_domain == Domain::PositiveInt) {
        for (Sign s : to_add) {
            if (s == Sign::Zero) {
                throw std::invalid_argument(
                    "Contradiction for symbol '" + symbol +
                    "': domain PositiveInt incompatible with Zero sign");
            }
            if (s == Sign::NonPositive) {
                throw std::invalid_argument(
                    "Contradiction for symbol '" + symbol +
                    "': domain PositiveInt incompatible with NonPositive sign");
            }
        }
    }

    // Store all signs
    props.signs.insert(to_add.begin(), to_add.end());

    // Zero implies Integer domain
    if (sign == Sign::Zero) {
        if (domain_specificity(Domain::Integer) > domain_specificity(props.most_specific_domain)) {
            props.most_specific_domain = Domain::Integer;
        }
    }
}

// ============================================================
// declare_parity
// ============================================================

void PropertyStore::declare_parity(const std::string& symbol, Parity parity) {
    auto& props = properties_[symbol];

    // Idempotent: same parity already set -> no-op
    if (props.parity == parity) {
        return;
    }

    // Unknown can always be set
    if (parity == Parity::Unknown) {
        props.parity = parity;
        return;
    }

    // Contradiction: Even vs Odd
    if (props.parity != Parity::Unknown && props.parity != parity) {
        throw std::invalid_argument(
            "Contradiction for symbol '" + symbol + "': parity conflict (Even vs Odd)");
    }

    // Auto-promote to Integer domain when Even or Odd is declared
    if (domain_specificity(Domain::Integer) > domain_specificity(props.most_specific_domain)) {
        props.most_specific_domain = Domain::Integer;
    }

    props.parity = parity;
}

// ============================================================
// declare_bounded
// ============================================================

void PropertyStore::declare_bounded(const std::string& symbol, Boundedness bounded,
                                    std::optional<Interval> bounds) {
    auto& props = properties_[symbol];

    // Idempotent: same boundedness already set -> no-op
    if (props.boundedness == bounded) {
        return;
    }

    // Unknown can always be set
    if (bounded == Boundedness::Unknown) {
        props.boundedness = bounded;
        props.bounds = std::nullopt;
        return;
    }

    // Contradiction: Bounded vs Unbounded
    if (props.boundedness != Boundedness::Unknown && props.boundedness != bounded) {
        throw std::invalid_argument(
            "Contradiction for symbol '" + symbol +
            "': boundedness conflict (Bounded vs Unbounded)");
    }

    props.boundedness = bounded;
    if (bounds.has_value()) {
        props.bounds = bounds;
    }
}

// ============================================================
// Query methods
// ============================================================

Domain PropertyStore::get_domain(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Domain::Complex;
    }
    return it->second.most_specific_domain;
}

std::unordered_set<Sign, SignHash> PropertyStore::get_signs(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return {};
    }
    return it->second.signs;
}

Parity PropertyStore::get_parity(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Parity::Unknown;
    }
    return it->second.parity;
}

Boundedness PropertyStore::get_boundedness(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Boundedness::Unknown;
    }
    return it->second.boundedness;
}

std::optional<Interval> PropertyStore::get_bounds(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return std::nullopt;
    }
    return it->second.bounds;
}

bool PropertyStore::has_sign(const std::string& symbol, Sign sign) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return false;
    }
    return it->second.signs.count(sign) > 0;
}

bool PropertyStore::has_domain(const std::string& symbol, Domain domain) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        // Default is Complex; only Complex query returns true
        return domain == Domain::Complex;
    }
    // The symbol has domain D. It "has" domain X if X is an ancestor of D (or equal to D).
    // i.e., the symbol's specificity is >= the queried domain's specificity.
    return domain_specificity(it->second.most_specific_domain) >= domain_specificity(domain);
}

} // namespace lamina
