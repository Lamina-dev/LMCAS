/**
 * @file property_store.cpp
 * @brief Implementation of the PropertyStore class.
 */

#include "property_store.hpp"
#include <algorithm>
#include <limits>
#include <sstream>

namespace lamina {

// String conversion helpers for diagnostic messages

namespace {

std::string domain_str(Domain d) {
    switch (d) {
        case Domain::Complex:     return "Complex";
        case Domain::Real:        return "Real";
        case Domain::Algebraic:   return "Algebraic";
        case Domain::Rational:    return "Rational";
        case Domain::Integer:     return "Integer";
        case Domain::Natural:     return "Natural";
        case Domain::PositiveInt: return "PositiveInt";
    }
    return "Complex";
}

std::string sign_str(Sign s) {
    switch (s) {
        case Sign::Positive:    return "Positive";
        case Sign::Negative:    return "Negative";
        case Sign::NonNegative: return "NonNegative";
        case Sign::NonPositive: return "NonPositive";
        case Sign::Zero:        return "Zero";
        case Sign::NonZero:     return "NonZero";
    }
    return "Unknown";
}

} // anonymous namespace

// Domain specificity ordering

int PropertyStore::domain_specificity(Domain domain) const {
    switch (domain) {
        case Domain::Complex:     return 0;
        case Domain::Real:        return 1;
        case Domain::Algebraic:   return 2;
        case Domain::Rational:    return 3;
        case Domain::Integer:     return 4;
        case Domain::Natural:     return 5;
        case Domain::PositiveInt: return 6;
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
            ancestors.push_back(Domain::Algebraic);
            ancestors.push_back(Domain::Real);
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Natural:
            ancestors.push_back(Domain::Integer);
            ancestors.push_back(Domain::Rational);
            ancestors.push_back(Domain::Algebraic);
            ancestors.push_back(Domain::Real);
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Integer:
            ancestors.push_back(Domain::Rational);
            ancestors.push_back(Domain::Algebraic);
            ancestors.push_back(Domain::Real);
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Rational:
            ancestors.push_back(Domain::Algebraic);
            ancestors.push_back(Domain::Real);
            ancestors.push_back(Domain::Complex);
            break;
        case Domain::Algebraic:
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

// Sign implication and contradiction

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
                "Contradiction for symbol '" + symbol + "': cannot declare sign " +
                sign_str(new_sign) + " because it conflicts with existing sign " +
                sign_str(existing_sign));
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
                "': domain " + domain_str(new_domain) +
                " is incompatible with existing sign Negative");
        }
    }

    // PositiveInt is additionally incompatible with Zero and NonPositive
    if (new_domain == Domain::PositiveInt) {
        if (props.signs.count(Sign::Zero)) {
            throw std::invalid_argument(
                "Contradiction for symbol '" + symbol +
                "': domain PositiveInt is incompatible with existing sign Zero");
        }
        if (props.signs.count(Sign::NonPositive)) {
            throw std::invalid_argument(
                "Contradiction for symbol '" + symbol +
                "': domain PositiveInt is incompatible with existing sign NonPositive");
        }
    }
}

// declare_domain

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

    // Transcendental guard: transcendental symbols cannot be Algebraic or more specific
    if (props.transcendental && domain_specificity(domain) > domain_specificity(Domain::Real)) {
        throw std::invalid_argument(
            "Contradiction for symbol '" + symbol +
            "': cannot declare domain " + domain_str(domain) +
            " because symbol is Transcendental (existing domain: " +
            domain_str(props.most_specific_domain) + ")");
    }

    // Check domain-sign cross-constraints before applying
    check_domain_sign_consistency(symbol, props, domain);

    // Set the new most-specific domain
    props.most_specific_domain = domain;
}

// declare_sign

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
                    "': cannot declare sign Negative because it conflicts with domain " +
                    domain_str(props.most_specific_domain));
            }
        }
    }
    if (props.most_specific_domain == Domain::PositiveInt) {
        for (Sign s : to_add) {
            if (s == Sign::Zero) {
                throw std::invalid_argument(
                    "Contradiction for symbol '" + symbol +
                    "': cannot declare sign Zero because it conflicts with domain PositiveInt");
            }
            if (s == Sign::NonPositive) {
                throw std::invalid_argument(
                    "Contradiction for symbol '" + symbol +
                    "': cannot declare sign NonPositive because it conflicts with domain PositiveInt");
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

// declare_parity

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

// declare_bounded

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

// Query methods

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

// Interval overlap detection (private helper)

/**
 * Helper: returns true if two intervals have a non-empty intersection.
 * Uses numeric comparison of endpoints.
 */
static bool intervals_overlap_impl(const Interval& a, const Interval& b) {
    // If either is empty, no overlap.
    if (a.is_empty() || b.is_empty()) return false;

    // If both are entire line, they overlap.
    if (a.is_entire_line() || b.is_entire_line()) return true;

    // Get numeric bounds for comparison.
    auto get_lower = [](const Interval& iv) -> double {
        if (iv.lower.is_neg_infinity) return -std::numeric_limits<double>::infinity();
        if (iv.lower.value) return iv.lower.value->to_numeric();
        return 0.0;
    };
    auto get_upper = [](const Interval& iv) -> double {
        if (iv.upper.is_pos_infinity) return std::numeric_limits<double>::infinity();
        if (iv.upper.value) return iv.upper.value->to_numeric();
        return 0.0;
    };

    double a_lo = get_lower(a), a_hi = get_upper(a);
    double b_lo = get_lower(b), b_hi = get_upper(b);

    // No overlap if one ends before the other starts.
    if (a_hi < b_lo || b_hi < a_lo) return false;

    // Boundary cases: if they touch at a single point, check open/closed.
    if (a_hi == b_lo) {
        // Overlap only if a's upper is closed AND b's lower is closed.
        return !a.upper.is_open && !b.lower.is_open;
    }
    if (b_hi == a_lo) {
        return !b.upper.is_open && !a.lower.is_open;
    }

    return true;
}

bool PropertyStore::intervals_overlap(const Interval& a, const Interval& b) const {
    return intervals_overlap_impl(a, b);
}

/**
 * Helper: returns true if outer interval fully contains inner interval.
 * outer covers inner iff outer.lower <= inner.lower AND inner.upper <= outer.upper
 * (with appropriate open/closed boundary handling).
 */
static bool interval_covers_impl(const Interval& outer, const Interval& inner) {
    if (inner.is_empty()) return true;
    if (outer.is_entire_line()) return true;
    if (outer.is_empty()) return false;

    auto get_lower = [](const Interval& iv) -> double {
        if (iv.lower.is_neg_infinity) return -std::numeric_limits<double>::infinity();
        if (iv.lower.value) return iv.lower.value->to_numeric();
        return 0.0;
    };
    auto get_upper = [](const Interval& iv) -> double {
        if (iv.upper.is_pos_infinity) return std::numeric_limits<double>::infinity();
        if (iv.upper.value) return iv.upper.value->to_numeric();
        return 0.0;
    };

    double o_lo = get_lower(outer), o_hi = get_upper(outer);
    double i_lo = get_lower(inner), i_hi = get_upper(inner);

    // Check lower bound: outer.lower must be <= inner.lower
    if (o_lo > i_lo) return false;
    if (o_lo == i_lo && outer.lower.is_open && !inner.lower.is_open) return false;

    // Check upper bound: outer.upper must be >= inner.upper
    if (o_hi < i_hi) return false;
    if (o_hi == i_hi && outer.upper.is_open && !inner.upper.is_open) return false;

    return true;
}

bool PropertyStore::interval_covers(const Interval& outer, const Interval& inner) const {
    return interval_covers_impl(outer, inner);
}

// Continuity and differentiability declarations

void PropertyStore::declare_continuous(const std::string& symbol, const Interval& interval) {
    auto& props = properties_[symbol];

    // Check for contradictions with existing declarations.
    // Declaring continuous-only on an interval that overlaps with an existing
    // differentiable declaration is a downgrade contradiction.
    for (const auto& decl : props.continuity_decls) {
        if (intervals_overlap_impl(decl.interval, interval)) {
            if (decl.is_differentiable) {
                // Existing differentiable declaration on overlapping interval.
                // Declaring continuous-only would be a downgrade — contradiction.
                throw std::invalid_argument(
                    "Contradiction for symbol '" + symbol +
                    "': cannot declare continuous-only on interval overlapping "
                    "an existing differentiable declaration");
            }
            // Existing continuous-only on overlapping interval — idempotent, no conflict.
        }
    }

    // Store the declaration.
    props.continuity_decls.push_back({interval, false});
}

void PropertyStore::declare_differentiable(const std::string& symbol, const Interval& interval) {
    auto& props = properties_[symbol];

    // Check for contradictions with existing declarations.
    // Declaring differentiable on an interval that overlaps with an existing
    // continuous-only declaration is an upgrade — that's fine.
    // Declaring differentiable on an interval that overlaps with an existing
    // differentiable declaration is idempotent — also fine.
    // No contradiction cases for upgrading to differentiable.

    // Store the differentiable declaration (implies continuous).
    props.continuity_decls.push_back({interval, true});
}

bool PropertyStore::is_continuous(const std::string& symbol, const Interval& interval) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return false;
    }

    // A symbol is continuous on the queried interval if any stored declaration
    // (continuous or differentiable) covers the queried interval.
    for (const auto& decl : it->second.continuity_decls) {
        if (interval_covers_impl(decl.interval, interval)) {
            return true;
        }
    }

    return false;
}

bool PropertyStore::is_differentiable(const std::string& symbol, const Interval& interval) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return false;
    }

    // A symbol is differentiable on the queried interval only if a differentiable
    // declaration covers the queried interval.
    for (const auto& decl : it->second.continuity_decls) {
        if (decl.is_differentiable && interval_covers_impl(decl.interval, interval)) {
            return true;
        }
    }

    return false;
}

// Monotonicity declarations

void PropertyStore::declare_monotonicity(const std::string& symbol, const std::string& variable,
                                         const Interval& interval, Monotonicity mono) {
    auto& props = properties_[symbol];
    props.monotonicity_decls.push_back({variable, interval, mono});
}

Monotonicity PropertyStore::get_monotonicity(const std::string& symbol, const std::string& variable,
                                             const Interval& interval) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Monotonicity::Unknown;
    }

    for (const auto& decl : it->second.monotonicity_decls) {
        if (decl.variable == variable && interval_covers_impl(decl.interval, interval)) {
            return decl.type;
        }
    }

    return Monotonicity::Unknown;
}

// Transcendental classification

void PropertyStore::declare_transcendental(const std::string& symbol) {
    auto& props = properties_[symbol];

    // Idempotent: already transcendental → no-op
    if (props.transcendental) {
        return;
    }

    // Transcendental numbers are Real but NOT Algebraic (or any sub-domain).
    // If the symbol already has a domain more specific than Real, that contradicts.
    if (domain_specificity(props.most_specific_domain) > domain_specificity(Domain::Real)) {
        throw std::invalid_argument(
            "Contradiction for symbol '" + symbol +
            "': cannot declare Transcendental because existing domain is " +
            domain_str(props.most_specific_domain) +
            " (Transcendental contradicts Algebraic or more specific domains)");
    }

    // Set domain to Real (if currently less specific)
    if (domain_specificity(props.most_specific_domain) < domain_specificity(Domain::Real)) {
        props.most_specific_domain = Domain::Real;
    }

    props.transcendental = true;
}

bool PropertyStore::is_transcendental(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return false;
    }
    return it->second.transcendental;
}

// Finiteness classification

void PropertyStore::declare_finiteness(const std::string& symbol, Finiteness f) {
    auto& props = properties_[symbol];

    // Idempotent: same value → no-op
    if (props.finiteness == f) {
        return;
    }

    // Unknown can always be set
    if (f == Finiteness::Unknown) {
        props.finiteness = f;
        return;
    }

    // Contradiction: Finite + Divergent
    if (props.finiteness != Finiteness::Unknown && props.finiteness != f) {
        throw std::invalid_argument(
            "Contradiction for symbol '" + symbol +
            "': cannot be both Finite and Divergent");
    }

    props.finiteness = f;

    // Implication: Finite → Bounded
    if (f == Finiteness::Finite) {
        declare_bounded(symbol, Boundedness::Bounded);
    }
}

Finiteness PropertyStore::get_finiteness(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Finiteness::Unknown;
    }
    return it->second.finiteness;
}

// Matrix definiteness

void PropertyStore::declare_definiteness(const std::string& symbol, Definiteness d) {
    auto& props = properties_[symbol];

    // Idempotent: same value → no-op
    if (props.definiteness == d) {
        return;
    }

    // Unknown can always be set
    if (d == Definiteness::Unknown) {
        props.definiteness = d;
        return;
    }

    // Check contradictions against existing definiteness
    if (props.definiteness != Definiteness::Unknown) {
        auto conflicts = [](Definiteness existing, Definiteness incoming) -> bool {
            // PositiveDefinite contradicts NegativeDefinite, NegativeSemiDefinite, Indefinite
            if (existing == Definiteness::PositiveDefinite &&
                (incoming == Definiteness::NegativeDefinite ||
                 incoming == Definiteness::NegativeSemiDefinite ||
                 incoming == Definiteness::Indefinite))
                return true;
            // NegativeDefinite contradicts PositiveDefinite, PositiveSemiDefinite, Indefinite
            if (existing == Definiteness::NegativeDefinite &&
                (incoming == Definiteness::PositiveDefinite ||
                 incoming == Definiteness::PositiveSemiDefinite ||
                 incoming == Definiteness::Indefinite))
                return true;
            // PositiveSemiDefinite contradicts NegativeDefinite, Indefinite
            if (existing == Definiteness::PositiveSemiDefinite &&
                (incoming == Definiteness::NegativeDefinite ||
                 incoming == Definiteness::Indefinite))
                return true;
            // NegativeSemiDefinite contradicts PositiveDefinite, Indefinite
            if (existing == Definiteness::NegativeSemiDefinite &&
                (incoming == Definiteness::PositiveDefinite ||
                 incoming == Definiteness::Indefinite))
                return true;
            // Indefinite contradicts PositiveDefinite, NegativeDefinite,
            // PositiveSemiDefinite, NegativeSemiDefinite
            if (existing == Definiteness::Indefinite &&
                (incoming == Definiteness::PositiveDefinite ||
                 incoming == Definiteness::NegativeDefinite ||
                 incoming == Definiteness::PositiveSemiDefinite ||
                 incoming == Definiteness::NegativeSemiDefinite))
                return true;
            return false;
        };

        if (conflicts(props.definiteness, d)) {
            throw std::invalid_argument(
                "Contradiction for symbol '" + symbol +
                "': definiteness conflict between existing and new declaration");
        }
        if (props.definiteness == Definiteness::PositiveSemiDefinite &&
            d == Definiteness::PositiveDefinite) {
            props.definiteness = d;
            return;
        }
        if (props.definiteness == Definiteness::NegativeSemiDefinite &&
            d == Definiteness::NegativeDefinite) {
            props.definiteness = d;
            return;
        }
        if (props.definiteness == Definiteness::PositiveDefinite &&
            d == Definiteness::PositiveSemiDefinite) {
            return; // already stronger
        }
        if (props.definiteness == Definiteness::NegativeDefinite &&
            d == Definiteness::NegativeSemiDefinite) {
            return; // already stronger
        }
    }

    props.definiteness = d;
}

Definiteness PropertyStore::get_definiteness(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Definiteness::Unknown;
    }
    return it->second.definiteness;
}

// Periodicity

void PropertyStore::declare_periodic(const std::string& symbol,
                                     const std::shared_ptr<SymbolicExpr>& period) {
    if (!period) {
        throw std::invalid_argument(
            "declare_periodic: period must not be null for symbol '" + symbol + "'");
    }
    auto& props = properties_[symbol];
    props.period = period;
}

std::optional<std::shared_ptr<SymbolicExpr>> PropertyStore::get_period(
    const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return std::nullopt;
    }
    return it->second.period;
}

bool PropertyStore::is_periodic(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return false;
    }
    return it->second.period.has_value();
}

// get_all_symbols

std::vector<std::string> PropertyStore::get_all_symbols() const {
    std::vector<std::string> symbols;
    symbols.reserve(properties_.size());
    for (const auto& [name, _] : properties_) {
        symbols.push_back(name);
    }
    // Sort for deterministic serialization order
    std::sort(symbols.begin(), symbols.end());
    return symbols;
}

std::vector<PropertyStore::ContinuityInfo> PropertyStore::get_continuity_decls(
    const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) return {};
    std::vector<ContinuityInfo> result;
    for (const auto& decl : it->second.continuity_decls) {
        result.push_back({decl.interval, decl.is_differentiable});
    }
    return result;
}

std::vector<PropertyStore::MonotonicityInfo> PropertyStore::get_monotonicity_decls(
    const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) return {};
    std::vector<MonotonicityInfo> result;
    for (const auto& decl : it->second.monotonicity_decls) {
        result.push_back({decl.variable, decl.interval, decl.type});
    }
    return result;
}

} // namespace lamina
