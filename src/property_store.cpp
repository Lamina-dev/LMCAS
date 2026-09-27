#include "internal/property_store_support.hpp"
#include "internal/symbolic_ast.hpp"
#include <algorithm>

namespace LMCAS {

using namespace property_detail;

PropertyStore::PropertyStore(PropertyStore&& other) noexcept
    : properties_(std::move(other.properties_)) {
    ++other.revision_;
}

PropertyStore& PropertyStore::operator=(const PropertyStore& other) {
    if (this != &other) {
        auto candidate = other.properties_;
        properties_.swap(candidate);
        ++revision_;
    }
    return *this;
}

PropertyStore& PropertyStore::operator=(PropertyStore&& other) noexcept {
    if (this != &other) {
        properties_ = std::move(other.properties_);
        ++revision_;
        ++other.revision_;
    }
    return *this;
}

namespace {

constexpr Domain domain_order[] = {
    Domain::Complex, Domain::Real, Domain::Algebraic, Domain::Rational,
    Domain::Integer, Domain::Natural, Domain::PositiveInt
};

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
bool signs_conflict(Sign a, Sign b) {
    struct Conflict {
        Sign incoming;
        std::initializer_list<Sign> existing;
    };
    const Conflict conflicts[] = {
        {Sign::Positive, {Sign::Negative, Sign::Zero, Sign::NonPositive}},
        {Sign::Negative, {Sign::Positive, Sign::Zero, Sign::NonNegative}},
        {Sign::NonNegative, {Sign::Negative}},
        {Sign::NonPositive, {Sign::Positive}},
        {Sign::Zero, {Sign::Positive, Sign::Negative, Sign::NonZero}},
        {Sign::NonZero, {Sign::Zero}}
    };
    for (const auto& conflict : conflicts) {
        if (conflict.incoming != a) {
            continue;
        }
        return std::find(conflict.existing.begin(), conflict.existing.end(), b)
            != conflict.existing.end();
    }
    return false;
}

}

int PropertyStore::domain_specificity(Domain domain) const {
    for (int i = 0; i < 7; ++i) {
        if (domain_order[i] == domain) return i;
    }
    return 0;
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
    for (Sign existing_sign : existing) {
        if (signs_conflict(new_sign, existing_sign)) {
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

void PropertyStore::declare_domain_unchecked(const std::string& symbol, Domain domain) {
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

    /// 超越域约束将符号域限制为 Real 或 Complex.
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

PropertyStoreResult PropertyStore::declare_domain(
    const std::string& symbol, Domain domain) {
    return declare_domain_checked(symbol, domain);
}

PropertyStoreResult PropertyStore::declare_domain_checked(
    const std::string& symbol,
    Domain domain) {
    return checked_property_update(*this, symbol, "declare_domain",
        [&](PropertyStore& candidate) {
            candidate.declare_domain_unchecked(symbol, domain);
        });
}

// declare_sign

void PropertyStore::declare_sign_unchecked(const std::string& symbol, Sign sign) {
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

PropertyStoreResult PropertyStore::declare_sign(
    const std::string& symbol, Sign sign) {
    return declare_sign_checked(symbol, sign);
}

PropertyStoreResult PropertyStore::declare_sign_checked(
    const std::string& symbol,
    Sign sign) {
    return checked_property_update(*this, symbol, "declare_sign",
        [&](PropertyStore& candidate) {
            candidate.declare_sign_unchecked(symbol, sign);
        });
}

// declare_parity

void PropertyStore::declare_parity_unchecked(const std::string& symbol, Parity parity) {
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

PropertyStoreResult PropertyStore::declare_parity(
    const std::string& symbol, Parity parity) {
    return declare_parity_checked(symbol, parity);
}

PropertyStoreResult PropertyStore::declare_parity_checked(
    const std::string& symbol,
    Parity parity) {
    return checked_property_update(*this, symbol, "declare_parity",
        [&](PropertyStore& candidate) {
            candidate.declare_parity_unchecked(symbol, parity);
        });
}

// declare_bounded

void PropertyStore::declare_bounded_unchecked(const std::string& symbol,
                                              Boundedness bounded,
                                              std::optional<Interval> bounds) {
    auto& props = properties_[symbol];

    if (props.boundedness == bounded) {
        if (bounded == Boundedness::Bounded && bounds) {
            props.bounds = std::move(bounds);
        }
        return;
    }

    if (bounded == Boundedness::Unknown) {
        props.boundedness = bounded;
        props.bounds = std::nullopt;
        return;
    }

    if (props.boundedness != Boundedness::Unknown && props.boundedness != bounded) {
        throw std::invalid_argument(
            "Contradiction for symbol '" + symbol +
            "': boundedness conflict (Bounded vs Unbounded)");
    }

    props.boundedness = bounded;
    props.bounds = std::move(bounds);
}

namespace {

using BoundsCandidateResult = Result<std::optional<Interval>>;

BoundsCandidateResult bounded_interval_candidate(
    const Interval& bounds, const std::optional<Interval>& existing,
    ComputationContext& context) {
    constexpr const char* operation = "declare_bounded";
    auto incoming = IntervalUnion::from_intervals_checked(
        {bounds}, context);
    if (!incoming) {
        if (incoming.error().code != CasErrc::Inconclusive || existing) {
            return BoundsCandidateResult::failure(incoming.error());
        }
        return BoundsCandidateResult::success(bounds);
    }
    if (incoming.value().is_empty()) {
        return BoundsCandidateResult::failure(
            CasErrc::InvalidArgument,
            "bounded interval must not be empty", operation);
    }

    std::optional<Interval> candidate =
        incoming.value().intervals().front();
    if (!existing) {
        return BoundsCandidateResult::success(std::move(candidate));
    }

    auto current = IntervalUnion::from_intervals_checked(
        {*existing}, context);
    if (!current) {
        return BoundsCandidateResult::failure(current.error());
    }
    auto intersection = current.value().intersect_checked(
        incoming.value(), context);
    if (!intersection) {
        return BoundsCandidateResult::failure(intersection.error());
    }
    if (intersection.value().is_empty()) {
        return BoundsCandidateResult::failure(
            CasErrc::InvalidArgument,
            "bounded interval contradicts existing bounds", operation);
    }
    candidate = intersection.value().intervals().front();
    return BoundsCandidateResult::success(std::move(candidate));
}

}


PropertyStoreResult PropertyStore::declare_bounded(
    const std::string& symbol,
    Boundedness bounded,
    std::optional<Interval> bounds) {
    return declare_bounded_checked(symbol, bounded, std::move(bounds));
}

PropertyStoreResult PropertyStore::declare_bounded_checked(
    const std::string& symbol,
    Boundedness bounded,
    std::optional<Interval> bounds) {
    ComputationContext context;
    return declare_bounded_checked(
        symbol, bounded, std::move(bounds), context);
}

PropertyStoreResult PropertyStore::declare_bounded_checked(
    const std::string& symbol,
    Boundedness bounded,
    std::optional<Interval> bounds,
    ComputationContext& context) {
    constexpr const char* operation = "declare_bounded";
    if (symbol.empty()) {
        return invalid_empty_symbol(operation);
    }
    if (bounds && bounded != Boundedness::Bounded) {
        return PropertyStoreResult::failure(
            CasErrc::InvalidArgument,
            "interval bounds require Bounded classification", operation);
    }

    try {
        std::optional<Interval> committed_bounds;
        if (bounds) {
            auto candidate = bounded_interval_candidate(
                *bounds, get_bounds(symbol), context);
            if (!candidate) {
                return PropertyStoreResult::failure(candidate.error());
            }
            committed_bounds = std::move(candidate.value());
        }

        PropertyStore candidate = *this;
        candidate.declare_bounded_unchecked(
            symbol, bounded, std::move(committed_bounds));
        *this = std::move(candidate);
        return PropertyStoreResult::success();
    } catch (const std::bad_alloc&) {
        return PropertyStoreResult::failure(
            CasErrc::ResourceLimit,
            "property-store allocation failed", operation);
    } catch (const std::invalid_argument& error) {
        return PropertyStoreResult::failure(
            CasErrc::InvalidArgument, error.what(), operation);
    } catch (const std::exception& error) {
        return PropertyStoreResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
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

}
