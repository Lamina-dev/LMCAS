/**
 * @file property_store.cpp
 * @brief Implementation of the PropertyStore class.
 */

#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <algorithm>
#include <limits>
#include <sstream>

namespace LMCAS {

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

PropertyStoreResult invalid_empty_symbol(const std::string& operation) {
    return PropertyStoreResult::failure(
        CasErrc::InvalidArgument, "symbol name must not be empty", operation);
}


template <typename F>
PropertyStoreResult checked_property_update(PropertyStore& store,
                                            const std::string& symbol,
                                            const std::string& operation,
                                            F&& update) {
    if (symbol.empty()) {
        return invalid_empty_symbol(operation);
    }
    try {
        PropertyStore candidate = store;
        update(candidate);
        store = std::move(candidate);
    } catch (const std::bad_alloc&) {
        return PropertyStoreResult::failure(
            CasErrc::ResourceLimit, "property-store allocation failed", operation);
    } catch (const std::invalid_argument& ex) {
        return PropertyStoreResult::failure(CasErrc::InvalidArgument, ex.what(), operation);
    } catch (const std::exception& ex) {
        return PropertyStoreResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
    return PropertyStoreResult::success();
}

PropertyStoreResult validate_property_interval(const Interval& interval,
                                               ComputationContext& context,
                                               const std::string& operation) {
    auto normalized = normalize_intervals_checked({interval}, context);
    if (!normalized) {
        return PropertyStoreResult::failure(
            normalized.error().code, normalized.error().message, operation);
    }
    if (normalized.value().empty()) {
        return PropertyStoreResult::failure(
            CasErrc::InvalidArgument,
            "property declaration interval must not be empty",
            operation);
    }
    return PropertyStoreResult::success();
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
    return checked_property_update(*this, symbol, "declare_bounded",
        [&](PropertyStore& candidate) {
            candidate.declare_bounded_unchecked(symbol, bounded, bounds);
        });
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

static bool endpoints_equivalent(const Endpoint& left, const Endpoint& right) {
    if (left.is_open != right.is_open ||
        left.is_neg_infinity != right.is_neg_infinity ||
        left.is_pos_infinity != right.is_pos_infinity) {
        return false;
    }
    if (!left.value || !right.value) {
        return left.value == right.value;
    }
    if (!LMCAS::detail::node(left.value) || !LMCAS::detail::node(right.value)) {
        return LMCAS::detail::node(left.value) == LMCAS::detail::node(right.value);
    }
    return LMCAS::detail::node(left.value)->compare(*LMCAS::detail::node(right.value)) == 0;
}

static bool intervals_equivalent(const Interval& left, const Interval& right) {
    return endpoints_equivalent(left.lower, right.lower) &&
           endpoints_equivalent(left.upper, right.upper);
}

/**
 * 辅助函数:判断两个区间是否具有已证明的交集.
 * 端点排序委托给受检区间运算;当前支持域之外的端点保持"交集未证明"状态.
 */
static Result<bool> intervals_overlap_checked_impl(
    const Interval& a,
    const Interval& b,
    ComputationContext& context) {
    auto left = IntervalUnion::from_intervals_checked({a}, context);
    if (!left) return Result<bool>::failure(left.error());
    auto right = IntervalUnion::from_intervals_checked({b}, context);
    if (!right) return Result<bool>::failure(right.error());
    auto intersection = left.value().intersect_checked(right.value(), context);
    if (!intersection) return Result<bool>::failure(intersection.error());
    return Result<bool>::success(!intersection.value().is_empty());
}

/**
 * Helper: returns true if outer interval fully contains inner interval.
 * outer covers inner iff outer.lower <= inner.lower AND inner.upper <= outer.upper
 * (with appropriate open/closed boundary handling).
 */
static Result<bool> interval_covers_checked_impl(
    const Interval& outer,
    const Interval& inner,
    ComputationContext& context) {
    auto normalized_outer = IntervalUnion::from_intervals_checked({outer}, context);
    if (!normalized_outer) {
        return Result<bool>::failure(normalized_outer.error());
    }
    auto normalized_inner = IntervalUnion::from_intervals_checked({inner}, context);
    if (!normalized_inner) {
        return Result<bool>::failure(normalized_inner.error());
    }
    if (normalized_inner.value().is_empty()) return Result<bool>::success(true);
    if (normalized_outer.value().is_empty()) return Result<bool>::success(false);

    auto intersection = normalized_outer.value().intersect_checked(
        normalized_inner.value(), context);
    if (!intersection) {
        return Result<bool>::failure(intersection.error());
    }
    const auto& intersection_intervals = intersection.value().intervals();
    const auto& inner_intervals = normalized_inner.value().intervals();
    if (intersection_intervals.size() != inner_intervals.size()) {
        return Result<bool>::success(false);
    }
    for (std::size_t i = 0; i < inner_intervals.size(); ++i) {
        if (!intervals_equivalent(intersection_intervals[i], inner_intervals[i])) {
            return Result<bool>::success(false);
        }
    }
    return Result<bool>::success(true);
}

// Continuity and differentiability declarations

PropertyStoreResult PropertyStore::declare_continuous(
    const std::string& symbol, const Interval& interval) {
    return declare_continuous_checked(symbol, interval);
}

PropertyStoreResult PropertyStore::declare_continuous_checked(
    const std::string& symbol,
    const Interval& interval,
    ComputationContext& context) {
    constexpr const char* operation = "declare_continuous";
    if (symbol.empty()) return invalid_empty_symbol(operation);
    auto valid = validate_property_interval(interval, context, operation);
    if (!valid) return valid;

    try {
        PropertyStore candidate = *this;
        auto& declarations = candidate.properties_[symbol].continuity_decls;
        for (const auto& declaration : declarations) {
            if (!declaration.is_differentiable &&
                intervals_equivalent(declaration.interval, interval)) {
                return PropertyStoreResult::success();
            }
            auto overlap = intervals_overlap_checked_impl(
                declaration.interval, interval, context);
            if (!overlap) {
                return PropertyStoreResult::failure(
                    overlap.error().code, overlap.error().message, operation);
            }
            if (overlap.value() && declaration.is_differentiable) {
                return PropertyStoreResult::failure(
                    CasErrc::InvalidArgument,
                    "Contradiction for symbol '" + symbol +
                        "': continuous-only declaration overlaps an existing "
                        "differentiable declaration",
                    operation);
            }
        }
        declarations.push_back({interval, false});
        *this = std::move(candidate);
    } catch (const std::bad_alloc&) {
        return PropertyStoreResult::failure(
            CasErrc::ResourceLimit, "property-store allocation failed", operation);
    } catch (const std::exception& ex) {
        return PropertyStoreResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
    return PropertyStoreResult::success();
}

PropertyStoreResult PropertyStore::declare_continuous_checked(
    const std::string& symbol,
    const Interval& interval) {
    ComputationContext context;
    return declare_continuous_checked(symbol, interval, context);
}

PropertyStoreResult PropertyStore::declare_differentiable(
    const std::string& symbol, const Interval& interval) {
    return declare_differentiable_checked(symbol, interval);
}

PropertyStoreResult PropertyStore::declare_differentiable_checked(
    const std::string& symbol,
    const Interval& interval,
    ComputationContext& context) {
    constexpr const char* operation = "declare_differentiable";
    if (symbol.empty()) return invalid_empty_symbol(operation);
    auto valid = validate_property_interval(interval, context, operation);
    if (!valid) return valid;

    try {
        PropertyStore candidate = *this;
        auto& declarations = candidate.properties_[symbol].continuity_decls;
        for (const auto& declaration : declarations) {
            if (declaration.is_differentiable &&
                intervals_equivalent(declaration.interval, interval)) {
                return PropertyStoreResult::success();
            }
        }
        declarations.push_back({interval, true});
        *this = std::move(candidate);
    } catch (const std::bad_alloc&) {
        return PropertyStoreResult::failure(
            CasErrc::ResourceLimit, "property-store allocation failed", operation);
    } catch (const std::exception& ex) {
        return PropertyStoreResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
    return PropertyStoreResult::success();
}

PropertyStoreResult PropertyStore::declare_differentiable_checked(
    const std::string& symbol,
    const Interval& interval) {
    ComputationContext context;
    return declare_differentiable_checked(symbol, interval, context);
}

Result<bool> PropertyStore::is_continuous(
    const std::string& symbol, const Interval& interval) const {
    return is_continuous_checked(symbol, interval);
}

Result<bool> PropertyStore::is_continuous_checked(
    const std::string& symbol,
    const Interval& interval,
    ComputationContext& context) const {
    constexpr const char* operation = "is_continuous";
    if (symbol.empty()) {
        return Result<bool>::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    auto valid = validate_property_interval(interval, context, operation);
    if (!valid) return Result<bool>::failure(valid.error());

    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Result<bool>::success(false);
    }

    for (const auto& decl : it->second.continuity_decls) {
        auto covers = interval_covers_checked_impl(decl.interval, interval, context);
        if (!covers) {
            return Result<bool>::failure(
                covers.error().code, covers.error().message, operation);
        }
        if (covers.value()) return Result<bool>::success(true);
    }
    return Result<bool>::success(false);
}

Result<bool> PropertyStore::is_continuous_checked(
    const std::string& symbol,
    const Interval& interval) const {
    ComputationContext context;
    return is_continuous_checked(symbol, interval, context);
}

Result<bool> PropertyStore::is_differentiable(
    const std::string& symbol, const Interval& interval) const {
    return is_differentiable_checked(symbol, interval);
}

Result<bool> PropertyStore::is_differentiable_checked(
    const std::string& symbol,
    const Interval& interval,
    ComputationContext& context) const {
    constexpr const char* operation = "is_differentiable";
    if (symbol.empty()) {
        return Result<bool>::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    auto valid = validate_property_interval(interval, context, operation);
    if (!valid) return Result<bool>::failure(valid.error());

    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Result<bool>::success(false);
    }

    for (const auto& decl : it->second.continuity_decls) {
        if (!decl.is_differentiable) continue;
        auto covers = interval_covers_checked_impl(decl.interval, interval, context);
        if (!covers) {
            return Result<bool>::failure(
                covers.error().code, covers.error().message, operation);
        }
        if (covers.value()) return Result<bool>::success(true);
    }
    return Result<bool>::success(false);
}

Result<bool> PropertyStore::is_differentiable_checked(
    const std::string& symbol,
    const Interval& interval) const {
    ComputationContext context;
    return is_differentiable_checked(symbol, interval, context);
}

// Monotonicity declarations

PropertyStoreResult PropertyStore::declare_monotonicity(
    const std::string& symbol,
    const std::string& variable,
    const Interval& interval,
    Monotonicity mono) {
    return declare_monotonicity_checked(symbol, variable, interval, mono);
}

PropertyStoreResult PropertyStore::declare_monotonicity_checked(
    const std::string& symbol,
    const std::string& variable,
    const Interval& interval,
    Monotonicity mono,
    ComputationContext& context) {
    constexpr const char* operation = "declare_monotonicity";
    if (symbol.empty()) return invalid_empty_symbol(operation);
    if (variable.empty()) {
        return PropertyStoreResult::failure(
            CasErrc::InvalidArgument, "variable name must not be empty", operation);
    }
    auto valid = validate_property_interval(interval, context, operation);
    if (!valid) return valid;

    try {
        PropertyStore candidate = *this;
        auto& declarations = candidate.properties_[symbol].monotonicity_decls;
        for (const auto& declaration : declarations) {
            if (declaration.variable == variable && declaration.type == mono &&
                intervals_equivalent(declaration.interval, interval)) {
                return PropertyStoreResult::success();
            }
        }
        declarations.push_back({variable, interval, mono});
        *this = std::move(candidate);
    } catch (const std::bad_alloc&) {
        return PropertyStoreResult::failure(
            CasErrc::ResourceLimit, "property-store allocation failed", operation);
    } catch (const std::exception& ex) {
        return PropertyStoreResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
    return PropertyStoreResult::success();
}

PropertyStoreResult PropertyStore::declare_monotonicity_checked(
    const std::string& symbol,
    const std::string& variable,
    const Interval& interval,
    Monotonicity mono) {
    ComputationContext context;
    return declare_monotonicity_checked(symbol, variable, interval, mono, context);
}

Result<Monotonicity> PropertyStore::get_monotonicity(
    const std::string& symbol,
    const std::string& variable,
    const Interval& interval) const {
    return get_monotonicity_checked(symbol, variable, interval);
}

Result<Monotonicity> PropertyStore::get_monotonicity_checked(
    const std::string& symbol,
    const std::string& variable,
    const Interval& interval,
    ComputationContext& context) const {
    constexpr const char* operation = "get_monotonicity";
    if (symbol.empty() || variable.empty()) {
        return Result<Monotonicity>::failure(
            CasErrc::InvalidArgument,
            symbol.empty() ? "symbol name must not be empty" : "variable name must not be empty",
            operation);
    }
    auto valid = validate_property_interval(interval, context, operation);
    if (!valid) return Result<Monotonicity>::failure(valid.error());

    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Result<Monotonicity>::success(Monotonicity::Unknown);
    }

    for (const auto& decl : it->second.monotonicity_decls) {
        if (decl.variable != variable) continue;
        auto covers = interval_covers_checked_impl(decl.interval, interval, context);
        if (!covers) {
            return Result<Monotonicity>::failure(
                covers.error().code, covers.error().message, operation);
        }
        if (covers.value()) return Result<Monotonicity>::success(decl.type);
    }
    return Result<Monotonicity>::success(Monotonicity::Unknown);
}

Result<Monotonicity> PropertyStore::get_monotonicity_checked(
    const std::string& symbol,
    const std::string& variable,
    const Interval& interval) const {
    ComputationContext context;
    return get_monotonicity_checked(symbol, variable, interval, context);
}

// Transcendental classification

void PropertyStore::declare_transcendental_unchecked(const std::string& symbol) {
    auto& props = properties_[symbol];

    // Idempotent: already transcendental -> no-op
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

PropertyStoreResult PropertyStore::declare_transcendental(
    const std::string& symbol) {
    return declare_transcendental_checked(symbol);
}

PropertyStoreResult PropertyStore::declare_transcendental_checked(
    const std::string& symbol) {
    return checked_property_update(*this, symbol, "declare_transcendental",
        [&](PropertyStore& candidate) {
            candidate.declare_transcendental_unchecked(symbol);
        });
}

bool PropertyStore::is_transcendental(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return false;
    }
    return it->second.transcendental;
}

// Finiteness classification

void PropertyStore::declare_finiteness_unchecked(const std::string& symbol, Finiteness f) {
    auto& props = properties_[symbol];

    // Idempotent: same value -> no-op
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

    // Implication: Finite -> Bounded
    if (f == Finiteness::Finite) {
        declare_bounded_unchecked(symbol, Boundedness::Bounded, std::nullopt);
    }
}

PropertyStoreResult PropertyStore::declare_finiteness(
    const std::string& symbol, Finiteness f) {
    return declare_finiteness_checked(symbol, f);
}

PropertyStoreResult PropertyStore::declare_finiteness_checked(
    const std::string& symbol,
    Finiteness f) {
    return checked_property_update(*this, symbol, "declare_finiteness",
        [&](PropertyStore& candidate) {
            candidate.declare_finiteness_unchecked(symbol, f);
        });
}

Finiteness PropertyStore::get_finiteness(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Finiteness::Unknown;
    }
    return it->second.finiteness;
}

// Matrix definiteness

void PropertyStore::declare_definiteness_unchecked(const std::string& symbol, Definiteness d) {
    auto& props = properties_[symbol];

    // Idempotent: same value -> no-op
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

PropertyStoreResult PropertyStore::declare_definiteness(
    const std::string& symbol, Definiteness d) {
    return declare_definiteness_checked(symbol, d);
}

PropertyStoreResult PropertyStore::declare_definiteness_checked(
    const std::string& symbol,
    Definiteness d) {
    return checked_property_update(*this, symbol, "declare_definiteness",
        [&](PropertyStore& candidate) {
            candidate.declare_definiteness_unchecked(symbol, d);
        });
}

Definiteness PropertyStore::get_definiteness(const std::string& symbol) const {
    auto it = properties_.find(symbol);
    if (it == properties_.end()) {
        return Definiteness::Unknown;
    }
    return it->second.definiteness;
}

// Periodicity

void PropertyStore::declare_periodic_unchecked(
    const std::string& symbol,
    const std::shared_ptr<SymbolicExpr>& period) {
    if (!period) {
        throw std::invalid_argument(
            "declare_periodic: period must not be null for symbol '" + symbol + "'");
    }
    auto& props = properties_[symbol];
    props.period = period;
}

PropertyStoreResult PropertyStore::declare_periodic(
    const std::string& symbol,
    const std::shared_ptr<SymbolicExpr>& period) {
    return declare_periodic_checked(symbol, period);
}

PropertyStoreResult PropertyStore::declare_periodic_checked(
    const std::string& symbol,
    const std::shared_ptr<SymbolicExpr>& period) {
    if (!period) {
        return PropertyStoreResult::failure(
            CasErrc::InvalidArgument, "period must not be null", "declare_periodic");
    }
    return checked_property_update(*this, symbol, "declare_periodic",
        [&](PropertyStore& candidate) {
            candidate.declare_periodic_unchecked(symbol, period);
        });
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

} // namespace LMCAS
