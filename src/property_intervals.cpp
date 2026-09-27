#include "internal/property_store_support.hpp"
#include "internal/symbolic_ast.hpp"
#include <algorithm>

namespace LMCAS {

using namespace property_detail;

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
 * @brief 判断两个区间是否存在已证明的交集。
 * 端点排序采用受检区间运算；超出支持域时，交集判定仍为未决。
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
 * @brief 判断 outer 是否完全包含 inner。
 * 包含条件为 outer.lower <= inner.lower 且 inner.upper <= outer.upper，
 * 并结合端点的开闭属性判断边界。
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


Result<bool> PropertyStore::check_continuous_declarations(
    const std::vector<SymbolProperties::ContinuityDecl>& declarations,
    const std::string& symbol, const Interval& interval,
    ComputationContext& context) {
    constexpr const char* operation = "declare_continuous";
    for (const auto& declaration : declarations) {
        if (!declaration.is_differentiable &&
            intervals_equivalent(declaration.interval, interval)) {
            return Result<bool>::success(false);
        }
        auto overlap = intervals_overlap_checked_impl(
            declaration.interval, interval, context);
        if (!overlap) {
            return Result<bool>::failure(
                overlap.error().code, overlap.error().message, operation);
        }
        if (overlap.value() && declaration.is_differentiable) {
            return Result<bool>::failure(
                CasErrc::InvalidArgument,
                "Contradiction for symbol '" + symbol +
                    "': continuous-only declaration overlaps an existing "
                    "differentiable declaration",
                operation);
        }
    }
    return Result<bool>::success(true);
}

bool PropertyStore::has_differentiable_declaration(
    const std::vector<SymbolProperties::ContinuityDecl>& declarations,
    const Interval& interval) {
    for (const auto& declaration : declarations) {
        if (declaration.is_differentiable &&
            intervals_equivalent(declaration.interval, interval)) {
            return true;
        }
    }
    return false;
}

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
        auto checked = check_continuous_declarations(declarations, symbol, interval, context);
        if (!checked) return PropertyStoreResult::failure(checked.error());
        if (!checked.value()) return PropertyStoreResult::success();
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
        if (has_differentiable_declaration(declarations, interval)) {
            return PropertyStoreResult::success();
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

}
