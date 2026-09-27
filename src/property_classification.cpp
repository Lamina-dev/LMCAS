#include "internal/property_store_support.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/facts_query.hpp"
#include <algorithm>

namespace LMCAS {

using namespace property_detail;

namespace {

bool definiteness_conflicts(Definiteness existing, Definiteness incoming) {
    using D = Definiteness;
    struct Conflict {
        D existing;
        std::initializer_list<D> incoming;
    };
    const Conflict conflicts[] = {
        {D::PositiveDefinite, {D::NegativeDefinite, D::NegativeSemiDefinite, D::Indefinite}},
        {D::NegativeDefinite, {D::PositiveDefinite, D::PositiveSemiDefinite, D::Indefinite}},
        {D::PositiveSemiDefinite, {D::NegativeDefinite, D::Indefinite}},
        {D::NegativeSemiDefinite, {D::PositiveDefinite, D::Indefinite}},
        {D::Indefinite, {D::PositiveDefinite, D::NegativeDefinite,
                         D::PositiveSemiDefinite, D::NegativeSemiDefinite}}
    };
    for (const auto& conflict : conflicts) {
        if (existing != conflict.existing) {
            continue;
        }
        return std::find(conflict.incoming.begin(), conflict.incoming.end(), incoming)
            != conflict.incoming.end();
    }
    return false;
}

}


void PropertyStore::declare_transcendental_unchecked(const std::string& symbol) {
    auto& props = properties_[symbol];

    if (props.transcendental) {
        return;
    }
    if (domain_specificity(props.most_specific_domain) > domain_specificity(Domain::Real)) {
        throw std::invalid_argument(
            "Contradiction for symbol '" + symbol +
            "': cannot declare Transcendental because existing domain is " +
            domain_str(props.most_specific_domain) +
            " (Transcendental contradicts Algebraic or more specific domains)");
    }

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


void PropertyStore::declare_finiteness_unchecked(const std::string& symbol, Finiteness f) {
    auto& props = properties_[symbol];

    if (props.finiteness == f) {
        return;
    }
    if (f == Finiteness::Unknown) {
        props.finiteness = f;
        return;
    }
    if (props.finiteness != Finiteness::Unknown && props.finiteness != f) {
        throw std::invalid_argument(
            "Contradiction for symbol '" + symbol +
            "': cannot be both Finite and Divergent");
    }

    props.finiteness = f;
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


void PropertyStore::declare_definiteness_unchecked(const std::string& symbol, Definiteness d) {
    auto& props = properties_[symbol];

    if (props.definiteness == d) {
        return;
    }
    if (d == Definiteness::Unknown) {
        props.definiteness = d;
        return;
    }

    if (props.definiteness != Definiteness::Unknown) {
        if (definiteness_conflicts(props.definiteness, d)) {
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
            return; /** @brief 保留更强的正定声明。 */
        }
        if (props.definiteness == Definiteness::NegativeDefinite &&
            d == Definiteness::NegativeSemiDefinite) {
            return; /** @brief 保留更强的负定声明。 */
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


void PropertyStore::declare_periodic_unchecked(
    const std::string& symbol, const std::string& variable,
    const SymbolicExpr& period) {
    properties_[symbol].periods.insert_or_assign(variable, period);
}

PropertyStoreResult PropertyStore::declare_periodic(
    const std::string& symbol, const std::string& variable,
    const SymbolicExpr& period) {
    return declare_periodic_checked(symbol, variable, period);
}

PropertyStoreResult PropertyStore::declare_periodic_checked(
    const std::string& symbol, const std::string& variable,
    const SymbolicExpr& period) {
    if (symbol.empty() || variable.empty() || !detail::node(period)) {
        return PropertyStoreResult::failure(CasErrc::InvalidArgument,
            "symbol, independent variable, and period must be nonempty", "declare_periodic");
    }
    ComputationContext context;
    auto defined = detail::query_definedness(
        detail::node(period), detail::no_facts(), Domain::Real, context);
    if (!defined) return PropertyStoreResult::failure(defined.error());
    auto positive = detail::query_positive_value(
        detail::node(period), detail::no_facts(), context);
    if (!positive) return PropertyStoreResult::failure(positive.error());
    if (defined.value() == Tribool::False || positive.value() == Tribool::False) {
        return PropertyStoreResult::failure(CasErrc::InvalidArgument,
            "period must be a positive defined real value", "declare_periodic");
    }
    return checked_property_update(*this, symbol, "declare_periodic",
        [&](PropertyStore& candidate) {
            candidate.declare_periodic_unchecked(symbol, variable, period);
        });
}

const PropertyStore::PeriodDeclarations& PropertyStore::get_period_decls(
    const std::string& symbol) const {
    static const PeriodDeclarations empty;
    const auto it = properties_.find(symbol);
    return it == properties_.end() ? empty : it->second.periods;
}

std::optional<SymbolicExpr> PropertyStore::get_period(
    const std::string& symbol, const std::string& variable) const {
    const auto& declarations = get_period_decls(symbol);
    const auto it = declarations.find(variable);
    if (it == declarations.end()) return std::nullopt;
    return it->second;
}

bool PropertyStore::is_periodic(const std::string& symbol, const std::string& variable) const {
    return get_period(symbol, variable).has_value();
}


std::vector<std::string> PropertyStore::get_all_symbols() const {
    std::vector<std::string> symbols;
    symbols.reserve(properties_.size());
    for (const auto& [name, _] : properties_) {
        symbols.push_back(name);
    }
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

}
