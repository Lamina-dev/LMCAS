/**
 * @file property_store.hpp
 * @brief PropertyStore class for mapping symbol names to declared properties.
 *
 * Provides storage and validation for domain, sign, parity, and boundedness
 * properties with implication chains and contradiction detection.
 */
#pragma once

#include "assumption.hpp"
#include "interval.hpp"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <vector>
#include <stdexcept>

namespace lamina {

/// Hash functor for Sign enum to allow use in unordered_set.
struct SignHash {
    std::size_t operator()(Sign s) const noexcept {
        return std::hash<int>{}(static_cast<int>(s));
    }
};

/**
 * @brief Stores and validates mathematical properties for named symbols.
 *
 * Supports domain hierarchy implication, sign implication, contradiction
 * detection, and cross-constraint checking between domains and signs.
 */
class LAMINA_API PropertyStore {
public:
    /// Declare domain for a symbol. Throws std::invalid_argument on contradiction.
    void declare_domain(const std::string& symbol, Domain domain);

    /// Declare sign for a symbol. Throws std::invalid_argument on contradiction.
    void declare_sign(const std::string& symbol, Sign sign);

    /// Declare parity for a symbol. Throws std::invalid_argument on contradiction.
    void declare_parity(const std::string& symbol, Parity parity);

    /// Declare boundedness for a symbol, optionally with interval bounds.
    void declare_bounded(const std::string& symbol, Boundedness bounded,
                         std::optional<Interval> bounds = std::nullopt);

    /// Query the most specific domain for a symbol (default: Complex).
    Domain get_domain(const std::string& symbol) const;

    /// Query all signs stored for a symbol.
    std::unordered_set<Sign, SignHash> get_signs(const std::string& symbol) const;

    /// Query parity for a symbol (default: Unknown).
    Parity get_parity(const std::string& symbol) const;

    /// Query boundedness for a symbol (default: Unknown).
    Boundedness get_boundedness(const std::string& symbol) const;

    /// Get stored interval bounds (if any).
    std::optional<Interval> get_bounds(const std::string& symbol) const;

    /// Check if symbol has a specific sign (explicit or implied).
    bool has_sign(const std::string& symbol, Sign sign) const;

    /// Check if symbol has at least the given domain specificity.
    bool has_domain(const std::string& symbol, Domain domain) const;

private:
    struct SymbolProperties {
        Domain most_specific_domain = Domain::Complex;
        std::unordered_set<Sign, SignHash> signs;
        Parity parity = Parity::Unknown;
        Boundedness boundedness = Boundedness::Unknown;
        std::optional<Interval> bounds;
    };

    std::unordered_map<std::string, SymbolProperties> properties_;

    /// Check domain-sign cross-constraints. Throws on contradiction.
    void check_domain_sign_consistency(const std::string& symbol,
                                       const SymbolProperties& props,
                                       Domain new_domain);

    /// Check sign contradiction pairs. Throws on contradiction.
    void check_sign_contradiction(const std::string& symbol,
                                  const std::unordered_set<Sign, SignHash>& existing,
                                  Sign new_sign);

    /// Get all signs implied by a given sign.
    std::unordered_set<Sign, SignHash> get_implied_signs(Sign sign) const;

    /// Get ancestor domains (less specific) for a given domain.
    std::vector<Domain> get_ancestor_domains(Domain domain) const;

    /// Get the specificity level of a domain (higher = more specific).
    int domain_specificity(Domain domain) const;
};

} // namespace lamina
