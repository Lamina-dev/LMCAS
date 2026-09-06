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
#include "result.hpp"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <vector>
#include <stdexcept>

namespace LMCAS {

using PropertyStoreResult = Result<void>;

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
class LMCAS_API PropertyStore {
public:
    /// Declare domain for a symbol and report contradictions.
    PropertyStoreResult declare_domain(const std::string& symbol, Domain domain);

    /** @brief Declares a domain and reports contradictions. */
    PropertyStoreResult declare_domain_checked(const std::string& symbol, Domain domain);

    /// Declare sign for a symbol and report contradictions.
    PropertyStoreResult declare_sign(const std::string& symbol, Sign sign);

    /** @brief Declares a sign and reports contradictions. */
    PropertyStoreResult declare_sign_checked(const std::string& symbol, Sign sign);

    /// Declare parity for a symbol and report contradictions.
    PropertyStoreResult declare_parity(const std::string& symbol, Parity parity);

    /** @brief Declares parity and reports contradictions. */
    PropertyStoreResult declare_parity_checked(const std::string& symbol, Parity parity);

    /// Declare boundedness for a symbol, optionally with interval bounds.
    PropertyStoreResult declare_bounded(const std::string& symbol, Boundedness bounded,
                                        std::optional<Interval> bounds = std::nullopt);

    /** @brief Declares boundedness and reports contradictions. */
    PropertyStoreResult declare_bounded_checked(
        const std::string& symbol,
        Boundedness bounded,
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

    /**
     * @brief Declare a symbol as continuous on a specified interval.
     * @param symbol Symbol name.
     * @param interval Interval on which continuity holds.
     * Contradictory overlapping declarations are returned as `CasError`.
     */
    PropertyStoreResult declare_continuous(
        const std::string& symbol, const Interval& interval);

    PropertyStoreResult declare_continuous_checked(
        const std::string& symbol,
        const Interval& interval,
        ComputationContext& context);
    PropertyStoreResult declare_continuous_checked(
        const std::string& symbol,
        const Interval& interval);

    /**
     * @brief Declare a symbol as differentiable on a specified interval.
     *
     * Differentiability implies continuity, so this also ensures the symbol
     * is recorded as continuous on the same interval.
     *
     * @param symbol Symbol name.
     * @param interval Interval on which differentiability holds.
     * Contradictory overlapping declarations are returned as `CasError`.
     */
    PropertyStoreResult declare_differentiable(
        const std::string& symbol, const Interval& interval);

    PropertyStoreResult declare_differentiable_checked(
        const std::string& symbol,
        const Interval& interval,
        ComputationContext& context);
    PropertyStoreResult declare_differentiable_checked(
        const std::string& symbol,
        const Interval& interval);

    /**
     * @brief Query whether a symbol is continuous on a given interval.
     * @param symbol Symbol name.
     * @param interval Interval to check.
     * @return true if the symbol has been declared continuous (or differentiable)
     *         on an interval that covers the queried interval.
     */
    Result<bool> is_continuous(const std::string& symbol, const Interval& interval) const;

    Result<bool> is_continuous_checked(
        const std::string& symbol,
        const Interval& interval,
        ComputationContext& context) const;
    Result<bool> is_continuous_checked(
        const std::string& symbol,
        const Interval& interval) const;

    /**
     * @brief Query whether a symbol is differentiable on a given interval.
     * @param symbol Symbol name.
     * @param interval Interval to check.
     * @return true if the symbol has been declared differentiable on an interval
     *         that covers the queried interval.
     */
    Result<bool> is_differentiable(const std::string& symbol, const Interval& interval) const;

    Result<bool> is_differentiable_checked(
        const std::string& symbol,
        const Interval& interval,
        ComputationContext& context) const;
    Result<bool> is_differentiable_checked(
        const std::string& symbol,
        const Interval& interval) const;


    /**
     * @brief Declare a symbol as transcendental.
     *
     * Sets the symbol's domain to Real and marks it transcendental.
     * Returns an error if the symbol already has an algebraic sub-domain.
     *
     * @param symbol Symbol name
     */
    PropertyStoreResult declare_transcendental(const std::string& symbol);

    /** @brief Declares transcendence and reports domain contradictions. */
    PropertyStoreResult declare_transcendental_checked(const std::string& symbol);

    /**
     * @brief Query whether a symbol has been declared transcendental.
     * @param symbol Symbol name
     * @return true if the symbol is marked transcendental
     */
    bool is_transcendental(const std::string& symbol) const;


    /**
     * @brief Declare finiteness for a symbol.
     *
     * Finite implies Bounded. Contradictions are returned as `CasError`.
     *
     * @param symbol Symbol name
     * @param f Finiteness value (Finite, Divergent, or Unknown)
     */
    PropertyStoreResult declare_finiteness(const std::string& symbol, Finiteness f);

    /** @brief Declares finiteness and reports contradictions. */
    PropertyStoreResult declare_finiteness_checked(const std::string& symbol, Finiteness f);

    /**
     * @brief Query the finiteness classification of a symbol.
     * @param symbol Symbol name
     * @return Finiteness value (default: Unknown)
     */
    Finiteness get_finiteness(const std::string& symbol) const;


    /**
     * @brief Declare matrix definiteness for a symbol.
     *
     * PositiveDefinite implies PositiveSemiDefinite. NegativeDefinite implies
     * NegativeSemiDefinite. Contradictions are returned as `CasError`.
     *
     * @param symbol Symbol name
     * @param d Definiteness value
     */
    PropertyStoreResult declare_definiteness(const std::string& symbol, Definiteness d);

    /** @brief Declares matrix definiteness and reports contradictions. */
    PropertyStoreResult declare_definiteness_checked(const std::string& symbol, Definiteness d);

    /**
     * @brief Query the definiteness classification of a symbol.
     * @param symbol Symbol name
     * @return Definiteness value (default: Unknown)
     */
    Definiteness get_definiteness(const std::string& symbol) const;


    /**
     * @brief Declare a symbol as periodic with a given period expression.
     * @param symbol Symbol name
     * @param period The period as a symbolic expression (must be non-null)
     */
    PropertyStoreResult declare_periodic(
        const std::string& symbol, const std::shared_ptr<SymbolicExpr>& period);

    /** @brief Declares periodicity and reports invalid periods. */
    PropertyStoreResult declare_periodic_checked(
        const std::string& symbol,
        const std::shared_ptr<SymbolicExpr>& period);

    /**
     * @brief Get the period expression for a symbol, if declared.
     * @param symbol Symbol name
     * @return The period expression, or std::nullopt if not periodic
     */
    std::optional<std::shared_ptr<SymbolicExpr>> get_period(const std::string& symbol) const;

    /**
     * @brief Query whether a symbol has been declared periodic.
     * @param symbol Symbol name
     * @return true if the symbol has a declared period
     */
    bool is_periodic(const std::string& symbol) const;

    /**
     * @brief Get all symbol names that have any declared properties.
     * @return Vector of symbol names with at least one non-default property.
     */
    std::vector<std::string> get_all_symbols() const;

    /// Continuity declaration (public type for serialization access).
    struct ContinuityInfo {
        Interval interval;
        bool is_differentiable;
    };

    /// Monotonicity declaration (public type for serialization access).
    struct MonotonicityInfo {
        std::string variable;
        Interval interval;
        Monotonicity type;
    };

    /**
     * @brief Get all continuity/differentiability declarations for a symbol.
     * @param symbol Symbol name.
     * @return Vector of continuity declarations.
     */
    std::vector<ContinuityInfo> get_continuity_decls(const std::string& symbol) const;

    /**
     * @brief Get all monotonicity declarations for a symbol.
     * @param symbol Symbol name.
     * @return Vector of monotonicity declarations.
     */
    std::vector<MonotonicityInfo> get_monotonicity_decls(const std::string& symbol) const;

    /**
     * @brief Declare monotonicity of a symbol with respect to a variable on an interval.
     *
     * Stores a MonotonicityDecl recording that the symbol is monotonically
     * increasing/decreasing (or non-decreasing/non-increasing) with respect to
     * the given variable over the specified interval.
     *
     * @param symbol   The symbol (expression name) being declared monotone.
     * @param variable The variable with respect to which monotonicity holds.
     * @param interval The interval on which the monotonicity property holds.
     * @param mono     The type of monotonicity (Increasing, Decreasing, etc.).
     */
    PropertyStoreResult declare_monotonicity(
        const std::string& symbol, const std::string& variable,
        const Interval& interval, Monotonicity mono);

    PropertyStoreResult declare_monotonicity_checked(
        const std::string& symbol,
        const std::string& variable,
        const Interval& interval,
        Monotonicity mono,
        ComputationContext& context);
    PropertyStoreResult declare_monotonicity_checked(
        const std::string& symbol,
        const std::string& variable,
        const Interval& interval,
        Monotonicity mono);

    /**
     * @brief Query the monotonicity of a symbol with respect to a variable on an interval.
     *
     * Searches stored MonotonicityDecl entries for the symbol. Returns the
     * monotonicity type if the queried interval is covered by a stored
     * declaration (i.e., the stored interval contains the queried interval),
     * otherwise returns Monotonicity::Unknown.
     *
     * @param symbol   The symbol to query.
     * @param variable The variable with respect to which monotonicity is queried.
     * @param interval The interval to check coverage for.
     * @return The monotonicity type if covered, or Monotonicity::Unknown.
     */
    Result<Monotonicity> get_monotonicity(
        const std::string& symbol, const std::string& variable,
        const Interval& interval) const;

    Result<Monotonicity> get_monotonicity_checked(
        const std::string& symbol,
        const std::string& variable,
        const Interval& interval,
        ComputationContext& context) const;
    Result<Monotonicity> get_monotonicity_checked(
        const std::string& symbol,
        const std::string& variable,
        const Interval& interval) const;

private:
    struct SymbolProperties {
        Domain most_specific_domain = Domain::Complex;
        std::unordered_set<Sign, SignHash> signs;
        Parity parity = Parity::Unknown;
        Boundedness boundedness = Boundedness::Unknown;
        std::optional<Interval> bounds;

        /// Whether the symbol is classified as transcendental.
        bool transcendental = false;

        /// Finiteness classification (Finite, Divergent, or Unknown).
        Finiteness finiteness = Finiteness::Unknown;

        /// Matrix definiteness classification.
        Definiteness definiteness = Definiteness::Unknown;

        /// Period expression for periodic symbols (nullopt if not periodic).
        std::optional<std::shared_ptr<SymbolicExpr>> period;

        /// Continuity declaration over an interval.
        struct ContinuityDecl {
            Interval interval;          ///< Interval on which continuity/differentiability holds
            bool is_differentiable;     ///< True if differentiable (implies continuous)
        };

        /// Declared continuity/differentiability intervals.
        std::vector<ContinuityDecl> continuity_decls;

        /// Monotonicity declaration with respect to a variable over an interval.
        struct MonotonicityDecl {
            std::string variable;       ///< Variable with respect to which monotonicity holds
            Interval interval;          ///< Interval on which monotonicity holds
            Monotonicity type;          ///< Type of monotonicity
        };

        /// Declared monotonicity intervals.
        std::vector<MonotonicityDecl> monotonicity_decls;
    };

    std::unordered_map<std::string, SymbolProperties> properties_;

    void declare_domain_unchecked(const std::string& symbol, Domain domain);
    void declare_sign_unchecked(const std::string& symbol, Sign sign);
    void declare_parity_unchecked(const std::string& symbol, Parity parity);
    void declare_bounded_unchecked(const std::string& symbol,
                                   Boundedness bounded,
                                   std::optional<Interval> bounds);
    void declare_transcendental_unchecked(const std::string& symbol);
    void declare_finiteness_unchecked(const std::string& symbol, Finiteness f);
    void declare_definiteness_unchecked(const std::string& symbol, Definiteness d);
    void declare_periodic_unchecked(const std::string& symbol,
                                    const std::shared_ptr<SymbolicExpr>& period);

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

} // namespace LMCAS
