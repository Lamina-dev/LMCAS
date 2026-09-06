/**
 * @file assumption.hpp
 * @brief Core enumerations and forward declarations for the Assumption System.
 *
 * Defines Domain, Sign, Parity, Boundedness, and Tribool enumerations used
 * throughout the assumption-based reasoning layers (PropertyStore, RelationStore,
 * InferenceEngine, AssumptionContext, QueryInterface).
 */
#pragma once

namespace LMCAS {

/// Mathematical domain hierarchy (least specific → most specific).
/// Complex ⊃ Real ⊃ Algebraic ⊃ Rational ⊃ Integer ⊃ Natural ⊃ PositiveInt
enum class Domain {
    Complex,      ///< Complex numbers (default, least specific)
    Real,         ///< Real numbers
    Algebraic,    ///< Algebraic numbers (roots of polynomials with rational coefficients)
    Rational,     ///< Rational numbers
    Integer,      ///< Integers
    Natural,      ///< Non-negative integers (0, 1, 2, ...)
    PositiveInt   ///< Positive integers (1, 2, 3, ...)
};

/// Sign classification for symbols and expressions.
enum class Sign {
    Positive,     ///< > 0
    Negative,     ///< < 0
    NonNegative,  ///< >= 0
    NonPositive,  ///< <= 0
    Zero,         ///< == 0
    NonZero       ///< != 0
};

/// Integer parity classification.
enum class Parity {
    Even,         ///< Even integer
    Odd,          ///< Odd integer
    Unknown       ///< Parity not determined
};

/// Boundedness classification for symbols.
enum class Boundedness {
    Bounded,      ///< Symbol has finite bounds
    Unbounded,    ///< Symbol is unbounded
    Unknown       ///< Boundedness not determined
};

/// Tri-state logic value for property queries.
enum class Tribool {
    True,         ///< Property definitely holds
    False,        ///< Property definitely does not hold
    Unknown       ///< Property cannot be determined
};

/// Monotonicity classification for functions over intervals.
enum class Monotonicity {
    Increasing,      ///< Strictly increasing
    Decreasing,      ///< Strictly decreasing
    NonDecreasing,   ///< Non-decreasing (weakly increasing)
    NonIncreasing,   ///< Non-increasing (weakly decreasing)
    Unknown          ///< Monotonicity not determined
};

/// Matrix definiteness classification.
enum class Definiteness {
    PositiveDefinite,      ///< All eigenvalues positive
    PositiveSemiDefinite,  ///< All eigenvalues non-negative
    NegativeDefinite,      ///< All eigenvalues negative
    NegativeSemiDefinite,  ///< All eigenvalues non-positive
    Indefinite,            ///< Mixed eigenvalue signs
    Unknown                ///< Definiteness not determined
};

/// Finiteness classification for limits and sequences.
enum class Finiteness {
    Finite,    ///< Has a finite value or limit
    Divergent, ///< Diverges to infinity or oscillates
    Unknown    ///< Finiteness not determined
};

// Forward declarations for assumption system components
class PropertyStore;
class RelationStore;
class InferenceEngine;
class AssumptionContext;
class QueryInterface;

} // namespace LMCAS
