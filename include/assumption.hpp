/**
 * @file assumption.hpp
 * @brief Core enumerations and forward declarations for the Assumption System.
 *
 * Defines Domain, Sign, Parity, Boundedness, and Tribool enumerations used
 * throughout the assumption-based reasoning layers (PropertyStore, RelationStore,
 * InferenceEngine, AssumptionContext, QueryInterface).
 */
#pragma once

namespace lamina {

/// Mathematical domain hierarchy (least specific → most specific).
/// Complex ⊃ Real ⊃ Rational ⊃ Integer ⊃ Natural ⊃ PositiveInt
enum class Domain {
    Complex,      ///< Complex numbers (default, least specific)
    Real,         ///< Real numbers
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

// Forward declarations for assumption system components
class PropertyStore;
class RelationStore;
class InferenceEngine;
class AssumptionContext;
class QueryInterface;

} // namespace lamina
