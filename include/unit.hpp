#pragma once

#include <map>
#include <string>
#include <unordered_map>

#include "lamina_export.hpp"
#include "rational.hpp"
#include "result.hpp"

namespace lamina {

/** A canonical product of base dimensions raised to exact rational powers. */
class LAMINA_API DimensionSignature {
public:
    using Exponents = std::map<std::string, Rational>;

    DimensionSignature() = default;
    explicit DimensionSignature(Exponents exponents);

    static DimensionSignature base(std::string name);

    bool is_dimensionless() const noexcept { return exponents_.empty(); }
    const Exponents& exponents() const noexcept { return exponents_; }

    DimensionSignature multiplied_by(const DimensionSignature& other) const;
    DimensionSignature divided_by(const DimensionSignature& other) const;
    DimensionSignature raised_to(const Rational& exponent) const;
    std::string to_string() const;

    bool operator==(const DimensionSignature& other) const noexcept {
        return exponents_ == other.exponents_;
    }
    bool operator!=(const DimensionSignature& other) const noexcept {
        return !(*this == other);
    }

private:
    Exponents exponents_;
};

struct UnitDefinition {
    DimensionSignature dimension;
    Rational scale_to_base = Rational(1);
};

/** Owns unit declarations for one computation context. */
class LAMINA_API UnitSystem {
public:
    UnitSystem();

    Result<void> declare_base_unit(const std::string& name);
    Result<void> declare_derived_unit(const std::string& name,
                                      UnitDefinition definition);
    Result<UnitDefinition> resolve(const std::string& name) const;
    Result<Rational> conversion_factor(const std::string& source,
                                       const std::string& target) const;

private:
    std::unordered_map<std::string, UnitDefinition> units_;
};

} // namespace lamina
