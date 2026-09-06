#include "unit.hpp"

#include <sstream>
#include <utility>

namespace LMCAS {
namespace {

constexpr const char* kUnitOperation = "unit";

bool valid_unit_name(const std::string& name) {
    if (name.empty()) return false;
    for (char character : name) {
        const bool valid =
            (character >= 'a' && character <= 'z') ||
            (character >= 'A' && character <= 'Z') ||
            (character >= '0' && character <= '9') || character == '_';
        if (!valid) return false;
    }
    return true;
}

} // namespace

DimensionSignature::DimensionSignature(Exponents exponents)
    : exponents_(std::move(exponents)) {
    for (auto iterator = exponents_.begin(); iterator != exponents_.end();) {
        if (iterator->second.is_zero()) {
            iterator = exponents_.erase(iterator);
        } else {
            ++iterator;
        }
    }
}

DimensionSignature DimensionSignature::base(std::string name) {
    Exponents exponents;
    exponents.emplace(std::move(name), Rational(1));
    return DimensionSignature(std::move(exponents));
}

DimensionSignature DimensionSignature::multiplied_by(
    const DimensionSignature& other) const {
    Exponents result = exponents_;
    for (const auto& entry : other.exponents_) result[entry.first] = result[entry.first] + entry.second;
    return DimensionSignature(std::move(result));
}

DimensionSignature DimensionSignature::divided_by(
    const DimensionSignature& other) const {
    Exponents result = exponents_;
    for (const auto& entry : other.exponents_) result[entry.first] = result[entry.first] - entry.second;
    return DimensionSignature(std::move(result));
}

DimensionSignature DimensionSignature::raised_to(const Rational& exponent) const {
    Exponents result;
    for (const auto& entry : exponents_) result.emplace(entry.first, entry.second * exponent);
    return DimensionSignature(std::move(result));
}

std::string DimensionSignature::to_string() const {
    if (exponents_.empty()) return "1";
    std::ostringstream output;
    bool first = true;
    for (const auto& entry : exponents_) {
        if (!first) output << '*';
        output << entry.first;
        if (entry.second != Rational(1)) output << '^' << entry.second.to_string();
        first = false;
    }
    return output.str();
}

UnitSystem::UnitSystem() {
    const auto length = DimensionSignature::base("m");
    const auto mass = DimensionSignature::base("kg");
    const auto time = DimensionSignature::base("s");
    const auto current = DimensionSignature::base("A");
    const auto temperature = DimensionSignature::base("K");
    const auto amount = DimensionSignature::base("mol");
    const auto luminous_intensity = DimensionSignature::base("cd");

    units_.emplace("m", UnitDefinition{length, Rational(1)});
    units_.emplace("kg", UnitDefinition{mass, Rational(1)});
    units_.emplace("s", UnitDefinition{time, Rational(1)});
    units_.emplace("A", UnitDefinition{current, Rational(1)});
    units_.emplace("K", UnitDefinition{temperature, Rational(1)});
    units_.emplace("mol", UnitDefinition{amount, Rational(1)});
    units_.emplace("cd", UnitDefinition{luminous_intensity, Rational(1)});
    units_.emplace("km", UnitDefinition{length, Rational(1000)});
    units_.emplace("cm", UnitDefinition{length, Rational(1, 100)});
    units_.emplace("mm", UnitDefinition{length, Rational(1, 1000)});
    units_.emplace("min", UnitDefinition{time, Rational(60)});
    units_.emplace("h", UnitDefinition{time, Rational(3600)});
    units_.emplace("g", UnitDefinition{mass, Rational(1, 1000)});
}

Result<void> UnitSystem::declare_base_unit(const std::string& name) {
    if (!valid_unit_name(name)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "unit name is invalid", kUnitOperation);
    }
    if (units_.find(name) != units_.end()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "unit name is already declared", kUnitOperation);
    }
    units_.emplace(name, UnitDefinition{DimensionSignature::base(name), Rational(1)});
    return Result<void>::success();
}

Result<void> UnitSystem::declare_derived_unit(
    const std::string& name, UnitDefinition definition) {
    if (!valid_unit_name(name) || definition.scale_to_base.is_zero() ||
        definition.scale_to_base < Rational(0)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "derived unit definition is invalid", kUnitOperation);
    }
    if (units_.find(name) != units_.end()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "unit name is already declared", kUnitOperation);
    }
    units_.emplace(name, std::move(definition));
    return Result<void>::success();
}

Result<UnitDefinition> UnitSystem::resolve(const std::string& name) const {
    const auto iterator = units_.find(name);
    if (iterator == units_.end()) {
        return Result<UnitDefinition>::failure(
            CasErrc::UnitInvalid, "unit is not declared: " + name, kUnitOperation);
    }
    return Result<UnitDefinition>::success(iterator->second);
}

Result<Rational> UnitSystem::conversion_factor(
    const std::string& source, const std::string& target) const {
    auto source_unit = resolve(source);
    if (!source_unit) return Result<Rational>::failure(source_unit.error());
    auto target_unit = resolve(target);
    if (!target_unit) return Result<Rational>::failure(target_unit.error());
    if (source_unit.value().dimension != target_unit.value().dimension) {
        return Result<Rational>::failure(
            CasErrc::DimensionMismatch,
            "unit conversion requires equal dimensions", kUnitOperation);
    }
    return Result<Rational>::success(
        source_unit.value().scale_to_base / target_unit.value().scale_to_base);
}

} // namespace LMCAS
