#include "internal/equivalence_support.hpp"

#include <utility>

namespace LMCAS {
namespace equivalence_detail {

Result<void> validate_eqv_options(const EqvOptions& options) {
    if (options.budget.max_rewrite_steps == 0 ||
        options.budget.max_rewrite_depth == 0 ||
        options.budget.max_node_growth_factor == 0) {
        return Result<void>::failure(
            CasErrc::ResourceLimit,
            "equivalence rewrite budget exhausted before normalization",
            kEquivalentOperation);
    }
    return Result<void>::success();
}

bool exact_integer_node(const std::shared_ptr<const SymbolicNode>& node,
                        int expected) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return false;
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()) == BigInt(expected);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()) == Rational(expected);
    }
    return false;
}

}
namespace {

Result<EqvProfile> eqv_profile_failure(std::string message) {
    return Result<EqvProfile>::failure(
        CasErrc::UnsupportedExpression,
        std::move(message),
        kEquivalentProfileOperation);
}

}

Result<EqvProfile> eqv_profile_from_name(const std::string& name) {
    if (name == "Core") {
        return Result<EqvProfile>::success(EqvProfile::Core);
    }
    if (name == "Trig-Basic") {
        return Result<EqvProfile>::success(EqvProfile::TrigBasic);
    }
    if (name == "ExpLog-Basic") {
        return Result<EqvProfile>::success(EqvProfile::ExpLogBasic);
    }
    return eqv_profile_failure("unsupported equivalence profile: " + name);
}

Result<void> set_eqv_profile(EqvOptions& options,
                             const std::string& name) {
    auto profile = eqv_profile_from_name(name);
    if (!profile) {
        return Result<void>::failure(profile.error());
    }
    options.profile = profile.value();
    return Result<void>::success();
}

Result<void> set_eqv_budget(EqvOptions& options,
                            std::size_t steps,
                            std::size_t depth,
                            std::size_t growth) {
    EqvOptions candidate = options;
    candidate.budget.max_rewrite_steps = steps;
    candidate.budget.max_rewrite_depth = depth;
    candidate.budget.max_node_growth_factor = growth;
    auto valid = equivalence_detail::validate_eqv_options(candidate);
    if (!valid) {
        return valid;
    }
    options = candidate;
    return Result<void>::success();
}

}
