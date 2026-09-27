#pragma once

#include "interval.hpp"
#include "internal/exact_algebraic.hpp"
#include <optional>
#include <utility>

namespace LMCAS::detail {

struct ComparableEndpoint {
    int infinity = 0;
    Rational rational{};
    Rational radical_coefficient{};
    Rational radicand{};
    std::optional<LMCAS::detail::ExactRealAlgebraic> algebraic;

    ComparableEndpoint(int infinity_value = 0,
                       Rational rational_value = {},
                       Rational coefficient_value = {},
                       Rational radicand_value = {})
        : infinity(infinity_value),
          rational(std::move(rational_value)),
          radical_coefficient(std::move(coefficient_value)),
          radicand(std::move(radicand_value)) {}
};

struct CheckedInterval {
    Interval interval;
    ComparableEndpoint lower;
    ComparableEndpoint upper;
};

Result<Rational> exact_double_rational(double value, ComputationContext& context,
                                       const std::string& operation);
Result<ComparableEndpoint> comparable_endpoint(
    const Endpoint& endpoint, ComputationContext& context,
    const std::string& operation);
Result<int> compare_comparable(const ComparableEndpoint& left,
                               const ComparableEndpoint& right,
                               ComputationContext& context);

}
