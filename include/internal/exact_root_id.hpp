#pragma once

#include "polynomial.hpp"
#include <algorithm>

#include <cstddef>

namespace LMCAS::detail {

struct ExactRootId {
    Polynomial<Rational> polynomial;
    std::size_t index = 0;

    std::size_t hash() const {
        std::size_t seed = index;
        for (const auto& coefficient : polynomial.coeffs) {
            seed ^= coefficient.hash() + 0x9e3779b9U +
                    (seed << 6U) + (seed >> 2U);
        }
        return seed;
    }

    int compare(const ExactRootId& other) const {
        const auto common = std::min(
            polynomial.coeffs.size(), other.polynomial.coeffs.size());
        for (std::size_t position = common; position-- > 0;) {
            const auto& left = polynomial.coeffs[position];
            const auto& right = other.polynomial.coeffs[position];
            if (left < right) return -1;
            if (right < left) return 1;
        }
        if (polynomial.coeffs.size() != other.polynomial.coeffs.size()) {
            return polynomial.coeffs.size() < other.polynomial.coeffs.size()
                ? -1 : 1;
        }
        if (index == other.index) return 0;
        return index < other.index ? -1 : 1;
    }
};

} // namespace LMCAS::detail
