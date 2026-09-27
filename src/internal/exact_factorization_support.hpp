#pragma once

#include "bigint.hpp"

namespace LMCAS::detail {

inline BigInt symmetric_remainder(const BigInt& value,
                                  const BigInt& modulus) {
    if (modulus == BigInt(0)) return value;
    auto remainder = value % modulus;
    if (remainder < BigInt(0)) remainder = remainder + modulus;
    if (remainder > modulus / BigInt(2)) remainder = remainder - modulus;
    return remainder;
}

}
