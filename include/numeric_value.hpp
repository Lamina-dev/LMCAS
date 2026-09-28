#pragma once

namespace LMCAS {

enum class NumericStatus { Finite, PositiveInfinity, NegativeInfinity };

struct ApproxReal {
    double value = 0.0;
    // Certified bound on the true value's distance from value; +infinity means unknown.
    double absolute_error = 0.0;
    NumericStatus status = NumericStatus::Finite;

    bool is_finite() const noexcept { return status == NumericStatus::Finite; }
};

struct ApproxComplex {
    ApproxReal real;
    ApproxReal imag;

    bool is_finite() const noexcept {
        return real.is_finite() && imag.is_finite();
    }
};

}
