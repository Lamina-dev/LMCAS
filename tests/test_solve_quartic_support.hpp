#pragma once
#include "test_common.hpp"
#include "solve_polynomial.hpp"
#include <cmath>
#include <algorithm>
#include <random>
#include <sstream>
#include <string>
#include <stdexcept>

using namespace LMCAS;

inline double eval_quartic(double a, double b, double c, double d, double e, double x) {
    return a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
}
