#pragma once
#include "test_common.hpp"
#include "interval.hpp"
#include <random>
#include <sstream>
#include <cmath>
#include <limits>

using namespace LMCAS;

inline Interval make_closed_interval(double lo, double hi) {
    auto lower_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto upper_val = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::closed(lower_val);
    iv.upper = Endpoint::closed(upper_val);
    return iv;
}

inline Interval make_open_interval(double lo, double hi) {
    auto lo_expr = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(lo)));
    auto hi_expr = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(hi)));
    Interval iv;
    iv.lower = Endpoint::open(lo_expr);
    iv.upper = Endpoint::open(hi_expr);
    return iv;
}
