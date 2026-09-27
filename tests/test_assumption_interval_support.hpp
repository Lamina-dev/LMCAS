#pragma once

#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "interval.hpp"
#include "internal/symbolic_ast.hpp"
#include <memory>

using namespace LMCAS;

inline std::shared_ptr<const SymbolicNode> make_num(double value) {
    return LMCAS::detail::make_node<NumberNode>(
        static_cast<lmmc_real_t>(value));
}
