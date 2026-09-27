#pragma once

#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "internal/symbolic_ast.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace LMCAS;

static std::shared_ptr<const SymbolicNode> make_multiply(
    std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    return LMCAS::detail::make_node<MultiplyNode>(std::move(ops));
}
