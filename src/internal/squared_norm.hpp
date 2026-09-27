#pragma once

#include "internal/symbolic_ast.hpp"
#include <array>

namespace LMCAS::detail {

/** Recognize real-numeric norm forms without evaluating their squares. */
std::array<const PowerNode*, 2> squared_norm_terms(const SymbolicNode& root);

}
