#pragma once

#include "symbolic_ast.hpp"

#include <memory>
#include <string>

namespace lamina {

/** Returns whether a variable occurs free in an expression tree. */
bool expression_depends_on_variable(
    const std::shared_ptr<const SymbolicNode>& expression,
    const std::string& variable);

} // namespace lamina
