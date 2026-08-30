#pragma once

#include "symbolic_ast.hpp"

#include <set>
#include <string>

namespace lamina {

/** Returns every variable name that occurs free in an expression tree. */
LAMINA_API std::set<std::string> free_variables(
    const detail::SymbolicNodePtr& expression);

/** Returns whether a variable occurs free in an expression tree. */
LAMINA_API bool expression_depends_on_variable(
    const detail::SymbolicNodePtr& expression,
    const std::string& variable);

/** Capture-avoiding substitution of free occurrences of a variable. */
LAMINA_API detail::SymbolicNodePtr substitute_free(
    const detail::SymbolicNodePtr& expression,
    const std::string& variable,
    const detail::SymbolicNodePtr& replacement);

} // namespace lamina
