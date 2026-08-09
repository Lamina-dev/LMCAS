#pragma once

#include "symbolic_ast.hpp"

#include <memory>
#include <string>

namespace lamina {

// Internal node-level traversal used by algorithms that already operate on AST nodes.
bool depends_on_var(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable);

} // namespace lamina
