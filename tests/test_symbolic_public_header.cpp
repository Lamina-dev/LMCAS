#include "symbolic.hpp"

#ifdef LAMINA_INTERNAL_AST_INCLUDED
#error "symbolic.hpp must not include the AST implementation"
#endif

static_assert(sizeof(SymbolicExpr) > 0, "SymbolicExpr must be a complete value type");

int main() {
    auto one = SymbolicExpr::number(1);
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(x, one);

    if (!one || !x || !expression) return 1;
    if (!one->is_number() || one->get_int() != 1) return 1;
    return expression->to_string().empty() ? 1 : 0;
}
