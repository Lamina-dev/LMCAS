#pragma once
#include "symbolic.hpp"
#include <string>
#include <memory>

namespace lamina {

std::shared_ptr<SymbolicExpr> volume_of_revolution_x(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

std::shared_ptr<SymbolicExpr> arc_length_x(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

std::shared_ptr<SymbolicExpr> volume_of_revolution_y(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

std::shared_ptr<SymbolicExpr> arc_length_y(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

}
