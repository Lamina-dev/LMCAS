#pragma once
#include "symbolic.hpp"
#include <string>
#include <memory>

namespace lamina {

// 旋转体体积（绕x轴，y=f(x)在[a,b]）
std::shared_ptr<SymbolicExpr> volume_of_revolution_x(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

// 弧长（y=f(x)在[a,b]）
std::shared_ptr<SymbolicExpr> arc_length_x(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

// 旋转体体积（绕y轴，x=f(y)在[a,b]）
std::shared_ptr<SymbolicExpr> volume_of_revolution_y(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

// 弧长（x=f(y)在[a,b]）
std::shared_ptr<SymbolicExpr> arc_length_y(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
);

}
