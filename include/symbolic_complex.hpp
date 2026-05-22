#pragma once
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {

struct ComplexSymbolic {
    std::shared_ptr<SymbolicExpr> real;
    std::shared_ptr<SymbolicExpr> imag;
};

inline ComplexSymbolic make_complex(
    std::shared_ptr<SymbolicExpr> real,
    std::shared_ptr<SymbolicExpr> imag
) {
    return ComplexSymbolic{real, imag};
}

ComplexSymbolic complex_add(const ComplexSymbolic& a, const ComplexSymbolic& b);
ComplexSymbolic complex_sub(const ComplexSymbolic& a, const ComplexSymbolic& b);
ComplexSymbolic complex_mul(const ComplexSymbolic& a, const ComplexSymbolic& b);
ComplexSymbolic complex_div(const ComplexSymbolic& a, const ComplexSymbolic& b);

ComplexSymbolic complex_conj(const ComplexSymbolic& z);

std::shared_ptr<SymbolicExpr> complex_abs(const ComplexSymbolic& z);

std::shared_ptr<SymbolicExpr> complex_arg(const ComplexSymbolic& z);

ComplexSymbolic complex_exp_form(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta
);

ComplexSymbolic complex_trig_form(
    std::shared_ptr<SymbolicExpr> r,
    std::shared_ptr<SymbolicExpr> theta
);

std::vector<ComplexSymbolic> solve_complex_nth_root(
    std::shared_ptr<SymbolicExpr> c,
    int n
);

std::vector<ComplexSymbolic> solve_complex_quadratic(
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    std::shared_ptr<SymbolicExpr> c
);

std::shared_ptr<SymbolicExpr> complex_locus_circle(
    const ComplexSymbolic& a,
    std::shared_ptr<SymbolicExpr> r,
    const std::string& z_var = "z"
);
std::shared_ptr<SymbolicExpr> complex_locus_perpendicular_bisector(
    const ComplexSymbolic& a,
    const ComplexSymbolic& b,
    const std::string& z_var = "z"
);

}
