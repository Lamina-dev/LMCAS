#define _USE_MATH_DEFINES
#include "../include/symbolic_complex.hpp"
#include "../include/symbolic.hpp"
#include <cmath>
#include <vector>

namespace lamina {

ComplexSymbolic complex_add(const ComplexSymbolic& a, const ComplexSymbolic& b) {
    return make_complex(SymbolicExpr::add(a.real, b.real), SymbolicExpr::add(a.imag, b.imag));
}
ComplexSymbolic complex_sub(const ComplexSymbolic& a, const ComplexSymbolic& b) {
    return make_complex(SymbolicExpr::add(a.real, SymbolicExpr::multiply(SymbolicExpr::number(-1), b.real)), SymbolicExpr::add(a.imag, SymbolicExpr::multiply(SymbolicExpr::number(-1), b.imag)));
}
ComplexSymbolic complex_mul(const ComplexSymbolic& a, const ComplexSymbolic& b) {
    auto real = SymbolicExpr::add(SymbolicExpr::multiply(a.real, b.real), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a.imag, b.imag)));
    auto imag = SymbolicExpr::add(SymbolicExpr::multiply(a.real, b.imag), SymbolicExpr::multiply(a.imag, b.real));
    return make_complex(real, imag);
}
ComplexSymbolic complex_div(const ComplexSymbolic& a, const ComplexSymbolic& b) {
    auto denom = SymbolicExpr::add(SymbolicExpr::power(b.real, SymbolicExpr::number(2)), SymbolicExpr::power(b.imag, SymbolicExpr::number(2)));
    auto real = SymbolicExpr::divide(SymbolicExpr::add(SymbolicExpr::multiply(a.real, b.real), SymbolicExpr::multiply(a.imag, b.imag)), denom);
    auto imag = SymbolicExpr::divide(SymbolicExpr::add(SymbolicExpr::multiply(a.imag, b.real), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a.real, b.imag))), denom);
    return make_complex(real, imag);
}

ComplexSymbolic complex_conj(const ComplexSymbolic& z) {
    return make_complex(z.real, SymbolicExpr::multiply(SymbolicExpr::number(-1), z.imag));
}

std::shared_ptr<SymbolicExpr> complex_abs(const ComplexSymbolic& z) {
    return SymbolicExpr::sqrt(SymbolicExpr::add(SymbolicExpr::power(z.real, SymbolicExpr::number(2)), SymbolicExpr::power(z.imag, SymbolicExpr::number(2))));
}

std::shared_ptr<SymbolicExpr> complex_arg(const ComplexSymbolic& z) {
    return SymbolicExpr::atan2(z.imag, z.real);
}

ComplexSymbolic complex_exp_form(std::shared_ptr<SymbolicExpr> r, std::shared_ptr<SymbolicExpr> theta) {
    auto real = SymbolicExpr::multiply(r, SymbolicExpr::cos(theta));
    auto imag = SymbolicExpr::multiply(r, SymbolicExpr::sin(theta));
    return make_complex(real, imag);
}

ComplexSymbolic complex_trig_form(std::shared_ptr<SymbolicExpr> r, std::shared_ptr<SymbolicExpr> theta) {
    return complex_exp_form(r, theta);
}

std::vector<ComplexSymbolic> solve_complex_nth_root(std::shared_ptr<SymbolicExpr> c, int n) {

    ComplexSymbolic c_sym{SymbolicExpr::number(0), SymbolicExpr::number(0)};

    if (c->is_number()) {
        double val = c->to_numeric();
        double r = std::abs(val);
        double theta = std::arg(std::complex<double>(val, 0.0));
        std::vector<ComplexSymbolic> roots;
        for (int k = 0; k < n; ++k) {
            double root_r = std::pow(r, 1.0 / n);
            double root_theta = (theta + 2 * M_PI * k) / n;
            roots.push_back(complex_exp_form(SymbolicExpr::number(root_r), SymbolicExpr::number(root_theta)));
        }
        return roots;
    }
    return {};
}

std::vector<ComplexSymbolic> solve_complex_quadratic(std::shared_ptr<SymbolicExpr> a, std::shared_ptr<SymbolicExpr> b, std::shared_ptr<SymbolicExpr> c) {

    auto four_ac = SymbolicExpr::multiply(SymbolicExpr::number(4), SymbolicExpr::multiply(a, c));
    auto delta = SymbolicExpr::add(SymbolicExpr::power(b, SymbolicExpr::number(2)), SymbolicExpr::multiply(SymbolicExpr::number(-1), four_ac));
    auto sqrt_delta = SymbolicExpr::sqrt(delta);
    auto denom = SymbolicExpr::multiply(SymbolicExpr::number(2), a);
    auto z1 = complex_div(make_complex(SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(-1), b), sqrt_delta), SymbolicExpr::number(0)), make_complex(denom, SymbolicExpr::number(0)));
    auto z2 = complex_div(make_complex(SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(-1), b), SymbolicExpr::multiply(SymbolicExpr::number(-1), sqrt_delta)), SymbolicExpr::number(0)), make_complex(denom, SymbolicExpr::number(0)));
    return {z1, z2};
}

std::shared_ptr<SymbolicExpr> complex_locus_circle(const ComplexSymbolic& a, std::shared_ptr<SymbolicExpr> r, const std::string& z_var) {

    auto z = SymbolicExpr::variable(z_var);
    auto z_minus_a = complex_sub(make_complex(z, SymbolicExpr::number(0)), a);
    return SymbolicExpr::eq(complex_abs(z_minus_a), r);
}
std::shared_ptr<SymbolicExpr> complex_locus_perpendicular_bisector(const ComplexSymbolic& a, const ComplexSymbolic& b, const std::string& z_var) {

    auto z = SymbolicExpr::variable(z_var);
    auto z_minus_a = complex_sub(make_complex(z, SymbolicExpr::number(0)), a);
    auto z_minus_b = complex_sub(make_complex(z, SymbolicExpr::number(0)), b);
    return SymbolicExpr::eq(complex_abs(z_minus_a), complex_abs(z_minus_b));
}

}
