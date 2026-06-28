/**
 * @file numerical_integration.cpp
 * @brief 数值积分实现。
 */
#include "numerical_integration.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <cmath>
#include <functional>

namespace lamina {

std::shared_ptr<SymbolicExpr> quadrature_simpson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n) {
    
    if (n <= 0) return SymbolicExpr::number(0);
    if (n % 2 != 0) n++; // Simpson's rule requires even n
    
    auto h = SymbolicExpr::divide(
        SymbolicExpr::add(b, SymbolicExpr::multiply(a, SymbolicExpr::number(-1))),
        SymbolicExpr::number(n)
    );
    
    auto sum = SymbolicExpr::add(f->substitute(var, a), f->substitute(var, b));
    
    for (int i = 1; i < n; i++) {
        auto xi = SymbolicExpr::add(a, SymbolicExpr::multiply(h, SymbolicExpr::number(i)));
        auto f_xi = f->substitute(var, xi);
        
        int coeff = (i % 2 == 1) ? 4 : 2;
        sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(SymbolicExpr::number(coeff), f_xi));
    }
    
    auto res = SymbolicExpr::multiply(SymbolicExpr::divide(h, SymbolicExpr::number(3)), sum);
    return res->simplify();
}

std::shared_ptr<SymbolicExpr> quadrature_gaussian(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n) {
    
    // Convert integral from [a, b] to [-1, 1]
    // x = 0.5 * (b - a) * t + 0.5 * (a + b)
    // dx = 0.5 * (b - a) dt
    
    auto diff = SymbolicExpr::add(b, SymbolicExpr::multiply(a, SymbolicExpr::number(-1)));
    auto sum_ab = SymbolicExpr::add(a, b);
    auto half_diff = SymbolicExpr::multiply(SymbolicExpr::number(0.5), diff);
    auto half_sum = SymbolicExpr::multiply(SymbolicExpr::number(0.5), sum_ab);
    
    auto evaluate_at_t = [&](const std::shared_ptr<SymbolicExpr>& t) {
        auto x = SymbolicExpr::add(SymbolicExpr::multiply(half_diff, t), half_sum);
        return f->substitute(var, x);
    };
    
    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    std::vector<std::shared_ptr<SymbolicExpr>> weights;
    
    if (n <= 1) {
        roots = { SymbolicExpr::number(0) };
        weights = { SymbolicExpr::number(2) };
    } else if (n == 2) {
        auto r = SymbolicExpr::divide(SymbolicExpr::number(1), SymbolicExpr::sqrt(SymbolicExpr::number(3)));
        auto neg_r = SymbolicExpr::multiply(r, SymbolicExpr::number(-1));
        roots = { neg_r, r };
        weights = { SymbolicExpr::number(1), SymbolicExpr::number(1) };
    } else if (n == 3) {
        auto r = SymbolicExpr::sqrt(SymbolicExpr::divide(SymbolicExpr::number(3), SymbolicExpr::number(5)));
        auto neg_r = SymbolicExpr::multiply(r, SymbolicExpr::number(-1));
        roots = { neg_r, SymbolicExpr::number(0), r };
        weights = { 
            SymbolicExpr::divide(SymbolicExpr::number(5), SymbolicExpr::number(9)),
            SymbolicExpr::divide(SymbolicExpr::number(8), SymbolicExpr::number(9)),
            SymbolicExpr::divide(SymbolicExpr::number(5), SymbolicExpr::number(9))
        };
    } else {
        // Fallback to Simpson for higher symbolic degrees due to analytical complexity of Legendre roots
        return quadrature_simpson(f, var, a, b, n * 2);
    }
    
    auto integral_sum = SymbolicExpr::number(0);
    for (size_t i = 0; i < roots.size(); i++) {
        auto term = SymbolicExpr::multiply(weights[i], evaluate_at_t(roots[i]));
        integral_sum = SymbolicExpr::add(integral_sum, term);
    }
    
    auto res = SymbolicExpr::multiply(half_diff, integral_sum);
    return res->simplify();
}

namespace {

// 在数值点处求 f 的 double 值；失败抛出。
double eval_f(const std::shared_ptr<SymbolicExpr>& f, const std::string& var, double xval) {
    auto x = SymbolicExpr::number(xval);
    auto fx = f->substitute(var, x);
    if (!fx) throw std::runtime_error("eval_f: substitution failed");
    fx = fx->simplify();
    return fx->to_numeric();
}

// 单段辛普森
double simpson_seg(const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
                   double a, double b, double fa, double fb, double fm) {
    double h = b - a;
    return (h / 6.0) * (fa + 4.0 * fm + fb);
}

double adaptive_rec(const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
                    double a, double b, double fa, double fb, double fm,
                    double whole, double tol, int depth) {
    double m = 0.5 * (a + b);
    double lm = 0.5 * (a + m);
    double rm = 0.5 * (m + b);
    double flm = eval_f(f, var, lm);
    double frm = eval_f(f, var, rm);
    double left = simpson_seg(f, var, a, m, fa, fm, flm);
    double right = simpson_seg(f, var, m, b, fm, fb, frm);
    if (depth <= 0 || std::abs(left + right - whole) <= 15.0 * tol) {
        return left + right + (left + right - whole) / 15.0;
    }
    return adaptive_rec(f, var, a, m, fa, fm, flm, left, tol / 2.0, depth - 1)
         + adaptive_rec(f, var, m, b, fm, fb, frm, right, tol / 2.0, depth - 1);
}

} // anonymous namespace

std::shared_ptr<SymbolicExpr> adaptive_simpson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    double tol) {
    if (!f || !a || !b) return nullptr;
    double av, bv;
    try {
        av = a->simplify()->to_numeric();
        bv = b->simplify()->to_numeric();
    } catch (...) {
        return nullptr;
    }
    try {
        double fa = eval_f(f, var, av);
        double fb = eval_f(f, var, bv);
        double m = 0.5 * (av + bv);
        double fm = eval_f(f, var, m);
        double whole = simpson_seg(f, var, av, bv, fa, fb, fm);
        double result = adaptive_rec(f, var, av, bv, fa, fb, fm, whole, tol, 50);
        return SymbolicExpr::number(result);
    } catch (...) {
        return nullptr;
    }
}

std::shared_ptr<SymbolicExpr> numerical_integrate(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n) {
    return quadrature_simpson(f, var, a, b, n);
}

} // namespace lamina
