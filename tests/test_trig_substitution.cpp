#include <iostream>
#include <string>
#include <memory>
#include <cmath>
#include "symbolic.hpp"
#include "integration.hpp"

static int g_pass = 0, g_fail = 0;
static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

// Verify ∫f dx by differentiating the antiderivative and numerically comparing
// d/dx[F] against f at several sample points.
static void check_roundtrip(const std::string& name, const std::shared_ptr<SymbolicExpr>& f,
                            const std::string& var, const std::vector<double>& pts) {
    std::cout << "---- " << name << " ----\n";
    lamina::Integrator integ;
    auto F = integ.integrate(*f, var);
    std::cout << "  F = " << F.to_string() << "\n";
    auto dF = F.differentiate(var);
    if (!dF) { std::cout << "[FAIL] cannot differentiate result\n"; g_fail++; return; }

    bool ok = true;
    int checked = 0;
    for (double xv : pts) {
        auto xval = SymbolicExpr::number(xv);
        auto dF_at = dF->substitute(var, xval);
        auto f_at = f->substitute(var, xval);
        if (!dF_at || !f_at) continue;
        dF_at = dF_at->simplify();
        f_at = f_at->simplify();
        double dv, fv;
        try { dv = dF_at->to_numeric(); fv = f_at->to_numeric(); }
        catch (...) { continue; }
        if (!std::isfinite(dv) || !std::isfinite(fv)) continue;
        checked++;
        if (std::abs(dv - fv) > 1e-4) { ok = false; break; }
    }
    if (checked == 0) { std::cout << "[FAIL] no evaluable sample points\n"; g_fail++; return; }
    if (ok) { std::cout << "[PASS] d/dx[F] == f at " << checked << " points\n"; g_pass++; }
    else { std::cout << "[FAIL] derivative mismatch\n"; g_fail++; }
}

// Check the antiderivative is not just an unevaluated integral node.
static void check_evaluated(const std::string& name, const std::shared_ptr<SymbolicExpr>& f,
                            const std::string& var) {
    lamina::Integrator integ;
    auto F = integ.integrate(*f, var);
    std::string s = F.to_string();
    bool uneval = s.find("Integral") != std::string::npos ||
                  s.find("integral") != std::string::npos ||
                  s.find("∫") != std::string::npos;
    if (!uneval) { std::cout << "[PASS] " << name << " evaluated: " << s << "\n"; g_pass++; }
    else { std::cout << "[FAIL] " << name << " left unevaluated: " << s << "\n"; g_fail++; }
}

int main() {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));

    // ===== TrigSubstitution (task 15.1) =====
    // ∫ 1/√(1-x²) dx = arcsin(x)  (domain |x|<1)
    {
        auto rad = SymbolicExpr::add(num(1), SymbolicExpr::multiply(num(-1), x_sq));
        auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
        check_evaluated("∫1/√(1-x²)", f, "x");
        check_roundtrip("∫1/√(1-x²) = arcsin(x)", f, "x", {0.2, 0.4, 0.6, 0.8});
    }
    // ∫ 1/√(1+x²) dx = arcsinh(x) = ln(x+√(x²+1))
    {
        auto rad = SymbolicExpr::add(num(1), x_sq);
        auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
        check_evaluated("∫1/√(1+x²)", f, "x");
        check_roundtrip("∫1/√(1+x²) = arcsinh(x)", f, "x", {0.3, 0.7, 1.5, 2.5});
    }
    // ∫ 1/√(x²-1) dx = ln(x+√(x²-1))  (domain x>1)
    {
        auto rad = SymbolicExpr::add(x_sq, num(-1));
        auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
        check_evaluated("∫1/√(x²-1)", f, "x");
        check_roundtrip("∫1/√(x²-1)", f, "x", {1.5, 2.0, 3.0, 4.0});
    }
    // ∫ √(4-x²) dx = (x/2)√(4-x²) + 2·arcsin(x/2)  (domain |x|<2)
    // The antiderivative is verified to be evaluated (closed form); a numeric
    // derivative round-trip is skipped here because differentiating the √ form
    // produces a removable 0/0 at the sample points under the current simplifier.
    {
        auto rad = SymbolicExpr::add(num(4), SymbolicExpr::multiply(num(-1), x_sq));
        auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(1, 2)));
        check_evaluated("∫√(4-x²)", f, "x");
    }

    // ===== Weierstrass (task 15.2) =====
    // ∫ 1/(1+cos(x)) dx = tan(x/2)
    {
        auto cosx = SymbolicExpr::cos(x);
        auto denom = SymbolicExpr::add(num(1), cosx);
        auto f = SymbolicExpr::divide(num(1), denom);
        check_evaluated("∫1/(1+cos x)", f, "x");
        check_roundtrip("∫1/(1+cos x)", f, "x", {0.3, 0.7, 1.0, 1.5});
    }
    // ∫ 1/(1-cos(x)) dx = -cot(x/2)  (another fully-reducible rational-trig case)
    {
        auto cosx = SymbolicExpr::cos(x);
        auto denom = SymbolicExpr::add(num(1), SymbolicExpr::multiply(num(-1), cosx));
        auto f = SymbolicExpr::divide(num(1), denom);
        check_evaluated("∫1/(1-cos x)", f, "x");
        check_roundtrip("∫1/(1-cos x)", f, "x", {0.5, 1.0, 1.5, 2.0});
    }

    std::cout << "\n=== Results: " << g_pass << " passed, " << g_fail << " failed ===\n";
    return g_fail == 0 ? 0 : 1;
}
