#include "internal/integration_support.hpp"

namespace lamina {

namespace {

// Build sin(c*var) / cos(c*var) where c is a non-zero integer scale.
std::shared_ptr<SymbolicExpr> make_sin_scaled(long long c, const std::string& var) {
    auto v = SymbolicExpr::variable(var);
    if (c == 1) return SymbolicExpr::sin(v);
    auto cx = SymbolicExpr::multiply(SymbolicExpr::number(static_cast<long long>(c)), v);
    return SymbolicExpr::sin(cx);
}

std::shared_ptr<SymbolicExpr> make_cos_scaled(long long c, const std::string& var) {
    auto v = SymbolicExpr::variable(var);
    if (c == 1) return SymbolicExpr::cos(v);
    auto cx = SymbolicExpr::multiply(SymbolicExpr::number(static_cast<long long>(c)), v);
    return SymbolicExpr::cos(cx);
}

std::shared_ptr<SymbolicExpr> trig_rational(long long num, long long den) {
    return SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
}

// Binomial coefficient C(n,k) for small n; safe for n up to ~62.
long long binomial_ll(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    long long c = 1;
    for (int i = 0; i < k; ++i) {
        c = c * static_cast<long long>(n - i) / static_cast<long long>(i + 1);
    }
    return c;
}

// Match a single-argument FunctionNode whose argument is the integration
// variable itself. Returns 0 (sin), 1 (cos), 2 (tan), 3 (sec), or -1 on
// no match for this strategy.
int trig_match_of_var(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!fn) return -1;
    if (fn->arguments().size() != 1) return -1;
    auto v = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0]);
    if (!v || v->name() != var) return -1;
    using FT = FunctionNode::FuncType;
    switch (fn->type()) {
        case FT::Sin: return 0;
        case FT::Cos: return 1;
        case FT::Tan: return 2;
        case FT::Sec: return 3;
        default: return -1;
    }
}

// Try to extract a non-negative integer exponent.
bool trig_extract_nonneg_int(const std::shared_ptr<const SymbolicNode>& exp_node, int& n_out) {
    auto e = lamina::detail::expression_from_node(exp_node);
    auto simp = e.simplify();
    if (!simp || !simp->is_int()) return false;
    int n = simp->get_int();
    if (n < 0) return false;
    n_out = n;
    return true;
}

// Match a single factor of the form trig(var) or trig(var)^k where trig is
// sin/cos/tan/sec and k is a non-negative integer. Returns true on success
// with the kind (0..3) and the integer power.
bool trig_extract_factor(const std::shared_ptr<const SymbolicNode>& node,
                         const std::string& var,
                         int& kind_out, int& power_out) {
    int kind = trig_match_of_var(node, var);
    if (kind >= 0) {
        kind_out = kind;
        power_out = 1;
        return true;
    }
    auto pn = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pn) return false;
    int base_kind = trig_match_of_var(pn->base(), var);
    if (base_kind < 0) return false;
    int p = 0;
    if (!trig_extract_nonneg_int(pn->exponent(), p)) return false;
    kind_out = base_kind;
    power_out = p;
    return true;
}

} // anonymous namespace

bool TrigCombinationStrategy::extract_sin_cos_powers(
    const SymbolicExpr& expr, const std::string& var, int& m_out, int& n_out) {

    int m = 0, n = 0;
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        factors = mul->operands();
    } else {
        factors.push_back(lamina::detail::node(expr));
    }

    for (const auto& f : factors) {
        int kind = -1, p = 0;
        if (!trig_extract_factor(f, var, kind, p)) return false;
        // Only sin/cos contribute to (m,n); tan/sec disqualify this form.
        if (kind == 0) m += p;
        else if (kind == 1) n += p;
        else return false;
    }

    if (m == 0 && n == 0) return false; // not a non-trivial sin/cos product
    m_out = m;
    n_out = n;
    return true;
}

bool TrigCombinationStrategy::extract_tan_power(
    const SymbolicExpr& expr, const std::string& var, int& n_out) {
    int kind = -1, p = 0;
    if (!trig_extract_factor(lamina::detail::node(expr), var, kind, p)) return false;
    if (kind != 2) return false;
    n_out = p;
    return true;
}

bool TrigCombinationStrategy::extract_sec_power(
    const SymbolicExpr& expr, const std::string& var, int& n_out) {
    int kind = -1, p = 0;
    if (!trig_extract_factor(lamina::detail::node(expr), var, kind, p)) return false;
    if (kind != 3) return false;
    n_out = p;
    return true;
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext&, int depth) {

    // sin^m(x) * cos^n(x)
    int m = 0, n = 0;
    if (extract_sin_cos_powers(expr, var, m, n)) {
        if (m < 0 || n < 0) return nullptr;
        if (m + n > 8) return nullptr;
        return integrate_sin_m_cos_n(m, n, 1, var, ctx, depth);
    }

    // tan^n(x)
    int tn = 0;
    if (extract_tan_power(expr, var, tn)) {
        if (tn < 2 || tn > 8) return nullptr;
        return integrate_tan_power(tn, var, ctx, depth);
    }

    // sec^n(x) -- only even powers per design.
    int sn = 0;
    if (extract_sec_power(expr, var, sn)) {
        if (sn < 2 || sn > 8) return nullptr;
        if (sn % 2 != 0) return nullptr;
        return integrate_sec_power(sn, var, ctx, depth);
    }

    return nullptr;
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_sin_m_cos_n(
    int m, int n, long long scale,
    const std::string& var, Integrator& ctx, int depth) {

    if (m < 0 || n < 0) return nullptr;

    if (m == 0 && n == 0) {
        // Constant 1 with respect to var: integrand is essentially 1.
        return SymbolicExpr::variable(var);
    }
    if ((m % 2 == 1) || (n % 2 == 1)) {
        return integrate_odd_case(m, n, scale, var, ctx, depth);
    }
    return integrate_even_case(m, n, scale, var, ctx, depth);
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_odd_case(
    int m, int n, long long scale,
    const std::string& var, Integrator&, int) {

    if (scale == 0) return nullptr;

    // Pick which factor to peel:
    //   m odd  -> peel one sin, set u = cos(scale*x), du = -scale*sin(scale*x)dx
    //             so sin(scale*x)dx = -du/scale.
    //             remaining: sin^(m-1) cos^n = (1-u^2)^k * u^n, k=(m-1)/2.
    //   n odd  -> peel one cos, set u = sin(scale*x), du = +scale*cos(scale*x)dx
    //             so cos(scale*x)dx = du/scale.
    //             remaining: sin^m cos^(n-1) = u^m * (1-u^2)^k, k=(n-1)/2.
    bool peel_sin = (m % 2 == 1);
    int k = 0;
    int other_pow = 0;
    long long sign_factor = 1; // includes the sign from du sign.
    bool u_is_cos = false;

    if (peel_sin) {
        k = (m - 1) / 2;
        other_pow = n;
        sign_factor = -1;
        u_is_cos = true;
    } else {
        k = (n - 1) / 2;
        other_pow = m;
        sign_factor = 1;
        u_is_cos = false;
    }

    auto u_expr = u_is_cos ? make_cos_scaled(scale, var)
                            : make_sin_scaled(scale, var);

    // Result = sum_{i=0..k} C(k,i)*(-1)^i * u^(2i+other_pow+1) / (scale*(2i+other_pow+1)) * sign_factor
    std::vector<std::shared_ptr<const SymbolicNode>> add_terms;
    for (int i = 0; i <= k; ++i) {
        long long bin = binomial_ll(k, i);
        long long alt_sign = ((i % 2) == 0) ? 1 : -1;
        long long num = sign_factor * alt_sign * bin;
        int u_pow = 2 * i + other_pow + 1; // always >= 1
        long long den = scale * static_cast<long long>(u_pow);
        if (den == 0) return nullptr;

        std::shared_ptr<SymbolicExpr> u_to_pow;
        if (u_pow == 1) {
            u_to_pow = u_expr;
        } else {
            u_to_pow = SymbolicExpr::power(u_expr, SymbolicExpr::number(u_pow));
        }
        // Normalize sign so the rational denominator is positive.
        if (den < 0) { num = -num; den = -den; }
        auto coeff = SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
        auto term = SymbolicExpr::multiply(coeff, u_to_pow);
        add_terms.push_back(lamina::detail::node(term));
    }

    if (add_terms.empty()) return SymbolicExpr::number(0);
    if (add_terms.size() == 1) {
        return lamina::detail::make_expression_ptr(add_terms[0]);
    }
    return lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(add_terms));
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_even_case(
    int m, int n, long long scale,
    const std::string& var, Integrator& ctx, int depth) {

    if ((m % 2 != 0) || (n % 2 != 0)) return nullptr;
    int p = m / 2;
    int q = n / 2;
    int K = p + q;

    // sin^(2p)(c*x)*cos^(2q)(c*x)
    //   = (1/2^(p+q)) * sum_{i,j} C(p,i)(-1)^i C(q,j) cos^(i+j)(2c*x)
    // Compute coefficients a[k] = sum_{i+j=k} C(p,i)*C(q,j)*(-1)^i for k=0..K.
    std::vector<long long> a(K + 1, 0);
    for (int i = 0; i <= p; ++i) {
        long long bp = binomial_ll(p, i);
        long long si = ((i % 2) == 0) ? 1 : -1;
        for (int j = 0; j <= q; ++j) {
            long long bq = binomial_ll(q, j);
            a[i + j] += si * bp * bq;
        }
    }

    long long pow2 = 1;
    for (int t = 0; t < K; ++t) pow2 *= 2;
    if (pow2 == 0) pow2 = 1;

    std::vector<std::shared_ptr<const SymbolicNode>> add_terms;
    for (int k = 0; k <= K; ++k) {
        if (a[k] == 0) continue;
        // Recurse: integrate cos^k(2c*x) at the new scale 2c.
        auto inner = integrate_sin_m_cos_n(0, k, scale * 2, var, ctx, depth + 1);
        if (!inner) return nullptr;
        long long num = a[k];
        long long den = pow2;
        if (den < 0) { num = -num; den = -den; }
        auto coeff = SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
        auto term = SymbolicExpr::multiply(coeff, inner);
        add_terms.push_back(lamina::detail::node(term));
    }

    if (add_terms.empty()) return SymbolicExpr::number(0);
    if (add_terms.size() == 1) {
        return lamina::detail::make_expression_ptr(add_terms[0]);
    }
    return lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(add_terms));
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_tan_power(
    int n, const std::string& var, Integrator& ctx, int depth) {

    if (n < 0) return nullptr;
    auto v = SymbolicExpr::variable(var);

    if (n == 0) {
        // int 1 dx = x
        return v;
    }
    if (n == 1) {
        // int tan(x) dx = -ln(cos(x))
        auto cos_x = SymbolicExpr::cos(v);
        auto ln_cos = SymbolicExpr::ln(cos_x);
        return SymbolicExpr::multiply(SymbolicExpr::number(-1), ln_cos);
    }

    // n >= 2: tan^n = tan^(n-2)*(sec^2 - 1)
    //   int tan^n = tan^(n-1)/(n-1) - int tan^(n-2)
    auto tan_x = SymbolicExpr::tan(v);
    std::shared_ptr<SymbolicExpr> tan_pow_term;
    if (n - 1 == 1) {
        tan_pow_term = tan_x;
    } else {
        tan_pow_term = SymbolicExpr::power(tan_x, SymbolicExpr::number(n - 1));
    }
    auto first = SymbolicExpr::multiply(trig_rational(1, n - 1), tan_pow_term);

    auto rest = integrate_tan_power(n - 2, var, ctx, depth + 1);
    if (!rest) return nullptr;
    auto neg_rest = SymbolicExpr::multiply(SymbolicExpr::number(-1), rest);
    return SymbolicExpr::add(first, neg_rest);
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_sec_power(
    int n, const std::string& var, Integrator& ctx, int depth) {

    if (n < 2 || (n % 2) != 0) return nullptr;

    using FT = FunctionNode::FuncType;
    auto v = SymbolicExpr::variable(var);
    auto tan_x = SymbolicExpr::tan(v);
    auto make_sec_x = [&]() -> std::shared_ptr<SymbolicExpr> {
        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(FT::Sec,
                std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(v)}));
    };

    if (n == 2) {
        return tan_x;
    }

    auto sec_x = make_sec_x();
    std::shared_ptr<SymbolicExpr> sec_pow;
    if (n - 2 == 1) {
        sec_pow = sec_x;
    } else {
        sec_pow = SymbolicExpr::power(sec_x, SymbolicExpr::number(n - 2));
    }

    // First term: sec^(n-2)(x) * tan(x) / (n-1)
    auto inner = SymbolicExpr::multiply(sec_pow, tan_x);
    auto first = SymbolicExpr::multiply(trig_rational(1, n - 1), inner);

    // Second term: ((n-2)/(n-1)) * int sec^(n-2)
    auto rest = integrate_sec_power(n - 2, var, ctx, depth + 1);
    if (!rest) return nullptr;
    auto coeff = SymbolicExpr::number(Rational(BigInt(static_cast<long long>(n - 2)),
                                               BigInt(static_cast<long long>(n - 1))));
    auto second = SymbolicExpr::multiply(coeff, rest);

    return SymbolicExpr::add(first, second);
}

} // namespace lamina
