#include "internal/integration_support.hpp"

#include <limits>

namespace LMCAS {

namespace {

// Build sin(c*var) / cos(c*var) where c is a non-zero integer scale.
std::shared_ptr<SymbolicExpr> make_sin_scaled(const BigInt& scale, const std::string& var) {
    auto v = SymbolicExpr::variable(var);
    if (scale == BigInt(1)) { return SymbolicExpr::sin(v); }
    auto cx = SymbolicExpr::multiply(SymbolicExpr::number(scale), v);
    return SymbolicExpr::sin(cx);
}

std::shared_ptr<SymbolicExpr> make_cos_scaled(const BigInt& scale, const std::string& var) {
    auto v = SymbolicExpr::variable(var);
    if (scale == BigInt(1)) { return SymbolicExpr::cos(v); }
    auto cx = SymbolicExpr::multiply(SymbolicExpr::number(scale), v);
    return SymbolicExpr::cos(cx);
}
BigInt trig_binomial(int degree, int selected_degree) {
    return BigInt::ncr(BigInt(degree), BigInt(selected_degree));
}

// Match a single-argument FunctionNode whose argument is the integration
// variable itself. Returns 0 (sin), 1 (cos), 2 (tan), 3 (sec), or -1 on
// no match for this strategy.
int trig_match_of_var(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!fn) { return -1; }
    if (fn->arguments().size() != 1) { return -1; }
    auto v = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0]);
    if (!v || v->is_constant() || v->name() != var) { return -1; }
    using FT = FunctionNode::FuncType;
    switch (fn->type()) {
        case FT::Sin: return 0;
        case FT::Cos: return 1;
        case FT::Tan: return 2;
        case FT::Sec: return 3;
        default: return -1;
    }
}
bool trig_extract_nonneg_int(const std::shared_ptr<const SymbolicNode>& exp_node, int& degree_out) {
    auto e = LMCAS::detail::expression_from_node(exp_node);
    auto simp = e.simplify();
    if (!simp) { return false; }
    const auto number = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(simp));
    if (!number) { return false; }
    const auto degree = integration_bounded_exponent(*number, 0, std::numeric_limits<int>::max());
    if (!degree) { return false; }
    degree_out = *degree;
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
    if (!pn) { return false; }
    int base_kind = trig_match_of_var(pn->base(), var);
    if (base_kind < 0) { return false; }
    int p = 0;
    if (!trig_extract_nonneg_int(pn->exponent(), p)) { return false; }
    kind_out = base_kind;
    power_out = p;
    return true;
}

} // anonymous namespace

bool TrigCombinationStrategy::extract_sin_cos_powers(
    const SymbolicExpr& expr, const std::string& var, int& m_out, int& n_out) {

    int m = 0, n = 0;
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        factors = mul->operands();
    } else {
        factors.push_back(LMCAS::detail::node(expr));
    }

    for (const auto& f : factors) {
        int kind = -1, p = 0;
        if (!trig_extract_factor(f, var, kind, p)) { return false; }
        // Only sin/cos contribute to (m,n); tan/sec disqualify this form.
        if (kind == 0) {
            if (p > std::numeric_limits<int>::max() - m) { return false; }
            m += p;
        } else if (kind == 1) {
            if (p > std::numeric_limits<int>::max() - n) { return false; }
            n += p;
        }
        else return false;
    }

    if (m == 0 && n == 0) { return false; } /** @brief 要求 sin/cos 乘积至少含一个非零次幂。 */
    m_out = m;
    n_out = n;
    return true;
}

bool TrigCombinationStrategy::extract_tan_power(
    const SymbolicExpr& expr, const std::string& var, int& n_out) {
    int kind = -1, p = 0;
    if (!trig_extract_factor(LMCAS::detail::node(expr), var, kind, p)) { return false; }
    if (kind != 2) { return false; }
    n_out = p;
    return true;
}

bool TrigCombinationStrategy::extract_sec_power(
    const SymbolicExpr& expr, const std::string& var, int& n_out) {
    int kind = -1, p = 0;
    if (!trig_extract_factor(LMCAS::detail::node(expr), var, kind, p)) { return false; }
    if (kind != 3) { return false; }
    n_out = p;
    return true;
}

Result<std::shared_ptr<SymbolicExpr>> TrigCombinationStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext&, int depth) {

    // sin^m(x) * cos^n(x)
    int m = 0, n = 0;
    if (extract_sin_cos_powers(expr, var, m, n)) {
        if (m < 0 || n < 0) { return nullptr; }
        if (n > 8 || m > 8 - n) { return nullptr; }
        return integrate_sin_m_cos_n(m, n, BigInt(1), var, ctx, depth);
    }

    // tan^n(x)
    int tn = 0;
    if (extract_tan_power(expr, var, tn)) {
        if (tn < 2 || tn > 8) { return nullptr; }
        return integrate_tan_power(tn, var, ctx, depth);
    }

    // sec^n(x) -- only even powers per design.
    int sn = 0;
    if (extract_sec_power(expr, var, sn)) {
        if (sn < 2 || sn > 8) { return nullptr; }
        if (sn % 2 != 0) { return nullptr; }
        return integrate_sec_power(sn, var, ctx, depth);
    }

    return nullptr;
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_sin_m_cos_n(
    int m, int n, const BigInt& scale,
    const std::string& var, Integrator& ctx, int depth) {

    if (m < 0 || n < 0) { return nullptr; }

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
    int m, int n, const BigInt& scale,
    const std::string& var, Integrator&, int) {

    if (scale == BigInt(0)) { return nullptr; }

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
    BigInt sign_factor(1); /** @brief 包含换元微分 du 的符号。 */
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
        const BigInt binomial = trig_binomial(k, i);
        const BigInt alternating_sign((i % 2) == 0 ? 1 : -1);
        BigInt numerator = sign_factor * alternating_sign * binomial;
        int u_pow = 2 * i + other_pow + 1; /** @brief 上界为输入总次数加一。 */
        BigInt denominator = scale * BigInt(u_pow);
        if (denominator == BigInt(0)) { return nullptr; }

        std::shared_ptr<SymbolicExpr> u_to_pow;
        if (u_pow == 1) {
            u_to_pow = u_expr;
        } else {
            u_to_pow = SymbolicExpr::power(u_expr, SymbolicExpr::number(u_pow));
        }
        // Normalize sign so the rational denominator is positive.
        if (denominator.is_negative()) { numerator = -numerator; denominator = -denominator; }
        auto coeff = SymbolicExpr::number(Rational(numerator, denominator));
        auto term = SymbolicExpr::multiply(coeff, u_to_pow);
        add_terms.push_back(LMCAS::detail::node(term));
    }

    if (add_terms.empty()) { return SymbolicExpr::number(0); }
    if (add_terms.size() == 1) {
        return LMCAS::detail::make_expression_ptr(add_terms[0]);
    }
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<AddNode>(add_terms));
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_even_case(
    int m, int n, const BigInt& scale,
    const std::string& var, Integrator& ctx, int depth) {

    if ((m % 2 != 0) || (n % 2 != 0)) { return nullptr; }
    int p = m / 2;
    int q = n / 2;
    int K = p + q;

    // sin^(2p)(c*x)*cos^(2q)(c*x)
    //   = (1/2^(p+q)) * sum_{i,j} C(p,i)(-1)^i C(q,j) cos^(i+j)(2c*x)
    // Compute coefficients a[k] = sum_{i+j=k} C(p,i)*C(q,j)*(-1)^i for k=0..K.
    std::vector<BigInt> a(K + 1, BigInt(0));
    for (int i = 0; i <= p; ++i) {
        const BigInt bp = trig_binomial(p, i);
        const BigInt alternating_sign((i % 2) == 0 ? 1 : -1);
        for (int j = 0; j <= q; ++j) {
            const BigInt bq = trig_binomial(q, j);
            a[i + j] += alternating_sign * bp * bq;
        }
    }

    BigInt denominator(1);
    for (int t = 0; t < K; ++t) denominator *= BigInt(2);

    std::vector<std::shared_ptr<const SymbolicNode>> add_terms;
    for (int k = 0; k <= K; ++k) {
        if (a[k] == 0) { continue; }
        // Recurse: integrate cos^k(2c*x) at the new scale 2c.
        auto inner = integrate_sin_m_cos_n(0, k, scale * BigInt(2), var, ctx, depth + 1);
        if (!inner) { return nullptr; }
        auto coeff = SymbolicExpr::number(Rational(a[k], denominator));
        auto term = SymbolicExpr::multiply(coeff, inner);
        add_terms.push_back(LMCAS::detail::node(term));
    }

    if (add_terms.empty()) { return SymbolicExpr::number(0); }
    if (add_terms.size() == 1) {
        return LMCAS::detail::make_expression_ptr(add_terms[0]);
    }
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<AddNode>(add_terms));
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_tan_power(
    int n, const std::string& var, Integrator& ctx, int depth) {

    if (n < 0) { return nullptr; }
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
    auto first = SymbolicExpr::multiply(sym_rational(1, n - 1), tan_pow_term);

    auto rest = integrate_tan_power(n - 2, var, ctx, depth + 1);
    if (!rest) { return nullptr; }
    auto neg_rest = SymbolicExpr::multiply(SymbolicExpr::number(-1), rest);
    return SymbolicExpr::add(first, neg_rest);
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_sec_power(
    int n, const std::string& var, Integrator& ctx, int depth) {

    if (n < 2 || (n % 2) != 0) { return nullptr; }

    using FT = FunctionNode::FuncType;
    auto v = SymbolicExpr::variable(var);
    auto tan_x = SymbolicExpr::tan(v);
    auto make_sec_x = [&]() -> std::shared_ptr<SymbolicExpr> {
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<FunctionNode>(FT::Sec,
                std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(v)}));
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
    auto first = SymbolicExpr::multiply(sym_rational(1, n - 1), inner);

    // Second term: ((n-2)/(n-1)) * int sec^(n-2)
    auto rest = integrate_sec_power(n - 2, var, ctx, depth + 1);
    if (!rest) { return nullptr; }
    auto coeff = SymbolicExpr::number(Rational(BigInt(n - 2), BigInt(n - 1)));
    auto second = SymbolicExpr::multiply(coeff, rest);

    return SymbolicExpr::add(first, second);
}

} // namespace LMCAS
