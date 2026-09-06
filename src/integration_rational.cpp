#include "internal/integration_support.hpp"
#include "internal/exact_matrix.hpp"

namespace LMCAS {

namespace {

// Canonical builder utilities used throughout the rational strategy.
inline std::shared_ptr<SymbolicExpr> rd_num_rat(const Rational& r) {
    return SymbolicExpr::number(r);
}

inline std::shared_ptr<SymbolicExpr> rd_num_int(long long v) {
    return SymbolicExpr::number(BigInt(v));
}

// Build (var - r) as a SymbolicExpr.
inline std::shared_ptr<SymbolicExpr> rd_var_minus(const std::string& var, const Rational& r) {
    auto v = SymbolicExpr::variable(var);
    if (r == Rational(0)) return v;
    auto neg_r = SymbolicExpr::number(Rational(0) - r);
    return SymbolicExpr::add(v, neg_r);
}

// Test whether a rational-coefficient poly is the zero polynomial.
inline bool rd_is_zero_poly(const Polynomial<Rational>& p) {
    if (p.coeffs.empty()) return true;
    for (const auto& c : p.coeffs) if (!(c == Rational(0))) return false;
    return true;
}

// Convert a Polynomial<Rational> into a SymbolicExpr in the given variable.
inline std::shared_ptr<SymbolicExpr> rd_poly_to_sym(
    const Polynomial<Rational>& p, const std::string& var) {
    if (rd_is_zero_poly(p)) return SymbolicExpr::number(0);
    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    auto v = SymbolicExpr::variable(var);
    for (size_t i = 0; i < p.coeffs.size(); ++i) {
        if (p.coeffs[i] == Rational(0)) continue;
        std::shared_ptr<SymbolicExpr> t;
        if (i == 0) {
            t = rd_num_rat(p.coeffs[i]);
        } else if (i == 1) {
            if (p.coeffs[i] == Rational(1)) {
                t = v;
            } else {
                t = SymbolicExpr::multiply(rd_num_rat(p.coeffs[i]), v);
            }
        } else {
            auto pw = SymbolicExpr::power(v, rd_num_int(static_cast<long long>(i)));
            if (p.coeffs[i] == Rational(1)) {
                t = pw;
            } else {
                t = SymbolicExpr::multiply(rd_num_rat(p.coeffs[i]), pw);
            }
        }
        terms.push_back(t);
    }
    if (terms.empty()) return SymbolicExpr::number(0);
    if (terms.size() == 1) return terms[0];
    auto res = terms[0];
    for (size_t i = 1; i < terms.size(); ++i) res = SymbolicExpr::add(res, terms[i]);
    return res;
}

// Recursively decompose a SymbolicNode into (numerator-polys, denominator-polys)
// over the rationals. Each polynomial is assumed to be a polynomial in `var`
// with rational coefficients. Returns false on any non-rational sub-expression
// (functions, irrational/symbolic exponents, variables other than the integration
// variable). Multiplication aggregates by appending; division comes from
// PowerNode with negative integer exponent. Numbers and the integration
// variable are converted directly via symbolic_to_poly.
bool rd_collect_rational(const std::shared_ptr<const SymbolicNode>& node,
                         const std::string& var,
                         std::vector<Polynomial<Rational>>& num,
                         std::vector<Polynomial<Rational>>& den) {
    if (!node) return false;

    // A NumberNode is a constant; convert directly.
    if (auto n = std::dynamic_pointer_cast<const NumberNode>(node)) {
        Polynomial<Rational> p =
            symbolic_to_poly<Rational>(LMCAS::detail::make_expression_ptr(n), var);
        num.push_back(p);
        return true;
    }

    // The integration variable, or any other variable that's actually a constant.
    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (v->name() == var) {
            num.push_back(Polynomial<Rational>({Rational(0), Rational(1)}, var));
            return true;
        }
        /// 当前符号常量保持在符号域中,不提升为有理系数.
        return false;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        // Addition: combine all summands over a common denominator.
        // We collect each operand as P_i/Q_i and assemble
        //   sum P_i * (prod_{j != i} Q_j)  /  prod_j Q_j
        std::vector<std::pair<std::vector<Polynomial<Rational>>,
                              std::vector<Polynomial<Rational>>>> parts;
        parts.reserve(add->operands().size());
        for (const auto& op : add->operands()) {
            std::vector<Polynomial<Rational>> sub_num, sub_den;
            if (!rd_collect_rational(op, var, sub_num, sub_den)) return false;
            parts.push_back({std::move(sub_num), std::move(sub_den)});
        }

        // Combine: numerator = sum_i (prod_k num_k_i) * prod_{j != i} (prod_k den_k_j)
        //          denominator = prod_i (prod_k den_k_i)
        Polynomial<Rational> total_num(var);    // start as 0
        Polynomial<Rational> total_den({Rational(1)}, var);

        // Pre-compute Q_i (denominator product per term) and total denominator.
        std::vector<Polynomial<Rational>> Q_per_term;
        Q_per_term.reserve(parts.size());
        for (const auto& [pn, pd] : parts) {
            Polynomial<Rational> q({Rational(1)}, var);
            for (const auto& d : pd) q = q * d;
            Q_per_term.push_back(q);
            total_den = total_den * q;
        }

        for (size_t i = 0; i < parts.size(); ++i) {
            Polynomial<Rational> p({Rational(1)}, var);
            for (const auto& nn : parts[i].first) p = p * nn;
            // multiply by prod_{j != i} Q_j
            for (size_t j = 0; j < parts.size(); ++j) {
                if (j == i) continue;
                p = p * Q_per_term[j];
            }
            total_num = total_num + p;
        }

        num.push_back(total_num);
        den.push_back(total_den);
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (!rd_collect_rational(op, var, num, den)) return false;
        }
        return true;
    }

    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        // Only integer exponents (negative or non-negative) are allowed.
        auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
        if (!en) return false;

        long long exp_v = 0;
        bool ok = false;
        if (std::holds_alternative<BigInt>(en->value())) {
            const auto& bi = std::get<BigInt>(en->value());
            // Bound to keep matrix sizes sane.
            int v = bi.to_int();
            if (v >= -64 && v <= 64) { ok = true; exp_v = v; }
        } else if (std::holds_alternative<Rational>(en->value())) {
            const auto& r = std::get<Rational>(en->value());
            if (r.is_integer()) {
                int v = r.to_BigInt().to_int();
                if (v >= -64 && v <= 64) { ok = true; exp_v = v; }
            }
        } else if (std::holds_alternative<lmmc_real_t>(en->value())) {
            lmmc_real_t d = std::get<lmmc_real_t>(en->value());
            if (std::isfinite(d) && d == std::floor(d) && d >= -64.0 && d <= 64.0) {
                ok = true; exp_v = static_cast<long long>(d);
            }
        }
        if (!ok) return false;

        // Build the base polynomial.
        std::vector<Polynomial<Rational>> bn, bd;
        if (!rd_collect_rational(pw->base(), var, bn, bd)) return false;

        // base_num = product of bn ; base_den = product of bd
        Polynomial<Rational> base_num({Rational(1)}, var);
        for (const auto& p : bn) base_num = base_num * p;
        Polynomial<Rational> base_den({Rational(1)}, var);
        for (const auto& p : bd) base_den = base_den * p;

        if (exp_v == 0) {
            num.push_back(Polynomial<Rational>({Rational(1)}, var));
            return true;
        }
        long long e_abs = std::llabs(exp_v);
        Polynomial<Rational> n_pow({Rational(1)}, var);
        Polynomial<Rational> d_pow({Rational(1)}, var);
        for (long long i = 0; i < e_abs; ++i) {
            n_pow = n_pow * base_num;
            d_pow = d_pow * base_den;
        }
        if (exp_v > 0) {
            num.push_back(n_pow);
            den.push_back(d_pow);
        } else {
            // (P/Q)^(-k) = (Q/P)^k
            num.push_back(d_pow);
            den.push_back(n_pow);
        }
        return true;
    }

    // FunctionNode etc. -> not a rational function in `var`.
    return false;
}

} // anonymous namespace

bool RationalDecompositionStrategy::extract_rational(
    const SymbolicExpr& expr, const std::string& var,
    Polynomial<Rational>& P_out, Polynomial<Rational>& Q_out) {

    std::vector<Polynomial<Rational>> nums, dens;
    if (!rd_collect_rational(LMCAS::detail::node(expr), var, nums, dens)) return false;

    Polynomial<Rational> P({Rational(1)}, var);
    for (const auto& p : nums) P = P * p;
    Polynomial<Rational> Q({Rational(1)}, var);
    for (const auto& p : dens) Q = Q * p;

    if (rd_is_zero_poly(Q)) return false; // 0 in denominator, refuse

    // Reduce by GCD to keep things small (and keep variable name consistent).
    if (!rd_is_zero_poly(P)) {
        Polynomial<Rational> g = Polynomial<Rational>::gcd(P, Q);
        if (!rd_is_zero_poly(g) && g.degree() >= 1) {
            auto pq = P.div_mod(g);
            auto qq = Q.div_mod(g);
            if (rd_is_zero_poly(pq.second) && rd_is_zero_poly(qq.second)) {
                P = pq.first;
                Q = qq.first;
            }
        }
    }

    // Make Q monic (rescale numerator accordingly).
    if (Q.degree() >= 0) {
        Rational lc = Q.lead_coeff();
        if (!(lc == Rational(1)) && !(lc == Rational(0))) {
            Polynomial<Rational> P_scaled(var);
            P_scaled.coeffs.reserve(P.coeffs.size());
            for (const auto& c : P.coeffs) P_scaled.coeffs.push_back(c / lc);
            Polynomial<Rational> Q_scaled(var);
            Q_scaled.coeffs.reserve(Q.coeffs.size());
            for (const auto& c : Q.coeffs) Q_scaled.coeffs.push_back(c / lc);
            P_scaled.trim();
            Q_scaled.trim();
            P = P_scaled;
            Q = Q_scaled;
        }
    }

    P_out = P;
    Q_out = Q;
    return true;
}

void RationalDecompositionStrategy::poly_divide(
    const Polynomial<Rational>& P, const Polynomial<Rational>& Q,
    Polynomial<Rational>& quotient_out,
    Polynomial<Rational>& remainder_out) {
    auto qr = P.div_mod(Q);
    quotient_out = qr.first;
    remainder_out = qr.second;
}

bool RationalDecompositionStrategy::factor_denominator(
    const Polynomial<Rational>& Q,
    std::vector<std::pair<Polynomial<Rational>, int>>& factors_out) {
    factors_out.clear();
    if (Q.degree() <= 0) return false;

    // Step 1: square-free factorization to peel off multiplicities.
    auto sqfree = square_free_factorization(Q);
    if (sqfree.empty()) return false;

    /// 步骤 2:对每个无平方因子逐次提取有理根,剩余整体因子次数至多为 2;
    /// 高次剩余因子使该策略返回未匹配.
    for (auto& [piece, mult] : sqfree) {
        Polynomial<Rational> current = piece.make_monic();

        if (current.degree() == 0) {
            // pure constant factor; nothing to integrate, ignore
            continue;
        }
        if (current.degree() == 1 || current.degree() == 2) {
            // Irreducible quadratics can stay as-is; quadratic discriminant
            // determines splitting later in solve_coefficients/integrate_term.
            // For deg 2 we still try rational roots to split into two linears
            // if possible.
            if (current.degree() == 2) {
                auto roots = find_rational_roots(current);
                if (!roots.empty()) {
                    for (const auto& r : roots) {
                        Polynomial<Rational> linear({Rational(0) - r, Rational(1)}, current.variable_name);
                        auto qr = current.div_mod(linear);
                        if (!rd_is_zero_poly(qr.second)) {
                            return false;
                        }
                        factors_out.push_back({linear, mult});
                        current = qr.first;
                        if (current.degree() == 0) break;
                    }
                }
                if (current.degree() == 2) {
                    // Confirmed irreducible quadratic over Q.
                    factors_out.push_back({current.make_monic(), mult});
                    current = Polynomial<Rational>({Rational(1)}, current.variable_name);
                } else if (current.degree() == 1) {
                    factors_out.push_back({current.make_monic(), mult});
                    current = Polynomial<Rational>({Rational(1)}, current.variable_name);
                }
                continue;
            }
            // degree 1
            factors_out.push_back({current, mult});
            continue;
        }

        // degree >= 3: try to split off rational linear factors.
        auto roots = find_rational_roots(current);
        for (const auto& r : roots) {
            Polynomial<Rational> linear({Rational(0) - r, Rational(1)}, current.variable_name);
            auto qr = current.div_mod(linear);
            if (rd_is_zero_poly(qr.second)) {
                factors_out.push_back({linear, mult});
                current = qr.first;
                if (current.degree() <= 0) break;
            }
        }

        if (current.degree() <= 0) continue;
        if (current.degree() == 1) {
            factors_out.push_back({current.make_monic(), mult});
            continue;
        }
        if (current.degree() == 2) {
            factors_out.push_back({current.make_monic(), mult});
            continue;
        }
        /// Q 上的高次整体因子由调用方映射为未求值积分节点.
        return false;
    }
    return true;
}

Result<bool> RationalDecompositionStrategy::solve_coefficients(
    const Polynomial<Rational>& P,
    const Polynomial<Rational>& Q,
    const std::vector<std::pair<Polynomial<Rational>, int>>& factors,
    std::vector<Polynomial<Rational>>& numerators_out,
    ComputationContext& context) {

    numerators_out.clear();

    // Map each (factor index, power) to a list of unknowns.
    // Linear factor (x - r)^k contributes k unknowns A_1..A_k (constant numerators).
    // Quadratic factor q(x)^k contributes 2k unknowns: (B_l, C_l) per power l.
    //
    // We also remember, for each (factor index, power l), the partial-fraction
    // term's *full* numerator polynomial (constant or linear) as a function of
    // the unknowns. Then we multiply each term by the missing portion of Q to
    // get a polynomial whose coefficients (in x) are linear combinations of the
    // unknowns; equating to P gives the linear system.

    // Total number of unknowns.
    size_t num_unknowns = 0;
    std::vector<std::pair<int, int>> ifac_ipow; // (factor_index, power index l in [1..mult])
    std::vector<int> kind_per_unknown_block;    // 1 -> linear, 2 -> quadratic
    std::vector<size_t> unknown_offset_per_term;
    unknown_offset_per_term.reserve(factors.size());

    for (size_t i = 0; i < factors.size(); ++i) {
        unknown_offset_per_term.push_back(num_unknowns);
        const auto& [fpoly, mult] = factors[i];
        if (fpoly.degree() == 1) {
            num_unknowns += static_cast<size_t>(mult);
        } else if (fpoly.degree() == 2) {
            num_unknowns += static_cast<size_t>(2 * mult);
        } else {
            return false;
        }
    }
    if (num_unknowns == 0) return false;

    const std::string& var = Q.variable_name;
    int N = Q.degree();
    if (N < 0) return false;

    // Build the matrix: each column is the contribution (coefficient vector
    // of x^0..x^{N-1}) of one unknown to the LHS polynomial; each row is one
    // equation indexed by the power of x.
    // Number of equations = N (size of polynomial of degree < N for the RHS P
    // when properly partial-fractioned: deg(P) < N).
    size_t rows = static_cast<size_t>(N);
    size_t cols = num_unknowns + 1; // augmented column for P

    // Initialize augmented matrix as Rational.
    std::vector<std::vector<Rational>> M(rows, std::vector<Rational>(cols, Rational(0)));

    // RHS column: coefficients of P (degree < N expected; if deg(P) >= N
    // then we have a problem - caller must run poly_divide first).
    if (P.degree() >= N) return false;
    for (size_t k = 0; k < rows; ++k) {
        Rational c = (k < P.coeffs.size()) ? P.coeffs[k] : Rational(0);
        M[k][num_unknowns] = c;
    }

    // For each factor, for each multiplicity power l, build the contribution
    // poly = numerator_term * (Q / fpoly^l). Each unknown in this term
    // contributes one column.
    for (size_t i = 0; i < factors.size(); ++i) {
        const auto& [fpoly, mult] = factors[i];
        // Q / fpoly^l for l = 1..mult; build incrementally: start with Q/fpoly^mult
        // by dividing Q by fpoly mult times.
        Polynomial<Rational> remaining = Q;
        for (int t = 0; t < mult; ++t) {
            auto qr = remaining.div_mod(fpoly);
            if (!rd_is_zero_poly(qr.second)) {
                /// Q 与当前因子幂的余式应为零;残余项表示内部不变量失效.
                return false;
            }
            remaining = qr.first;
        }
        // remaining now equals Q / fpoly^mult.

        // Now build for each power l = 1..mult:
        //   term coefficient pattern is multiplied by Q / fpoly^l
        //   = remaining * fpoly^(mult - l)
        Polynomial<Rational> fp_pow({Rational(1)}, var);
        // fp_pow starts as fpoly^0 = 1; we will iterate l from mult down to 1
        // so we increment fp_pow by * fpoly each time after using it.
        // Pre-build the array of coefficient polynomials for each l in 1..mult.
        std::vector<Polynomial<Rational>> q_over_fpow_l(mult + 1, Polynomial<Rational>(var));
        // index by l from 1..mult: q_over_fpow_l[l] = remaining * fpoly^(mult - l)
        Polynomial<Rational> cur({Rational(1)}, var);
        for (int l = mult; l >= 1; --l) {
            // when l = mult: cur = fpoly^0 = 1, factor = remaining * 1 = remaining
            // when decreasing l: factor = remaining * fpoly^(mult - l)
            Polynomial<Rational> factor_poly = remaining * cur;
            q_over_fpow_l[l] = factor_poly;
            cur = cur * fpoly;
        }

        size_t col_off = unknown_offset_per_term[i];

        if (fpoly.degree() == 1) {
            // Linear factor. Each l contributes one unknown A_l (constant
            // numerator). Coefficient column = q_over_fpow_l[l].
            for (int l = 1; l <= mult; ++l) {
                const auto& qpl = q_over_fpow_l[l];
                size_t this_col = col_off + static_cast<size_t>(l - 1);
                for (size_t k = 0; k < rows; ++k) {
                    Rational c = (k < qpl.coeffs.size()) ? qpl.coeffs[k] : Rational(0);
                    M[k][this_col] = c;
                }
            }
        } else if (fpoly.degree() == 2) {
            // Quadratic factor. Each l contributes two unknowns (B_l, C_l).
            // Numerator at power l = B_l * x + C_l ; the coefficient column for
            // B_l is q_over_fpow_l[l] shifted by 1 (multiplied by x), for C_l
            // is q_over_fpow_l[l] (multiplied by 1).
            for (int l = 1; l <= mult; ++l) {
                const auto& qpl = q_over_fpow_l[l];
                size_t base_col = col_off + static_cast<size_t>(2 * (l - 1));
                size_t b_col = base_col;     // for B_l
                size_t c_col = base_col + 1; // for C_l
                // C_l contributes qpl directly.
                for (size_t k = 0; k < rows; ++k) {
                    Rational c = (k < qpl.coeffs.size()) ? qpl.coeffs[k] : Rational(0);
                    M[k][c_col] = c;
                }
                // B_l contributes qpl * x = qpl with index shifted by 1.
                for (size_t k = 0; k < rows; ++k) {
                    if (k == 0) {
                        M[k][b_col] = Rational(0);
                    } else {
                        size_t src = k - 1;
                        Rational c = (src < qpl.coeffs.size()) ? qpl.coeffs[src] : Rational(0);
                        M[k][b_col] = c;
                    }
                }
            }
        }
    }

    std::vector<Rational> augmented;
    augmented.reserve(rows * cols);
    for (const auto& row : M) {
        augmented.insert(augmented.end(), row.begin(), row.end());
    }
    auto solved = detail::solve_rational_unique(
        rows, num_unknowns, std::move(augmented), context,
        "integrate.rational.coefficients");
    if (!solved) {
        if (solved.error().code == CasErrc::Inconclusive ||
            solved.error().code == CasErrc::DomainError) {
            return false;
        }
        return Result<bool>::failure(solved.error());
    }
    auto sol = std::move(solved.value());

    // Re-package solution into per-(factor, power) numerator polynomials.
    numerators_out.reserve(num_unknowns);
    for (size_t i = 0; i < factors.size(); ++i) {
        const auto& [fpoly, mult] = factors[i];
        size_t col_off = unknown_offset_per_term[i];
        if (fpoly.degree() == 1) {
            for (int l = 1; l <= mult; ++l) {
                Rational A = sol[col_off + static_cast<size_t>(l - 1)];
                numerators_out.push_back(Polynomial<Rational>({A}, Q.variable_name));
            }
        } else { // degree 2
            for (int l = 1; l <= mult; ++l) {
                size_t base_col = col_off + static_cast<size_t>(2 * (l - 1));
                Rational B = sol[base_col];
                Rational C = sol[base_col + 1];
                numerators_out.push_back(
                    Polynomial<Rational>({C, B}, Q.variable_name));
            }
        }
    }
    return true;
}

std::shared_ptr<SymbolicExpr> RationalDecompositionStrategy::integrate_term(
    const Polynomial<Rational>& numerator,
    const Polynomial<Rational>& factor,
    int power, const std::string& var) {

    if (rd_is_zero_poly(numerator)) return SymbolicExpr::number(0);

    if (factor.degree() == 1) {
        // factor = x - r (monic linear). r = -coeffs[0].
        Rational r = Rational(0) - factor.coeffs[0];
        Rational A = (numerator.coeffs.size() > 0) ? numerator.coeffs[0] : Rational(0);
        if (A == Rational(0)) return SymbolicExpr::number(0);

        if (power == 1) {
            // A * ln|x - r|. We use ln(x - r) here (matching existing strategies
            // which omit the absolute value).
            auto inner = rd_var_minus(var, r);
            return SymbolicExpr::multiply(rd_num_rat(A), SymbolicExpr::ln(inner));
        }
        // power >= 2: -A / ((power-1) * (x - r)^(power-1))
        Rational coeff = Rational(0) - A / Rational(BigInt(power - 1));
        auto base = rd_var_minus(var, r);
        auto exp = rd_num_int(power - 1);
        auto pw = SymbolicExpr::power(base, exp);
        auto inv_pw = SymbolicExpr::power(pw, rd_num_int(-1));
        return SymbolicExpr::multiply(rd_num_rat(coeff), inv_pw);
    }

    if (factor.degree() == 2) {
        // factor = x^2 + p*x + q (monic quadratic). discriminant = p^2 - 4q.
        Rational q = (factor.coeffs.size() > 0) ? factor.coeffs[0] : Rational(0);
        Rational p = (factor.coeffs.size() > 1) ? factor.coeffs[1] : Rational(0);

        // Verify it's irreducible: discriminant must be < 0. (If it's >= 0
        // we should have factored it; but be defensive.)
        Rational disc = p * p - Rational(4) * q;

        // numerator = B*x + C
        Rational B = (numerator.coeffs.size() > 1) ? numerator.coeffs[1] : Rational(0);
        Rational C = (numerator.coeffs.size() > 0) ? numerator.coeffs[0] : Rational(0);

        if (power >= 2) {
            // Higher powers of irreducible quadratics: leave as unevaluated
            // integral; this triggers the safe fallback in try_integrate.
            // Build the symbolic integrand for the unevaluated node.
            auto num_sym = rd_poly_to_sym(numerator, var);
            auto den_sym_base = rd_poly_to_sym(factor, var);
            auto den_pw = SymbolicExpr::power(den_sym_base, rd_num_int(power));
            auto inv_den = SymbolicExpr::power(den_pw, rd_num_int(-1));
            auto integrand = SymbolicExpr::multiply(num_sym, inv_den);
            return LMCAS::detail::make_expression_ptr(
                LMCAS::detail::make_node<IntegralNode>(
                    LMCAS::detail::node(integrand), var));
        }

        // power == 1: integral (B x + C) / (x^2 + p x + q) dx
        //   Split numerator: (B x + C) = (B/2) * (2x + p) + (C - B p / 2)
        //   integral (B/2)(2x+p)/(x^2+px+q) dx = (B/2) * ln(x^2+px+q)
        //   integral (C - B p / 2) / (x^2 + p x + q) dx
        //     = (C - B p / 2) * (2 / sqrt(4q - p^2)) * arctan( (2x + p) / sqrt(4q - p^2) )
        std::shared_ptr<SymbolicExpr> result = SymbolicExpr::number(0);

        Rational B_half = B / Rational(2);
        if (!(B_half == Rational(0))) {
            auto den_sym = rd_poly_to_sym(factor, var);
            auto ln_part = SymbolicExpr::ln(den_sym);
            auto term_log = SymbolicExpr::multiply(rd_num_rat(B_half), ln_part);
            result = SymbolicExpr::add(result, term_log);
        }

        Rational atan_coeff = C - (B * p) / Rational(2);
        if (!(atan_coeff == Rational(0))) {
            // s2 = -disc = 4q - p^2 (positive for irreducible factor).
            Rational s2 = Rational(0) - disc;
            // Build sqrt(s2) symbolically (works even when s2 isn't a perfect square).
            auto s2_sym = rd_num_rat(s2);
            auto sqrt_s2 = SymbolicExpr::sqrt(s2_sym);
            // 2x + p
            auto two_x = SymbolicExpr::multiply(rd_num_int(2), SymbolicExpr::variable(var));
            std::shared_ptr<SymbolicExpr> two_x_plus_p;
            if (p == Rational(0)) {
                two_x_plus_p = two_x;
            } else {
                two_x_plus_p = SymbolicExpr::add(two_x, rd_num_rat(p));
            }
            auto inv_sqrt = SymbolicExpr::power(sqrt_s2, rd_num_int(-1));
            auto arctan_arg = SymbolicExpr::multiply(two_x_plus_p, inv_sqrt);
            auto atan_part = make_arctan(arctan_arg);
            auto two = rd_num_int(2);
            auto coeff_part = SymbolicExpr::multiply(
                SymbolicExpr::multiply(rd_num_rat(atan_coeff), two), inv_sqrt);
            auto term_atan = SymbolicExpr::multiply(coeff_part, atan_part);
            result = SymbolicExpr::add(result, term_atan);
        }

        return result;
    }

    // Should not reach here (factor degree > 2 is rejected by factor_denominator).
    return SymbolicExpr::number(0);
}

Result<std::shared_ptr<SymbolicExpr>> RationalDecompositionStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int depth) {
    (void)ctx;
    (void)depth;

    Polynomial<Rational> P, Q;
    try {
        if (!extract_rational(expr, var, P, Q)) {
            return nullptr; // not rational -> let next strategy try
        }
    } catch (const std::invalid_argument&) {
        return nullptr;
    } catch (const std::out_of_range&) {
        return nullptr;
    } catch (const std::runtime_error&) {
        return nullptr;
    }
    if (!P.is_zero() && P.degree() > 0 && P.degree() < Q.degree()) {
        Polynomial<Rational> reduced_denominator;
        Polynomial<Rational> remainder;
        poly_divide(Q, P, reduced_denominator, remainder);
        if (remainder.is_zero()) {
            P = Polynomial<Rational>(
                std::vector<Rational>{Rational(1)}, P.variable_name);
            Q = std::move(reduced_denominator);
        }
    }

    // Defensive: if Q is zero or constant, this is not the right strategy.
    if (rd_is_zero_poly(Q) || Q.degree() < 1) return nullptr;

    /// PartialFractionStrategy 处理 deg<=1 的简单情形;degree>=2 交给本策略,
    /// 因为它的 integrate_term 对不可约二次因子(-> arctan + ln)是完整的,
    /// 而 PartialFraction 对不可约二次式会失败(留下未求值积分).
    if (Q.degree() == 1 && P.degree() <= 0 &&
        !P.coeffs.empty() && Q.coeffs.size() > 1 &&
        Q.coeffs[1] != Rational(0)) {
        auto coefficient = SymbolicExpr::number(
            P.coeffs[0] / Q.coeffs[1]);
        auto logarithm = SymbolicExpr::ln(rd_poly_to_sym(Q, var));
        return SymbolicExpr::multiply(
            coefficient, logarithm)->simplify();
    }
    if (Q.degree() < 2) return nullptr;
    try {
        // Long division if needed.
        Polynomial<Rational> quot, rem;
        if (P.degree() >= Q.degree()) {
            poly_divide(P, Q, quot, rem);
        } else {
            quot = Polynomial<Rational>(Q.variable_name);
            rem = P;
        }

        // Factor denominator.
        std::vector<std::pair<Polynomial<Rational>, int>> factors;
        if (!factor_denominator(Q, factors)) {
            /// Q 上因式分解未决时保留未求值积分节点.
            return Integrator::depends_on(expr, var)
                ? LMCAS::detail::make_expression_ptr(
                      LMCAS::detail::make_node<IntegralNode>(
                          LMCAS::detail::node(expr), var))
                : nullptr;
        }
        if (factors.empty()) {
            // Means Q is constant after factoring, so really there's nothing left;
            // integrate quot only.
            if (rd_is_zero_poly(quot)) return SymbolicExpr::number(0);
            // integral quot(x) dx = poly_integral(quot)
            // fall through; handled below
        }

        // Solve coefficients on the proper part rem / Q.
        std::vector<Polynomial<Rational>> numerators;
        if (!rd_is_zero_poly(rem)) {
            auto coefficients = solve_coefficients(
                rem, Q, factors, numerators, computation);
            if (!coefficients) {
                return Result<std::shared_ptr<SymbolicExpr>>::failure(
                    coefficients.error());
            }
            if (!coefficients.value()) {
                return LMCAS::detail::make_expression_ptr(
                    LMCAS::detail::make_node<IntegralNode>(
                        LMCAS::detail::node(expr), var));
            }
        }

        // Build the result piece by piece.
        std::shared_ptr<SymbolicExpr> result = SymbolicExpr::number(0);

        // Polynomial part from long division.
        if (!rd_is_zero_poly(quot)) {
            // Antiderivative of x^k is x^(k+1) / (k+1).
            for (size_t k = 0; k < quot.coeffs.size(); ++k) {
                if (quot.coeffs[k] == Rational(0)) continue;
                Rational coeff = quot.coeffs[k] / Rational(BigInt(static_cast<long long>(k + 1)));
                std::shared_ptr<SymbolicExpr> term;
                auto v = SymbolicExpr::variable(var);
                auto pw = SymbolicExpr::power(v, rd_num_int(static_cast<long long>(k + 1)));
                if (coeff == Rational(1)) {
                    term = pw;
                } else {
                    term = SymbolicExpr::multiply(rd_num_rat(coeff), pw);
                }
                result = SymbolicExpr::add(result, term);
            }
        }

        // Partial-fraction terms.
        if (!numerators.empty()) {
            size_t idx = 0;
            for (size_t i = 0; i < factors.size(); ++i) {
                const auto& [fpoly, mult] = factors[i];
                for (int l = 1; l <= mult; ++l, ++idx) {
                    if (idx >= numerators.size()) break;
                    auto term = integrate_term(numerators[idx], fpoly, l, var);
                    if (!term) {
                        // Fallback: unevaluated integral over the original.
                        return LMCAS::detail::make_expression_ptr(
                            LMCAS::detail::make_node<IntegralNode>(
                                LMCAS::detail::node(expr), var));
                    }
                    result = SymbolicExpr::add(result, term);
                }
            }
        }

        auto simplified = result->simplify();
        if (simplified) {
            // If the simplified form contains an unevaluated integral marker
            // (e.g. high-power irreducible quadratic), prefer keeping the
            // un-simplified result so downstream consumers can still read the
            // partial-fraction form.
            return simplified;
        }
        return result;
    } catch (const std::invalid_argument&) {
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<IntegralNode>(
                LMCAS::detail::node(expr), var));
    } catch (const std::out_of_range&) {
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<IntegralNode>(
                LMCAS::detail::node(expr), var));
    } catch (const std::runtime_error&) {
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<IntegralNode>(
                LMCAS::detail::node(expr), var));
    }
}

} // namespace LMCAS
