#include "internal/solver_groebner_internal.hpp"
#include <algorithm>

namespace LMCAS::groebner_detail {

        bool MonomialLess::operator()(const Monomial& a, const Monomial& b) const {

            for (size_t i = 0; i < a.size(); ++i) {
                if (a[i] != b[i]) {
                    return a[i] > b[i];
                }
            }
            return false;
        }

        Monomial Poly::lead_monomial() const {
            if (terms.empty()) return std::vector<int>(num_vars, 0);
            return terms.begin()->first;
        }

        Rational Poly::lead_coeff() const {
            if (terms.empty()) return Rational(0);
            return terms.begin()->second;
        }

        Term Poly::lead_term() const {
             if (terms.empty()) return {Rational(0), std::vector<int>(num_vars, 0)};
             return {terms.begin()->second, terms.begin()->first};
        }

        void Poly::add_term(const Monomial& m, const Rational& c) {
            if (c == Rational(0)) return;
            auto it = terms.find(m);
            if (it != terms.end()) {

                Rational new_c = it->second + c;
                if (new_c == Rational(0)) {
                    terms.erase(it);
                } else {
                    it->second = new_c;
                }
            } else {
                terms[m] = c;
            }
        }

    Monomial mul_mono(const Monomial& a, const Monomial& b) {
        Monomial res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = a[i] + b[i];
        return res;
    }

    Monomial lcm_mono(const Monomial& a, const Monomial& b) {
        Monomial res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = std::max(a[i], b[i]);
        return res;
    }

    bool divides_mono(const Monomial& a, const Monomial& b) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] > b[i]) return false;
        }
        return true;
    }

    Monomial div_mono(const Monomial& num, const Monomial& den) {
        Monomial res(num.size());
        for (size_t i = 0; i < num.size(); ++i) res[i] = num[i] - den[i];
        return res;
    }

    Poly add_poly(const Poly& a, const Poly& b) {
        Poly res = a;
        for (auto const& [m, c] : b.terms) {
            res.add_term(m, c);
        }
        return res;
    }

    Poly sub_poly(const Poly& a, const Poly& b) {
        Poly res = a;
        for (auto const& [m, c] : b.terms) {

            Rational c_neg = c * Rational(-1);
            res.add_term(m, c_neg);
        }
        return res;
    }

    Poly mul_poly_term(const Poly& p, const Term& t) {
        Poly res(p.num_vars);
        if (t.coeff == Rational(0)) return res;
        for (auto const& [m, c] : p.terms) {
            res.add_term(mul_mono(m, t.mono), c * t.coeff);
        }
        return res;
    }

    Poly mul_poly(const Poly& a, const Poly& b) {
        Poly res(a.num_vars);
        for (auto const& [ma, ca] : a.terms) {
            for (auto const& [mb, cb] : b.terms) {
                res.add_term(mul_mono(ma, mb), ca * cb);
            }
        }
        return res;
    }

    Poly reduce(Poly p, const std::vector<Poly>& G) {
        Poly r(p.num_vars);

        while (!p.is_zero()) {
            bool reduced = false;
            Monomial p_lm = p.lead_monomial();

            for (const auto& g : G) {
                if (divides_mono(g.lead_monomial(), p_lm)) {

                    Term factor;
                    factor.coeff = p.lead_coeff() / g.lead_coeff();
                    factor.mono = div_mono(p_lm, g.lead_monomial());

                    Poly term_g = mul_poly_term(g, factor);
                    p = sub_poly(p, term_g);
                    reduced = true;
                    break;
                }
            }

            if (!reduced) {
                Term lt = p.lead_term();
                r.add_term(lt.mono, lt.coeff);

                p.terms.erase(p.terms.begin());
            }
        }
        return r;
    }

    Poly s_poly(const Poly& f, const Poly& g) {
        Monomial L = lcm_mono(f.lead_monomial(), g.lead_monomial());

        Term t1;
        t1.mono = div_mono(L, f.lead_monomial());
        t1.coeff = Rational(1) / f.lead_coeff();

        Term t2;
        t2.mono = div_mono(L, g.lead_monomial());
        t2.coeff = Rational(1) / g.lead_coeff();

        return sub_poly(mul_poly_term(f, t1), mul_poly_term(g, t2));
    }

    bool coprime_leading_monomials(const Poly& f, const Poly& g) {
        if (f.is_zero() || g.is_zero()) {
            return false;
        }
        const Monomial& lm_f = f.lead_monomial();
        const Monomial& lm_g = g.lead_monomial();
        for (size_t i = 0; i < lm_f.size() && i < lm_g.size(); ++i) {
            if (lm_f[i] > 0 && lm_g[i] > 0) {
                return false;
            }
        }
        return true;
    }

    bool chain_criterion(const std::vector<Poly>& G, size_t i, size_t j,
                         const std::set<std::pair<size_t,size_t>>& processed_pairs) {
        Monomial L = lcm_mono(G[i].lead_monomial(), G[j].lead_monomial());
        for (size_t k = 0; k < G.size(); ++k) {
            if (k == i || k == j) {
                continue;
            }
            if (G[k].is_zero()) {
                continue;
            }
            if (divides_mono(G[k].lead_monomial(), L)) {
                auto pair_ik = (i < k) ? std::make_pair(i, k) : std::make_pair(k, i);
                auto pair_kj = (k < j) ? std::make_pair(k, j) : std::make_pair(j, k);
                if (processed_pairs.count(pair_ik) && processed_pairs.count(pair_kj)) {
                    return true;
                }
            }
        }
        return false;
    }

    int compute_spoly_sugar(const Poly& f, const Poly& g) {
        Monomial L = lcm_mono(f.lead_monomial(), g.lead_monomial());

        int deg_L = 0;
        for (int e : L) deg_L += e;

        int deg_lm_f = 0;
        for (int e : f.lead_monomial()) deg_lm_f += e;

        int deg_lm_g = 0;
        for (int e : g.lead_monomial()) deg_lm_g += e;

        int sugar_f = f.sugar + deg_L - deg_lm_f;
        int sugar_g = g.sugar + deg_L - deg_lm_g;

        return std::max(sugar_f, sugar_g);
    }

void extend_variables(Poly& polynomial, size_t count) {
    if (polynomial.num_vars >= count) return;
    Poly extended(count);
    for (const auto& [monomial, coefficient] : polynomial.terms) {
        Monomial resized = monomial;
        resized.resize(count, 0);
        extended.add_term(resized, coefficient);
    }
    polynomial = std::move(extended);
}

void align_variables(Poly& a, Poly& b) {
    if (a.num_vars < b.num_vars) {
        extend_variables(a, b.num_vars);
    } else if (b.num_vars < a.num_vars) {
        extend_variables(b, a.num_vars);
    }
}
}
