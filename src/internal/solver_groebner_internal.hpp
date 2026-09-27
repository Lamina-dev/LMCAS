#pragma once

#include "rational.hpp"
#include "symbolic.hpp"
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

namespace LMCAS {
class SymbolicNode;
}

namespace LMCAS::groebner_detail {
using Monomial = std::vector<int>;

struct MonomialLess {
    bool operator()(const Monomial& a, const Monomial& b) const;
};

struct Term {
    Rational coeff;
    Monomial mono;
};

using PolyTerms = std::map<Monomial, Rational, MonomialLess>;
struct Poly {
    PolyTerms terms;
    size_t num_vars;
    int sugar = 0;

    Poly() : num_vars(0) {}
    Poly(size_t count) : num_vars(count) {}
    bool is_zero() const { return terms.empty(); }
    Monomial lead_monomial() const;
    Rational lead_coeff() const;
    Term lead_term() const;
    void add_term(const Monomial& monomial, const Rational& coefficient);
};

struct PolyContext {
    std::vector<std::string> ext_vars;
    std::unordered_map<std::string, size_t> transcendental_map;
    std::unordered_map<size_t, std::shared_ptr<const SymbolicNode>> aux_to_node;
    size_t num_original_vars;

    PolyContext(const std::vector<std::string>& variables)
        : ext_vars(variables), num_original_vars(variables.size()) {}
};

struct SugarPair {
    size_t i, j;
    int sugar_degree;
    bool operator>(const SugarPair& other) const {
        return sugar_degree > other.sugar_degree;
    }
};

Monomial mul_mono(const Monomial& a, const Monomial& b);
Monomial lcm_mono(const Monomial& a, const Monomial& b);
bool divides_mono(const Monomial& a, const Monomial& b);
Monomial div_mono(const Monomial& numerator, const Monomial& denominator);
Poly add_poly(const Poly& a, const Poly& b);
Poly sub_poly(const Poly& a, const Poly& b);
Poly mul_poly_term(const Poly& polynomial, const Term& term);
Poly mul_poly(const Poly& a, const Poly& b);
void extend_variables(Poly& polynomial, size_t count);
void align_variables(Poly& a, Poly& b);
Poly reduce(Poly polynomial, const std::vector<Poly>& basis);
Poly s_poly(const Poly& f, const Poly& g);
bool coprime_leading_monomials(const Poly& f, const Poly& g);
bool chain_criterion(const std::vector<Poly>& basis, size_t i, size_t j,
                     const std::set<std::pair<size_t, size_t>>& processed_pairs);
int compute_spoly_sugar(const Poly& f, const Poly& g);
void buchberger(std::vector<Poly>& basis, size_t num_vars);
void reduce_basis(std::vector<Poly>& basis);
Poly to_poly(const SymbolicExpr& expression, PolyContext& context);
SymbolicExpr from_poly_ext(
    const Poly& polynomial, const PolyContext& context);
}
