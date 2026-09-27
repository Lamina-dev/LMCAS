#include "internal/solver_groebner_internal.hpp"
#include <algorithm>
#include <functional>
#include <queue>

namespace LMCAS::groebner_detail {
namespace {
using PairQueue = std::priority_queue<SugarPair, std::vector<SugarPair>, std::greater<SugarPair>>;

void initialize_sugar(std::vector<Poly>& basis) {
    for (auto& polynomial : basis) {
        if (polynomial.is_zero()) continue;
        int maximum = 0;
        for (const auto& [monomial, coefficient] : polynomial.terms) {
            (void)coefficient;
            int degree = 0;
            for (int exponent : monomial) degree += exponent;
            maximum = std::max(maximum, degree);
        }
        polynomial.sugar = maximum;
    }
}

PairQueue initial_pairs(const std::vector<Poly>& basis) {
    PairQueue pairs;
    for (size_t i = 0; i < basis.size(); ++i) {
        for (size_t j = i + 1; j < basis.size(); ++j) {
            pairs.push({i, j, compute_spoly_sugar(basis[i], basis[j])});
        }
    }
    return pairs;
}

std::vector<Poly> minimal_basis(const std::vector<Poly>& basis) {
    std::vector<bool> marked(basis.size(), false);
    for (size_t i = 0; i < basis.size(); ++i) {
        if (marked[i]) continue;
        for (size_t j = 0; j < basis.size(); ++j) {
            if (i == j || marked[j]) continue;
            if (divides_mono(basis[j].lead_monomial(), basis[i].lead_monomial())) {
                marked[i] = true;
                break;
            }
        }
    }
    std::vector<Poly> minimal;
    for (size_t i = 0; i < basis.size(); ++i) {
        if (!marked[i]) minimal.push_back(basis[i]);
    }
    return minimal;
}

void make_monic(Poly& polynomial) {
    Rational leading = polynomial.lead_coeff();
    if (leading == Rational(0) || leading == Rational(1)) return;
    Rational inverse = Rational(1) / leading;
    Poly monic(polynomial.num_vars);
    for (const auto& [monomial, coefficient] : polynomial.terms) {
        monic.add_term(monomial, coefficient * inverse);
    }
    polynomial = std::move(monic);
}
}

void buchberger(std::vector<Poly>& basis, size_t num_vars) {
    initialize_sugar(basis);
    auto pairs = initial_pairs(basis);
    std::set<std::pair<size_t, size_t>> processed_pairs;
    while (!pairs.empty()) {
        auto [i, j, sugar_degree] = pairs.top();
        pairs.pop();
        auto canonical_pair = i < j ? std::make_pair(i, j) : std::make_pair(j, i);
        processed_pairs.insert(canonical_pair);
        if (coprime_leading_monomials(basis[i], basis[j])) continue;
        if (chain_criterion(basis, i, j, processed_pairs)) continue;
        Poly polynomial = s_poly(basis[i], basis[j]);
        extend_variables(polynomial, num_vars);
        Poly remainder = reduce(polynomial, basis);
        if (remainder.is_zero()) continue;
        remainder.sugar = sugar_degree;
        size_t new_index = basis.size();
        basis.push_back(remainder);
        for (size_t k = 0; k < new_index; ++k) {
            pairs.push({k, new_index, compute_spoly_sugar(basis[k], basis[new_index])});
        }
    }
}

void reduce_basis(std::vector<Poly>& basis) {
    auto minimal = minimal_basis(basis);
    for (size_t i = 0; i < minimal.size(); ++i) {
        std::vector<Poly> others;
        for (size_t j = 0; j < minimal.size(); ++j) {
            if (j != i) others.push_back(minimal[j]);
        }
        minimal[i] = reduce(minimal[i], others);
    }
    std::vector<Poly> reduced;
    for (const auto& polynomial : minimal) {
        if (!polynomial.is_zero()) reduced.push_back(polynomial);
    }
    for (auto& polynomial : reduced) make_monic(polynomial);
    basis = std::move(reduced);
}
}
