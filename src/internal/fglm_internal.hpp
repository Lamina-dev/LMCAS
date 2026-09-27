#pragma once

#include "fglm.hpp"
#include <set>
#include <vector>

namespace LMCAS::fglm_detail {
class MonomialEnumerator {
public:
    MonomialEnumerator(size_t num_vars, const MonomialOrder& order, int max_degree);
    bool next(Monomial& out);

private:
    struct HeapCmp {
        const MonomialOrder* order;
        bool operator()(const Monomial& a, const Monomial& b) const {
            return (*order)(a, b);
        }
    };

    size_t num_vars_;
    MonomialOrder order_;
    int max_degree_;
    std::vector<Monomial> heap_;
    std::set<Monomial> visited_;
    HeapCmp cmp_{&order_};

    void push(const Monomial& monomial);
};
class GaussianEliminator {
public:
    bool add_vector(const std::vector<Rational>& vector,
                    std::vector<Rational>& combination);

private:
    std::vector<std::vector<Rational>> rows_;
    std::vector<size_t> pivots_;
    std::vector<std::vector<Rational>> combinations_;
    size_t basis_count_ = 0;
};

bool is_standard(const Monomial& monomial, const std::vector<Monomial>& leading);
bool has_pure_generators(const std::vector<Monomial>& leading, size_t num_vars);
std::vector<Monomial> leading_monomials(const std::vector<FGLMPoly>& basis);
void generate_degree(Monomial& current, size_t variable, int remaining,
                     std::vector<Monomial>& monomials);

}
