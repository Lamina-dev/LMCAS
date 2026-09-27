#pragma once
#include "exact_factorization.hpp"

namespace LMCAS::berlekamp_detail {
std::vector<int64_t> bk_reduce_to_mod_coeffs(const Polynomial<Rational>& poly, int64_t p);
int64_t bk_mod_inverse(int64_t value, int64_t prime);
void bk_div_mod(const std::vector<int64_t>& a, const std::vector<int64_t>& b,
                int64_t p, std::vector<int64_t>& quotient, std::vector<int64_t>& remainder);
std::vector<int64_t> bk_gcd_mod(std::vector<int64_t> a, std::vector<int64_t> b, int64_t p);
bool bk_is_suitable_prime(const Polynomial<Rational>& poly, int64_t prime);
int64_t bk_select_prime(const Polynomial<Rational>& poly);
std::vector<std::vector<int64_t>> bk_build_q_matrix(const std::vector<int64_t>& f, int64_t p);
std::vector<std::vector<int64_t>> bk_null_space(
    const std::vector<std::vector<int64_t>>& matrix, int64_t p);
std::vector<std::vector<int64_t>> bk_split_factors(
    const std::vector<int64_t>& polynomial,
    const std::vector<std::vector<int64_t>>& basis, int dimension, int64_t prime);
}
