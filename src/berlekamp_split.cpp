#include "internal/berlekamp_support.hpp"
#include <algorithm>
#include <limits>

namespace LMCAS::berlekamp_detail {

static bool bk_split_one_factor(std::vector<std::vector<int64_t>>& factor_list,
                                size_t fi, const std::vector<int64_t>& v_poly,
                                int null_dim, int64_t p) {
    for (int64_t c = 0; c < p; ++c) {
        if (static_cast<int>(factor_list.size()) >= null_dim) return false;
        std::vector<int64_t> v_minus_c = v_poly;
        if (v_minus_c.empty()) {
            v_minus_c.push_back((p - c) % p);
        } else {
            v_minus_c[0] = (v_minus_c[0] - c % p + p) % p;
        }
        while (!v_minus_c.empty() && v_minus_c.back() == 0) {
            v_minus_c.pop_back();
        }

        if (v_minus_c.empty()) continue;
        std::vector<int64_t> h = bk_gcd_mod(factor_list[fi], v_minus_c, p);
        int deg_h = static_cast<int>(h.size()) - 1;
        int deg_g = static_cast<int>(factor_list[fi].size()) - 1;

        if (deg_h > 0 && deg_h < deg_g) {
            std::vector<int64_t> quotient, remainder;
            bk_div_mod(factor_list[fi], h, p, quotient, remainder);

            factor_list[fi] = h;
            factor_list.push_back(quotient);

            return true;
        }
    }
    return false;
}

std::vector<std::vector<int64_t>> bk_split_factors(
    const std::vector<int64_t>& f_coeffs,
    const std::vector<std::vector<int64_t>>& null_basis,
    int null_dim,
    int64_t p) {

    std::vector<std::vector<int64_t>> factor_list;
    factor_list.push_back(f_coeffs);

    for (const auto& basis_vec : null_basis) {
        bool is_trivial = true;
        if (!basis_vec.empty() && basis_vec[0] == 1) {
            for (size_t i = 1; i < basis_vec.size(); ++i) {
                if (basis_vec[i] != 0) {
                    is_trivial = false;
                    break;
                }
            }
        } else {
            is_trivial = false;
        }
        if (is_trivial) continue;

        if (static_cast<int>(factor_list.size()) >= null_dim) break;
        std::vector<int64_t> v_poly = basis_vec;
        while (!v_poly.empty() && v_poly.back() == 0) {
            v_poly.pop_back();
        }
        for (size_t fi = 0; fi < factor_list.size(); ++fi) {
            if (static_cast<int>(factor_list.size()) >= null_dim) break;
            if (factor_list[fi].size() <= 2) continue;
            bk_split_one_factor(factor_list, fi, v_poly, null_dim, p);
        }
    }
    return factor_list;
}

}
