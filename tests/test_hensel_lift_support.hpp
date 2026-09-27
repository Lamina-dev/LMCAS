#pragma once
/**
 * @file test_hensel_lift_support.hpp
 * @brief Hensel 提升测试支持，涵盖 Mignotte 界、提升高度与二因子二次提升。
 */

#include "test_common.hpp"
#include "exact_factorization.hpp"
#include "transcendental_factor.hpp"
#include "polynomial.hpp"
#include "bigint.hpp"
#include "modular_arithmetic.hpp"

#include <vector>
#include <cmath>

using namespace LMCAS;

inline std::vector<BigInt> test_poly_mul(const std::vector<BigInt>& a,
                                  const std::vector<BigInt>& b) {
    if (a.empty() || b.empty()) { return {}; }
    size_t n = a.size() + b.size() - 1;
    std::vector<BigInt> result(n, BigInt(0));
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            result[i + j] = result[i + j] + a[i] * b[j];
        }
    }
    return result;
}

inline BigInt test_sym_mod(const BigInt& c, const BigInt& m) {
    if (m.is_zero()) { return c; }
    BigInt r = c % m;
    if (r.is_negative()) { r = r + m; }
    BigInt half_m = m / BigInt(2);
    if (r > half_m) { r = r - m; }
    return r;
}

inline std::vector<BigInt> test_reduce(const std::vector<BigInt>& poly, const BigInt& m) {
    std::vector<BigInt> result = poly;
    for (auto& c : result) {
        c = test_sym_mod(c, m);
    }
    while (!result.empty() && result.back().is_zero()) {
        result.pop_back();
    }
    return result;
}

inline bool test_verify_factorization(const std::vector<BigInt>& f,
                               const std::vector<BigInt>& g,
                               const std::vector<BigInt>& h,
                               const BigInt& m) {
    auto product = test_poly_mul(g, h);
    auto reduced_product = test_reduce(product, m);
    auto reduced_f = test_reduce(f, m);

    size_t n = std::max(reduced_product.size(), reduced_f.size());
    reduced_product.resize(n, BigInt(0));
    reduced_f.resize(n, BigInt(0));

    for (size_t i = 0; i < n; ++i) {
        if (reduced_product[i] != reduced_f[i]) { return false; }
    }
    return true;
}

inline bool test_verify_bezout(const std::vector<BigInt>& s,
                        const std::vector<BigInt>& g,
                        const std::vector<BigInt>& t,
                        const std::vector<BigInt>& h,
                        const BigInt& m) {
    auto sg = test_poly_mul(s, g);
    auto th = test_poly_mul(t, h);

    size_t n = std::max(sg.size(), th.size());
    std::vector<BigInt> sum(n, BigInt(0));
    for (size_t i = 0; i < n; ++i) {
        BigInt a = (i < sg.size()) ? sg[i] : BigInt(0);
        BigInt b = (i < th.size()) ? th[i] : BigInt(0);
        sum[i] = a + b;
    }

    auto reduced = test_reduce(sum, m);

    if (reduced.empty()) { return false; }
    if (reduced[0] != BigInt(1)) { return false; }
    for (size_t i = 1; i < reduced.size(); ++i) {
        if (!reduced[i].is_zero()) { return false; }
    }
    return true;
}
