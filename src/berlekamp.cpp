#include "internal/berlekamp_support.hpp"

namespace LMCAS {
using namespace berlekamp_detail;

/**
 * @brief 将模 p 系数向量转换为 Polynomial<ModInt>.
 *
 * @param[in] coeffs 系数向量(升幂排列)
 * @param[in] p      素数模数
 * @param[in] var    变量名
 * @return 对应的 Polynomial<ModInt>
 * @internal
 */
static Polynomial<ModInt> bk_to_poly_modint(
    const std::vector<int64_t>& coeffs,
    int64_t p,
    const std::string& var) {

    std::vector<ModInt> mod_coeffs;
    mod_coeffs.reserve(coeffs.size());
    for (int64_t c : coeffs) {
        mod_coeffs.emplace_back(c, p);
    }
    /// 直接构造,不调用 trim()(已保证无高次零系数)
    Polynomial<ModInt> result(var);
    result.coeffs = std::move(mod_coeffs);
    return result;
}

BerlekampResult berlekamp_factor(
    const Polynomial<Rational>& poly,
    int64_t prime) {

    BerlekampResult result;

    /// 零多项式或常数多项式:无因子
    if (poly.is_zero() || poly.degree() <= 0) {
        return result;
    }

    /// 显式和自动模数都必须保持次数且使约化多项式无平方。
    int64_t p = prime;
    if (p > 0) {
        if (!bk_is_suitable_prime(poly, p)) {
            return result;
        }
    } else {
        p = bk_select_prime(poly);
        if (p < 0) {
            return result;
        }
    }

    result.prime = p;

    /// 将多项式约化到 F_p
    std::vector<int64_t> f_coeffs = bk_reduce_to_mod_coeffs(poly, p);

    /// 首一化
    if (!f_coeffs.empty()) {
        int64_t lc = f_coeffs.back();
        if (lc != 1) {
            int64_t lc_inv = bk_mod_inverse(lc, p);
            for (auto& c : f_coeffs) {
                c = (c * lc_inv) % p;
            }
        }
    }

    /// 线性多项式:本身不可约
    if (f_coeffs.size() <= 2) {
        result.factors.push_back(bk_to_poly_modint(f_coeffs, p, poly.variable_name));
        return result;
    }

    /// 构造 Berlekamp 矩阵 Q
    auto Q = bk_build_q_matrix(f_coeffs, p);

    /// 计算 (Q - I)^T 的零空间
    auto null_basis = bk_null_space(Q, p);
    int null_dim = static_cast<int>(null_basis.size());

    /// 存储零空间信息
    result.null_space_dim = null_dim;
    result.null_space_basis = null_basis;

    /// 零空间维度 = 1:多项式在 F_p 上不可约
    if (null_dim <= 1) {
        result.factors.push_back(bk_to_poly_modint(f_coeffs, p, poly.variable_name));
        return result;
    }

    /// 零空间维度 > 1:多项式可分解,使用基向量分裂
    auto split = bk_split_factors(f_coeffs, null_basis, null_dim, p);
    for (const auto& factor_coeffs : split) {
        result.factors.push_back(bk_to_poly_modint(factor_coeffs, p, poly.variable_name));
    }
    return result;
}

}
