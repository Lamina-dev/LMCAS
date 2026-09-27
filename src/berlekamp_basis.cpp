#include "internal/berlekamp_support.hpp"
#include <algorithm>
#include <limits>

namespace LMCAS::berlekamp_detail {

/**
 * @brief 在 F_p 上计算 a(x) * b(x) mod f(x)。
 *
 * @param[in] a 第一个多项式系数向量(升幂排列)
 * @param[in] b 第二个多项式系数向量(升幂排列)
 * @param[in] f 模多项式系数向量(首一,升幂排列)
 * @param[in] p 素数模数
 * @return a * b mod f 的系数向量
 * @internal
 */
static std::vector<int64_t> bk_poly_mul_mod(
    const std::vector<int64_t>& a,
    const std::vector<int64_t>& b,
    const std::vector<int64_t>& f,
    int64_t p) {

    if (a.empty() || b.empty()) {
        return {};
    }
    size_t prod_size = a.size() + b.size() - 1;
    std::vector<int64_t> product(prod_size, 0);

    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == 0) continue;
        for (size_t j = 0; j < b.size(); ++j) {
            product[i + j] = (product[i + j] + a[i] * b[j]) % p;
            if (product[i + j] < 0) product[i + j] += p;
        }
    }

    std::vector<int64_t> q, r;
    bk_div_mod(product, f, p, q, r);

    return r;
}

/**
 * @brief 在 F_p 上用平方乘法计算 base^exp mod f(x)。
 *
 * @param[in] base 底多项式系数向量(升幂排列)
 * @param[in] exp  指数(非负整数)
 * @param[in] f    模多项式系数向量(首一,升幂排列)
 * @param[in] p    素数模数
 * @return base^exp mod f 的系数向量
 *
 * @see Berlekamp, E.R. "Factoring polynomials over finite fields."
 * @internal
 */
static std::vector<int64_t> bk_poly_pow_mod(
    const std::vector<int64_t>& base,
    int64_t exp,
    const std::vector<int64_t>& f,
    int64_t p) {

    if (exp == 0) {
        return {1};
    }

    if (exp == 1) {
        std::vector<int64_t> q, r;
        bk_div_mod(base, f, p, q, r);
        return r;
    }
    std::vector<int64_t> result = {1};
    std::vector<int64_t> cur = base;
    {
        std::vector<int64_t> q, r;
        bk_div_mod(cur, f, p, q, r);
        cur = std::move(r);
    }

    int64_t e = exp;
    while (e > 0) {
        if (e & 1) {
            result = bk_poly_mul_mod(result, cur, f, p);
        }
        e >>= 1;
        if (e > 0) {
            cur = bk_poly_mul_mod(cur, cur, f, p);
        }
    }

    return result;
}

/**
 * @brief 构造 F_p 上的 Berlekamp 矩阵 Q。
 *
 * 对次数为 n 的首一多项式 f，Q 的第 i 行为 x^(i*p) mod f(x) 的系数向量。
 * 先用平方乘法求 h = x^p mod f，再由第 0 行 [1, 0, ..., 0] 逐行乘 h 并取余。
 *
 * @param[in] f 首一多项式系数向量(升幂排列,deg(f) = n)
 * @param[in] p 素数模数
 * @return nxn Berlekamp 矩阵(行优先存储)
 *
 * @see Berlekamp, E.R. "Factoring polynomials over finite fields."
 *      Bell System Technical Journal, 46(8), 1967.
 * @internal
 */
std::vector<std::vector<int64_t>> bk_build_q_matrix(
    const std::vector<int64_t>& f,
    int64_t p) {

    int n = static_cast<int>(f.size()) - 1;  /**< f 的次数。 */

    std::vector<std::vector<int64_t>> Q(n, std::vector<int64_t>(n, 0));

    if (n <= 0) {
        return Q;
    }
    Q[0][0] = 1;

    if (n == 1) {
        return Q;
    }
    std::vector<int64_t> x_poly = {0, 1};
    std::vector<int64_t> h = bk_poly_pow_mod(x_poly, p, f, p);
    for (size_t j = 0; j < h.size() && j < static_cast<size_t>(n); ++j) {
        Q[1][j] = h[j];
    }
    std::vector<int64_t> current = h;
    for (int i = 2; i < n; ++i) {
        current = bk_poly_mul_mod(current, h, f, p);
        for (size_t j = 0; j < current.size() && j < static_cast<size_t>(n); ++j) {
            Q[i][j] = current[j];
        }
    }

    return Q;
}

static void bk_eliminate_column(std::vector<std::vector<int64_t>>& matrix,
                                int current_row, int col, int64_t p) {
    const int n = static_cast<int>(matrix.size());
    int64_t inverse = bk_mod_inverse(matrix[current_row][col], p);
    for (int j = 0; j < n; ++j) {
        matrix[current_row][j] = (matrix[current_row][j] * inverse) % p;
    }
    for (int row = 0; row < n; ++row) {
        if (row == current_row) continue;
        int64_t factor = matrix[row][col];
        if (factor == 0) continue;
        for (int j = 0; j < n; ++j) {
            matrix[row][j] =
                (matrix[row][j] - factor * matrix[current_row][j] % p + p) % p;
        }
    }
}

static void bk_reduce_matrix(std::vector<std::vector<int64_t>>& M, int64_t p,
                             std::vector<int>& pivot_col, std::vector<int>& pivot_row) {
    const int n = static_cast<int>(M.size());
    int current_row = 0;
    for (int col = 0; col < n && current_row < n; ++col) {
        int pivot = -1;
        for (int row = current_row; row < n; ++row) {
            if (M[row][col] != 0) {
                pivot = row;
                break;
            }
        }

        if (pivot == -1) {
            continue;
        }

        if (pivot != current_row) {
            std::swap(M[pivot], M[current_row]);
        }

        pivot_col[current_row] = col;
        pivot_row[col] = current_row;
        bk_eliminate_column(M, current_row, col, p);
        ++current_row;
    }
}

std::vector<std::vector<int64_t>> bk_null_space(
    const std::vector<std::vector<int64_t>>& Q,
    int64_t p) {

    int n = static_cast<int>(Q.size());
    if (n == 0) {
        return {};
    }
    std::vector<std::vector<int64_t>> M(n, std::vector<int64_t>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int64_t val = Q[j][i];
            if (i == j) {
                val = (val - 1 + p) % p;
            }
            M[i][j] = val;
        }
    }
    std::vector<int> pivot_col(n, -1);
    std::vector<int> pivot_row(n, -1);

    bk_reduce_matrix(M, p, pivot_col, pivot_row);
    std::vector<std::vector<int64_t>> basis;

    for (int col = 0; col < n; ++col) {
        if (pivot_row[col] != -1) {
            continue;
        }

        std::vector<int64_t> vec(n, 0);
        vec[col] = 1;  /**< 当前自由变量取 1，其余自由变量取 0。 */

        for (int row = 0; row < n; ++row) {
            int pc = pivot_col[row];
            if (pc == -1) continue;
            /**
             * @brief 由行简化阶梯形回代约束分量。
             * M[row][pc] = 1，故 x[pc] + sum_{自由列 j} M[row][j] * x[j] = 0，
             * 当前基向量满足 x[pc] = -M[row][col] mod p。
             */
            vec[pc] = (p - M[row][col]) % p;
        }

        basis.push_back(std::move(vec));
    }

    return basis;
}

}
