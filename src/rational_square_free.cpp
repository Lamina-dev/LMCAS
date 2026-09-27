#include "rational_polynomial.hpp"
#include "exact_factorization.hpp"

namespace LMCAS {

std::vector<std::pair<Polynomial<Rational>, int>> square_free_factorization(
    const Polynomial<Rational>& poly) {

    std::vector<std::pair<Polynomial<Rational>, int>> factors;

    if (poly.is_zero() || poly.degree() <= 0) {
        return factors;
    }

    Polynomial<Rational> f = poly.make_monic();

    if (f.degree() == 1) {
        factors.push_back({f, 1});
        return factors;
    }

    Polynomial<Rational> f_prime = f.differentiate();
    Polynomial<Rational> g = Polynomial<Rational>::gcd(f, f_prime);

    if (g.degree() <= 0) {
        factors.push_back({f, 1});
        return factors;
    }

    auto [w, w_rem] = f.div_mod(g);

    int multiplicity = 1;

    while (w.degree() > 0) {

        Polynomial<Rational> y = Polynomial<Rational>::gcd(g, w);

        auto [z, z_rem] = w.div_mod(y);

        if (z.degree() > 0) {

            factors.push_back({z.make_monic(), multiplicity});
        }

        auto [g_next, g_rem] = g.div_mod(y);
        g = g_next;
        w = y;
        multiplicity++;
    }

    if (g.degree() > 0) {
        factors.push_back({g.make_monic(), multiplicity});
    }

    return factors;
}

/**
 * @brief 计算多项式的无平方因子部分，用于 Berlekamp 模分解预处理。
 *
 * 令 f' = derivative(f)、g = gcd(f,f')；g 为常数时 f 已无重因子，
 * 否则取 f/g，并将结果首一化。
 *
 * @param[in] poly 输入的有理系数多项式
 * @return 无平方因子预处理结果
 */
TfSquareFreeResult tf_square_free(const Polynomial<Rational>& poly) {
    TfSquareFreeResult result;
    if (poly.degree() <= 0) {
        result.square_free = poly;
        result.repeated_factor = Polynomial<Rational>(Rational(1), poly.variable_name);
        result.had_repeated_factors = false;
        return result;
    }
    if (poly.degree() == 1) {
        result.square_free = poly;
        result.repeated_factor = Polynomial<Rational>(Rational(1), poly.variable_name);
        result.had_repeated_factors = false;
        return result;
    }

    Polynomial<Rational> deriv = poly.differentiate();

    Polynomial<Rational> g = Polynomial<Rational>::gcd(poly, deriv);

    if (g.degree() <= 0) {
        result.square_free = poly.make_monic();
        result.repeated_factor = Polynomial<Rational>(Rational(1), poly.variable_name);
        result.had_repeated_factors = false;
        return result;
    }

    auto [quotient, remainder] = poly.div_mod(g);
    result.square_free = quotient.make_monic();
    result.repeated_factor = g.make_monic();
    result.had_repeated_factors = true;
    return result;
}

}
