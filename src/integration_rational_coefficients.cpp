#include "internal/integration_rational_support.hpp"
#include "internal/exact_rational_matrix.hpp"

namespace LMCAS {

namespace {
using RationalFactors = std::vector<std::pair<Polynomial<Rational>, int>>;

bool rd_unknown_offsets(const RationalFactors& factors,
    std::vector<size_t>& unknown_offset_per_term, size_t& num_unknowns) {
    for (size_t i = 0; i < factors.size(); ++i) {
        unknown_offset_per_term.push_back(num_unknowns);
        const auto& [fpoly, mult] = factors[i];
        if (fpoly.degree() == 1) {
            num_unknowns += static_cast<size_t>(mult);
        } else if (fpoly.degree() == 2) {
            num_unknowns += static_cast<size_t>(2 * mult);
        } else {
            return false;
        }
    }
    return num_unknowns != 0;
}

bool rd_factor_complements(const Polynomial<Rational>& Q,
    const Polynomial<Rational>& fpoly, int mult,
    std::vector<Polynomial<Rational>>& q_over_fpow_l) {
    const auto& var = Q.variable_name;
        Polynomial<Rational> remaining = Q;
        for (int t = 0; t < mult; ++t) {
            auto qr = remaining.div_mod(fpoly);
            if (!qr.second.is_zero()) {
                return false;
            }
            remaining = qr.first;
        }
        Polynomial<Rational> cur({Rational(1)}, var);
        for (int l = mult; l >= 1; --l) {
            Polynomial<Rational> factor_poly = remaining * cur;
            q_over_fpow_l[l] = factor_poly;
            cur = cur * fpoly;
        }
    return true;
}

void rd_fill_factor_columns(const Polynomial<Rational>& fpoly, int mult,
    const std::vector<Polynomial<Rational>>& q_over_fpow_l,
    size_t col_off, std::vector<std::vector<Rational>>& M) {
    const size_t rows = M.size();
        if (fpoly.degree() == 1) {
            for (int l = 1; l <= mult; ++l) {
                const auto& qpl = q_over_fpow_l[l];
                size_t this_col = col_off + static_cast<size_t>(l - 1);
                for (size_t k = 0; k < rows; ++k) {
                    Rational c = (k < qpl.coeffs.size()) ? qpl.coeffs[k] : Rational(0);
                    M[k][this_col] = c;
                }
            }
        } else if (fpoly.degree() == 2) {
            for (int l = 1; l <= mult; ++l) {
                const auto& qpl = q_over_fpow_l[l];
                size_t base_col = col_off + static_cast<size_t>(2 * (l - 1));
                size_t b_col = base_col;
                size_t c_col = base_col + 1;
                for (size_t k = 0; k < rows; ++k) {
                    Rational c = (k < qpl.coeffs.size()) ? qpl.coeffs[k] : Rational(0);
                    M[k][c_col] = c;
                }
                for (size_t k = 0; k < rows; ++k) {
                    if (k == 0) {
                        M[k][b_col] = Rational(0);
                    } else {
                        size_t src = k - 1;
                        Rational c = (src < qpl.coeffs.size()) ? qpl.coeffs[src] : Rational(0);
                        M[k][b_col] = c;
                    }
                }
            }
        }
}

bool rd_fill_coefficient_matrix(const Polynomial<Rational>& Q,
    const RationalFactors& factors, const std::vector<size_t>& unknown_offset_per_term,
    std::vector<std::vector<Rational>>& M) {
    for (size_t i = 0; i < factors.size(); ++i) {
        const auto& [fpoly, mult] = factors[i];
        std::vector<Polynomial<Rational>> complements(mult + 1, Polynomial<Rational>(Q.variable_name));
        if (!rd_factor_complements(Q, fpoly, mult, complements)) { return false; }
        rd_fill_factor_columns(fpoly, mult, complements, unknown_offset_per_term[i], M);
    }
    return true;
}

void rd_unpack_numerators(const RationalFactors& factors,
    const std::vector<size_t>& unknown_offset_per_term, const std::vector<Rational>& sol,
    const std::string& var, std::vector<Polynomial<Rational>>& numerators_out) {
    for (size_t i = 0; i < factors.size(); ++i) {
        const auto& [fpoly, mult] = factors[i];
        size_t col_off = unknown_offset_per_term[i];
        if (fpoly.degree() == 1) {
            for (int l = 1; l <= mult; ++l) {
                Rational A = sol[col_off + static_cast<size_t>(l - 1)];
                numerators_out.push_back(Polynomial<Rational>({A}, var));
            }
        } else { /** @brief 二次因子。 */
            for (int l = 1; l <= mult; ++l) {
                size_t base_col = col_off + static_cast<size_t>(2 * (l - 1));
                Rational B = sol[base_col];
                Rational C = sol[base_col + 1];
                numerators_out.push_back(
                    Polynomial<Rational>({C, B}, var));
            }
        }
    }
}
}

Result<bool> RationalDecompositionStrategy::solve_coefficients(
    const Polynomial<Rational>& P,
    const Polynomial<Rational>& Q,
    const std::vector<std::pair<Polynomial<Rational>, int>>& factors,
    std::vector<Polynomial<Rational>>& numerators_out,
    ComputationContext& context) {

    numerators_out.clear();

    size_t num_unknowns = 0;
    std::vector<size_t> unknown_offset_per_term;
    unknown_offset_per_term.reserve(factors.size());
    if (!rd_unknown_offsets(factors, unknown_offset_per_term, num_unknowns)) { return false; }
    int N = Q.degree();
    if (N < 0) { return false; }
    size_t rows = static_cast<size_t>(N);
    size_t cols = num_unknowns + 1; /**< P 的增广列。 */

    std::vector<std::vector<Rational>> M(rows, std::vector<Rational>(cols, Rational(0)));
    if (P.degree() >= N) { return false; }
    for (size_t k = 0; k < rows; ++k) {
        Rational c = (k < P.coeffs.size()) ? P.coeffs[k] : Rational(0);
        M[k][num_unknowns] = c;
    }

    if (!rd_fill_coefficient_matrix(Q, factors, unknown_offset_per_term, M)) { return false; }
    std::vector<Rational> augmented;
    augmented.reserve(rows * cols);
    for (const auto& row : M) {
        augmented.insert(augmented.end(), row.begin(), row.end());
    }
    auto solved = detail::solve_rational_unique(
        rows, num_unknowns, std::move(augmented), context,
        "integrate.rational.coefficients");
    if (!solved) {
        if (solved.error().code == CasErrc::Inconclusive ||
            solved.error().code == CasErrc::DomainError) {
            return false;
        }
        return Result<bool>::failure(solved.error());
    }
    auto sol = std::move(solved.value());
    numerators_out.reserve(num_unknowns);
    rd_unpack_numerators(factors, unknown_offset_per_term, sol, Q.variable_name, numerators_out);
    return true;
}

}
