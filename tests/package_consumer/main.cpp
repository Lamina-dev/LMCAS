#include "assumption_context.hpp"
#include "poly_utils.hpp"
#include "property_store.hpp"
#include "query_interface.hpp"
#include "symbolic.hpp"
#include "lsr_expr.hpp"
#include "solve_strategies.hpp"
#include "lmmc/lsr_stdlib.h"

#include <iostream>
#include <string>

int main() {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(x, SymbolicExpr::number(1));
    if (!expr) {
        std::cerr << "failed to construct expression\n";
        return 1;
    }

    lamina::ComputationContext context;
    auto solved = lamina::solve_dispatch_checked(
        expr, "x", context, lamina::SolveOptions{});
    if (!solved || solved.value().kind() != lamina::SolutionSet::Kind::Finite ||
        solved.value().finite_solutions().size() != 1) {
        std::cerr << "failed to solve expression\n";
        return 2;
    }

    auto polynomial = lamina::symbolic_to_poly<Rational>(expr, "x");
    if (polynomial.coeffs.size() != 2 ||
        polynomial.coeffs[0] != Rational(1) ||
        polynomial.coeffs[1] != Rational(1)) {
        std::cerr << "failed to convert expression to polynomial\n";
        return 3;
    }

    auto interval = lamina::Interval::point(SymbolicExpr::number(0));
    lamina::ComputationContext interval_context;
    auto interval_union = lamina::IntervalUnion::from_intervals_checked(
        {interval}, interval_context);
    if (!interval_union || interval_union.value().intervals().size() != 1) {
        std::cerr << "failed to construct interval union\n";
        return 4;
    }

    lamina::PropertyStore properties;
    lamina::ComputationContext property_context;
    auto declared = properties.declare_continuous_checked(
        "f", interval, property_context);
    if (!declared) {
        std::cerr << "failed to declare property\n";
        return 5;
    }
    lamina::ComputationContext query_context;
    auto continuous = properties.is_continuous_checked(
        "f", interval, query_context);
    if (!continuous || !continuous.value()) {
        std::cerr << "failed to query property\n";
        return 6;
    }

    lamina::AssumptionContext assumptions;
    lamina::QueryInterface queries(assumptions);
    auto one = SymbolicExpr::number(1);
    auto positive = queries.query_positive_checked(*one);
    if (!positive || positive.value() != lamina::Tribool::True) {
        std::cerr << "failed to query positivity\n";
        return 7;
    }

    auto assumed_continuous = assumptions.current_properties().declare_continuous_checked(
        "g", interval);
    if (!assumed_continuous) {
        std::cerr << "failed to declare assumption property\n";
        return 8;
    }
    auto context_continuous = assumptions.is_continuous_checked("g", interval);
    if (!context_continuous || context_continuous.value() != lamina::Tribool::True) {
        std::cerr << "failed to query assumption property\n";
        return 9;
    }

    auto i = lamina::lsr::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LSR imaginary unit\n";
        return 10;
    }
    auto reserved_i = lamina::lsr::sym("i");
    if (reserved_i ||
        std::string(lamina::lsr::error_name(reserved_i.error())) !=
            "ImaginaryUnitReserved") {
        std::cerr << "failed to expose LSR diagnostic names\n";
        return 10;
    }
    auto i_squared = SymbolicExpr::multiply(i.value(), i.value());
    lamina::ComputationContext lsr_context;
    auto i_rule = lamina::lsr::equivalent_core(
        *i_squared, *SymbolicExpr::number(-1), lsr_context);
    if (!i_rule || !i_rule.value()) {
        std::cerr << "failed to prove LSR i*i == -1\n";
        return 11;
    }

    auto lsr_roots = lamina::lsr::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    if (!lsr_roots || lsr_roots.value().size() != 2) {
        std::cerr << "failed to lower LSR solve result to set<Expr>\n";
        return 12;
    }

    auto four_i = lamina::lsr::complex(SymbolicExpr::number(0),
                                       SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    auto lowered_complex = lamina::lsr::eval_complex(*three_plus_four_i);
    if (!lowered_complex || !lowered_complex.value().is_finite() ||
        lowered_complex.value().real.value != 3.0 ||
        lowered_complex.value().imag.value != 4.0) {
        std::cerr << "failed to explicitly lower LSR Expr to complex\n";
        return 13;
    }
    auto lsr_real = lamina::lsr::real(three_plus_four_i);
    auto lsr_imag = lamina::lsr::imag(three_plus_four_i);
    auto lsr_abs = lamina::lsr::abs(three_plus_four_i);
    if (!lsr_real || !lsr_imag || !lsr_abs ||
        !lamina::lsr::structurally_equal(*lsr_real.value(),
                                         *SymbolicExpr::number(3)) ||
        !lamina::lsr::structurally_equal(*lsr_imag.value(),
                                         *SymbolicExpr::number(4))) {
        std::cerr << "failed to call LSR complex part facade\n";
        return 14;
    }

    lmmc_complex_t numeric_i;
    lmmc_real_t numeric_abs = 0;
    if (lmmc_lsr_math_i(&numeric_i) != LMMC_STATUS_OK ||
        lmmc_lsr_math_complex_abs(&numeric_i, &numeric_abs) != LMMC_STATUS_OK ||
        numeric_abs != 1.0 ||
        std::string(lmmc_lsr_error_name(LMMC_STATUS_OUT_OF_RANGE)) !=
            "DomainError") {
        std::cerr << "failed to call LMMC LSR std.math adapter\n";
        return 15;
    }

    lmmc_mat_t numeric_a = {0};
    lmmc_mat_t numeric_inv = {0};
    lmmc_mat_t numeric_rhs = {0};
    lmmc_mat_t numeric_solution = {0};
    lmmc_eigen_gen_full_result_t numeric_eig = {0};
    lmmc_svd_result_t numeric_svd = {0};
    lmmc_lsr_eig_table_t numeric_eig_table = {0};
    lmmc_lsr_svd_table_t numeric_svd_table = {0};
    const lmmc_mat_t* numeric_named = nullptr;
    size_t numeric_rows = 0;
    size_t numeric_cols = 0;
    size_t numeric_rank = 0;
    lmmc_real_t numeric_det = 0;
    if (lmmc_mat_create(2, 2, &numeric_a) != LMMC_STATUS_OK ||
        lmmc_mat_create(2, 1, &numeric_rhs) != LMMC_STATUS_OK) {
        std::cerr << "failed to allocate LMMC matrices\n";
        return 16;
    }
    numeric_a.data[0] = 1;
    numeric_a.data[1] = 2;
    numeric_a.data[numeric_a.stride] = 3;
    numeric_a.data[numeric_a.stride + 1] = 4;
    numeric_rhs.data[0] = 5;
    numeric_rhs.data[numeric_rhs.stride] = 11;
    if (lmmc_lsr_linalg_shape(&numeric_a, &numeric_rows, &numeric_cols) !=
            LMMC_STATUS_OK ||
        numeric_rows != 2 || numeric_cols != 2 ||
        lmmc_lsr_linalg_det(&numeric_a, &numeric_det) != LMMC_STATUS_OK ||
        numeric_det != -2.0 ||
        lmmc_lsr_linalg_rank(&numeric_a, &numeric_rank) != LMMC_STATUS_OK ||
        numeric_rank != 2 ||
        lmmc_lsr_linalg_inv(&numeric_a, &numeric_inv) != LMMC_STATUS_OK ||
        lmmc_lsr_linalg_solve_left(&numeric_a, &numeric_rhs,
                                   &numeric_solution) != LMMC_STATUS_OK ||
        numeric_solution.data[0] != 1.0 ||
        numeric_solution.data[numeric_solution.stride] != 2.0 ||
        lmmc_lsr_linalg_eig(&numeric_a, &numeric_eig) != LMMC_STATUS_OK ||
        numeric_eig.real_parts.size != 2 ||
        lmmc_lsr_linalg_svd(&numeric_a, &numeric_svd) != LMMC_STATUS_OK ||
        numeric_svd.sigma.size != 2 ||
        lmmc_lsr_linalg_eig_table(&numeric_a, &numeric_eig_table) !=
            LMMC_STATUS_OK ||
        !(numeric_named = lmmc_lsr_eig_table_get(&numeric_eig_table,
                                                 "values_real")) ||
        numeric_named->rows != 2 || numeric_named->cols != 1 ||
        lmmc_lsr_linalg_svd_table(&numeric_a, &numeric_svd_table) !=
            LMMC_STATUS_OK ||
        !(numeric_named = lmmc_lsr_svd_table_get(&numeric_svd_table, "S")) ||
        numeric_named->rows != 2 || numeric_named->cols != 2) {
        std::cerr << "failed to call LMMC LSR std.linalg adapter\n";
        lmmc_lsr_svd_table_destroy(&numeric_svd_table);
        lmmc_lsr_eig_table_destroy(&numeric_eig_table);
        lmmc_svd_result_destroy(&numeric_svd);
        lmmc_eigen_gen_full_result_destroy(&numeric_eig);
        lmmc_mat_destroy(&numeric_solution);
        lmmc_mat_destroy(&numeric_inv);
        lmmc_mat_destroy(&numeric_rhs);
        lmmc_mat_destroy(&numeric_a);
        return 17;
    }
    lmmc_lsr_svd_table_destroy(&numeric_svd_table);
    lmmc_lsr_eig_table_destroy(&numeric_eig_table);
    lmmc_svd_result_destroy(&numeric_svd);
    lmmc_eigen_gen_full_result_destroy(&numeric_eig);
    lmmc_mat_destroy(&numeric_solution);
    lmmc_mat_destroy(&numeric_inv);
    lmmc_mat_destroy(&numeric_rhs);
    lmmc_mat_destroy(&numeric_a);

    std::cout << expr->to_string() << '\n';
    if (expr->to_string().empty()) {
        std::cerr << "expression string is empty\n";
        return 18;
    }
    return 0;
}
