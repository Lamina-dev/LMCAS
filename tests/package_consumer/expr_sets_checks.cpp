#include "assumption_context.hpp"
#include "expr.hpp"
#include "value.hpp"
#include "symbolic.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

using namespace LMCAS;

static int check_polynomial_roots() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_roots = LMCAS::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    if (!lsr_roots || lsr_roots.value().size() != 2) {
        std::cerr << "failed to lower LMCAS solve result to set<Expr>\n";
        return 12;
    }
    return 0;
}

static int check_named_solvers() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_named_roots = LMCAS::roots(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    auto lsr_named_solve = LMCAS::solve(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(-1)),
        "x");
    if (!lsr_named_roots || lsr_named_roots.value().size() != 2 ||
        !lsr_named_solve || lsr_named_solve.value().size() != 2) {
        std::cerr << "failed to call LMCAS roots/solve set<Expr> aliases\n";
        return 12;
    }
    return 0;
}

static int check_repeated_roots() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_repeated_roots = LMCAS::roots(
        SymbolicExpr::power(x, SymbolicExpr::number(2)), "x");
    if (!lsr_repeated_roots || lsr_repeated_roots.value().size() != 1 ||
        !lsr_repeated_roots.value().contains(*SymbolicExpr::number(0))) {
        std::cerr << "failed to lower repeated LMCAS roots to set<Expr>\n";
        return 12;
    }
    return 0;
}

static int check_algebraic_roots() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_cubic_roots = LMCAS::roots(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(3)),
                          SymbolicExpr::number(-2)),
        "x");
    if (!lsr_cubic_roots || lsr_cubic_roots.value().size() != 3) {
        std::cerr << "failed to lower RootOf LMCAS roots to set<Expr>\n";
        return 12;
    }
    return 0;
}

static int check_complex_roots() {
    auto x = SymbolicExpr::variable("x");
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto lsr_complex_roots = LMCAS::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(1)),
        "x");
    auto negative_i = LMCAS::complex(SymbolicExpr::number(0),
                                           SymbolicExpr::number(-1));
    if (!lsr_complex_roots || lsr_complex_roots.value().size() != 2 ||
        !lsr_complex_roots.value().contains(*i.value()) ||
        !negative_i ||
        !lsr_complex_roots.value().contains(*negative_i.value())) {
        std::cerr << "failed to lower LMCAS complex roots to set<Expr>\n";
        return 12;
    }
    return 0;
}

static int check_complex_root_domain() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_complex_roots = LMCAS::solve_expr_set(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)),
                          SymbolicExpr::number(1)),
        "x");
    auto lsr_complex_roots_subset_c = LMCAS::expr_set_subset_domain(
        lsr_complex_roots.value(), LMCAS::complexes());
    if (!lsr_complex_roots_subset_c ||
        !lsr_complex_roots_subset_c.value()) {
        std::cerr << "failed to prove installed LMCAS complex roots subset C\n";
        return 12;
    }
    return 0;
}

static int check_set_construction() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = LMCAS::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    if (!lsr_set_left || !lsr_set_right ||
        lsr_set_left.value().size() != 2 ||
        !lsr_set_left.value().contains(*SymbolicExpr::number(1))) {
        std::cerr << "failed to construct LMCAS set<Expr>\n";
        return 12;
    }
    return 0;
}

static int check_empty_set() {
    auto lsr_empty_set = LMCAS::expr_set({});
    if (!lsr_empty_set || !lsr_empty_set.value().empty()) {
        std::cerr << "failed to construct empty LMCAS set<Expr>\n";
        return 12;
    }
    return 0;
}

static int check_set_union_intersection() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = LMCAS::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto lsr_set_union = lsr_set_left.value().set_union(lsr_set_right.value());
    auto lsr_set_intersection =
        lsr_set_left.value().intersection(lsr_set_right.value());
    if (lsr_set_union.size() != 3 ||
        lsr_set_intersection.size() != 1 ||
        !lsr_set_intersection.contains(*SymbolicExpr::number(2))) {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    return 0;
}

static int check_set_differences() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = LMCAS::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto lsr_set_difference =
        lsr_set_left.value().difference(lsr_set_right.value());
    auto lsr_set_symmetric =
        lsr_set_left.value().symmetric_difference(lsr_set_right.value());
    if (lsr_set_difference.size() != 1 ||
        !lsr_set_difference.contains(*SymbolicExpr::number(1)) ||
        lsr_set_symmetric.size() != 2 ||
        !lsr_set_symmetric.contains(*SymbolicExpr::number(1)) ||
        !lsr_set_symmetric.contains(*SymbolicExpr::number(3))) {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    return 0;
}

static int check_subset_and_empty_identities() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = LMCAS::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto lsr_empty_set = LMCAS::expr_set({});
    auto lsr_set_union = lsr_set_left.value().set_union(lsr_set_right.value());
    auto lsr_set_intersection =
        lsr_set_left.value().intersection(lsr_set_right.value());
    if (!lsr_set_intersection.subset_of(lsr_set_union) ||
        !lsr_empty_set.value().subset_of(lsr_set_left.value()) ||
        lsr_set_left.value().set_union(lsr_empty_set.value()).size() !=
            lsr_set_left.value().size() ||
        !lsr_set_left.value().intersection(lsr_empty_set.value()).empty() ||
        lsr_set_left.value().difference(lsr_empty_set.value()).size() !=
            lsr_set_left.value().size()) {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    return 0;
}

static int check_facade_membership() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_facade_contains = LMCAS::expr_set_contains(
        lsr_set_left.value(), SymbolicExpr::number(1));
    auto lsr_facade_not_contains = LMCAS::expr_set_not_contains(
        lsr_set_left.value(), SymbolicExpr::number(3));
    if (!lsr_facade_contains ||
        !lsr_facade_contains.value() ||
        !lsr_facade_not_contains ||
        !lsr_facade_not_contains.value()) {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    return 0;
}

static int check_facade_union_intersection() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = LMCAS::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto lsr_facade_union = LMCAS::expr_set_union(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_intersection = LMCAS::expr_set_intersection(
        lsr_set_left.value(), lsr_set_right.value());
    if (!lsr_facade_union ||
        lsr_facade_union.value().size() != 3 ||
        !lsr_facade_intersection ||
        lsr_facade_intersection.value().size() != 1) {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    return 0;
}

static int check_facade_differences() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = LMCAS::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto lsr_facade_difference = LMCAS::expr_set_difference(
        lsr_set_left.value(), lsr_set_right.value());
    auto lsr_facade_symmetric = LMCAS::expr_set_symmetric_difference(
        lsr_set_left.value(), lsr_set_right.value());
    if (!lsr_facade_difference ||
        lsr_facade_difference.value().size() != 1 ||
        !lsr_facade_symmetric ||
        lsr_facade_symmetric.value().size() != 2) {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    return 0;
}

static int check_facade_subset() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_set_right = LMCAS::expr_set({
        SymbolicExpr::number(2), SymbolicExpr::number(3)});
    auto lsr_set_union = lsr_set_left.value().set_union(lsr_set_right.value());
    auto lsr_set_intersection =
        lsr_set_left.value().intersection(lsr_set_right.value());
    auto lsr_facade_subset = LMCAS::expr_set_subset(
        lsr_set_intersection, lsr_set_union);
    if (!lsr_facade_subset ||
        !lsr_facade_subset.value()) {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    return 0;
}

static int check_invalid_membership() {
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_facade_null_membership = LMCAS::expr_set_contains(
        lsr_set_left.value(), nullptr);
    if (lsr_facade_null_membership ||
        std::string(LMCAS::error_name(
            lsr_facade_null_membership.error())) != "SetElementTypeMismatch") {
        std::cerr << "failed to call LMCAS set<Expr> operations\n";
        return 12;
    }
    return 0;
}

static int check_domain_names() {
    auto lsr_domain_z = LMCAS::integers();
    auto lsr_domain_q = LMCAS::rationals();
    auto lsr_domain_r = LMCAS::reals();
    auto lsr_domain_c = LMCAS::complexes();
    auto lsr_domain_expr = LMCAS::expressions();
    if (std::string(lsr_domain_z.name()) != "Z" ||
        std::string(lsr_domain_q.name()) != "Q" ||
        std::string(lsr_domain_r.name()) != "R" ||
        std::string(lsr_domain_c.name()) != "C" ||
        std::string(lsr_domain_expr.name()) != "Expr") {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    return 0;
}

static int check_domain_inclusion() {
    auto lsr_domain_z = LMCAS::integers();
    auto lsr_domain_q = LMCAS::rationals();
    auto lsr_domain_r = LMCAS::reals();
    auto lsr_domain_c = LMCAS::complexes();
    auto lsr_domain_expr = LMCAS::expressions();
    auto lsr_domain_chain_left =
        LMCAS::domain_subset(lsr_domain_z, lsr_domain_q);
    auto lsr_domain_chain_right =
        LMCAS::domain_subset(lsr_domain_r, lsr_domain_c);
    auto lsr_domain_chain_expr =
        LMCAS::domain_subset(lsr_domain_c, lsr_domain_expr);
    auto lsr_domain_reverse =
        LMCAS::domain_subset(lsr_domain_c, lsr_domain_r);
    if (!lsr_domain_chain_left ||
        !lsr_domain_chain_left.value() ||
        !lsr_domain_chain_right ||
        !lsr_domain_chain_right.value() ||
        !lsr_domain_chain_expr ||
        !lsr_domain_chain_expr.value() ||
        !lsr_domain_reverse ||
        lsr_domain_reverse.value()) {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    return 0;
}

static int check_numeric_membership() {
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto lsr_domain_z = LMCAS::integers();
    auto lsr_domain_r = LMCAS::reals();
    auto lsr_domain_c = LMCAS::complexes();
    auto lsr_domain_exact = LMCAS::domain_contains(
        lsr_domain_z, SymbolicExpr::number(2));
    auto lsr_domain_real = LMCAS::domain_contains(
        lsr_domain_r, SymbolicExpr::number(0.25));
    auto lsr_domain_complex = LMCAS::domain_contains(
        lsr_domain_c, i.value());
    if (!lsr_domain_exact ||
        !lsr_domain_exact.value() ||
        !lsr_domain_real ||
        !lsr_domain_real.value() ||
        !lsr_domain_complex ||
        !lsr_domain_complex.value()) {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    return 0;
}

static int check_ordinary_i_membership() {
    auto lsr_domain_r = LMCAS::reals();
    auto lsr_domain_c = LMCAS::complexes();
    auto lsr_domain_legacy_i = LMCAS::domain_contains(
        lsr_domain_c, SymbolicExpr::variable("i"));
    auto lsr_domain_legacy_i_not_real = LMCAS::domain_contains(
        lsr_domain_r, SymbolicExpr::variable("i"));
    auto lsr_legacy_four_i = SymbolicExpr::multiply(
        SymbolicExpr::number(4), SymbolicExpr::variable("i"));
    auto lsr_legacy_three_plus_four_i = SymbolicExpr::add(
        SymbolicExpr::number(3), lsr_legacy_four_i);
    auto lsr_domain_legacy_complex_arithmetic =
        LMCAS::domain_contains(lsr_domain_c,
                                     lsr_legacy_three_plus_four_i);
    if (lsr_domain_legacy_i ||
        std::string(LMCAS::error_name(
            lsr_domain_legacy_i.error())) != "Inconclusive" ||
        lsr_domain_legacy_i_not_real ||
        std::string(LMCAS::error_name(
            lsr_domain_legacy_i_not_real.error())) != "Inconclusive" ||
        lsr_domain_legacy_complex_arithmetic ||
        std::string(LMCAS::error_name(
            lsr_domain_legacy_complex_arithmetic.error())) != "Inconclusive") {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    return 0;
}

static int check_symbol_and_real_sets() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_set_left = LMCAS::expr_set({
        SymbolicExpr::number(1), SymbolicExpr::number(1),
        SymbolicExpr::number(2)});
    auto lsr_domain_r = LMCAS::reals();
    auto lsr_domain_c = LMCAS::complexes();
    auto lsr_domain_expr = LMCAS::expressions();
    auto lsr_domain_expr_symbol = LMCAS::domain_contains(
        lsr_domain_expr, x);
    auto lsr_set_subset_r = LMCAS::expr_set_subset_domain(
        lsr_set_left.value(), lsr_domain_r);
    auto lsr_set_subset_c = LMCAS::expr_set_subset_domain(
        lsr_set_left.value(), lsr_domain_c);
    if (!lsr_domain_expr_symbol ||
        !lsr_domain_expr_symbol.value() ||
        !lsr_set_subset_r ||
        !lsr_set_subset_r.value() ||
        !lsr_set_subset_c ||
        !lsr_set_subset_c.value()) {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    return 0;
}

static int check_complex_set_domains() {
    auto i = LMCAS::imaginary_unit();
    if (!i) {
        std::cerr << "failed to construct LMCAS imaginary unit\n";
        return 10;
    }
    auto lsr_domain_r = LMCAS::reals();
    auto lsr_domain_c = LMCAS::complexes();
    auto lsr_complex_set = LMCAS::expr_set({i.value()});
    auto lsr_complex_set_subset_r =
        lsr_complex_set ? LMCAS::expr_set_subset_domain(
                              lsr_complex_set.value(), lsr_domain_r)
                        : LMCAS::Result<bool>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "complex set construction failed", "consumer");
    auto lsr_complex_set_subset_c =
        lsr_complex_set ? LMCAS::expr_set_subset_domain(
                              lsr_complex_set.value(), lsr_domain_c)
                        : LMCAS::Result<bool>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "complex set construction failed", "consumer");
    if (!lsr_complex_set_subset_r ||
        lsr_complex_set_subset_r.value() ||
        !lsr_complex_set_subset_c ||
        !lsr_complex_set_subset_c.value()) {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    return 0;
}

static int check_ordinary_i_set_domains() {
    auto lsr_domain_c = LMCAS::complexes();
    auto lsr_legacy_four_i = SymbolicExpr::multiply(
        SymbolicExpr::number(4), SymbolicExpr::variable("i"));
    auto lsr_legacy_three_plus_four_i = SymbolicExpr::add(
        SymbolicExpr::number(3), lsr_legacy_four_i);
    auto lsr_legacy_i_set =
        LMCAS::expr_set({SymbolicExpr::variable("i")});
    auto lsr_legacy_i_set_subset_c =
        lsr_legacy_i_set ? LMCAS::expr_set_subset_domain(
                               lsr_legacy_i_set.value(), lsr_domain_c)
                         : LMCAS::Result<bool>::failure(
                               LMCAS::CasErrc::InternalInvariant,
                               "legacy i set construction failed", "consumer");
    auto lsr_legacy_complex_arithmetic_set =
        LMCAS::expr_set({lsr_legacy_three_plus_four_i});
    auto lsr_legacy_complex_arithmetic_set_subset_c =
        lsr_legacy_complex_arithmetic_set
            ? LMCAS::expr_set_subset_domain(
                  lsr_legacy_complex_arithmetic_set.value(), lsr_domain_c)
            : LMCAS::Result<bool>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "legacy complex arithmetic set construction failed",
                  "consumer");
    if (lsr_legacy_i_set_subset_c ||
        std::string(LMCAS::error_name(
            lsr_legacy_i_set_subset_c.error())) != "Inconclusive" ||
        lsr_legacy_complex_arithmetic_set_subset_c ||
        std::string(LMCAS::error_name(
            lsr_legacy_complex_arithmetic_set_subset_c.error())) !=
            "Inconclusive") {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    return 0;
}

static int check_symbol_set_domains() {
    auto x = SymbolicExpr::variable("x");
    auto lsr_domain_r = LMCAS::reals();
    auto lsr_domain_expr = LMCAS::expressions();
    auto lsr_domain_unknown = LMCAS::domain_contains(
        lsr_domain_r, x);
    auto lsr_unknown_set = LMCAS::expr_set({x});
    auto lsr_unknown_set_subset_r =
        lsr_unknown_set ? LMCAS::expr_set_subset_domain(
                              lsr_unknown_set.value(), lsr_domain_r)
                        : LMCAS::Result<bool>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "unknown set construction failed", "consumer");
    auto lsr_unknown_set_subset_expr =
        lsr_unknown_set ? LMCAS::expr_set_subset_domain(
                              lsr_unknown_set.value(), lsr_domain_expr)
                        : LMCAS::Result<bool>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "unknown set construction failed", "consumer");
    if (!lsr_unknown_set_subset_expr ||
        !lsr_unknown_set_subset_expr.value() ||
        lsr_unknown_set_subset_r ||
        std::string(LMCAS::error_name(
            lsr_unknown_set_subset_r.error())) != "Inconclusive" ||
        lsr_domain_unknown ||
        std::string(LMCAS::error_name(lsr_domain_unknown.error())) !=
            "Inconclusive") {
        std::cerr << "failed to call LMCAS predefined number domain sets\n";
        return 12;
    }
    return 0;
}

int run_expr_sets_consumer_checks() {
    const auto checks = {
        check_polynomial_roots,
        check_named_solvers,
        check_repeated_roots,
        check_algebraic_roots,
        check_complex_roots,
        check_complex_root_domain,
        check_set_construction,
        check_empty_set,
        check_set_union_intersection,
        check_set_differences,
        check_subset_and_empty_identities,
        check_facade_membership,
        check_facade_union_intersection,
        check_facade_differences,
        check_facade_subset,
        check_invalid_membership,
        check_domain_names,
        check_domain_inclusion,
        check_numeric_membership,
        check_ordinary_i_membership,
        check_symbol_and_real_sets,
        check_complex_set_domains,
        check_ordinary_i_set_domains,
        check_symbol_set_domains,
    };
    for (const auto check : checks) {
        if (const int status = check(); status != 0) {
            return status;
        }
    }
    return 0;
}
