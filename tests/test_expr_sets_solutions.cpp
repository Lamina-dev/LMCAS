#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>
#include <iostream>
#include "assumption_context.hpp"
#include "residual_verification.hpp"

using namespace LMCAS;

namespace {

template <class T>
void report_error(const Result<T> &result) {
    if (!result)
        std::cerr << "CasError " << static_cast<int>(result.error().code)
                  << " [" << result.error().operation << "]: " << result.error().message << '\n';
}

TEST(ExprSetsSolutions, FiniteAndEmptySolutionVariants) {
    auto x = LMCAS::sym("x");
    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));
    LMCAS::ComputationContext context;
    auto solved = LMCAS::solve_set(equation, "x", context);
    const auto *finite_solutions = solved
                                       ? std::get_if<LMCAS::FiniteSolutions>(&solved.value())
                                       : nullptr;
    EXPECT_TRUE((finite_solutions && finite_solutions->values.size() == 2)) << "solve_set preserves both finite roots";

    auto empty = LMCAS::solve_set(SymbolicExpr::number(1), "x");
    EXPECT_TRUE((empty &&
                 std::holds_alternative<LMCAS::EmptySolutions>(
                     empty.value())))
        << "solve_set represents mathematical no-solution as Empty";
}

TEST(ExprSetsSolutions, FiniteSolutionCollection) {
    auto x = LMCAS::sym("x");
    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));
    auto finite_solved = LMCAS::solve_expr_set(equation, "x");
    EXPECT_TRUE((finite_solved && finite_solved.value().size() == 2)) << "solve_expr_set lowers finite solutions to set<Expr>";

    auto repeated_root_equation =
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2));
    auto repeated_roots = LMCAS::roots(repeated_root_equation, "x");
    EXPECT_TRUE((repeated_roots && repeated_roots.value().size() == 1 &&
                 repeated_roots.value().contains(*SymbolicExpr::number(0))))
        << "roots lowers repeated roots to one set<Expr> member";
}

TEST(ExprSetsSolutions, HigherDegreeExactRoots) {
    auto x = LMCAS::sym("x");
    auto quintic_equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(5)),
        SymbolicExpr::number(-2));
    auto quintic_roots = LMCAS::roots(quintic_equation, "x");
    EXPECT_TRUE((quintic_roots && quintic_roots.value().size() == 5)) << "roots lowers exact degree-five roots to finite set<Expr>";
}

TEST(ExprSetsSolutions, RootsAndSolveCollections) {
    auto x = LMCAS::sym("x");
    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));
    auto roots_solved = LMCAS::roots(equation, "x");
    bool roots_match = false;
    if (roots_solved) {
        roots_match = roots_solved.value().size() == 2 &&
                      roots_solved.value().contains(*SymbolicExpr::number(-1)) &&
                      roots_solved.value().contains(*SymbolicExpr::number(1));
    }
    EXPECT_TRUE((roots_match)) << "roots returns the LMCAS set<Expr> finite root collection";

    auto solve_solved = LMCAS::solve(equation, "x");
    bool solutions_match = false;
    if (solve_solved) {
        solutions_match = solve_solved.value().size() == 2 &&
                          solve_solved.value().contains(*SymbolicExpr::number(-1)) &&
                          solve_solved.value().contains(*SymbolicExpr::number(1));
    }
    EXPECT_TRUE((solutions_match)) << "solve returns the LMCAS set<Expr> finite solution collection";
}

TEST(ExprSetsSolutions, ComplexSolutionMembership) {
    auto x = LMCAS::sym("x");
    auto i = LMCAS::imaginary_unit();
    auto domain_c = LMCAS::complexes();
    auto complex_equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    auto complex_solved = LMCAS::solve_expr_set(complex_equation, "x");
    auto negative_i = LMCAS::complex(SymbolicExpr::number(0),
                                     SymbolicExpr::number(-1));
    EXPECT_TRUE((complex_solved && complex_solved.value().size() == 2)) << "solve_expr_set returns both complex roots for x^2 + 1";
    auto complex_roots_subset_c =
        complex_solved ? LMCAS::expr_set_subset_domain(
                             complex_solved.value(), domain_c)
                       : LMCAS::Result<bool>::failure(
                             LMCAS::CasErrc::InternalInvariant,
                             "complex solve set construction failed",
                             "test_expr");
    EXPECT_TRUE((complex_roots_subset_c && complex_roots_subset_c.value())) << "solve(x^2 + 1, x) subset C is directly checkable";
    EXPECT_TRUE((complex_solved && i && negative_i &&
                 complex_solved.value().contains(*i.value()) &&
                 complex_solved.value().contains(*negative_i.value())))
        << "solve_expr_set lowers x^2 + 1 roots to explicit LMCAS complex expressions";
}

TEST(ExprSetsSolutions, ComplexSolutionValues) {
    auto complex_equation = LMCAS::parse_expr("x^2 + 1");
    ASSERT_TRUE(complex_equation.has_value());
    auto complex_solved = LMCAS::solve_expr_set(complex_equation.value(), "x");
    ASSERT_TRUE(complex_solved.has_value());
    ASSERT_EQ(complex_solved.value().size(), 2u);
    bool saw_positive_i = false;
    bool saw_negative_i = false;
    for (const auto &root : complex_solved.value().elements()) {
        auto lowered_root = LMCAS::eval_complex(*root);
        ASSERT_TRUE(lowered_root.has_value()) << "complex solve roots explicitly lower to complex values";
        const auto &value = lowered_root.value();
        ASSERT_TRUE(value.real.is_finite());
        ASSERT_TRUE(value.imag.is_finite());
        if (value.real.value == 0.0 && value.imag.value == 1.0) {
            saw_positive_i = true;
        }
        if (value.real.value == 0.0 && value.imag.value == -1.0) {
            saw_negative_i = true;
        }
    }
    EXPECT_TRUE(saw_positive_i) << "complex solve roots evaluate to i";
    EXPECT_TRUE(saw_negative_i) << "complex solve roots evaluate to -i";
}

TEST(ExprSetsSolutions, ShiftedComplexSolutionValues) {
    auto x = LMCAS::sym("x");
    auto shifted_complex_equation = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
                          SymbolicExpr::multiply(SymbolicExpr::number(2), x.value())),
        SymbolicExpr::number(2));
    auto shifted_complex_solved =
        LMCAS::solve_expr_set(shifted_complex_equation, "x");
    EXPECT_TRUE((shifted_complex_solved &&
                 shifted_complex_solved.value().size() == 2))
        << "solve_expr_set returns both complex roots for x^2 + 2x + 2";
    if (shifted_complex_solved) {
        bool saw_negative_one_plus_i = false;
        bool saw_negative_one_minus_i = false;
        for (const auto &root : shifted_complex_solved.value().elements()) {
            auto lowered_root = LMCAS::eval_complex(*root);
            EXPECT_TRUE((lowered_root && lowered_root.value().is_finite())) << "shifted complex roots explicitly lower to complex values";
            if (!lowered_root) {
                continue;
            }
            const auto &value = lowered_root.value();
            if (value.real.value == -1.0 && value.imag.value == 1.0) {
                saw_negative_one_plus_i = true;
            }
            if (value.real.value == -1.0 && value.imag.value == -1.0) {
                saw_negative_one_minus_i = true;
            }
        }
        EXPECT_TRUE((saw_negative_one_plus_i && saw_negative_one_minus_i)) << "complex solve roots preserve nonzero real components";
    }
}

TEST(ExprSetsSolutions, NonfiniteSolutionCollections) {
    auto empty_solved = LMCAS::solve_expr_set(SymbolicExpr::number(1), "x");
    EXPECT_TRUE((empty_solved && empty_solved.value().empty())) << "solve_expr_set lowers mathematical no-solution to empty set<Expr>";

    auto universal_solved = LMCAS::solve_expr_set(SymbolicExpr::number(0), "x");
    EXPECT_TRUE((!universal_solved &&
                 universal_solved.error().code == LMCAS::CasErrc::Inconclusive))
        << "solve_expr_set does not pretend universal solutions are finite set<Expr>";
}

TEST(ExprSetsSolutions, NullSolutionInputs) {
    auto null_solve_input = LMCAS::solve_expr_set(nullptr, "x");
    auto null_solve_set_input = LMCAS::solve_set(nullptr, "x");
    auto null_roots_input = LMCAS::roots(nullptr, "x");
    auto null_solve_alias_input = LMCAS::solve(nullptr, "x");
    EXPECT_TRUE((!null_solve_input &&
                 null_solve_input.error().code == LMCAS::CasErrc::InvalidArgument))
        << "solve_expr_set rejects null equations before lowering";
    EXPECT_TRUE((!null_solve_set_input &&
                 null_solve_set_input.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "solve_set rejects null equations";
    EXPECT_TRUE((!null_roots_input &&
                 null_roots_input.error().code == LMCAS::CasErrc::InvalidArgument))
        << "roots rejects null expressions";
    EXPECT_TRUE((!null_solve_alias_input &&
                 null_solve_alias_input.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "solve rejects null equations";
}

TEST(ExprSetsSolutions, EmptySolutionVariables) {
    auto x = LMCAS::sym("x");
    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x.value(), SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));
    auto empty_solve_set_variable = LMCAS::solve_set(equation, "");
    auto empty_solve_variable = LMCAS::solve_expr_set(equation, "");
    auto empty_roots_variable = LMCAS::roots(equation, "");
    auto empty_solve_alias_variable = LMCAS::solve(equation, "");
    EXPECT_TRUE((!empty_solve_set_variable &&
                 empty_solve_set_variable.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "solve_set rejects empty variable names";
    EXPECT_TRUE((!empty_solve_variable &&
                 empty_solve_variable.error().code == LMCAS::CasErrc::InvalidArgument))
        << "solve_expr_set rejects empty variable names";
    EXPECT_TRUE((!empty_roots_variable &&
                 empty_roots_variable.error().code == LMCAS::CasErrc::InvalidArgument))
        << "roots rejects empty variable names";
    EXPECT_TRUE((!empty_solve_alias_variable &&
                 empty_solve_alias_variable.error().code == LMCAS::CasErrc::InvalidArgument))
        << "solve rejects empty variable names";
}

TEST(ExprSetsSolutions, ComplexQuadraticConsumerContracts) {
    auto polynomial = parse_expr("x^2 + 1");
    ASSERT_TRUE((polynomial.has_value())) << "quadratic polynomial parses";
    if (!polynomial) {
        return;
    }
    auto solved = roots(polynomial.value(), "x");
    EXPECT_TRUE((solved && solved.value().size() == 2)) << "finite quadratic set has both roots";
    if (!solved) {
        return;
    }
    for (const auto &root : solved.value().elements()) {
        auto operand = std::make_shared<SymbolicExpr>(*root);
        auto residual = polynomial.value()->substitute("x", operand)->simplify();
        auto residual_value = LMCAS::eval_complex(*residual);
        EXPECT_TRUE((residual_value && residual_value.value().is_finite() &&
                     residual_value.value().real.value == 0.0 &&
                     residual_value.value().imag.value == 0.0))
            << "complex evaluation verifies each finite root's polynomial residual";
        auto modulus = abs(operand);
        EXPECT_TRUE((modulus && modulus.value()->simplify()->is_one())) << "finite quadratic root has exact modulus one";
        auto conjugate = conj(operand);
        EXPECT_TRUE((conjugate && solved.value().contains(*conjugate.value()))) << "finite quadratic root set is closed under conjugation";
    }
}

void logarithm_domain_exclusions() {
    for (const char *source : {"ln(x)/(x-1)", "x*ln(x)"}) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << "domain-sensitive equation parses";
        if (!expression) {
            continue;
        }
        auto result = solve_expr_set(expression.value(), "x");
        ASSERT_TRUE((result.has_value())) << "complete real roots project without unresolved conditions";
        if (!result) {
            continue;
        }
        if (std::string(source) == "ln(x)/(x-1)") {
            EXPECT_TRUE((result.value().empty())) << "the only logarithm zero is a denominator pole";
        } else {
            EXPECT_TRUE((result.value().size() == 1 && result.value().contains(*SymbolicExpr::number(1)))) << "zero factor x=0 is excluded by the other factor's logarithm domain";
        }
    }
}

void logarithm_product_union() {
    auto product = parse_expr("(x-2)*ln(x)");
    auto roots = product ? solve_expr_set(product.value(), "x") : ExprSetResult::failure(product.error());
    EXPECT_TRUE((roots && roots.value().size() == 2 &&
                 roots.value().contains(*SymbolicExpr::number(1)) &&
                 roots.value().contains(*SymbolicExpr::number(2))))
        << "all zero-producing factors contribute their roots";
}

void logarithm_legal_denominator() {
    auto allowed = parse_expr("ln(x)/(x-2)");
    auto finite = allowed ? solve_expr_set(allowed.value(), "x") : ExprSetResult::failure(allowed.error());
    EXPECT_TRUE((finite && finite.value().size() == 1 && finite.value().contains(*SymbolicExpr::number(1)))) << "a legal denominator does not exclude the logarithm root";
}

void parameter_denominator_conditions() {
    auto parameter = parse_expr("ln(x)/(x-a)");
    ASSERT_TRUE((parameter.has_value())) << "parameter denominator parses";
    if (!parameter) {
        return;
    }
    auto result = solve_set(parameter.value(), "x");
    report_error(result);
    const auto *values = result ? std::get_if<FiniteSolutions>(&result.value()) : nullptr;
    EXPECT_TRUE((values && values->values.size() == 1 && values->values[0].value->is_one())) << "parameter denominator retains the exact root";
    if (values && values->values.size() == 1) {
        EXPECT_TRUE((!values->values[0].conditions.empty())) << "remaining denominator restriction is explicit";
        AssumptionContext facts;
        bool rejects_pole = false;
        for (const auto &condition : values->values[0].conditions) {
            auto bound = condition->substitute("a", SymbolicExpr::number(1));
            rejects_pole = rejects_pole || facts.evaluate_condition(*bound) == Tribool::False;
        }
        EXPECT_TRUE((rejects_pole)) << "the retained condition excludes a=1";
    }
    auto projection = solve_expr_set(parameter.value(), "x");
    EXPECT_TRUE((!projection && projection.error().code == CasErrc::Inconclusive)) << "conditional finite roots are not unconditional finite sets";
    auto at_pole = solve_set(parameter.value()->substitute("a", SymbolicExpr::number(1)), "x");
    EXPECT_TRUE((at_pole && std::holds_alternative<EmptySolutions>(at_pole.value()))) << "specializing the pole yields mathematical Empty";
}

void equality_original_domains() {
    auto equality = parse_expr("ln(x)/(x-1)==0");
    auto result = equality ? solve_expr_set(equality.value(), "x") : ExprSetResult::failure(equality.error());
    report_error(result);
    EXPECT_TRUE((result && result.value().empty())) << "the equality facade retains both original operand domains";
}

void periodic_denominator_exclusions(const ExprPtr &sine_ratio) {
    auto sine_result = solve_set(sine_ratio, "x");
    report_error(sine_result);
    const auto *families = sine_result ? std::get_if<ParametricSolutions>(&sine_result.value()) : nullptr;
    EXPECT_TRUE((families && families->values.size() == 1)) << "sine divided by its argument retains an infinite integer family";
    if (families && families->values.size() == 1) {
        const auto &family = families->values[0];
        EXPECT_TRUE((family.integer_parameters.size() == 1)) << "the period parameter is explicitly integer";
        if (family.integer_parameters.size() == 1) {
            AssumptionContext facts;
            bool rejects_zero = false;
            for (const auto &condition : family.conditions) {
                auto bound = condition->substitute(family.integer_parameters[0], SymbolicExpr::number(0));
                rejects_zero = rejects_zero || facts.evaluate_condition(*bound) == Tribool::False;
            }
            EXPECT_TRUE((rejects_zero)) << "the denominator excludes the zero integer index";
        }
    }
}

void product_binder_capture_avoidance() {
    auto occupied = parse_expr("sin(x)/(x-_k)");
    ResourceLimits limits;
    limits.max_steps = 100000;
    ComputationContext bounded(limits);
    auto fresh_result = solve_set(occupied.value(), "x", bounded);
    report_error(fresh_result);
    const auto *fresh = fresh_result ? std::get_if<ParametricSolutions>(&fresh_result.value()) : nullptr;
    EXPECT_TRUE((fresh && fresh->values.size() == 1 &&
                 fresh->values[0].integer_parameters.size() == 1 &&
                 fresh->values[0].integer_parameters[0] != "_k"))
        << "integer binders cannot capture names occurring only in another product factor";
}

void caller_scope_binder_isolation(const ExprPtr &sine_ratio) {
    auto assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE((assumptions->assume_sign("_k", Sign::Zero).has_value())) << "unrelated parent zero fact declared";
    assumptions->push();
    EXPECT_TRUE((assumptions->assume_sign("_k1", Sign::Negative).has_value())) << "unrelated child negative fact declared";
    ComputationContext scoped;
    EXPECT_TRUE((scoped.set_assumptions(assumptions).has_value())) << "caller scope stack installed";
    auto isolated = solve_set(sine_ratio, "x", scoped);
    report_error(isolated);
    const auto *isolated_families = isolated ? std::get_if<ParametricSolutions>(&isolated.value()) : nullptr;
    EXPECT_TRUE((isolated_families && isolated_families->values.size() == 1)) << "unrelated caller signs cannot erase an infinite solution family";
    if (isolated_families && isolated_families->values.size() == 1) {
        const auto &parameters = isolated_families->values[0].integer_parameters;
        EXPECT_TRUE((parameters.size() == 1 && parameters[0] != "_k" && parameters[0] != "_k1")) << "integer binders avoid names declared in every caller scope";
    }
    EXPECT_TRUE((assumptions->has_sign("_k", Sign::Zero) && assumptions->has_sign("_k1", Sign::Negative))) << "fresh-name selection preserves the caller's assumptions";
}

TEST(ExprSetsSolutions, OriginalDomainsAndProductUnions) {
    logarithm_domain_exclusions();
    logarithm_product_union();
    logarithm_legal_denominator();
    parameter_denominator_conditions();
    equality_original_domains();
    auto sine_ratio = parse_expr("sin(x)/x");
    periodic_denominator_exclusions(sine_ratio.value());
    product_binder_capture_avoidance();
    caller_scope_binder_isolation(sine_ratio.value());
}

TEST(ExprSetsSolutions, CompleteRootRepresentationPolicy) {
    for (const char *source : {"x*(x^5-x-1)", "(x-1)^2*(x^5-x-1)"}) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << "mixed-degree polynomial parses";
        if (!expression) {
            continue;
        }
        SolveOptions options;
        options.return_rootof = false;
        auto forbidden = solve_set(expression.value(), "x", options);
        EXPECT_TRUE((!forbidden && forbidden.error().code == CasErrc::Inconclusive)) << "forbidding a required representation fails the whole solve";
        options.return_rootof = true;
        auto complete = solve_set(expression.value(), "x", options);
        report_error(complete);
        const auto *finite = complete ? std::get_if<FiniteSolutions>(&complete.value()) : nullptr;
        EXPECT_TRUE((finite && finite->values.size() == 6)) << "the rational root and all five algebraic roots survive";
        if (finite) {
            std::size_t total = 0;
            for (const auto &root : finite->values) {
                total += root.multiplicity;
                EXPECT_TRUE((root.conditions.empty())) << "exact polynomial roots have no artificial real-domain restrictions";
            }
            EXPECT_TRUE((total == (std::string(source)[0] == 'x' ? 6u : 7u))) << "the root ledger preserves the original polynomial degree";
        }
    }
    auto expression = parse_expr("x^2-1");
    SolveOptions options;
    options.return_rootof = false;
    auto represented = solve_expr_set(expression.value(), "x", options);
    EXPECT_TRUE((represented && represented.value().size() == 2)) << "a polynomial with representable complete roots still projects";
    auto cubic = parse_expr("x^3-2");
    auto forbidden_cubic = solve_set(cubic.value(), "x", options);
    EXPECT_TRUE((!forbidden_cubic && forbidden_cubic.error().code == CasErrc::Inconclusive)) << "the same representation policy applies to exact cubics";
    auto approximate = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(3)),
        SymbolicExpr::number(-2.0));
    auto disallowed = solve_set(approximate, "x");
    EXPECT_TRUE((!disallowed && disallowed.error().code == CasErrc::Inconclusive)) << "an explicitly approximate cubic cannot silently claim exact default roots";
}

TEST(ExprSetsSolutions, SymbolicCoefficientDegeneracy) {
    for (const char *source : {"a*x", "a*x-1", "a*x+b", "(a-b)*x+c", "a^2*x+b", "a*x^2+x-1"}) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << "symbolic coefficient polynomial parses";
        if (!expression) {
            continue;
        }
        auto result = solve_set(expression.value(), "x");
        report_error(result);
        EXPECT_TRUE((result && std::holds_alternative<ConditionalSolutions>(result.value()))) << "unproved leading coefficients retain the complete predicate";
        auto projected = solve_finite_checked(expression.value(), "x");
        EXPECT_TRUE((!projected && projected.error().code == CasErrc::Inconclusive)) << "generic parameter solutions cannot be projected unconditionally";
    }
    auto ax = parse_expr("a*x");
    auto zero = solve_set(ax.value()->substitute("a", SymbolicExpr::number(0)), "x");
    EXPECT_TRUE((zero && std::holds_alternative<UniversalSolutions>(zero.value()))) << "a=0 leaves every x";
    auto ax_minus_one = parse_expr("a*x-1");
    auto empty = solve_set(ax_minus_one.value()->substitute("a", SymbolicExpr::number(0)), "x");
    EXPECT_TRUE((empty && std::holds_alternative<EmptySolutions>(empty.value()))) << "zero times x cannot equal one";
    {
        auto quadratic = parse_expr("a*x^2+x-1");
        auto specialized = solve_expr_set(quadratic.value()->substitute("a", SymbolicExpr::number(0)), "x");
        EXPECT_TRUE((specialized && specialized.value().size() == 1 &&
                     specialized.value().contains(*SymbolicExpr::number(1))))
            << "a vanished quadratic leading term reduces to the complete linear solution";
    }
    auto assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE((assumptions->assume_sign_checked("a", Sign::NonZero).has_value())) << "nonzero coefficient assumption is accepted";
    ComputationContext context;
    EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "solve context retains assumptions";
    auto known = solve_set(ax_minus_one.value(), "x", context);
    report_error(known);
    const auto *finite = known ? std::get_if<FiniteSolutions>(&known.value()) : nullptr;
    EXPECT_TRUE((finite && finite->values.size() == 1 && finite->values[0].conditions.empty())) << "a proved nonzero linear coefficient allows an unconditional exact root";
}

} // namespace
