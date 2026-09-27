#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprEvaluation, NumericBindings) {
    auto x = LMCAS::sym("x");
    auto linear = SymbolicExpr::add(x.value(), SymbolicExpr::number(2));
    auto unbound = LMCAS::evalf(*linear);
    EXPECT_TRUE((!unbound &&
                 unbound.error().code == LMCAS::CasErrc::UnboundSymbol))
        << "evalf without a required binding reports UnboundSymbol";
    EXPECT_TRUE((!unbound &&
                 std::string(LMCAS::error_name(unbound.error())) ==
                     "UnboundSymbol"))
        << "evalf keeps the generic unbound symbol diagnostic";

    LMCAS::NumericBindings bindings{{"x", 3.0}};
    auto evaluated = LMCAS::evalf(*linear, bindings);
    EXPECT_TRUE((evaluated && evaluated.value().is_finite())) << "evalf with a binding succeeds";
    EXPECT_NEAR(evaluated.value().value, 5.0, 0.0) << "evalf computes the numeric value";
}

TEST(ExprEvaluation, NonfiniteNumericValues) {
    auto x = LMCAS::sym("x");
    auto nonfinite_binding = LMCAS::evalf(
        *x.value(), LMCAS::NumericBindings{{"x", INFINITY}});
    EXPECT_TRUE((!nonfinite_binding &&
                 nonfinite_binding.error().code ==
                     LMCAS::CasErrc::NumericFailure))
        << "evalf rejects non-finite numeric bindings";
    EXPECT_TRUE((!nonfinite_binding &&
                 std::string(LMCAS::error_name(
                     nonfinite_binding.error())) == "NumericFailure"))
        << "evalf reports the LMCAS numeric failure diagnostic for non-finite bindings";

    auto nonfinite_expression = LMCAS::evalf(*SymbolicExpr::infinity());
    EXPECT_TRUE((!nonfinite_expression &&
                 nonfinite_expression.error().code ==
                     LMCAS::CasErrc::NumericFailure))
        << "evalf rejects expressions that evaluate to infinity";
}

TEST(ExprEvaluation, EvaluationBudget) {
    auto x = LMCAS::sym("x");
    auto linear = SymbolicExpr::add(x.value(), SymbolicExpr::number(2));
    LMCAS::NumericBindings bindings{{"x", 3.0}};
    LMCAS::ResourceLimits exhausted_evalf_limits;
    exhausted_evalf_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_evalf_context(exhausted_evalf_limits);
    auto exhausted_evalf =
        LMCAS::evalf(*linear, bindings, exhausted_evalf_context);
    EXPECT_TRUE((!exhausted_evalf &&
                 exhausted_evalf.error().code ==
                     LMCAS::CasErrc::ResourceLimit))
        << "evalf reports ResourceLimit when the computation budget is exhausted";
    EXPECT_TRUE((!exhausted_evalf &&
                 std::string(LMCAS::error_name(
                     exhausted_evalf.error())) == "ResourceLimit"))
        << "evalf exposes the LMCAS resource-limit diagnostic";
}

TEST(ExprEvaluation, SingleAndTypedSubstitution) {
    auto x = LMCAS::sym("x");
    auto linear = SymbolicExpr::add(x.value(), SymbolicExpr::number(2));
    auto substituted = LMCAS::substitute(
        linear, "x", SymbolicExpr::number(3));
    EXPECT_TRUE((substituted.has_value())) << "substitute(x + 2, x => 3) succeeds";
    auto substituted_value =
        substituted ? LMCAS::evalf(*substituted.value())
                    : LMCAS::Result<LMCAS::ApproxReal>::failure(
                          LMCAS::CasErrc::InternalInvariant,
                          "substitution failed", "test");
    EXPECT_TRUE((substituted_value && substituted_value.value().is_finite())) << "substituted Expr can be explicitly evaluated";
    EXPECT_NEAR(substituted_value.value().value, 5.0, 0.0) << "substitute(x + 2, x => 3) evaluates to 5";

    auto typed_binding = LMCAS::binding(
        x.value(), SymbolicExpr::number(3));
    EXPECT_TRUE((typed_binding.has_value())) << "Expr Binding can be constructed from a symbol and Expr value";
    auto typed_substitution = typed_binding
                                  ? LMCAS::substitute(linear, typed_binding.value())
                                  : LMCAS::ExprResult::failure(
                                        LMCAS::CasErrc::InternalInvariant,
                                        "typed binding construction failed", "test");
    EXPECT_TRUE((typed_substitution && substituted &&
                 LMCAS::structurally_equal(
                     *typed_substitution.value(), *substituted.value())))
        << "substitute accepts an Expr Binding without a string variable name";
}

TEST(ExprEvaluation, BatchSubstitution) {
    auto x = LMCAS::sym("x");
    auto typed_binding = LMCAS::binding(x.value(), SymbolicExpr::number(3));
    auto y_symbol = LMCAS::sym("y");
    auto y_binding = LMCAS::binding(
        y_symbol.value(), SymbolicExpr::number(4));
    auto x_plus_y = SymbolicExpr::add(x.value(), y_symbol.value());
    auto batch_substitution = typed_binding && y_binding
                                  ? LMCAS::substitute(
                                        x_plus_y,
                                        std::vector<LMCAS::Binding>{typed_binding.value(),
                                                                    y_binding.value()})
                                  : LMCAS::ExprResult::failure(
                                        LMCAS::CasErrc::InternalInvariant,
                                        "binding list construction failed", "test");
    auto batch_value = batch_substitution
                           ? LMCAS::evalf(*batch_substitution.value())
                           : LMCAS::Result<LMCAS::ApproxReal>::failure(
                                 LMCAS::CasErrc::InternalInvariant,
                                 "binding list substitution failed", "test");
    EXPECT_TRUE((batch_value && batch_value.value().value == 7.0)) << "substitute accepts a deterministic list of Expr bindings";
}

TEST(ExprEvaluation, AbsentAndInvalidBindings) {
    auto x = LMCAS::sym("x");
    auto linear = SymbolicExpr::add(x.value(), SymbolicExpr::number(2));
    auto invalid_binding = LMCAS::binding(
        SymbolicExpr::number(1), SymbolicExpr::number(2));
    EXPECT_TRUE((!invalid_binding &&
                 invalid_binding.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "Expr Binding rejects a non-symbol left-hand side";

    auto unchanged_substitution = LMCAS::substitute(
        linear, "y", SymbolicExpr::number(9));
    EXPECT_TRUE((unchanged_substitution &&
                 LMCAS::structurally_equal(
                     *unchanged_substitution.value(), *linear)))
        << "substituting an absent symbol leaves the Expr unchanged";

    auto null_substitution_expr = LMCAS::substitute(
        nullptr, "x", SymbolicExpr::number(1));
    EXPECT_TRUE((!null_substitution_expr &&
                 null_substitution_expr.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "substitute rejects a null Expr";

    auto empty_substitution_var = LMCAS::substitute(
        linear, "", SymbolicExpr::number(1));
    EXPECT_TRUE((!empty_substitution_var &&
                 empty_substitution_var.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "substitute rejects an empty variable name";

    auto null_substitution_value = LMCAS::substitute(linear, "x", nullptr);
    EXPECT_TRUE((!null_substitution_value &&
                 null_substitution_value.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "substitute rejects a null replacement Expr";
}

TEST(ExprEvaluation, SubstitutionBudget) {
    auto x = LMCAS::sym("x");
    auto linear = SymbolicExpr::add(x.value(), SymbolicExpr::number(2));
    LMCAS::ResourceLimits exhausted_substitution_limits;
    exhausted_substitution_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_substitution_context(
        exhausted_substitution_limits);
    auto exhausted_substitution = LMCAS::substitute(
        linear, "x", SymbolicExpr::number(1), exhausted_substitution_context);
    EXPECT_TRUE((!exhausted_substitution &&
                 exhausted_substitution.error().code ==
                     LMCAS::CasErrc::ResourceLimit))
        << "substitute observes the computation budget";
}

} // namespace
