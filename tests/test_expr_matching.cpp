#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprMatching, CommutativeAddition) {
    auto wildcard_a = SymbolicExpr::variable("A");
    auto match_pattern =
        SymbolicExpr::add(wildcard_a, SymbolicExpr::number(1));
    auto match_target =
        SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::variable("x"));
    auto matched =
        LMCAS::expr_match(match_pattern, match_target, {"A"});
    EXPECT_TRUE((matched && matched.value().matched)) << "expr_match succeeds for a commutative additive pattern";
    EXPECT_TRUE((matched && matched.value().bindings.size() == 1)) << "expr_match returns one wildcard binding";
    EXPECT_TRUE((matched && matched.value().bindings[0].name == "A")) << "expr_match sorts bindings by wildcard name";
    EXPECT_TRUE((matched && matched.value().bindings[0].value &&
                 LMCAS::structurally_equal(
                     *matched.value().bindings[0].value,
                     *SymbolicExpr::variable("x"))))
        << "expr_match binds A to x";
}

TEST(ExprMatching, CommutativeMultiplication) {
    auto multiply_pattern = SymbolicExpr::multiply(SymbolicExpr::variable("B"),
                                                   SymbolicExpr::number(2));
    auto multiply_target = SymbolicExpr::multiply(SymbolicExpr::number(2),
                                                  SymbolicExpr::variable("y"));
    auto multiply_matched =
        LMCAS::expr_match(multiply_pattern, multiply_target, {"B"});
    EXPECT_TRUE((multiply_matched && multiply_matched.value().matched &&
                 multiply_matched.value().bindings.size() == 1 &&
                 LMCAS::structurally_equal(
                     *multiply_matched.value().bindings[0].value,
                     *SymbolicExpr::variable("y"))))
        << "expr_match succeeds for a commutative multiplicative pattern";
}

TEST(ExprMatching, PowerBindings) {
    auto power_pattern = SymbolicExpr::power(SymbolicExpr::variable("U"),
                                             SymbolicExpr::variable("N"));
    auto power_target = SymbolicExpr::power(
        SymbolicExpr::sin(SymbolicExpr::variable("theta")),
        SymbolicExpr::number(2));
    auto power_matched =
        LMCAS::expr_match(power_pattern, power_target, {"N", "U"});
    const char *message =
        "expr_match exposes LMCAS power pattern bindings deterministically";
    if (!power_matched) {
        ADD_FAILURE() << message;
        return;
    }
    const auto &match = power_matched.value();
    if (!match.matched || match.bindings.size() != 2) {
        ADD_FAILURE() << message;
        return;
    }
    EXPECT_TRUE((match.bindings[0].name == "N" &&
                 match.bindings[1].name == "U" &&
                 LMCAS::structurally_equal(
                     *match.bindings[0].value, *SymbolicExpr::number(2)) &&
                 LMCAS::structurally_equal(
                     *match.bindings[1].value,
                     *SymbolicExpr::sin(SymbolicExpr::variable("theta")))))
        << message;
}

TEST(ExprMatching, FunctionBindings) {
    auto function_matched = LMCAS::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("U")),
        SymbolicExpr::sin(SymbolicExpr::add(SymbolicExpr::variable("theta"),
                                            SymbolicExpr::number(1))),
        {"U"});
    EXPECT_TRUE((function_matched && function_matched.value().matched &&
                 function_matched.value().bindings.size() == 1 &&
                 LMCAS::structurally_equal(
                     *function_matched.value().bindings[0].value,
                     *SymbolicExpr::add(SymbolicExpr::variable("theta"),
                                        SymbolicExpr::number(1)))))
        << "expr_match supports function-node patterns";
}

TEST(ExprMatching, RepeatedWildcards) {
    auto repeated_wildcard_match = LMCAS::expr_match(
        SymbolicExpr::add(SymbolicExpr::variable("A"),
                          SymbolicExpr::variable("A")),
        SymbolicExpr::add(SymbolicExpr::variable("z"),
                          SymbolicExpr::variable("z")),
        {"A"});
    auto inconsistent_wildcard_match = LMCAS::expr_match(
        SymbolicExpr::add(SymbolicExpr::variable("A"),
                          SymbolicExpr::variable("A")),
        SymbolicExpr::add(SymbolicExpr::variable("z"),
                          SymbolicExpr::variable("w")),
        {"A"});
    EXPECT_TRUE((repeated_wildcard_match &&
                 repeated_wildcard_match.value().matched))
        << "expr_match accepts repeated wildcards with identical bindings";
    EXPECT_TRUE((inconsistent_wildcard_match &&
                 !inconsistent_wildcard_match.value().matched))
        << "expr_match rejects repeated wildcards with inconsistent bindings";
}

TEST(ExprMatching, StructuralNonmatch) {
    auto unmatched = LMCAS::expr_match(
        SymbolicExpr::sin(SymbolicExpr::variable("A")),
        SymbolicExpr::cos(SymbolicExpr::variable("x")), {"A"});
    EXPECT_TRUE((unmatched && !unmatched.value().matched &&
                 unmatched.value().bindings.empty()))
        << "expr_match reports structural non-matches without error";
}

TEST(ExprMatching, InvalidMatchInputs) {
    auto wildcard_a = SymbolicExpr::variable("A");
    auto match_pattern =
        SymbolicExpr::add(wildcard_a, SymbolicExpr::number(1));
    auto match_target =
        SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::variable("x"));
    auto null_pattern_match =
        LMCAS::expr_match(nullptr, match_target, {"A"});
    auto null_target_match =
        LMCAS::expr_match(match_pattern, nullptr, {"A"});
    auto empty_wildcard_match =
        LMCAS::expr_match(match_pattern, match_target, {""});
    auto duplicate_wildcard_match =
        LMCAS::expr_match(match_pattern, match_target, {"A", "A"});
    EXPECT_TRUE((!null_pattern_match &&
                 null_pattern_match.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "expr_match rejects a null pattern";
    EXPECT_TRUE((!null_target_match &&
                 null_target_match.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "expr_match rejects a null target";
    EXPECT_TRUE((!empty_wildcard_match &&
                 empty_wildcard_match.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "expr_match rejects an empty wildcard name";
    EXPECT_TRUE((!duplicate_wildcard_match &&
                 duplicate_wildcard_match.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "expr_match rejects duplicate wildcard names";
}

TEST(ExprMatching, MatchingBudget) {
    auto wildcard_a = SymbolicExpr::variable("A");
    auto match_pattern =
        SymbolicExpr::add(wildcard_a, SymbolicExpr::number(1));
    auto match_target =
        SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::variable("x"));
    LMCAS::ResourceLimits exhausted_match_limits;
    exhausted_match_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_match_context(exhausted_match_limits);
    auto exhausted_match = LMCAS::expr_match(
        match_pattern, match_target, {"A"}, exhausted_match_context);
    EXPECT_TRUE((!exhausted_match &&
                 exhausted_match.error().code ==
                     LMCAS::CasErrc::ResourceLimit))
        << "expr_match observes the computation budget";
}

} // namespace
