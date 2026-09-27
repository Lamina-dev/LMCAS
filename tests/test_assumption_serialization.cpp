
#include "test_common.hpp"
#include <rapidcheck.h>
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "computation_context.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "expr.hpp"
#include "query_interface.hpp"
#include "residual_verification.hpp"
#include <cstring>
#include <limits>
#include <string>
#include <vector>
#include <memory>
#include <variant>

using namespace LMCAS;

static AssumptionContext deserialize_success(const std::string &data) {
    auto result = AssumptionContext::deserialize(data);
    RC_ASSERT(result.has_value());
    return std::move(result.value());
}

/// Generate a random Domain (excluding Complex which is the default).
static Domain random_domain() {
    std::vector<Domain> domains = {
        Domain::Real, Domain::Algebraic, Domain::Rational,
        Domain::Integer, Domain::Natural, Domain::PositiveInt};
    return *rc::gen::elementOf(domains);
}

/// Generate a random Sign.
static Sign random_sign() {
    std::vector<Sign> signs = {
        Sign::Positive, Sign::Negative, Sign::NonNegative,
        Sign::NonPositive, Sign::NonZero};
    return *rc::gen::elementOf(signs);
}

/// Generate a random Parity.
static Parity random_parity() {
    std::vector<Parity> parities = {Parity::Even, Parity::Odd};
    return *rc::gen::elementOf(parities);
}

/// Generate a random Boundedness.
static Boundedness random_boundedness() {
    std::vector<Boundedness> values = {Boundedness::Bounded, Boundedness::Unbounded};
    return *rc::gen::elementOf(values);
}

/// Generate a random Finiteness (non-Unknown).
static Finiteness random_finiteness() {
    std::vector<Finiteness> values = {Finiteness::Finite, Finiteness::Divergent};
    return *rc::gen::elementOf(values);
}

/// Generate a random Definiteness (non-Unknown).
static Definiteness random_definiteness() {
    std::vector<Definiteness> values = {
        Definiteness::PositiveDefinite, Definiteness::PositiveSemiDefinite,
        Definiteness::NegativeDefinite, Definiteness::NegativeSemiDefinite,
        Definiteness::Indefinite};
    return *rc::gen::elementOf(values);
}

/// Generate a unique variable name.
static std::string random_var_name(int prefix_id) {
    return "s" + std::to_string(prefix_id) + "_" + std::to_string(*rc::gen::inRange(0, (99) + 1));
}

TEST(LmcasAssumptionSerialization, SerializationDomainRoundtrip) {
    EXPECT_TRUE(rc::check("Domain declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(0);
        Domain dom = random_domain();

        RC_ASSERT(ctx.assume_domain(var, dom).has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.get_domain(var) == dom);
        RC_ASSERT(restored.has_domain(var, dom));
    }));
}

TEST(LmcasAssumptionSerialization, SerializationSignRoundtrip) {
    EXPECT_TRUE(rc::check("Sign declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(1);
        Sign sign = random_sign();

        RC_ASSERT(ctx.assume_sign(var, sign).has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.has_sign(var, sign));
    }));
}

TEST(LmcasAssumptionSerialization, SerializationParityRoundtrip) {
    EXPECT_TRUE(rc::check("Parity declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(2);
        Parity par = random_parity();

        RC_ASSERT(
            ctx.current_properties().declare_parity(var, par).has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.get_parity(var) == par);
    }));
}

TEST(LmcasAssumptionSerialization, SerializationBoundednessRoundtrip) {
    EXPECT_TRUE(rc::check("Boundedness declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(3);
        Boundedness bnd = random_boundedness();

        RC_ASSERT(
            ctx.current_properties().declare_bounded(var, bnd).has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.get_boundedness(var) == bnd);
    }));
}

TEST(LmcasAssumptionSerialization, SerializationTranscendentalRoundtrip) {
    EXPECT_TRUE(rc::check("Transcendental declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(4);

        RC_ASSERT(
            ctx.current_properties().declare_transcendental(var).has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        // Transcendental implies Real domain
        RC_ASSERT(restored.get_domain(var) == Domain::Real);
        RC_ASSERT(restored.current_properties().is_transcendental(var));
    }));
}

TEST(LmcasAssumptionSerialization, SerializationFinitenessRoundtrip) {
    EXPECT_TRUE(rc::check("Finiteness declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(5);
        Finiteness fin = random_finiteness();

        RC_ASSERT(
            ctx.current_properties().declare_finiteness(var, fin).has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.current_properties().get_finiteness(var) == fin);

        // If Finite, boundedness should also be Bounded
        if (fin == Finiteness::Finite) {
            RC_ASSERT(restored.get_boundedness(var) == Boundedness::Bounded);
        }
    }));
}

TEST(LmcasAssumptionSerialization, SerializationDefinitenessRoundtrip) {
    EXPECT_TRUE(rc::check("Definiteness declarations survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(6);
        Definiteness def = random_definiteness();

        RC_ASSERT(
            ctx.current_properties().declare_definiteness(var, def).has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.current_properties().get_definiteness(var) == def);
    }));
}

TEST(LmcasAssumptionSerialization, SerializationMultiScopeRoundtrip) {
    EXPECT_TRUE(rc::check("Multi-scope contexts survive serialize/deserialize round-trip", []() {
        AssumptionContext ctx;

        std::string var_a = "a_" + std::to_string(*rc::gen::inRange(0, (99) + 1));
        Domain dom_a = random_domain();
        RC_ASSERT(ctx.assume_domain(var_a, dom_a).has_value());

        ctx.push();
        std::string var_b = "b_" + std::to_string(*rc::gen::inRange(0, (99) + 1));
        Sign sign_b = random_sign();
        RC_ASSERT(ctx.assume_sign(var_b, sign_b).has_value());

        int depth_before = ctx.depth();

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.depth() == depth_before);

        // Root scope properties should be visible (read-through)
        RC_ASSERT(restored.get_domain(var_a) == dom_a);

        RC_ASSERT(restored.has_sign(var_b, sign_b));
    }));
}

TEST(LmcasAssumptionSerialization, SerializationCombinedPropertiesRoundtrip) {
    EXPECT_TRUE(rc::check("Multiple properties on a single symbol survive round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(7);

        // Declare domain (must be compatible with sign)
        // Use Integer domain with Positive sign (compatible)
        RC_ASSERT(ctx.assume_domain(var, Domain::Integer).has_value());
        RC_ASSERT(ctx.assume_sign(var, Sign::Positive).has_value());
        RC_ASSERT(
            ctx.current_properties().declare_parity(var, Parity::Odd).has_value());
        RC_ASSERT(
            ctx.current_properties()
                .declare_bounded(var, Boundedness::Bounded)
                .has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.get_domain(var) == Domain::Integer);
        RC_ASSERT(restored.has_sign(var, Sign::Positive));
        RC_ASSERT(restored.get_parity(var) == Parity::Odd);
        RC_ASSERT(restored.get_boundedness(var) == Boundedness::Bounded);
    }));
}

TEST(LmcasAssumptionSerialization, SerializationSimpleRelationRoundtrip) {
    EXPECT_TRUE(rc::check("Simple variable-vs-zero relations survive round-trip", []() {
        AssumptionContext ctx;
        std::string var = random_var_name(8);

        auto var_node = LMCAS::detail::make_node<VariableNode>(var);
        auto zero_node = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        auto rel_node = LMCAS::detail::make_node<RelationalNode>(
            var_node, zero_node, RelationalNode::Op::GT);
        auto rel_expr = LMCAS::detail::expression_from_node(rel_node);
        RC_ASSERT(ctx.assume(rel_expr).has_value());

        AssumptionContext restored = deserialize_success(ctx.serialize());

        // The sign property derived from the relation should be preserved
        RC_ASSERT(restored.has_sign(var, Sign::Positive));
    }));
}

TEST(LmcasAssumptionSerialization, SerializationEmptyRoundtrip) {
    EXPECT_TRUE(rc::check("Empty context survives serialize/deserialize round-trip", []() {
        AssumptionContext ctx;

        AssumptionContext restored = deserialize_success(ctx.serialize());

        RC_ASSERT(restored.depth() == ctx.depth());
        // No properties should be set
        std::string var = "nonexistent";
        RC_ASSERT(restored.get_domain(var) == Domain::Complex);
    }));
}

TEST(LmcasAssumptionSerialization, CheckedDeserializationContractsRelationErrors) {
    auto malformed = AssumptionContext::deserialize_checked("SCOPE 0\nRELATION x 0\nEND\n");
    EXPECT_TRUE(!malformed.has_value()) << "checked deserialize rejects malformed relation";
    EXPECT_TRUE(malformed.error().code == CasErrc::ParseError) << "checked deserialize reports ParseError for malformed input";

    auto contradictory = AssumptionContext::deserialize_checked(
        "SCOPE 0\nSIGN x Positive\nRELATION x LT 0\nEND\n");
    EXPECT_TRUE(!contradictory.has_value()) << "checked deserialize rejects contradictory derived relation property";
    EXPECT_TRUE(contradictory.error().code == CasErrc::ParseError) << "checked deserialize maps relation contradiction to ParseError";

    auto property_contradiction = AssumptionContext::deserialize_checked(
        "SCOPE 0\nDOMAIN n Natural\nSIGN n Negative\nEND\n");
    EXPECT_TRUE(!property_contradiction.has_value()) << "checked deserialize rejects contradictory property declarations";
    EXPECT_TRUE(property_contradiction.error().code == CasErrc::ParseError) << "checked deserialize maps property contradiction to ParseError";

    auto canonical = AssumptionContext::deserialize(
        "SCOPE 0\nDOMAIN n Natural\nSIGN n Negative\nEND\n");
    EXPECT_TRUE(!canonical.has_value()) << "canonical deserialize returns the checked ParseError";
    EXPECT_TRUE(canonical.error().code == CasErrc::ParseError) << "canonical deserialize preserves ParseError";
}

TEST(LmcasAssumptionSerialization, CheckedDeserializationContractsIntervals) {
    auto exact_interval = AssumptionContext::deserialize_checked(
        "SCOPE 0\nCONTINUOUS f [1/3, 2/3)\nEND\n");
    ASSERT_TRUE(exact_interval) << "checked deserialize accepts exact rational interval endpoints";
    const auto declarations =
        exact_interval.value().current_properties().get_continuity_decls("f");
    ASSERT_EQ(declarations.size(), 1u);
    ASSERT_TRUE(declarations.front().interval.lower.value);
    ASSERT_TRUE(declarations.front().interval.upper.value);
    const auto *lower = dynamic_cast<const NumberNode *>(
        detail::node(*declarations.front().interval.lower.value).get());
    const auto *upper = dynamic_cast<const NumberNode *>(
        detail::node(*declarations.front().interval.upper.value).get());
    ASSERT_NE(lower, nullptr);
    ASSERT_NE(upper, nullptr);
    ASSERT_TRUE(std::holds_alternative<Rational>(lower->value()));
    ASSERT_TRUE(std::holds_alternative<Rational>(upper->value()));
    EXPECT_EQ(std::get<Rational>(lower->value()), Rational(1, 3));
    EXPECT_EQ(std::get<Rational>(upper->value()), Rational(2, 3));
    EXPECT_FALSE(declarations.front().interval.lower.is_open);
    EXPECT_TRUE(declarations.front().interval.upper.is_open);

    auto malformed_interval = AssumptionContext::deserialize_checked(
        "SCOPE 0\nCONTINUOUS f [0, 1}\nEND\n");
    EXPECT_TRUE(!malformed_interval.has_value()) << "checked deserialize rejects malformed interval delimiters";

    auto trailing_endpoint = AssumptionContext::deserialize_checked(
        "SCOPE 0\nCONTINUOUS f [0junk, 1]\nEND\n");
    EXPECT_TRUE(!trailing_endpoint.has_value()) << "checked deserialize rejects trailing endpoint text";

    auto nonfinite_endpoint = AssumptionContext::deserialize_checked(
        "SCOPE 0\nCONTINUOUS f [nan, 1]\nEND\n");
    EXPECT_TRUE(!nonfinite_endpoint.has_value()) << "checked deserialize rejects non-finite interval endpoints";

    auto closed_infinity = AssumptionContext::deserialize_checked(
        "SCOPE 0\nCONTINUOUS f [-inf, 1]\nEND\n");
    EXPECT_TRUE(!closed_infinity.has_value()) << "checked deserialize rejects a closed infinite endpoint";
}

TEST(LmcasAssumptionSerialization, CheckedDeserializationContractsExactAtoms) {
    auto exact_atoms = AssumptionContext::deserialize_checked(
        "SCOPE 0\n"
        "PERIODIC f x 1/3\n"
        "RELATION x GT 1/3\n"
        "CONDITIONAL (x GT 1/3) => (y LT 2/3)\n"
        "END\n");
    EXPECT_TRUE(exact_atoms.has_value()) << "checked deserialize accepts exact rational atoms";
    if (!exact_atoms) {
        return;
    }
    auto restored = AssumptionContext::deserialize_checked(exact_atoms.value().serialize());
    EXPECT_TRUE(restored.has_value()) << "exact atom serialization reparses";
    if (!restored) {
        return;
    }
    auto period = restored.value().get_period("f", "x");
    auto third = SymbolicExpr::number(Rational(1, 3));
    auto canonical_period = period ? period->simplify() : nullptr;
    EXPECT_TRUE(canonical_period && detail::node(canonical_period)->equals(*detail::node(third))) << "declared period retains its exact rational value";
    EXPECT_TRUE(restored.value().has_relation(*SymbolicExpr::variable("x"), *third, RelationOp::GT)) << "exact relational fact survives the roundtrip";
    const auto conditionals = restored.value().get_active_conditionals();
    auto condition = gt(SymbolicExpr::variable("x"), third);
    auto conclusion = lt(SymbolicExpr::variable("y"), SymbolicExpr::number(Rational(2, 3)));
    EXPECT_TRUE(conditionals.size() == 1 && condition && conclusion &&
                detail::node(conditionals.front().condition)->equals(*detail::node(*condition.value())) &&
                detail::node(conditionals.front().conclusion)->equals(*detail::node(*conclusion.value())))
        << "conditional predicates retain their exact numeric atoms";
}

TEST(LmcasAssumptionSerialization, CheckedDeserializationContractsScopeGrammar) {
    auto noncanonical_scope = AssumptionContext::deserialize_checked(
        "SCOPE 2\nEND\n");
    EXPECT_TRUE(!noncanonical_scope.has_value()) << "checked deserialize rejects non-sequential scopes";

    auto trailing_field = AssumptionContext::deserialize_checked(
        "SCOPE 0\nDOMAIN x Real ignored\nEND\n");
    EXPECT_TRUE(!trailing_field.has_value()) << "checked deserialize rejects trailing declaration fields";

    auto data_after_end = AssumptionContext::deserialize_checked(
        "SCOPE 0\nEND\nSIGN x Positive\n");
    EXPECT_TRUE(!data_after_end.has_value()) << "checked deserialize rejects data after END";
}

TEST(LmcasAssumptionSerialization, CheckedDeserializationContracts) {
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value()) << "test setup accepts positive sign";
    auto rel_node = LMCAS::detail::make_node<RelationalNode>(
        LMCAS::detail::make_node<VariableNode>("x"),
        LMCAS::detail::make_node<NumberNode>(BigInt(0)),
        RelationalNode::Op::GT);
    auto relation = LMCAS::detail::expression_from_node(rel_node);
    EXPECT_TRUE(ctx.assume(relation).has_value()) << "test setup accepts positive relation";

    auto ok = AssumptionContext::deserialize_checked(ctx.serialize());
    EXPECT_TRUE(ok.has_value()) << "checked deserialize accepts serialized context";
    if (ok) {
        EXPECT_TRUE(ok.value().has_sign("x", Sign::Positive)) << "checked deserialize preserves positive sign";
    }
}

TEST(LmcasAssumptionSerialization, SerializedVariablePeriods) {
    auto parsed = AssumptionContext::deserialize_checked(
        "SCOPE 0\nPERIODIC f x 4*pi\nPERIODIC f y 7\n"
        "SCOPE 1\nPERIODIC f x 5\nEND\n");
    EXPECT_TRUE(parsed.has_value()) << "composite periods parse";
    if (!parsed) {
        return;
    }
    auto restored = AssumptionContext::deserialize_checked(parsed.value().serialize());
    EXPECT_TRUE(restored.has_value()) << "scoped periods reparse";
    if (!restored) {
        return;
    }
    auto f = SymbolicExpr::variable("f");
    QueryInterface query(restored.value());
    const auto check = [&](const char *variable, const char *expected_source) {
        auto actual = query.get_period(*f, variable);
        auto expected = parse_expr(expected_source);
        EXPECT_TRUE(actual && actual.value() && expected) << "period and exact reference exist";
        if (!actual || !actual.value() || !expected) {
            return;
        }
        ComputationContext context;
        auto equal = check_equivalent(detail::make_expression_ptr(*actual.value()), expected.value(), context);
        EXPECT_TRUE(equal && std::holds_alternative<ProvedZeroResidual>(equal.value())) << "serialized period preserves its mathematical value";
    };
    check("x", "5");
    check("y", "7");
    EXPECT_TRUE(restored.value().pop().has_value()) << "serialized child scope pops";
    check("x", "4*pi");
    check("y", "7");
    auto legacy = AssumptionContext::deserialize_checked("SCOPE 0\nPERIODIC f 5\nEND\n");
    EXPECT_TRUE(!legacy && legacy.error().code == CasErrc::ParseError) << "legacy variableless period declarations fail rather than guessing x";
}
