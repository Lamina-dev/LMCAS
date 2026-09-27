
#include "test_common.hpp"
#include "relation_store.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "assumption_context.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "query_interface.hpp"
#include "computation_context.hpp"
#include <algorithm>

using namespace LMCAS;

/// Helper: create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var(const std::string &name) {
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
}

/// Helper: create a SymbolicExpr wrapping a NumberNode with value 0.
static SymbolicExpr make_zero() {
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(0)));
}

TEST(AssumptionRelationExt, Reversed0LtVarPositive) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("x");

    // 0 < x means x > 0 -> Positive
    EXPECT_TRUE((rs.add_relation(zero, var, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("x", Sign::Positive))) << "0 LT x should derive Positive for x";
    EXPECT_TRUE((ps.has_sign("x", Sign::NonNegative))) << "Positive implies NonNegative";
    EXPECT_TRUE((ps.has_sign("x", Sign::NonZero))) << "Positive implies NonZero";
}

TEST(AssumptionRelationExt, Reversed0GtVarNegative) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("y");

    // 0 > y means y < 0 -> Negative
    EXPECT_TRUE((rs.add_relation(zero, var, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("y", Sign::Negative))) << "0 GT y should derive Negative for y";
    EXPECT_TRUE((ps.has_sign("y", Sign::NonPositive))) << "Negative implies NonPositive";
    EXPECT_TRUE((ps.has_sign("y", Sign::NonZero))) << "Negative implies NonZero";
}

TEST(AssumptionRelationExt, Reversed0GeqVarNonpositive) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("z");

    // 0 >= z means z <= 0 -> NonPositive
    EXPECT_TRUE((rs.add_relation(zero, var, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("z", Sign::NonPositive))) << "0 GEQ z should derive NonPositive for z";
}

TEST(AssumptionRelationExt, Reversed0LeqVarNonnegative) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("w");

    // 0 <= w means w >= 0 -> NonNegative
    EXPECT_TRUE((rs.add_relation(zero, var, RelationalNode::Op::LEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("w", Sign::NonNegative))) << "0 LEQ w should derive NonNegative for w";
}

TEST(AssumptionRelationExt, Reversed0NeqVarNonzero) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr zero = make_zero();
    SymbolicExpr var = make_var("a");

    // 0 != a means a != 0 -> NonZero
    EXPECT_TRUE((rs.add_relation(zero, var, RelationalNode::Op::NEQ, ps).has_value())) << "relation insertion succeeds";

    EXPECT_TRUE((ps.has_sign("a", Sign::NonZero))) << "0 NEQ a should derive NonZero for a";
}

TEST(AssumptionRelationExt, ReversedAllOperatorsComprehensive) {
    struct TestCase {
        RelationalNode::Op op;
        Sign expected_sign;
        std::string desc;
    };

    std::vector<TestCase> cases = {
        {RelationalNode::Op::LT, Sign::Positive, "0 LT var → Positive"},
        {RelationalNode::Op::GT, Sign::Negative, "0 GT var → Negative"},
        {RelationalNode::Op::GEQ, Sign::NonPositive, "0 GEQ var → NonPositive"},
        {RelationalNode::Op::LEQ, Sign::NonNegative, "0 LEQ var → NonNegative"},
        {RelationalNode::Op::NEQ, Sign::NonZero, "0 NEQ var → NonZero"},
    };

    for (const auto &tc : cases) {
        RelationStore rs;
        PropertyStore ps;

        SymbolicExpr zero = make_zero();
        SymbolicExpr var = make_var("v");

        EXPECT_TRUE((rs.add_relation(zero, var, tc.op, ps).has_value())) << "relation insertion succeeds";

        EXPECT_TRUE((ps.has_sign("v", tc.expected_sign))) << tc.desc;
    }
}

TEST(AssumptionRelationExt, TransitiveChain3Gt) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    // Transitive closure should deduce x GT z
    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GT))) << "x GT y, y GT z → x GT z should be deduced";
}

TEST(AssumptionRelationExt, TransitiveChainGeqGt) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    // GEQ + GT -> GT
    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GT))) << "x GEQ y, y GT z → x GT z should be deduced";
}

TEST(AssumptionRelationExt, TransitiveChainGtGeq) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";

    // GT + GEQ -> GT
    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GT))) << "x GT y, y GEQ z → x GT z should be deduced";
}

TEST(AssumptionRelationExt, TransitiveChainGeqGeq) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::GEQ, ps).has_value())) << "relation insertion succeeds";

    // GEQ + GEQ -> GEQ
    EXPECT_TRUE((rs.has_relation(x, z, RelationalNode::Op::GEQ))) << "x GEQ y, y GEQ z → x GEQ z should be deduced";
}

TEST(AssumptionRelationExt, TransitiveChain4Variables) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr a = make_var("a");
    SymbolicExpr b = make_var("b");
    SymbolicExpr c = make_var("c");
    SymbolicExpr d = make_var("d");

    EXPECT_TRUE((rs.add_relation(a, b, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(b, c, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(c, d, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";

    // Should deduce a GT c, a GT d, b GT d
    EXPECT_TRUE((rs.has_relation(a, c, RelationalNode::Op::GT))) << "a GT b, b GT c → a GT c";
    EXPECT_TRUE((rs.has_relation(a, d, RelationalNode::Op::GT))) << "a GT c, c GT d → a GT d";
    EXPECT_TRUE((rs.has_relation(b, d, RelationalNode::Op::GT))) << "b GT c, c GT d → b GT d";
}

TEST(AssumptionRelationExt, TransitiveNoClosureForLt) {
    RelationStore rs;
    PropertyStore ps;

    SymbolicExpr x = make_var("x");
    SymbolicExpr y = make_var("y");
    SymbolicExpr z = make_var("z");

    // LT is not a transitive operator in the implementation (only GT/GEQ are)
    EXPECT_TRUE((rs.add_relation(x, y, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";
    EXPECT_TRUE((rs.add_relation(y, z, RelationalNode::Op::LT, ps).has_value())) << "relation insertion succeeds";

    // No transitive deduction for LT (the implementation only handles GT/GEQ)
    EXPECT_FALSE((rs.has_relation(x, z, RelationalNode::Op::LT))) << "LT does not participate in transitive closure";
}

TEST(AssumptionRelationExt, TransitiveCap64) {
    RelationStore rs;
    PropertyStore ps;

    // Create a long chain: v0 GT v1, v1 GT v2, ..., v(N-1) GT vN
    // With N variables in a chain, adding the last relation can trigger many
    // transitive deductions. We need enough variables so that the total
    // deductions from a single add_relation call would exceed 64.
    //
    // Strategy: build a chain of 12 variables first (v0..v11), then add one
    // more relation that connects to the chain. With 12 existing nodes,
    // adding v12 would only produce 12 new relations. Instead, we build
    // a large chain incrementally and check that the total stored relations
    // are bounded.
    //
    // Actually, the cap is per add_relation call. Let's create a scenario
    // where a single add_relation triggers many deductions:
    // First, add many independent chains that share an endpoint.

    // Build a star pattern: many variables all GT than a central variable "center"
    // Then add "center GT bottom" - this should trigger deductions for all
    // star variables GT bottom.

    // With 70 star variables + center GT bottom, the single add_relation(center, bottom, GT)
    // would try to deduce 70 new relations (star_i GT bottom), but cap at 64.

    SymbolicExpr center = make_var("center");
    SymbolicExpr bottom = make_var("bottom");

    // Add 70 relations: star_i GT center
    const int num_star = 70;
    for (int i = 0; i < num_star; ++i) {
        SymbolicExpr star_i = make_var("star_" + std::to_string(i));
        EXPECT_TRUE((rs.add_relation(star_i, center, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    }

    // Now add: center GT bottom
    // This should trigger backward chaining: for each star_i GT center,
    // deduce star_i GT bottom. But cap at 64.
    size_t relations_before = rs.get_relations().size();
    EXPECT_TRUE((rs.add_relation(center, bottom, RelationalNode::Op::GT, ps).has_value())) << "relation insertion succeeds";
    size_t relations_after = rs.get_relations().size();

    // The new relations added should be: 1 (center GT bottom) + at most 64 deduced
    size_t new_relations = relations_after - relations_before;
    EXPECT_TRUE((new_relations <= 65)) << "Should add at most 1 + 64 = 65 new relations (cap at 64 deductions)";

    // Verify that at least some deductions were made
    EXPECT_TRUE((new_relations > 1)) << "Should have deduced at least some transitive relations";

    // Verify the cap was actually hit (70 potential deductions, only 64 allowed)
    // The 1 is for the directly added relation, the rest are deductions
    size_t deductions = new_relations - 1;
    EXPECT_TRUE((deductions <= 64)) << "Deductions should be capped at 64";
    EXPECT_TRUE((deductions == 64)) << "With 70 potential deductions, exactly 64 should be produced";

    // Verify that some star variables DO have the deduced relation
    int found_count = 0;
    for (int i = 0; i < num_star; ++i) {
        SymbolicExpr star_i = make_var("star_" + std::to_string(i));
        if (rs.has_relation(star_i, bottom, RelationalNode::Op::GT)) {
            ++found_count;
        }
    }
    EXPECT_TRUE((found_count == 64)) << "Exactly 64 star variables should have deduced GT bottom relation";

    /// 闭包上限触发后,部分星形节点保持关系未知.
    EXPECT_TRUE((found_count < num_star)) << "Not all 70 star variables should have the deduced relation (cap hit)";
}

TEST(AssumptionRelationExt, ExistingConclusionGetsAlternativeAfterNewConclusionCap) {
    AssumptionContext ctx;
    const auto a = make_var("a");
    const auto b = make_var("b");
    const auto c = make_var("c");
    const auto d = make_var("d");
    auto& relations = ctx.current_relations();
    auto& properties = ctx.current_properties();
    ASSERT_TRUE(relations.add_relation(a, b, RelationOp::GT, properties));
    ASSERT_TRUE(relations.add_relation(b, c, RelationOp::GT, properties));
    ASSERT_TRUE(ctx.has_relation(a, c, RelationOp::GT));

    for (int i = 0; i < 65; ++i) {
        ASSERT_TRUE(relations.add_relation(d, make_var("target_" + std::to_string(i)),
                                           RelationOp::GT, properties));
    }
    // This route is scanned only after all 65 potential new conclusions.
    ASSERT_TRUE(relations.add_relation(d, c, RelationOp::GT, properties));
    const auto before = relations.get_relations().size();
    ASSERT_TRUE(relations.add_relation(a, d, RelationOp::GT, properties));
    EXPECT_EQ(relations.get_relations().size() - before, 65u);
    for (int i = 0; i < 64; ++i) {
        EXPECT_TRUE(ctx.has_relation(a, make_var("target_" + std::to_string(i)),
                                    RelationOp::GT));
    }
    EXPECT_FALSE(ctx.has_relation(a, make_var("target_64"), RelationOp::GT));

    ctx.push();
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative));
    EXPECT_FALSE(ctx.has_relation(a, b, RelationOp::GT));
    EXPECT_TRUE(ctx.has_relation(a, c, RelationOp::GT));
    ASSERT_TRUE(ctx.assume_sign("d", Sign::Negative));
    EXPECT_FALSE(ctx.has_relation(a, c, RelationOp::GT));
    ASSERT_TRUE(ctx.pop());
    EXPECT_TRUE(ctx.has_relation(a, c, RelationOp::GT));
}

TEST(AssumptionScopeShadowing, AlternativeProofsPreserveDiamondAndDownstream) {
    AssumptionContext ctx;
    const auto a = SymbolicExpr::add(
        SymbolicExpr::variable("x"), SymbolicExpr::variable("z"));
    const auto b = SymbolicExpr::variable("b");
    const auto c = SymbolicExpr::variable("c");
    const auto d = SymbolicExpr::variable("d");
    const auto downstream = SymbolicExpr::add(
        SymbolicExpr::variable("u"), SymbolicExpr::variable("v"));
    const auto one = SymbolicExpr::number(1);
    for (const auto& endpoints : {
             std::make_pair(a, b), std::make_pair(b, c),
             std::make_pair(a, d), std::make_pair(d, c),
             std::make_pair(downstream, a), std::make_pair(c, one)}) {
        ASSERT_TRUE(ctx.current_relations().add_relation(
            *endpoints.first, *endpoints.second, RelationOp::GT, ctx.current_properties()));
    }
    QueryInterface query(ctx);
    auto positive = query.query_positive(*downstream);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::True);
    ctx.push();
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative));
    EXPECT_TRUE(ctx.has_relation(*a, *c, RelationOp::GT));
    EXPECT_TRUE(ctx.has_relation(*downstream, *c, RelationOp::GT));
    positive = query.query_positive(*downstream);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::True);
    ctx.push();
    ASSERT_TRUE(ctx.assume_sign("d", Sign::Negative));
    EXPECT_FALSE(ctx.has_relation(*a, *c, RelationOp::GT));
    EXPECT_FALSE(ctx.has_relation(*downstream, *c, RelationOp::GT));
    const auto visible = ctx.get_visible_relations();
    EXPECT_FALSE(std::any_of(visible.begin(), visible.end(), [&](const Relation& relation) {
        return relation.op == RelationOp::GT &&
            detail::node(relation.lhs)->equals(*detail::node(*a)) &&
            detail::node(relation.rhs)->equals(*detail::node(*c));
    }));
    positive = query.query_positive(*downstream);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::Unknown);
    ASSERT_TRUE(ctx.pop());
    EXPECT_TRUE(ctx.has_relation(*a, *c, RelationOp::GT));
    positive = query.query_positive(*downstream);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::True);
    ASSERT_TRUE(ctx.pop());

    ASSERT_TRUE(ctx.current_relations().add_relation(
        *a, *c, RelationOp::GT, ctx.current_properties()));
    ctx.push();
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative));
    ASSERT_TRUE(ctx.assume_sign("d", Sign::Negative));
    EXPECT_TRUE(ctx.has_relation(*a, *c, RelationOp::GT));
    EXPECT_TRUE(ctx.has_relation(*downstream, *c, RelationOp::GT));
    positive = query.query_positive(*downstream);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::True);
}

TEST(AssumptionScopeShadowing, ProofCyclesCannotReplaceHiddenDeclaredRoots) {
    AssumptionContext ctx;
    const auto a = SymbolicExpr::add(
        SymbolicExpr::variable("x"), SymbolicExpr::variable("z"));
    const auto b = SymbolicExpr::variable("b");
    const auto c = SymbolicExpr::variable("c");
    const auto d = SymbolicExpr::variable("d");
    const auto one = SymbolicExpr::number(1);
    for (const auto& endpoints : {
             std::make_pair(a, b), std::make_pair(b, c),
             std::make_pair(c, d), std::make_pair(d, a)}) {
        ASSERT_TRUE(ctx.current_relations().add_relation(
            *endpoints.first, *endpoints.second, RelationOp::GEQ, ctx.current_properties()));
    }
    ASSERT_TRUE(ctx.current_relations().add_relation(
        *b, *one, RelationOp::GT, ctx.current_properties()));
    QueryInterface query(ctx);
    auto positive = query.query_positive(*a);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::True);
    ctx.push();
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative));
    ASSERT_TRUE(ctx.assume_sign("d", Sign::Negative));
    EXPECT_FALSE(ctx.has_relation(*a, *c, RelationOp::GEQ));
    EXPECT_FALSE(ctx.has_relation(*c, *a, RelationOp::GEQ));
    EXPECT_FALSE(ctx.has_relation(*a, *one, RelationOp::GT));
    EXPECT_TRUE(ctx.get_visible_relations().empty());
    positive = query.query_positive(*a);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::Unknown);
    ASSERT_TRUE(ctx.pop());
    positive = query.query_positive(*a);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::True);
}

template <typename T>
static void expect_relation_query_error(const Result<T>& result, CasErrc code) {
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, code);
}

TEST(AssumptionScopeShadowing, CheckedDiamondQueriesShareCallerBudget) {
    AssumptionContext ctx;
    const auto a = SymbolicExpr::add(
        SymbolicExpr::variable("x"), SymbolicExpr::variable("z"));
    const auto b = SymbolicExpr::variable("b");
    const auto c = SymbolicExpr::variable("c");
    const auto d = SymbolicExpr::variable("d");
    const auto one = SymbolicExpr::number(1);
    for (const auto& endpoints : {
             std::make_pair(a, b), std::make_pair(b, c),
             std::make_pair(a, d), std::make_pair(d, c),
             std::make_pair(c, one)}) {
        ASSERT_TRUE(ctx.current_relations().add_relation(
            *endpoints.first, *endpoints.second, RelationOp::GT, ctx.current_properties()));
    }
    ctx.push();
    ASSERT_TRUE(ctx.assume_sign("b", Sign::Negative));
    ResourceLimits limits;
    limits.max_steps = 20;
    ComputationContext limited_relation(limits);
    expect_relation_query_error(
        ctx.has_relation_checked(*a, *c, RelationOp::GT, limited_relation), CasErrc::ResourceLimit);
    ComputationContext limited_visible(limits);
    expect_relation_query_error(
        ctx.get_visible_relations_checked(limited_visible), CasErrc::ResourceLimit);
    ComputationContext limited_inference(limits);
    expect_relation_query_error(
        ctx.is_positive_checked(*a, limited_inference), CasErrc::ResourceLimit);

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    expect_relation_query_error(
        ctx.has_relation_checked(*a, *c, RelationOp::GT, cancelled), CasErrc::Cancelled);
    expect_relation_query_error(ctx.get_visible_relations_checked(cancelled), CasErrc::Cancelled);
    expect_relation_query_error(ctx.is_positive_checked(*a, cancelled), CasErrc::Cancelled);

    ResourceLimits read_only_limits;
    read_only_limits.max_ast_nodes = 0;
    ComputationContext ordinary(read_only_limits);
    auto relation = ctx.has_relation_checked(*a, *c, RelationOp::GT, ordinary);
    ASSERT_TRUE(relation);
    EXPECT_TRUE(relation.value());
    auto visible = ctx.get_visible_relations_checked(ordinary);
    ASSERT_TRUE(visible);
    EXPECT_TRUE(std::any_of(visible.value().begin(), visible.value().end(),
        [&](const Relation& item) {
            return item.op == RelationOp::GT &&
                detail::node(item.lhs)->equals(*detail::node(*a)) &&
                detail::node(item.rhs)->equals(*detail::node(*c));
        }));
    ComputationContext inference;
    auto positive = ctx.is_positive_checked(*a, inference);
    ASSERT_TRUE(positive);
    EXPECT_EQ(positive.value(), Tribool::True);
}
