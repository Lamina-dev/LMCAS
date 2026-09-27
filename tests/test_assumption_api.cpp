
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "query_interface.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>

using namespace LMCAS;

static SymbolicExpr make_var_expr(const std::string &name) {
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
    return expr;
}

static SymbolicExpr make_num_expr(int v) {
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<NumberNode>(BigInt(v)));
    return expr;
}

static SymbolicExpr make_relation_expr(const std::string &var_name,
                                       RelationalNode::Op op, int rhs_val) {
    auto lhs = LMCAS::detail::make_node<VariableNode>(var_name);
    auto rhs = LMCAS::detail::make_node<NumberNode>(BigInt(rhs_val));
    auto expr = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(lhs, rhs, op));
    return expr;
}

TEST(AssumptionApi, AssumeRelationDerivesSign) {
    struct RelTestCase {
        std::string var;
        RelationalNode::Op op;
        int rhs;
        std::string label;
        Sign expected_sign;
    };

    std::vector<RelTestCase> cases = {
        {"x", RelationalNode::Op::GT, 0, "x > 0", Sign::Positive},
        {"y", RelationalNode::Op::GEQ, 0, "y >= 0", Sign::NonNegative},
        {"z", RelationalNode::Op::LT, 0, "z < 0", Sign::Negative},
        {"w", RelationalNode::Op::LEQ, 0, "w <= 0", Sign::NonPositive},
        {"v", RelationalNode::Op::NEQ, 0, "v != 0", Sign::NonZero},
    };

    for (const auto &tc : cases) {
        // Method A: Use convenience API
        AssumptionContext ctx_a;
        auto rel_expr = make_relation_expr(tc.var, tc.op, tc.rhs);
        EXPECT_TRUE((ctx_a.assume(rel_expr).has_value())) << "relation assumption succeeds";

        SymbolicExpr lhs_expr = make_var_expr(tc.var);
        SymbolicExpr rhs_expr = make_num_expr(tc.rhs);

        bool stored_a = ctx_a.current_relations().has_relation(lhs_expr, rhs_expr, tc.op);
        EXPECT_TRUE((stored_a)) << tc.label + " — relation is stored via convenience API";

        bool sign_a = ctx_a.has_sign(tc.var, tc.expected_sign);
        EXPECT_TRUE((sign_a)) << tc.label + " — sign property derived via convenience API";
    }
}

TEST(AssumptionApi, AssumeDomainEmptyNameReturnsError) {
    AssumptionContext ctx;
    auto result = ctx.assume_domain("", Domain::Real);
    EXPECT_TRUE((!result.has_value()));
}

TEST(AssumptionApi, AssumeSignEmptyNameReturnsError) {
    AssumptionContext ctx;
    auto result = ctx.assume_sign("", Sign::Positive);
    EXPECT_TRUE((!result.has_value()));
}

TEST(AssumptionApi, AssumeNonRelationalExprReturnsError) {
    AssumptionContext ctx;
    auto x = LMCAS::detail::make_node<VariableNode>("x");
    auto y = LMCAS::detail::make_node<VariableNode>("y");
    auto add = LMCAS::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x, y});
    auto add_expr = LMCAS::detail::expression_from_node(add);
    auto result = ctx.assume(add_expr);
    EXPECT_TRUE((!result.has_value()));
}

TEST(AssumptionApi, IsPositiveUndeclaredReturnsUnknown) {
    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_var");
    Tribool result = ctx.is_positive(var_expr).value();
    EXPECT_TRUE((result == Tribool::Unknown)) << "is_positive(undeclared_var) should return Tribool::Unknown";
}

TEST(AssumptionApi, IsNegativeUndeclaredReturnsUnknown) {
    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_var");
    Tribool result = ctx.is_negative(var_expr).value();
    EXPECT_TRUE((result == Tribool::Unknown)) << "is_negative(undeclared_var) should return Tribool::Unknown";
}

TEST(AssumptionApi, IsNonnegativeUndeclaredReturnsUnknown) {
    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_z");
    Tribool result = ctx.is_nonnegative(var_expr).value();
    EXPECT_TRUE((result == Tribool::Unknown)) << "is_nonnegative(undeclared_z) should return Tribool::Unknown";
}

TEST(AssumptionApi, IsRealUndeclaredReturnsUnknown) {
    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_w");
    Tribool result = ctx.is_real(var_expr).value();
    EXPECT_TRUE((result == Tribool::Unknown)) << "is_real(undeclared_w) should return Tribool::Unknown";
}

TEST(AssumptionApi, IsIntegerUndeclaredReturnsUnknown) {
    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_alpha");
    Tribool result = ctx.is_integer(var_expr).value();
    EXPECT_TRUE((result == Tribool::Unknown)) << "is_integer(undeclared_alpha) should return Tribool::Unknown";
}

TEST(AssumptionApi, IsNonzeroUndeclaredReturnsUnknown) {
    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_beta");
    Tribool result = ctx.is_nonzero(var_expr).value();
    EXPECT_TRUE((result == Tribool::Unknown)) << "is_nonzero(undeclared_beta) should return Tribool::Unknown";
}

static void test_checked_assumption_context_contracts_domains_and_signs(AssumptionContext &ctx, uint64_t &generation) {
    auto bad_domain = ctx.assume_domain_checked("", Domain::Real);
    EXPECT_TRUE((!bad_domain.has_value())) << "checked assume_domain rejects empty variable";
    EXPECT_TRUE((bad_domain.error().code == CasErrc::InvalidArgument)) << "checked assume_domain reports InvalidArgument";
    EXPECT_TRUE((ctx.cache_generation() == generation)) << "failed checked assume_domain does not mutate generation";

    auto ok_domain = ctx.assume_domain_checked("x", Domain::Real);
    EXPECT_TRUE((ok_domain.has_value())) << "checked assume_domain succeeds";
    EXPECT_TRUE((ctx.has_domain("x", Domain::Real))) << "checked assume_domain updates domain facts";

    generation = ctx.cache_generation();
    auto bad_sign = ctx.assume_sign_checked("", Sign::Positive);
    EXPECT_TRUE((!bad_sign.has_value())) << "checked assume_sign rejects empty variable";
    EXPECT_TRUE((bad_sign.error().code == CasErrc::InvalidArgument)) << "checked assume_sign reports InvalidArgument";
    EXPECT_TRUE((ctx.cache_generation() == generation)) << "failed checked assume_sign does not mutate generation";

    auto ok_sign = ctx.assume_sign_checked("x", Sign::Positive);
    EXPECT_TRUE((ok_sign.has_value())) << "checked assume_sign succeeds";

    auto natural = ctx.assume_domain_checked("n", Domain::Natural);
    EXPECT_TRUE((natural.has_value())) << "checked assume_domain stores Natural domain";
    generation = ctx.cache_generation();
    auto sign_conflict = ctx.assume_sign_checked("n", Sign::Negative);
    EXPECT_TRUE((!sign_conflict.has_value())) << "checked assume_sign rejects a domain/sign contradiction";
    EXPECT_TRUE((sign_conflict.error().code == CasErrc::InvalidArgument)) << "checked assume_sign reports InvalidArgument for contradiction";
    EXPECT_TRUE((ctx.cache_generation() == generation)) << "failed checked assume_sign does not advance cache generation";
    EXPECT_FALSE((ctx.has_sign("n", Sign::Negative))) << "failed checked assume_sign does not commit a contradictory sign";
}

static void test_checked_assumption_context_contracts_queries_and_relations(AssumptionContext &ctx, uint64_t &generation) {
    SymbolicExpr x = make_var_expr("x");
    auto positive = ctx.is_positive_checked(x);
    ASSERT_TRUE((positive.has_value())) << "checked is_positive succeeds";
    if (positive) {
        EXPECT_TRUE((positive.value() == Tribool::True)) << "checked is_positive sees declared positive sign";
    }

    auto relation = make_relation_expr("x", RelationalNode::Op::GT, 0);
    auto relation_result = ctx.assume_checked(relation);
    EXPECT_TRUE((relation_result.has_value())) << "checked assume accepts relational expression";

    generation = ctx.cache_generation();
    const auto relation_count = ctx.current_relations().get_relations().size();
    auto conflict_relation = make_relation_expr("x", RelationalNode::Op::LT, 0);
    auto conflict_result = ctx.assume_checked(conflict_relation);
    EXPECT_TRUE((!conflict_result.has_value())) << "checked assume rejects relation with contradictory derived property";
    EXPECT_TRUE((conflict_result.error().code == CasErrc::InvalidArgument)) << "checked assume reports InvalidArgument for contradictory relation";
    EXPECT_TRUE((ctx.cache_generation() == generation)) << "failed checked assume does not mutate generation";
    EXPECT_TRUE((ctx.current_relations().get_relations().size() == relation_count)) << "failed checked assume does not store contradictory relation";
    EXPECT_FALSE((ctx.current_relations().has_relation(
        make_var_expr("x"), make_num_expr(0), RelationalNode::Op::LT)))
        << "failed checked assume leaves relation store unchanged";
    EXPECT_FALSE((ctx.has_sign("x", Sign::Negative))) << "failed checked assume does not apply contradictory sign";

    auto bad_relation = ctx.assume_checked(x);
    EXPECT_TRUE((!bad_relation.has_value())) << "checked assume rejects non-relational expression";
    EXPECT_TRUE((bad_relation.error().code == CasErrc::InvalidArgument)) << "checked assume reports InvalidArgument for non-relational expression";

    auto nonzero = ctx.is_nonzero_checked(x);
    EXPECT_TRUE((nonzero.has_value())) << "checked is_nonzero succeeds";
    auto integer = ctx.is_integer_checked(x);
    EXPECT_TRUE((integer.has_value())) << "checked is_integer succeeds";
    auto nonnegative = ctx.is_nonnegative_checked(x);
    EXPECT_TRUE((nonnegative.has_value())) << "checked is_nonnegative succeeds";
    auto negative = ctx.is_negative_checked(x);
    EXPECT_TRUE((negative.has_value())) << "checked is_negative succeeds";
}

TEST(AssumptionApi, CheckedAssumptionContextContracts) {
    AssumptionContext ctx;
    uint64_t generation = ctx.cache_generation();

    test_checked_assumption_context_contracts_domains_and_signs(ctx, generation);
    test_checked_assumption_context_contracts_queries_and_relations(ctx, generation);
}
