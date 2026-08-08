
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "query_interface.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>

using namespace lamina;


static SymbolicExpr make_var_expr(const std::string& name) {
    auto expr = lamina::detail::expression_from_node(lamina::detail::make_node<VariableNode>(name));
    return expr;
}

static SymbolicExpr make_num_expr(int v) {
    auto expr = lamina::detail::expression_from_node(lamina::detail::make_node<NumberNode>(BigInt(v)));
    return expr;
}

static SymbolicExpr make_relation_expr(const std::string& var_name,
                                       RelationalNode::Op op, int rhs_val) {
    auto lhs = lamina::detail::make_node<VariableNode>(var_name);
    auto rhs = lamina::detail::make_node<NumberNode>(BigInt(rhs_val));
    auto expr = lamina::detail::expression_from_node(lamina::detail::make_node<RelationalNode>(lhs, rhs, op));
    return expr;
}

static void EXPECT_TRIBOOL(Tribool actual, Tribool expected, const std::string& msg) {
    const char* names[] = {"True", "False", "Unknown"};
    int ai = static_cast<int>(actual);
    int ei = static_cast<int>(expected);
    if (ai == ei) {
        std::cout << "[PASS] " << msg << std::endl;
        g_passes++;
    } else {
        std::cerr << "[FAIL] " << msg
                  << "\n  Expected: " << names[ei]
                  << "\n  Got:      " << names[ai] << std::endl;
        g_failures++;
    }
}


void test_assume_domain_equivalence() {
    TEST_CASE("Property 32: assume_domain equivalence with PropertyStore.declare_domain");

    // Test all domain values
    std::vector<Domain> domains = {
        Domain::Complex, Domain::Real, Domain::Algebraic, Domain::Rational,
        Domain::Integer, Domain::Natural, Domain::PositiveInt
    };
    std::vector<std::string> domain_names = {
        "Complex", "Real", "Algebraic", "Rational", "Integer", "Natural", "PositiveInt"
    };

    std::vector<std::string> variables = {"x", "alpha", "var_123"};

    for (size_t di = 0; di < domains.size(); ++di) {
        for (const auto& var : variables) {
            // Method A: Use convenience API
            AssumptionContext ctx_a;
            ctx_a.assume_domain(var, domains[di]);

            // Method B: Use PropertyStore directly
            AssumptionContext ctx_b;
            ctx_b.current_properties().declare_domain(var, domains[di]);

            // Compare: get_domain should return the same result
            Domain result_a = ctx_a.get_domain(var);
            Domain result_b = ctx_b.get_domain(var);

            std::string label = "assume_domain(" + var + ", " + domain_names[di] + ")";
            EXPECT_TRUE(result_a == result_b,
                        label + " — get_domain matches direct PropertyStore");

            // Also verify has_domain for the declared domain
            bool has_a = ctx_a.has_domain(var, domains[di]);
            bool has_b = ctx_b.has_domain(var, domains[di]);
            EXPECT_TRUE(has_a == has_b,
                        label + " — has_domain matches direct PropertyStore");
        }
    }
}


void test_assume_sign_equivalence() {
    TEST_CASE("Property 32: assume_sign equivalence with PropertyStore.declare_sign");

    std::vector<Sign> signs = {
        Sign::Positive, Sign::Negative, Sign::NonNegative,
        Sign::NonPositive, Sign::Zero, Sign::NonZero
    };
    std::vector<std::string> sign_names = {
        "Positive", "Negative", "NonNegative",
        "NonPositive", "Zero", "NonZero"
    };

    std::vector<std::string> variables = {"x", "y", "beta"};

    for (size_t si = 0; si < signs.size(); ++si) {
        for (const auto& var : variables) {
            // Method A: Use convenience API
            AssumptionContext ctx_a;
            ctx_a.assume_sign(var, signs[si]);

            // Method B: Use PropertyStore directly
            AssumptionContext ctx_b;
            ctx_b.current_properties().declare_sign(var, signs[si]);

            // Compare: has_sign should return the same result for the declared sign
            bool has_a = ctx_a.has_sign(var, signs[si]);
            bool has_b = ctx_b.has_sign(var, signs[si]);

            std::string label = "assume_sign(" + var + ", " + sign_names[si] + ")";
            EXPECT_TRUE(has_a == has_b,
                        label + " — has_sign matches direct PropertyStore");

            // Also verify all implied signs match
            for (size_t si2 = 0; si2 < signs.size(); ++si2) {
                bool impl_a = ctx_a.has_sign(var, signs[si2]);
                bool impl_b = ctx_b.has_sign(var, signs[si2]);
                EXPECT_TRUE(impl_a == impl_b,
                            label + " — implied sign " + sign_names[si2] + " matches");
            }
        }
    }
}


void test_assume_relation_equivalence() {
    TEST_CASE("Property 32: assume(relation) equivalence with RelationStore.add_relation");

    struct RelTestCase {
        std::string var;
        RelationalNode::Op op;
        int rhs;
        std::string label;
        Sign expected_sign;
    };

    std::vector<RelTestCase> cases = {
        {"x", RelationalNode::Op::GT,  0, "x > 0",  Sign::Positive},
        {"y", RelationalNode::Op::GEQ, 0, "y >= 0", Sign::NonNegative},
        {"z", RelationalNode::Op::LT,  0, "z < 0",  Sign::Negative},
        {"w", RelationalNode::Op::LEQ, 0, "w <= 0", Sign::NonPositive},
        {"v", RelationalNode::Op::NEQ, 0, "v != 0", Sign::NonZero},
    };

    for (const auto& tc : cases) {
        // Method A: Use convenience API
        AssumptionContext ctx_a;
        auto rel_expr = make_relation_expr(tc.var, tc.op, tc.rhs);
        ctx_a.assume(rel_expr);

        // Method B: Use RelationStore directly
        AssumptionContext ctx_b;
        SymbolicExpr lhs_expr = make_var_expr(tc.var);
        SymbolicExpr rhs_expr = make_num_expr(tc.rhs);
        ctx_b.current_relations().add_relation(lhs_expr, rhs_expr, tc.op,
                                               ctx_b.current_properties());

        // Compare: the relation should be stored in both
        bool stored_a = ctx_a.current_relations().has_relation(lhs_expr, rhs_expr, tc.op);
        bool stored_b = ctx_b.current_relations().has_relation(lhs_expr, rhs_expr, tc.op);

        EXPECT_TRUE(stored_a == stored_b,
                    tc.label + " — relation stored equivalently");
        EXPECT_TRUE(stored_a, tc.label + " — relation is stored via convenience API");

        // Compare: sign property should be derived in both
        bool sign_a = ctx_a.has_sign(tc.var, tc.expected_sign);
        bool sign_b = ctx_b.has_sign(tc.var, tc.expected_sign);

        EXPECT_TRUE(sign_a == sign_b,
                    tc.label + " — derived sign property matches");
        EXPECT_TRUE(sign_a, tc.label + " — sign property derived via convenience API");
    }
}


void test_is_positive_equivalence() {
    TEST_CASE("Property 32: is_positive equivalence with QueryInterface.query_positive");

    // Set up context with some assumptions
    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Negative);
    ctx.assume_sign("z", Sign::NonNegative);

    // Test various expressions
    std::vector<std::pair<SymbolicExpr, std::string>> exprs = {
        {make_var_expr("x"), "x (Positive)"},
        {make_var_expr("y"), "y (Negative)"},
        {make_var_expr("z"), "z (NonNegative)"},
        {make_var_expr("undeclared"), "undeclared"},
        {make_num_expr(5), "5"},
        {make_num_expr(-3), "-3"},
        {make_num_expr(0), "0"},
    };

    QueryInterface qi(ctx);

    for (const auto& [expr, label] : exprs) {
        Tribool via_convenience = ctx.is_positive(expr);
        Tribool via_direct = qi.query_positive(expr);

        EXPECT_TRIBOOL(via_convenience, via_direct,
                       "is_positive(" + label + ") matches QueryInterface");
    }
}


void test_is_negative_equivalence() {
    TEST_CASE("Property 32: is_negative equivalence with QueryInterface.query_negative");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Negative);
    ctx.assume_sign("z", Sign::NonPositive);

    std::vector<std::pair<SymbolicExpr, std::string>> exprs = {
        {make_var_expr("x"), "x (Positive)"},
        {make_var_expr("y"), "y (Negative)"},
        {make_var_expr("z"), "z (NonPositive)"},
        {make_var_expr("undeclared"), "undeclared"},
        {make_num_expr(5), "5"},
        {make_num_expr(-3), "-3"},
        {make_num_expr(0), "0"},
    };

    QueryInterface qi(ctx);

    for (const auto& [expr, label] : exprs) {
        Tribool via_convenience = ctx.is_negative(expr);
        Tribool via_direct = qi.query_negative(expr);

        EXPECT_TRIBOOL(via_convenience, via_direct,
                       "is_negative(" + label + ") matches QueryInterface");
    }
}


void test_is_real_equivalence() {
    TEST_CASE("Property 32: is_real equivalence with QueryInterface.query_real");

    AssumptionContext ctx;
    ctx.assume_domain("x", Domain::Real);
    ctx.assume_domain("n", Domain::Integer);

    std::vector<std::pair<SymbolicExpr, std::string>> exprs = {
        {make_var_expr("x"), "x (Real)"},
        {make_var_expr("n"), "n (Integer)"},
        {make_var_expr("undeclared"), "undeclared"},
        {make_num_expr(42), "42"},
    };

    QueryInterface qi(ctx);

    for (const auto& [expr, label] : exprs) {
        Tribool via_convenience = ctx.is_real(expr);
        Tribool via_direct = qi.query_real(expr);

        EXPECT_TRIBOOL(via_convenience, via_direct,
                       "is_real(" + label + ") matches QueryInterface");
    }
}

void test_is_integer_equivalence() {
    TEST_CASE("Property 32: is_integer equivalence with QueryInterface.query_integer");

    AssumptionContext ctx;
    ctx.assume_domain("n", Domain::Integer);
    ctx.assume_domain("x", Domain::Real);

    std::vector<std::pair<SymbolicExpr, std::string>> exprs = {
        {make_var_expr("n"), "n (Integer)"},
        {make_var_expr("x"), "x (Real)"},
        {make_var_expr("undeclared"), "undeclared"},
        {make_num_expr(7), "7"},
    };

    QueryInterface qi(ctx);

    for (const auto& [expr, label] : exprs) {
        Tribool via_convenience = ctx.is_integer(expr);
        Tribool via_direct = qi.query_integer(expr);

        EXPECT_TRIBOOL(via_convenience, via_direct,
                       "is_integer(" + label + ") matches QueryInterface");
    }
}

void test_is_nonnegative_equivalence() {
    TEST_CASE("Property 32: is_nonnegative equivalence with QueryInterface.query_nonnegative");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Negative);
    ctx.assume_sign("z", Sign::NonNegative);

    std::vector<std::pair<SymbolicExpr, std::string>> exprs = {
        {make_var_expr("x"), "x (Positive)"},
        {make_var_expr("y"), "y (Negative)"},
        {make_var_expr("z"), "z (NonNegative)"},
        {make_var_expr("undeclared"), "undeclared"},
        {make_num_expr(0), "0"},
        {make_num_expr(-1), "-1"},
    };

    QueryInterface qi(ctx);

    for (const auto& [expr, label] : exprs) {
        Tribool via_convenience = ctx.is_nonnegative(expr);
        Tribool via_direct = qi.query_nonnegative(expr);

        EXPECT_TRIBOOL(via_convenience, via_direct,
                       "is_nonnegative(" + label + ") matches QueryInterface");
    }
}

void test_is_nonzero_equivalence() {
    TEST_CASE("Property 32: is_nonzero equivalence with QueryInterface.query_nonzero");

    AssumptionContext ctx;
    ctx.assume_sign("x", Sign::Positive);
    ctx.assume_sign("y", Sign::Zero);
    ctx.assume_sign("z", Sign::NonZero);

    std::vector<std::pair<SymbolicExpr, std::string>> exprs = {
        {make_var_expr("x"), "x (Positive)"},
        {make_var_expr("y"), "y (Zero)"},
        {make_var_expr("z"), "z (NonZero)"},
        {make_var_expr("undeclared"), "undeclared"},
        {make_num_expr(0), "0"},
        {make_num_expr(5), "5"},
    };

    QueryInterface qi(ctx);

    for (const auto& [expr, label] : exprs) {
        Tribool via_convenience = ctx.is_nonzero(expr);
        Tribool via_direct = qi.query_nonzero(expr);

        EXPECT_TRIBOOL(via_convenience, via_direct,
                       "is_nonzero(" + label + ") matches QueryInterface");
    }
}


void test_combined_assumptions_equivalence() {
    TEST_CASE("Property 32: Combined domain+sign assumptions — full equivalence");

    // Method A: Use convenience API for everything
    AssumptionContext ctx_a;
    ctx_a.assume_domain("x", Domain::Real);
    ctx_a.assume_sign("x", Sign::Positive);
    ctx_a.assume_domain("n", Domain::Integer);
    ctx_a.assume_sign("n", Sign::NonNegative);

    // Method B: Use PropertyStore directly
    AssumptionContext ctx_b;
    ctx_b.current_properties().declare_domain("x", Domain::Real);
    ctx_b.current_properties().declare_sign("x", Sign::Positive);
    ctx_b.current_properties().declare_domain("n", Domain::Integer);
    ctx_b.current_properties().declare_sign("n", Sign::NonNegative);

    // Compare domain queries
    EXPECT_TRUE(ctx_a.get_domain("x") == ctx_b.get_domain("x"),
                "Combined: get_domain(x) matches");
    EXPECT_TRUE(ctx_a.get_domain("n") == ctx_b.get_domain("n"),
                "Combined: get_domain(n) matches");

    // Compare sign queries
    EXPECT_TRUE(ctx_a.has_sign("x", Sign::Positive) == ctx_b.has_sign("x", Sign::Positive),
                "Combined: has_sign(x, Positive) matches");
    EXPECT_TRUE(ctx_a.has_sign("x", Sign::NonNegative) == ctx_b.has_sign("x", Sign::NonNegative),
                "Combined: has_sign(x, NonNegative) matches (implied)");
    EXPECT_TRUE(ctx_a.has_sign("x", Sign::NonZero) == ctx_b.has_sign("x", Sign::NonZero),
                "Combined: has_sign(x, NonZero) matches (implied)");
    EXPECT_TRUE(ctx_a.has_sign("n", Sign::NonNegative) == ctx_b.has_sign("n", Sign::NonNegative),
                "Combined: has_sign(n, NonNegative) matches");

    // Compare query results via QueryInterface
    QueryInterface qi_a(ctx_a);
    QueryInterface qi_b(ctx_b);

    auto x_expr = make_var_expr("x");
    auto n_expr = make_var_expr("n");

    EXPECT_TRIBOOL(qi_a.query_positive(x_expr), qi_b.query_positive(x_expr),
                   "Combined: query_positive(x) matches");
    EXPECT_TRIBOOL(qi_a.query_real(x_expr), qi_b.query_real(x_expr),
                   "Combined: query_real(x) matches");
    EXPECT_TRIBOOL(qi_a.query_integer(n_expr), qi_b.query_integer(n_expr),
                   "Combined: query_integer(n) matches");
    EXPECT_TRIBOOL(qi_a.query_nonnegative(n_expr), qi_b.query_nonnegative(n_expr),
                   "Combined: query_nonnegative(n) matches");
}


void test_scoped_convenience_equivalence() {
    TEST_CASE("Property 32: Scoped convenience API — push/pop equivalence");

    // Method A: Use convenience API within a pushed scope
    AssumptionContext ctx_a;
    ctx_a.assume_sign("x", Sign::Positive);
    ctx_a.push();
    ctx_a.assume_sign("x", Sign::Negative);

    // Method B: Use PropertyStore directly within a pushed scope
    AssumptionContext ctx_b;
    ctx_b.current_properties().declare_sign("x", Sign::Positive);
    ctx_b.push();
    ctx_b.current_properties().declare_sign("x", Sign::Negative);

    // In child scope, x should be Negative in both
    EXPECT_TRUE(ctx_a.has_sign("x", Sign::Negative) == ctx_b.has_sign("x", Sign::Negative),
                "Scoped: child scope has_sign(x, Negative) matches");

    // Pop both
    ctx_a.pop();
    ctx_b.pop();

    // After pop, x should be Positive in both
    EXPECT_TRUE(ctx_a.has_sign("x", Sign::Positive) == ctx_b.has_sign("x", Sign::Positive),
                "Scoped: after pop has_sign(x, Positive) matches");
    EXPECT_TRUE(ctx_a.has_sign("x", Sign::Negative) == ctx_b.has_sign("x", Sign::Negative),
                "Scoped: after pop has_sign(x, Negative) matches (should be false)");
}



void test_assume_domain_empty_name_throws() {
    TEST_CASE("Req 13.6: assume_domain with empty name throws std::invalid_argument");

    AssumptionContext ctx;
    bool threw = false;
    try {
        ctx.assume_domain("", Domain::Real);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "assume_domain(\"\", Domain::Real) should throw std::invalid_argument");
}

void test_assume_sign_empty_name_throws() {
    TEST_CASE("Req 13.6: assume_sign with empty name throws std::invalid_argument");

    AssumptionContext ctx;
    bool threw = false;
    try {
        ctx.assume_sign("", Sign::Positive);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "assume_sign(\"\", Sign::Positive) should throw std::invalid_argument");
}


void test_assume_non_relational_expr_throws() {
    TEST_CASE("Req 13.7: assume with AddNode root throws std::invalid_argument");

    AssumptionContext ctx;
    // Create an expression with AddNode root (x + y)
    auto x = lamina::detail::make_node<VariableNode>("x");
    auto y = lamina::detail::make_node<VariableNode>("y");
    auto add = lamina::detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{x, y});
    auto add_expr = lamina::detail::expression_from_node(add);
    bool threw = false;
    try {
        ctx.assume(add_expr);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "assume(expr_with_AddNode_root) should throw std::invalid_argument");
}


void test_is_positive_undeclared_returns_unknown() {
    TEST_CASE("Req 13.5: is_positive on undeclared variable returns Tribool::Unknown");

    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_var");
    Tribool result = ctx.is_positive(var_expr);
    EXPECT_TRUE(result == Tribool::Unknown,
                "is_positive(undeclared_var) should return Tribool::Unknown");
}

void test_is_negative_undeclared_returns_unknown() {
    TEST_CASE("Req 13.5: is_negative on undeclared variable returns Tribool::Unknown");

    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_var");
    Tribool result = ctx.is_negative(var_expr);
    EXPECT_TRUE(result == Tribool::Unknown,
                "is_negative(undeclared_var) should return Tribool::Unknown");
}

void test_is_nonnegative_undeclared_returns_unknown() {
    TEST_CASE("Req 13.5: is_nonnegative on undeclared variable returns Tribool::Unknown");

    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_z");
    Tribool result = ctx.is_nonnegative(var_expr);
    EXPECT_TRUE(result == Tribool::Unknown,
                "is_nonnegative(undeclared_z) should return Tribool::Unknown");
}

void test_is_real_undeclared_returns_unknown() {
    TEST_CASE("Req 13.5: is_real on undeclared variable returns Tribool::Unknown");

    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_w");
    Tribool result = ctx.is_real(var_expr);
    EXPECT_TRUE(result == Tribool::Unknown,
                "is_real(undeclared_w) should return Tribool::Unknown");
}

void test_is_integer_undeclared_returns_unknown() {
    TEST_CASE("Req 13.5: is_integer on undeclared variable returns Tribool::Unknown");

    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_alpha");
    Tribool result = ctx.is_integer(var_expr);
    EXPECT_TRUE(result == Tribool::Unknown,
                "is_integer(undeclared_alpha) should return Tribool::Unknown");
}

void test_is_nonzero_undeclared_returns_unknown() {
    TEST_CASE("Req 13.5: is_nonzero on undeclared variable returns Tribool::Unknown");

    AssumptionContext ctx;
    SymbolicExpr var_expr = make_var_expr("undeclared_beta");
    Tribool result = ctx.is_nonzero(var_expr);
    EXPECT_TRUE(result == Tribool::Unknown,
                "is_nonzero(undeclared_beta) should return Tribool::Unknown");
}

void test_checked_assumption_context_contracts() {
    TEST_CASE("AssumptionContext checked APIs: explicit errors and successful queries");

    AssumptionContext ctx;
    uint64_t generation = ctx.cache_generation();

    auto bad_domain = ctx.assume_domain_checked("", Domain::Real);
    EXPECT_TRUE(!bad_domain.has_value(), "checked assume_domain rejects empty variable");
    EXPECT_TRUE(bad_domain.error().code == CasErrc::InvalidArgument,
                "checked assume_domain reports InvalidArgument");
    EXPECT_TRUE(ctx.cache_generation() == generation,
                "failed checked assume_domain does not mutate generation");

    auto ok_domain = ctx.assume_domain_checked("x", Domain::Real);
    EXPECT_TRUE(ok_domain.has_value(), "checked assume_domain succeeds");
    EXPECT_TRUE(ctx.has_domain("x", Domain::Real),
                "checked assume_domain updates domain facts");

    generation = ctx.cache_generation();
    auto bad_sign = ctx.assume_sign_checked("", Sign::Positive);
    EXPECT_TRUE(!bad_sign.has_value(), "checked assume_sign rejects empty variable");
    EXPECT_TRUE(bad_sign.error().code == CasErrc::InvalidArgument,
                "checked assume_sign reports InvalidArgument");
    EXPECT_TRUE(ctx.cache_generation() == generation,
                "failed checked assume_sign does not mutate generation");

    auto ok_sign = ctx.assume_sign_checked("x", Sign::Positive);
    EXPECT_TRUE(ok_sign.has_value(), "checked assume_sign succeeds");

    auto natural = ctx.assume_domain_checked("n", Domain::Natural);
    EXPECT_TRUE(natural.has_value(), "checked assume_domain stores Natural domain");
    generation = ctx.cache_generation();
    auto sign_conflict = ctx.assume_sign_checked("n", Sign::Negative);
    EXPECT_TRUE(!sign_conflict.has_value(),
                "checked assume_sign rejects a domain/sign contradiction");
    EXPECT_TRUE(sign_conflict.error().code == CasErrc::InvalidArgument,
                "checked assume_sign reports InvalidArgument for contradiction");
    EXPECT_TRUE(ctx.cache_generation() == generation,
                "failed checked assume_sign does not advance cache generation");
    EXPECT_FALSE(ctx.has_sign("n", Sign::Negative),
                 "failed checked assume_sign does not commit a contradictory sign");

    SymbolicExpr x = make_var_expr("x");
    auto positive = ctx.is_positive_checked(x);
    EXPECT_TRUE(positive.has_value(), "checked is_positive succeeds");
    if (positive) {
        EXPECT_TRUE(positive.value() == Tribool::True,
                    "checked is_positive sees declared positive sign");
    }

    auto relation = make_relation_expr("x", RelationalNode::Op::GT, 0);
    auto relation_result = ctx.assume_checked(relation);
    EXPECT_TRUE(relation_result.has_value(), "checked assume accepts relational expression");

    generation = ctx.cache_generation();
    const auto relation_count = ctx.current_relations().get_relations().size();
    auto conflict_relation = make_relation_expr("x", RelationalNode::Op::LT, 0);
    auto conflict_result = ctx.assume_checked(conflict_relation);
    EXPECT_TRUE(!conflict_result.has_value(),
                "checked assume rejects relation with contradictory derived property");
    EXPECT_TRUE(conflict_result.error().code == CasErrc::InvalidArgument,
                "checked assume reports InvalidArgument for contradictory relation");
    EXPECT_TRUE(ctx.cache_generation() == generation,
                "failed checked assume does not mutate generation");
    EXPECT_TRUE(ctx.current_relations().get_relations().size() == relation_count,
                "failed checked assume does not store contradictory relation");
    EXPECT_FALSE(ctx.current_relations().has_relation(
                     make_var_expr("x"), make_num_expr(0), RelationalNode::Op::LT),
                 "failed checked assume leaves relation store unchanged");
    EXPECT_FALSE(ctx.has_sign("x", Sign::Negative),
                 "failed checked assume does not apply contradictory sign");

    auto bad_relation = ctx.assume_checked(x);
    EXPECT_TRUE(!bad_relation.has_value(), "checked assume rejects non-relational expression");
    EXPECT_TRUE(bad_relation.error().code == CasErrc::InvalidArgument,
                "checked assume reports InvalidArgument for non-relational expression");

    auto nonzero = ctx.is_nonzero_checked(x);
    EXPECT_TRUE(nonzero.has_value(), "checked is_nonzero succeeds");
    auto integer = ctx.is_integer_checked(x);
    EXPECT_TRUE(integer.has_value(), "checked is_integer succeeds");
    auto nonnegative = ctx.is_nonnegative_checked(x);
    EXPECT_TRUE(nonnegative.has_value(), "checked is_nonnegative succeeds");
    auto negative = ctx.is_negative_checked(x);
    EXPECT_TRUE(negative.has_value(), "checked is_negative succeeds");
}


int main() {
    // Test 1: assume_domain equivalence
    test_assume_domain_equivalence();

    // Test 2: assume_sign equivalence
    test_assume_sign_equivalence();

    // Test 3: assume(relation) equivalence
    test_assume_relation_equivalence();

    // Test 4: is_positive equivalence
    test_is_positive_equivalence();

    // Test 5: is_negative equivalence
    test_is_negative_equivalence();

    // Test 6: is_real, is_integer, is_nonnegative, is_nonzero equivalence
    test_is_real_equivalence();
    test_is_integer_equivalence();
    test_is_nonnegative_equivalence();
    test_is_nonzero_equivalence();

    // Combined scenario
    test_combined_assumptions_equivalence();

    // Scoped equivalence
    test_scoped_convenience_equivalence();


    test_assume_domain_empty_name_throws();
    test_assume_sign_empty_name_throws();

    test_assume_non_relational_expr_throws();

    test_is_positive_undeclared_returns_unknown();
    test_is_negative_undeclared_returns_unknown();
    test_is_nonnegative_undeclared_returns_unknown();
    test_is_real_undeclared_returns_unknown();
    test_is_integer_undeclared_returns_unknown();
    test_is_nonzero_undeclared_returns_unknown();
    test_checked_assumption_context_contracts();

    return TEST_REPORT();
}
