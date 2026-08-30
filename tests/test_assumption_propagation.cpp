
#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "property_store.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace lamina;


static std::shared_ptr<const SymbolicNode> make_var(const std::string& name) {
    return lamina::detail::make_node<VariableNode>(name);
}

static std::shared_ptr<const SymbolicNode> make_number(int val) {
    return lamina::detail::make_node<NumberNode>(BigInt(val));
}

static std::shared_ptr<const SymbolicNode> make_power(
    std::shared_ptr<const SymbolicNode> base,
    std::shared_ptr<const SymbolicNode> exp) {
    return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exp));
}

static std::shared_ptr<const SymbolicNode> make_abs(std::shared_ptr<const SymbolicNode> arg) {
    return lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{std::move(arg)});
}

static SymbolicExpr wrap_expr(std::shared_ptr<const SymbolicNode> node) {
    auto expr = lamina::detail::expression_from_node(std::move(node));
    return expr;
}

/// Build x^2 expression: PowerNode(VariableNode(name), NumberNode(2))
static SymbolicExpr make_x_squared(const std::string& var_name) {
    auto var = make_var(var_name);
    auto two = make_number(2);
    auto pow = make_power(var, two);
    return wrap_expr(pow);
}

/// Build |x| expression: FunctionNode(Abs, [VariableNode(name)])
static SymbolicExpr make_abs_x(const std::string& var_name) {
    auto var = make_var(var_name);
    auto abs_node = make_abs(var);
    return wrap_expr(abs_node);
}

/// Generate a random variable name for property tests
static std::string random_var_name() {
    static const std::vector<std::string> prefixes = {"x", "y", "z", "a", "b", "t", "u", "v", "w"};
    std::string prefix = rc::gen::elementOf(prefixes);
    return prefix + "_" + std::to_string(rc::gen::inRange(0, 999));
}

/// Generate a random Real-or-more-specific domain (Real, Algebraic, Rational, Integer, Natural, PositiveInt)
static Domain random_real_domain() {
    static const std::vector<Domain> real_domains = {
        Domain::Real, Domain::Algebraic, Domain::Rational,
        Domain::Integer, Domain::Natural, Domain::PositiveInt
    };
    return rc::gen::elementOf(real_domains);
}

/// Generate a random Integer-or-more-specific domain (Integer, Natural, PositiveInt)
static Domain random_integer_domain() {
    static const std::vector<Domain> int_domains = {
        Domain::Integer, Domain::Natural, Domain::PositiveInt
    };
    return rc::gen::elementOf(int_domains);
}


static void test_real_var_x_squared_nonnegative() {
    TEST_CASE("Real variable x² is non-negative");

    rc::check("For any variable declared Real (or more specific), "
              "query_nonnegative on x² returns True", []() {
        std::string var_name = random_var_name();
        Domain domain = random_real_domain();

        AssumptionContext ctx;
        ctx.assume_domain(var_name, domain);
        InferenceEngine engine(ctx);

        auto x_squared = make_x_squared(var_name);

        // x^2 should be non-negative for any Real variable
        RC_ASSERT(engine.query_nonnegative_checked(x_squared).value() == Tribool::True);
    });
}


static void test_integer_var_x_squared_integer() {
    TEST_CASE("Integer variable x² has Integer domain");

    rc::check("For any variable declared Integer (or more specific), "
              "query_integer on x² returns True", []() {
        std::string var_name = random_var_name();
        Domain domain = random_integer_domain();

        AssumptionContext ctx;
        ctx.assume_domain(var_name, domain);
        InferenceEngine engine(ctx);

        auto x_squared = make_x_squared(var_name);

        // x^2 should have Integer domain when x is Integer
        RC_ASSERT(engine.query_integer_checked(x_squared).value() == Tribool::True);
    });
}


static void test_positive_var_abs_positive() {
    TEST_CASE("Positive variable |x| is positive");

    rc::check("For any variable declared Positive and Real, "
              "query_positive on |x| returns True", []() {
        std::string var_name = random_var_name();
        Domain domain = random_real_domain();

        AssumptionContext ctx;
        ctx.assume_domain(var_name, domain);
        ctx.assume_sign(var_name, Sign::Positive);
        InferenceEngine engine(ctx);

        auto abs_x = make_abs_x(var_name);

        // |x| should be positive when x is positive and real
        RC_ASSERT(engine.query_positive_checked(abs_x).value() == Tribool::True);
        RC_ASSERT(engine.query_nonnegative_checked(abs_x).value() == Tribool::True);
    });
}


static void test_propagation_is_lazy() {
    TEST_CASE("Propagation is lazy (query-time only)");

    rc::check("After declaring a variable Real, the PropertyStore does not "
              "contain any entry for x² — propagation happens at query time only", []() {
        std::string var_name = random_var_name();
        Domain domain = random_real_domain();

        AssumptionContext ctx;
        ctx.assume_domain(var_name, domain);

        /// PropertyStore 仅保存声明变量;x^2 等复合表达式按需推导.
        const PropertyStore& store = ctx.current_properties();

        // The store should have the variable we declared
        RC_ASSERT(store.get_domain(var_name) == domain);

        // A synthetic name like "var_name^2" should NOT be in the store
        // (since propagation is lazy, no derived expressions are stored)
        std::string squared_name = var_name + "^2";
        RC_ASSERT(store.get_domain(squared_name) == Domain::Complex); // default = not stored

        // But querying through the InferenceEngine should still work
        InferenceEngine engine(ctx);
        auto x_squared = make_x_squared(var_name);
        RC_ASSERT(engine.query_nonnegative_checked(x_squared).value() == Tribool::True);
    });
}


int main() {
    test_real_var_x_squared_nonnegative();
    test_integer_var_x_squared_integer();
    test_positive_var_abs_positive();
    test_propagation_is_lazy();

    return TEST_REPORT();
}
