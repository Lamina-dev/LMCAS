
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "expr.hpp"
#include "numeric_evaluation.hpp"
#include "query_interface.hpp"
#include "residual_verification.hpp"
#include <utility>
#include <memory>
#include <string>
#include <cmath>
#include <optional>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace LMCAS;

/// Create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var(const std::string &name) {
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
}

/// Create a FunctionNode expression (e.g., sin(x), cos(x), tan(x)).
static SymbolicExpr make_func(FunctionNode::FuncType type, const std::string &var_name) {
    auto var_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto func_node = LMCAS::detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{var_node});
    return LMCAS::detail::expression_from_node(func_node);
}

static SymbolicExpr make_number_expr(double value) {
    return LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<NumberNode>(
            static_cast<lmmc_real_t>(value)));
}

static std::optional<double> extract_numeric(const SymbolicExpr &expr) {
    auto value = evaluate_numeric(expr);
    if (!value)
        return std::nullopt;
    return value.value().value;
}

static void expect_exact_period_value(
    const std::optional<SymbolicExpr> &period, const std::string &source) {
    auto expected = parse_expr(source);
    EXPECT_TRUE((period.has_value() && expected.has_value())) << "exact period and reference exist";
    if (!period || !expected)
        return;
    ComputationContext context;
    auto equal = check_equivalent(detail::make_expression_ptr(*period), expected.value(), context);
    EXPECT_TRUE((equal && std::holds_alternative<ProvedZeroResidual>(equal.value()))) << "period agrees with an exact symbolic reference";
}

TEST(AssumptionPropertiesExt, PeriodicityDeclaredSymbolRoundtrip) {
    AssumptionContext ctx;
    auto period_expr = make_number_expr(5.0);
    EXPECT_TRUE((ctx.current_properties().declare_periodic("f", "x", period_expr).has_value())) << "period declaration succeeds";

    InferenceEngine engine(ctx);

    SymbolicExpr f_expr = make_var("f");

    // query_periodic should return True
    Tribool is_periodic = engine.query_periodic_checked(f_expr, "x").value();
    EXPECT_TRUE((is_periodic == Tribool::True)) << "Declared periodic symbol: query_periodic returns True";

    // infer_period should return the declared period
    auto inferred = engine.infer_period_checked(f_expr, "x").value();
    ASSERT_TRUE((inferred.has_value())) << "Declared periodic symbol: infer_period returns a value";

    if (inferred.has_value()) {
        auto val = extract_numeric(*inferred);
        EXPECT_TRUE((val.has_value() && std::abs(*val - 5.0) < 1e-10)) << "Declared periodic symbol: infer_period returns period = 5.0";
    }
}

TEST(AssumptionPropertiesExt, PeriodicityDeclaredSymbolVariousPeriods) {
    AssumptionContext ctx;

    // Declare several symbols with different periods
    EXPECT_TRUE((ctx.current_properties().declare_periodic("g", "x", make_number_expr(2.0)).has_value())) << "period declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_periodic("h", "x", make_number_expr(M_PI)).has_value())) << "period declaration succeeds";
    EXPECT_TRUE((ctx.current_properties().declare_periodic("k", "x", make_number_expr(100.0)).has_value())) << "period declaration succeeds";

    InferenceEngine engine(ctx);

    // g: period 2.0
    {
        SymbolicExpr g_expr = make_var("g");
        EXPECT_TRUE((engine.query_periodic_checked(g_expr, "x").value() == Tribool::True)) << "g is periodic";
        auto period = engine.infer_period_checked(g_expr, "x").value();
        ASSERT_TRUE((period.has_value())) << "g has inferred period";
        if (period.has_value()) {
            auto val = extract_numeric(*period);
            EXPECT_TRUE((val.has_value() && std::abs(*val - 2.0) < 1e-10)) << "g period = 2.0";
        }
    }

    // h: period pi
    {
        SymbolicExpr h_expr = make_var("h");
        EXPECT_TRUE((engine.query_periodic_checked(h_expr, "x").value() == Tribool::True)) << "h is periodic";
        auto period = engine.infer_period_checked(h_expr, "x").value();
        ASSERT_TRUE((period.has_value())) << "h has inferred period";
        if (period.has_value()) {
            auto val = extract_numeric(*period);
            EXPECT_TRUE((val.has_value() && std::abs(*val - M_PI) < 1e-10)) << "h period = pi";
        }
    }

    // k: period 100.0
    {
        SymbolicExpr k_expr = make_var("k");
        EXPECT_TRUE((engine.query_periodic_checked(k_expr, "x").value() == Tribool::True)) << "k is periodic";
        auto period = engine.infer_period_checked(k_expr, "x").value();
        ASSERT_TRUE((period.has_value())) << "k has inferred period";
        if (period.has_value()) {
            auto val = extract_numeric(*period);
            EXPECT_TRUE((val.has_value() && std::abs(*val - 100.0) < 1e-10)) << "k period = 100.0";
        }
    }
}

TEST(AssumptionPropertiesExt, PeriodicitySinAutoInferred) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr sin_x = make_func(FunctionNode::FuncType::Sin, "x");

    // sin should be periodic
    Tribool is_periodic = engine.query_periodic_checked(sin_x, "x").value();
    EXPECT_TRUE((is_periodic == Tribool::True)) << "sin(x) is periodic";

    // Period should be 2*pi
    auto period = engine.infer_period_checked(sin_x, "x").value();
    ASSERT_TRUE((period.has_value())) << "sin(x) has inferred period";

    if (period.has_value()) {
        expect_exact_period_value(period, "2*pi");
    }
}

TEST(AssumptionPropertiesExt, PeriodicityCosAutoInferred) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr cos_x = make_func(FunctionNode::FuncType::Cos, "x");

    Tribool is_periodic = engine.query_periodic_checked(cos_x, "x").value();
    EXPECT_TRUE((is_periodic == Tribool::True)) << "cos(x) is periodic";

    auto period = engine.infer_period_checked(cos_x, "x").value();
    ASSERT_TRUE((period.has_value())) << "cos(x) has inferred period";

    if (period.has_value()) {
        expect_exact_period_value(period, "2*pi");
    }
}

TEST(AssumptionPropertiesExt, PeriodicityTanAutoInferred) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr tan_x = make_func(FunctionNode::FuncType::Tan, "x");

    Tribool is_periodic = engine.query_periodic_checked(tan_x, "x").value();
    EXPECT_TRUE((is_periodic == Tribool::True)) << "tan(x) is periodic";

    auto period = engine.infer_period_checked(tan_x, "x").value();
    ASSERT_TRUE((period.has_value())) << "tan(x) has inferred period";

    if (period.has_value()) {
        expect_exact_period_value(period, "pi");
    }
}

TEST(AssumptionPropertiesExt, PeriodicityNonPeriodicFunction) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr exp_x = make_func(FunctionNode::FuncType::Exp, "x");

    Tribool is_periodic = engine.query_periodic_checked(exp_x, "x").value();
    EXPECT_TRUE((is_periodic == Tribool::Unknown)) << "exp(x) is not known to be periodic (returns Unknown)";

    auto period = engine.infer_period_checked(exp_x, "x").value();
    EXPECT_FALSE((period.has_value())) << "exp(x) has no inferred period";
}

TEST(AssumptionPropertiesExt, PeriodicityNonPeriodicVariable) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var("x");

    Tribool is_periodic = engine.query_periodic_checked(x_expr, "x").value();
    EXPECT_TRUE((is_periodic == Tribool::Unknown)) << "Undeclared variable: query_periodic returns Unknown";

    auto period = engine.infer_period_checked(x_expr, "x").value();
    EXPECT_FALSE((period.has_value())) << "Undeclared variable: infer_period returns nullopt";
}

TEST(AssumptionPropertiesExt, PeriodicityPropertyStoreRoundtrip) {
    PropertyStore store;

    // Not periodic initially
    EXPECT_FALSE((store.is_periodic("f", "x"))) << "f not periodic initially";
    EXPECT_FALSE((store.get_period("f", "x").has_value())) << "f has no period initially";

    // Declare periodic
    auto period = make_number_expr(7.0);
    EXPECT_TRUE((store.declare_periodic("f", "x", period).has_value())) << "period declaration succeeds";

    // Now periodic
    EXPECT_TRUE((store.is_periodic("f", "x"))) << "f is periodic after declaration";
    EXPECT_TRUE((store.get_period("f", "x").has_value())) << "f has period after declaration";

    // Period value matches
    auto retrieved = store.get_period("f", "x");
    if (retrieved.has_value()) {
        auto val = extract_numeric(*retrieved);
        EXPECT_TRUE((val.has_value() && std::abs(*val - 7.0) < 1e-10)) << "Retrieved period = 7.0";
    }
}

TEST(AssumptionPropertiesExt, PeriodicityLnNotPeriodic) {
    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr ln_x = make_func(FunctionNode::FuncType::Ln, "x");

    Tribool is_periodic = engine.query_periodic_checked(ln_x, "x").value();
    EXPECT_TRUE((is_periodic == Tribool::Unknown)) << "ln(x) is not known to be periodic";

    auto period = engine.infer_period_checked(ln_x, "x").value();
    EXPECT_FALSE((period.has_value())) << "ln(x) has no inferred period";
}

static void expect_affine_period_inference(InferenceEngine &engine) {
    const std::pair<const char *, const char *> cases[] = {
        {"sin(x/2)", "4*pi"}, {"cos(x/2)", "4*pi"}, {"sin(2*x)", "pi"}, {"tan(-3*x+1)", "pi/3"}, {"cos(x/2+y)", "4*pi"}, {"atan(sin(x/2))", "4*pi"}};
    for (const auto &[source, expected] : cases) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << "periodic input parses";
        if (!expression) {
            continue;
        }
        auto periodic = engine.query_periodic_checked(*expression.value(), "x");
        EXPECT_TRUE((periodic && periodic.value() == Tribool::True)) << "affine period is proved";
        auto period = engine.infer_period_checked(*expression.value(), "x");
        ASSERT_TRUE((period.has_value())) << "period inference succeeds";
        if (period) {
            expect_exact_period_value(period.value(), expected);
        }
    }
    for (const char *source : {"sin(x*x)", "sin(sin(x))", "sin(x*y)", "sin(x^0)"}) {
        auto expression = parse_expr(source);
        ASSERT_TRUE((expression.has_value())) << "unproved input parses";
        if (!expression) {
            continue;
        }
        auto periodic = engine.query_periodic_checked(*expression.value(), "x");
        auto period = engine.infer_period_checked(*expression.value(), "x");
        EXPECT_TRUE((periodic && periodic.value() == Tribool::Unknown && period && !period.value())) << "nonlinear, degenerate-unknown, or punctured-domain input remains unproved";
    }
}

static void expect_variable_keyed_period_cache(
    AssumptionContext &ctx, QueryInterface &cached, const ExprPtr &scaled) {
    {
        auto along_x = cached.query_periodic(*scaled, "x");
        auto along_y = cached.query_periodic(*scaled, "y");
        EXPECT_TRUE((along_x && along_x.value() == Tribool::True &&
                     along_y && along_y.value() == Tribool::Unknown))
            << "one cache separates the independent variable from an unproved real parameter";
    }
    EXPECT_TRUE(ctx.assume_domain("x", Domain::Real).has_value()) << "parameter realness is declared";
    {
        auto along_y = cached.query_periodic(*scaled, "y");
        auto constant_period = cached.get_period(*scaled, "y");
        EXPECT_TRUE((along_y && along_y.value() == Tribool::True &&
                     constant_period && !constant_period.value()))
            << "a real constant has no minimum positive period";
    }
    auto monotone = parse_expr("atan(x)").value();
    auto nonperiodic = cached.query_periodic(*monotone, "x");
    auto independent = cached.query_periodic(*monotone, "y");
    EXPECT_TRUE((nonperiodic && nonperiodic.value() == Tribool::False &&
                 independent && independent.value() == Tribool::True))
        << "cached strict-monotonic and constant queries cannot alias";
}

TEST(AssumptionPropertiesExt, AffinePeriodsAndVariableCache) {
    AssumptionContext ctx;
    EXPECT_TRUE(ctx.assume_domain("y", Domain::Real).has_value()) << "parameter is real";
    InferenceEngine engine(ctx);
    QueryInterface cached(ctx);
    expect_affine_period_inference(engine);
    auto scaled = parse_expr("sin(x/2)").value();
    expect_variable_keyed_period_cache(ctx, cached, scaled);
    for (const char *source : {"0", "sin(0)", "exp(y)"}) {
        auto expression = parse_expr(source).value();
        auto periodic = engine.query_periodic_checked(*expression, "x");
        auto period = engine.infer_period_checked(*expression, "x");
        EXPECT_TRUE((periodic && periodic.value() == Tribool::True && period && !period.value())) << "proved constants are periodic without a zero or minimum period";
    }
    EXPECT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value()) << "nonzero slope is proved";
    auto symbolic = parse_expr("sin(x*y)").value();
    auto symbolic_period = engine.infer_period_checked(*symbolic, "x");
    ASSERT_TRUE((symbolic_period.has_value())) << "symbolic affine query succeeds";
    if (symbolic_period) {
        expect_exact_period_value(symbolic_period.value(), "2*pi/y");
    }
    auto invalid = cached.query_periodic(*scaled, "");
    EXPECT_TRUE((!invalid && invalid.error().code == CasErrc::InvalidArgument)) << "empty independent-variable names are rejected";
    const int proof_depth = engine.get_max_depth();
    engine.set_max_depth(1);
    auto deep = parse_expr("atan(sin(x))").value();
    auto unknown = engine.query_periodic_checked(*deep, "x");
    ASSERT_TRUE(unknown);
    EXPECT_EQ(unknown.value(), Tribool::Unknown);
    auto unproved_period = engine.infer_period_checked(*deep, "x");
    ASSERT_TRUE(unproved_period);
    EXPECT_FALSE(unproved_period.value());

    engine.set_max_depth(proof_depth);
    ResourceLimits limits;
    limits.max_recursion_depth = 1;
    ComputationContext context(limits);
    auto exhausted = engine.infer_period_checked(*deep, "x", context);
    ASSERT_FALSE(exhausted);
    EXPECT_EQ(exhausted.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(context.recursion_depth(), 0u);
}

TEST(AssumptionPropertiesExt, PeriodValuesOwnTheirWrappersAndRefreshCaches) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("f", Domain::Real));
    SymbolicExpr caller_period = *SymbolicExpr::number(5);
    ASSERT_TRUE(ctx.current_properties().declare_periodic(
        "f", "x", caller_period));
    const SymbolicExpr y_period = *SymbolicExpr::number(17);
    ASSERT_TRUE(ctx.current_properties().declare_periodic(
        "f", "y", y_period));

    QueryInterface cached(ctx);
    auto symbol = SymbolicExpr::variable("f");
    const auto expect_period = [&](const char *variable,
                                   const char *expected) {
        auto result = cached.get_period(*symbol, variable);
        ASSERT_TRUE(result);
        expect_exact_period_value(result.value(), expected);
    };
    expect_period("x", "5");
    expect_period("y", "17");

    caller_period = *SymbolicExpr::number(7);
    expect_period("x", "5");
    auto retrieved = ctx.get_period("f", "x");
    ASSERT_TRUE(retrieved);
    *retrieved = *SymbolicExpr::number(9);
    expect_period("x", "5");

    const SymbolicExpr updated = *SymbolicExpr::number(11);
    ASSERT_TRUE(ctx.current_properties().declare_periodic(
        "f", "x", updated));
    expect_period("x", "11");

    ctx.push();
    const SymbolicExpr child = *SymbolicExpr::number(13);
    ASSERT_TRUE(ctx.current_properties().declare_periodic(
        "f", "x", child));
    expect_period("x", "13");
    expect_period("y", "17");
    ASSERT_TRUE(ctx.pop());
    expect_period("x", "11");
    expect_period("y", "17");

    const SymbolicExpr symbolic_period = *SymbolicExpr::variable("a");
    ASSERT_TRUE(ctx.current_properties().declare_periodic(
        "g", "x", symbolic_period));
    auto g = SymbolicExpr::variable("g");
    auto unknown = cached.query_periodic(*g, "x");
    ASSERT_TRUE(unknown);
    EXPECT_EQ(unknown.value(), Tribool::Unknown);
    ASSERT_TRUE(ctx.assume_sign("a", Sign::Positive));
    auto proved = cached.query_periodic(*g, "x");
    ASSERT_TRUE(proved);
    EXPECT_EQ(proved.value(), Tribool::True);
}
