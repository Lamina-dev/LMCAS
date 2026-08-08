
#include "test_common.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "assumption.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <memory>
#include <string>
#include <cmath>
#include <optional>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace lamina;


/// Create a SymbolicExpr wrapping a VariableNode.
static SymbolicExpr make_var(const std::string& name) {
    return lamina::detail::expression_from_node(lamina::detail::make_node<VariableNode>(name));
}

/// Create a FunctionNode expression (e.g., sin(x), cos(x), tan(x)).
static SymbolicExpr make_func(FunctionNode::FuncType type, const std::string& var_name) {
    auto var_node = lamina::detail::make_node<VariableNode>(var_name);
    auto func_node = lamina::detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{var_node});
    return lamina::detail::expression_from_node(func_node);
}

/// Create a numeric SymbolicExpr from a double value.
static std::shared_ptr<SymbolicExpr> make_number_expr(double val) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(val)));
}

/// Extract a numeric value from a SymbolicExpr (if it's a simple number or product).
static std::optional<double> extract_numeric(const SymbolicExpr& expr) {
    if (!lamina::detail::node(expr)) return std::nullopt;

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
        if (std::holds_alternative<lmmc_real_t>(num->value()))
            return std::get<lmmc_real_t>(num->value());
        if (std::holds_alternative<BigInt>(num->value()))
            return std::get<BigInt>(num->value()).to_double();
        if (std::holds_alternative<Rational>(num->value()))
            return std::get<Rational>(num->value()).to_double();
        return std::nullopt;
    }

    // Handle MultiplyNode (e.g., 2*pi)
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        double product = 1.0;
        for (const auto& op : mul->operands()) {
            auto num = std::dynamic_pointer_cast<const NumberNode>(op);
            if (!num) return std::nullopt;
            if (std::holds_alternative<lmmc_real_t>(num->value()))
                product *= std::get<lmmc_real_t>(num->value());
            else if (std::holds_alternative<BigInt>(num->value()))
                product *= std::get<BigInt>(num->value()).to_double();
            else if (std::holds_alternative<Rational>(num->value()))
                product *= std::get<Rational>(num->value()).to_double();
            else
                return std::nullopt;
        }
        return product;
    }

    return std::nullopt;
}


static void test_periodicity_declared_symbol_roundtrip() {
    TEST_CASE("Declared periodic symbol — get_period returns declared period");

    AssumptionContext ctx;
    auto period_expr = make_number_expr(5.0);
    ctx.current_properties().declare_periodic("f", period_expr);

    InferenceEngine engine(ctx);

    SymbolicExpr f_expr = make_var("f");

    // query_periodic should return True
    Tribool is_periodic = engine.query_periodic(f_expr);
    EXPECT_TRUE(is_periodic == Tribool::True,
        "Declared periodic symbol: query_periodic returns True");

    // infer_period should return the declared period
    auto inferred = engine.infer_period(f_expr);
    EXPECT_TRUE(inferred.has_value(),
        "Declared periodic symbol: infer_period returns a value");

    if (inferred.has_value()) {
        auto val = extract_numeric(*inferred);
        EXPECT_TRUE(val.has_value() && std::abs(*val - 5.0) < 1e-10,
            "Declared periodic symbol: infer_period returns period = 5.0");
    }
}

static void test_periodicity_declared_symbol_various_periods() {
    TEST_CASE("Multiple symbols with different periods — round-trip");

    AssumptionContext ctx;

    // Declare several symbols with different periods
    ctx.current_properties().declare_periodic("g", make_number_expr(2.0));
    ctx.current_properties().declare_periodic("h", make_number_expr(M_PI));
    ctx.current_properties().declare_periodic("k", make_number_expr(100.0));

    InferenceEngine engine(ctx);

    // g: period 2.0
    {
        SymbolicExpr g_expr = make_var("g");
        EXPECT_TRUE(engine.query_periodic(g_expr) == Tribool::True,
            "g is periodic");
        auto period = engine.infer_period(g_expr);
        EXPECT_TRUE(period.has_value(), "g has inferred period");
        if (period.has_value()) {
            auto val = extract_numeric(*period);
            EXPECT_TRUE(val.has_value() && std::abs(*val - 2.0) < 1e-10,
                "g period = 2.0");
        }
    }

    // h: period pi
    {
        SymbolicExpr h_expr = make_var("h");
        EXPECT_TRUE(engine.query_periodic(h_expr) == Tribool::True,
            "h is periodic");
        auto period = engine.infer_period(h_expr);
        EXPECT_TRUE(period.has_value(), "h has inferred period");
        if (period.has_value()) {
            auto val = extract_numeric(*period);
            EXPECT_TRUE(val.has_value() && std::abs(*val - M_PI) < 1e-10,
                "h period = pi");
        }
    }

    // k: period 100.0
    {
        SymbolicExpr k_expr = make_var("k");
        EXPECT_TRUE(engine.query_periodic(k_expr) == Tribool::True,
            "k is periodic");
        auto period = engine.infer_period(k_expr);
        EXPECT_TRUE(period.has_value(), "k has inferred period");
        if (period.has_value()) {
            auto val = extract_numeric(*period);
            EXPECT_TRUE(val.has_value() && std::abs(*val - 100.0) < 1e-10,
                "k period = 100.0");
        }
    }
}

static void test_periodicity_sin_auto_inferred() {
    TEST_CASE("sin(x) has auto-inferred period 2*pi");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr sin_x = make_func(FunctionNode::FuncType::Sin, "x");

    // sin should be periodic
    Tribool is_periodic = engine.query_periodic(sin_x);
    EXPECT_TRUE(is_periodic == Tribool::True,
        "sin(x) is periodic");

    // Period should be 2*pi
    auto period = engine.infer_period(sin_x);
    EXPECT_TRUE(period.has_value(), "sin(x) has inferred period");

    if (period.has_value()) {
        auto val = extract_numeric(*period);
        EXPECT_TRUE(val.has_value() && std::abs(*val - 2.0 * M_PI) < 1e-10,
            "sin(x) period = 2*pi (approx " + std::to_string(val.value_or(0.0)) + ")");
    }
}

static void test_periodicity_cos_auto_inferred() {
    TEST_CASE("cos(x) has auto-inferred period 2*pi");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr cos_x = make_func(FunctionNode::FuncType::Cos, "x");

    Tribool is_periodic = engine.query_periodic(cos_x);
    EXPECT_TRUE(is_periodic == Tribool::True,
        "cos(x) is periodic");

    auto period = engine.infer_period(cos_x);
    EXPECT_TRUE(period.has_value(), "cos(x) has inferred period");

    if (period.has_value()) {
        auto val = extract_numeric(*period);
        EXPECT_TRUE(val.has_value() && std::abs(*val - 2.0 * M_PI) < 1e-10,
            "cos(x) period = 2*pi");
    }
}

static void test_periodicity_tan_auto_inferred() {
    TEST_CASE("tan(x) has auto-inferred period pi");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr tan_x = make_func(FunctionNode::FuncType::Tan, "x");

    Tribool is_periodic = engine.query_periodic(tan_x);
    EXPECT_TRUE(is_periodic == Tribool::True,
        "tan(x) is periodic");

    auto period = engine.infer_period(tan_x);
    EXPECT_TRUE(period.has_value(), "tan(x) has inferred period");

    if (period.has_value()) {
        auto val = extract_numeric(*period);
        EXPECT_TRUE(val.has_value() && std::abs(*val - M_PI) < 1e-10,
            "tan(x) period = pi");
    }
}

static void test_periodicity_non_periodic_function() {
    TEST_CASE("exp(x) is NOT periodic");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr exp_x = make_func(FunctionNode::FuncType::Exp, "x");

    Tribool is_periodic = engine.query_periodic(exp_x);
    EXPECT_TRUE(is_periodic == Tribool::Unknown,
        "exp(x) is not known to be periodic (returns Unknown)");

    auto period = engine.infer_period(exp_x);
    EXPECT_FALSE(period.has_value(),
        "exp(x) has no inferred period");
}

static void test_periodicity_non_periodic_variable() {
    TEST_CASE("Undeclared variable is not periodic");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var("x");

    Tribool is_periodic = engine.query_periodic(x_expr);
    EXPECT_TRUE(is_periodic == Tribool::Unknown,
        "Undeclared variable: query_periodic returns Unknown");

    auto period = engine.infer_period(x_expr);
    EXPECT_FALSE(period.has_value(),
        "Undeclared variable: infer_period returns nullopt");
}

static void test_periodicity_property_store_roundtrip() {
    TEST_CASE("PropertyStore declare_periodic / get_period / is_periodic round-trip");

    PropertyStore store;

    // Not periodic initially
    EXPECT_FALSE(store.is_periodic("f"), "f not periodic initially");
    EXPECT_FALSE(store.get_period("f").has_value(), "f has no period initially");

    // Declare periodic
    auto period = make_number_expr(7.0);
    store.declare_periodic("f", period);

    // Now periodic
    EXPECT_TRUE(store.is_periodic("f"), "f is periodic after declaration");
    EXPECT_TRUE(store.get_period("f").has_value(), "f has period after declaration");

    // Period value matches
    auto retrieved = store.get_period("f");
    if (retrieved.has_value()) {
        auto val = test_numeric_eval(*retrieved);
        EXPECT_TRUE(val.has_value() && std::abs(*val - 7.0) < 1e-10,
            "Retrieved period = 7.0");
    }
}

static void test_periodicity_ln_not_periodic() {
    TEST_CASE("ln(x) is NOT periodic");

    AssumptionContext ctx;
    InferenceEngine engine(ctx);

    SymbolicExpr ln_x = make_func(FunctionNode::FuncType::Ln, "x");

    Tribool is_periodic = engine.query_periodic(ln_x);
    EXPECT_TRUE(is_periodic == Tribool::Unknown,
        "ln(x) is not known to be periodic");

    auto period = engine.infer_period(ln_x);
    EXPECT_FALSE(period.has_value(),
        "ln(x) has no inferred period");
}


int main() {
    test_periodicity_declared_symbol_roundtrip();
    test_periodicity_declared_symbol_various_periods();
    test_periodicity_sin_auto_inferred();
    test_periodicity_cos_auto_inferred();
    test_periodicity_tan_auto_inferred();
    test_periodicity_non_periodic_function();
    test_periodicity_non_periodic_variable();
    test_periodicity_property_store_roundtrip();
    test_periodicity_ln_not_periodic();

    return TEST_REPORT();
}
