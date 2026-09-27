
#include "test_interval_support.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "assumption.hpp"
#include "interval.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace LMCAS;

static SymbolicExpr make_var_expr(const std::string &name) {
    return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<VariableNode>(name));
}

static SymbolicExpr make_func_expr(FunctionNode::FuncType type, const std::string &var_name) {
    auto var_node = LMCAS::detail::make_node<VariableNode>(var_name);
    auto func_node = LMCAS::detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{var_node});
    return LMCAS::detail::expression_from_node(func_node);
}

TEST(AssumptionMonotonicityIntervals, ExpIncreasingOnReals) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");

    Interval entire = Interval::entire_line();
    Monotonicity mono = engine.infer_monotonicity(exp_x, "x", entire);

    EXPECT_TRUE((mono == Monotonicity::Increasing)) << "exp(x) is Increasing on entire real line";
}

TEST(AssumptionMonotonicityIntervals, ExpIncreasingOnFiniteInterval) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");

    Interval iv = make_closed_interval(0.0, 10.0);
    Monotonicity mono = engine.infer_monotonicity(exp_x, "x", iv);

    EXPECT_TRUE((mono == Monotonicity::Increasing)) << "exp(x) is Increasing on [0, 10]";
}

TEST(AssumptionMonotonicityIntervals, LnIncreasingOnPositiveReals) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");

    Interval pos_reals = make_open_interval(0.0, 1000.0);
    Monotonicity mono = engine.infer_monotonicity(ln_x, "x", pos_reals);

    EXPECT_TRUE((mono == Monotonicity::Increasing)) << "ln(x) is Increasing on (0, 1000)";
}

TEST(AssumptionMonotonicityIntervals, LnIncreasingOnClosedPositive) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");

    Interval iv = make_closed_interval(1.0, 100.0);
    Monotonicity mono = engine.infer_monotonicity(ln_x, "x", iv);

    EXPECT_TRUE((mono == Monotonicity::Increasing)) << "ln(x) is Increasing on [1, 100]";
}

TEST(AssumptionMonotonicityIntervals, NegationReversesMonotonicity) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto exp_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Exp,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto neg_one = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
    auto neg_exp = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{neg_one, exp_node});
    auto neg_exp_x = LMCAS::detail::expression_from_node(neg_exp);
    Interval entire = Interval::entire_line();
    Monotonicity mono = engine.infer_monotonicity(neg_exp_x, "x", entire);

    EXPECT_TRUE((mono == Monotonicity::Decreasing)) << "-exp(x) is Decreasing (negation reverses Increasing)";
}

TEST(AssumptionMonotonicityIntervals, NegationReversesLn) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto ln_node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Ln,
        std::vector<std::shared_ptr<const SymbolicNode>>{x_node});
    auto neg_one = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
    auto neg_ln = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{neg_one, ln_node});
    auto neg_ln_x = LMCAS::detail::expression_from_node(neg_ln);
    Interval iv = make_closed_interval(1.0, 100.0);
    Monotonicity mono = engine.infer_monotonicity(neg_ln_x, "x", iv);

    EXPECT_TRUE((mono == Monotonicity::Decreasing)) << "-ln(x) is Decreasing on [1, 100]";
}

TEST(AssumptionMonotonicityIntervals, DeclaredMonotonicityDeduction) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_domain("y", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");
    SymbolicExpr exp_y = make_func_expr(FunctionNode::FuncType::Exp, "y");

    EXPECT_TRUE((ctx.current_relations().has_relation(exp_x, exp_y, RelationalNode::Op::GT))) << "x > y with exp increasing => exp(x) > exp(y)";
}

TEST(AssumptionMonotonicityIntervals, LnDeductionFromInequality) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());
    ASSERT_TRUE(ctx.assume_sign("y", Sign::Positive).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr x_expr = make_var_expr("x");
    SymbolicExpr y_expr = make_var_expr("y");

    Relation rel{x_expr, y_expr, RelationalNode::Op::GT};
    EXPECT_TRUE((ctx.current_relations().add_relation(x_expr, y_expr, RelationalNode::Op::GT, ctx.current_properties()).has_value())) << "relation insertion succeeds";

    engine.apply_monotonicity_rules(rel, ctx.current_relations(), ctx.current_properties());

    SymbolicExpr ln_x = make_func_expr(FunctionNode::FuncType::Ln, "x");
    SymbolicExpr ln_y = make_func_expr(FunctionNode::FuncType::Ln, "y");

    EXPECT_TRUE((ctx.current_relations().has_relation(ln_x, ln_y, RelationalNode::Op::GT))) << "x > y with both Positive => ln(x) > ln(y)";
}

TEST(AssumptionMonotonicityIntervals, UnknownForNonMonotoneFunction) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr sin_x = make_func_expr(FunctionNode::FuncType::Sin, "x");

    Interval iv = make_closed_interval(0.0, 2.0 * M_PI);
    Monotonicity mono = engine.infer_monotonicity(sin_x, "x", iv);

    EXPECT_TRUE((mono == Monotonicity::Unknown)) << "sin(x) on [0, 2*pi] has Unknown monotonicity (not monotone on full period)";
}

TEST(AssumptionMonotonicityIntervals, WrongVariableReturnsUnknown) {
    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    InferenceEngine engine(ctx);

    SymbolicExpr exp_x = make_func_expr(FunctionNode::FuncType::Exp, "x");

    Interval iv = make_closed_interval(0.0, 10.0);
    Monotonicity mono = engine.infer_monotonicity(exp_x, "y", iv);

    EXPECT_TRUE((mono == Monotonicity::Unknown)) << "exp(x) w.r.t. y returns Unknown (wrong variable)";
}

TEST(AssumptionMonotonicityIntervals, ArctanMonotonicityAndPeriodicity) {
    AssumptionContext context;
    EXPECT_TRUE((context.assume_domain("x", Domain::Real).has_value())) << "实数域假设应成功";
    InferenceEngine engine(context);
    SymbolicExpr atan_x =
        make_func_expr(FunctionNode::FuncType::ArcTan, "x");

    EXPECT_TRUE((engine.infer_monotonicity(
                     atan_x, "x", Interval::entire_line()) ==
                 Monotonicity::Increasing))
        << "atan(x) 在实数域严格递增";
    EXPECT_TRUE((engine.query_periodic_checked(atan_x, "x").value() ==
                 Tribool::False))
        << "未声明周期性的自变量经 atan 后非周期";

    auto sin_x = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sin,
        std::vector<std::shared_ptr<const SymbolicNode>>{
            LMCAS::detail::make_node<VariableNode>("x")});
    SymbolicExpr atan_sin = LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ArcTan,
            std::vector<std::shared_ptr<const SymbolicNode>>{sin_x}));
    EXPECT_TRUE((engine.query_periodic_checked(atan_sin, "x").value() ==
                 Tribool::True))
        << "atan 应保留已证明参数的周期";
}
