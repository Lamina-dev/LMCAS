/**
 * @file test_transform_engine.cpp
 * @brief 积分变换引擎单元测试：Fourier 变换、逆 Fourier 变换、卷积。
 */

#include "test_common.hpp"
#include "transform_engine.hpp"

using namespace LMCAS;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using Expr = std::shared_ptr<SymbolicExpr>;

static Expr num(int n) { return SymbolicExpr::number(n); }
static Expr var(const std::string& name) { return SymbolicExpr::variable(name); }

// Evaluate a single-variable symbolic result at a numeric point.
static double eval_at(const Expr& e, const std::string& v, double x) {
    auto sub = e->substitute(v, SymbolicExpr::number(x));
    return sub->simplify()->to_numeric();
}

// Evaluate a two-variable symbolic result at numeric points.
static double eval_at2(const Expr& e, const std::string& v1, double x1,
                       const std::string& v2, double x2) {
    auto sub = e->substitute(v1, SymbolicExpr::number(x1))
                ->substitute(v2, SymbolicExpr::number(x2));
    return sub->simplify()->to_numeric();
}


static void test_fourier_gaussian() {
    TEST_CASE("fourier_transform: Gaussian e^(-t^2)");

    auto t = var("t");
    auto t_sq = SymbolicExpr::power(t, num(2));
    auto neg_t_sq = SymbolicExpr::multiply(num(-1), t_sq);
    auto f = SymbolicExpr::exp(neg_t_sq);

    auto result = LMCAS::fourier_transform_checked(f, "t", "omega").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    // ℱ{e^(-t²)} = sqrt(π) * e^(-ω²/4)
    // Check it's not an unevaluated TransformNode
    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "Result is evaluated (not TransformNode)");

    std::string s = result->to_string();
    std::cout << "  Fourier(e^(-t^2)) = " << s << std::endl;

    // Verify value at omega=0 (should be sqrt(pi)) and omega=2 (sqrt(pi)*e^-1)
    EXPECT_NEAR(eval_at(result, "omega", 0.0), std::sqrt(M_PI), 1e-6,
                "F{e^-t^2}(0) = sqrt(pi)");
    EXPECT_NEAR(eval_at(result, "omega", 2.0), std::sqrt(M_PI) * std::exp(-1.0), 1e-6,
                "F{e^-t^2}(2) = sqrt(pi)*e^-1");
}

static void test_fourier_gaussian_with_coeff() {
    TEST_CASE("fourier_transform: Gaussian e^(-2*t^2)");

    auto t = var("t");
    auto t_sq = SymbolicExpr::power(t, num(2));
    auto neg_2t_sq = SymbolicExpr::multiply(num(-2), t_sq);
    auto f = SymbolicExpr::exp(neg_2t_sq);

    auto result = LMCAS::fourier_transform_checked(f, "t", "omega").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "Result is evaluated (not TransformNode)");

    std::string s = result->to_string();
    std::cout << "  Fourier(e^(-2t^2)) = " << s << std::endl;
}

static void test_fourier_exp_decay() {
    TEST_CASE("fourier_transform: exponential decay e^(-a|t|)");

    auto t = var("t");
    auto a = var("a");
    auto abs_t = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(t)}));
    auto neg_a_abs_t = SymbolicExpr::multiply(num(-1),
        SymbolicExpr::multiply(a, abs_t));
    auto f = SymbolicExpr::exp(neg_a_abs_t);

    auto result = LMCAS::fourier_transform_checked(f, "t", "omega").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "Result is evaluated (not TransformNode)");

    // ℱ{e^(-a|t|)} = 2a/(a² + ω²)
    std::string s = result->to_string();
    std::cout << "  Fourier(e^(-a|t|)) = " << s << std::endl;

    EXPECT_NEAR(eval_at2(result, "a", 1.0, "omega", 0.0), 2.0, 1e-6,
                "F{e^-|t|}(0) = 2");
    EXPECT_NEAR(eval_at2(result, "a", 2.0, "omega", 2.0), 4.0 / 8.0, 1e-6,
                "F{e^-2|t|}(2) = 2*2/(4+4) = 0.5");
}

static void test_fourier_causal_exp() {
    TEST_CASE("fourier_transform: causal exponential e^(-at)");

    auto t = var("t");
    auto a = var("a");
    auto neg_at = SymbolicExpr::multiply(num(-1), SymbolicExpr::multiply(a, t));
    auto f = SymbolicExpr::exp(neg_at);

    auto result = LMCAS::fourier_transform_checked(f, "t", "omega").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "Result is evaluated (not TransformNode)");

    // ℱ{e^(-at)} = 1/(a + iω)
    std::string s = result->to_string();
    std::cout << "  Fourier(e^(-at)) = " << s << std::endl;
}

static void test_fourier_constant_returns_unevaluated() {
    TEST_CASE("fourier_transform: unsupported constant is Inconclusive");
    auto result = LMCAS::fourier_transform_checked(
        num(5), "t", "omega");
    EXPECT_TRUE(!result &&
                    result.error().code == LMCAS::CasErrc::Inconclusive,
                "constant Fourier transform is explicitly Inconclusive");
}

static void test_fourier_linearity() {
    TEST_CASE("fourier_transform: linearity for sum of Gaussians");

    auto t = var("t");
    auto t_sq = SymbolicExpr::power(t, num(2));
    auto neg_t_sq = SymbolicExpr::multiply(num(-1), t_sq);
    auto g1 = SymbolicExpr::exp(neg_t_sq);

    auto neg_2t_sq = SymbolicExpr::multiply(num(-2), t_sq);
    auto g2 = SymbolicExpr::exp(neg_2t_sq);

    auto f = SymbolicExpr::add(g1, g2);

    auto result = LMCAS::fourier_transform_checked(f, "t", "omega").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "Sum of Gaussians is evaluated");

    std::string s = result->to_string();
    std::cout << "  Fourier(e^(-t^2) + e^(-2t^2)) = " << s << std::endl;
}

static void test_fourier_with_constant_factor() {
    TEST_CASE("fourier_transform: constant factor 3*e^(-t^2)");

    auto t = var("t");
    auto t_sq = SymbolicExpr::power(t, num(2));
    auto neg_t_sq = SymbolicExpr::multiply(num(-1), t_sq);
    auto g = SymbolicExpr::exp(neg_t_sq);
    auto f = SymbolicExpr::multiply(num(3), g);

    auto result = LMCAS::fourier_transform_checked(f, "t", "omega").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "3*Gaussian is evaluated");

    std::string s = result->to_string();
    std::cout << "  Fourier(3*e^(-t^2)) = " << s << std::endl;
}

static void test_fourier_unknown_returns_unevaluated() {
    TEST_CASE("fourier_transform: unknown function is Inconclusive");
    auto t = var("t");
    auto result = LMCAS::fourier_transform_checked(
        SymbolicExpr::ln(t), "t", "omega");
    EXPECT_TRUE(!result &&
                    result.error().code == LMCAS::CasErrc::Inconclusive,
                "unknown Fourier transform is explicitly Inconclusive");
}


static void test_inverse_fourier_gaussian() {
    TEST_CASE("inverse_fourier_transform: approximate Gaussian round trip is unproved");
    auto w = var("omega");
    auto F = SymbolicExpr::exp(
        SymbolicExpr::multiply(
            num(-1), SymbolicExpr::power(w, num(2))));
    auto result = LMCAS::inverse_fourier_transform_checked(
        F, "omega", "t");
    EXPECT_TRUE(!result &&
                    result.error().code == LMCAS::CasErrc::Inconclusive,
                "approximate Gaussian inverse requires a round-trip proof");
}

static void test_inverse_fourier_constant_returns_unevaluated() {
    TEST_CASE("inverse_fourier_transform: constant is Inconclusive");
    auto result = LMCAS::inverse_fourier_transform_checked(
        num(3), "omega", "t");
    EXPECT_TRUE(!result &&
                    result.error().code == LMCAS::CasErrc::Inconclusive,
                "constant inverse Fourier transform is explicitly Inconclusive");
}


static void test_convolve_returns_result() {
    TEST_CASE("convolve: 正高斯双边卷积返回受证闭式");
    auto x = var("x");
    auto x_squared = SymbolicExpr::power(x, num(2));
    auto first = SymbolicExpr::multiply(
        num(2), SymbolicExpr::exp(
            SymbolicExpr::multiply(num(-1), x_squared)));
    auto second = SymbolicExpr::multiply(
        num(3), SymbolicExpr::exp(
            SymbolicExpr::multiply(num(-2), x_squared)));

    auto result = LMCAS::convolve_checked(first, second, "x");
    EXPECT_TRUE(result.has_value(), "正高斯卷积应返回闭式");
    if (result) {
        auto expression = result.value().value.expression;
        EXPECT_TRUE(expression != nullptr, "卷积闭式不能为空");
        EXPECT_TRUE(
            std::holds_alternative<LMCAS::ByConstructionProof>(
                result.value().certificate),
            "高斯卷积应具有构造证明");
        const auto text = expression->to_string();
        EXPECT_TRUE(
            text.find("exp") != std::string::npos &&
                text.find("pi") != std::string::npos,
            "闭式应包含高斯指数与 pi 归一化因子");
    }

    auto zero = LMCAS::convolve_checked(num(0), first, "x");
    EXPECT_TRUE(
        zero && zero.value().value.expression->is_zero(),
        "零函数卷积应为零");

    auto sum = SymbolicExpr::add(first, second);
    auto linear = LMCAS::convolve_checked(sum, first, "x");
    EXPECT_TRUE(linear.has_value(), "卷积应对加法保持线性");
}

static void test_convolve_null_inputs() {
    TEST_CASE("convolve: null inputs return InvalidArgument");
    auto first = LMCAS::convolve_checked(nullptr, num(1), "x");
    EXPECT_TRUE(!first &&
                    first.error().code == LMCAS::CasErrc::InvalidArgument,
                "null first convolution input is InvalidArgument");
    auto second = LMCAS::convolve_checked(num(1), nullptr, "x");
    EXPECT_TRUE(!second &&
                    second.error().code == LMCAS::CasErrc::InvalidArgument,
                "null second convolution input is InvalidArgument");
}


static void test_laplace_constant() {
    TEST_CASE("laplace_transform: constant c -> c/s");

    auto result = LMCAS::laplace_transform_checked(num(5), "t", "s").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "Constant Laplace is evaluated");

    std::string s = result->to_string();
    std::cout << "  Laplace(5) = " << s << std::endl;
    // ℒ{5} = 5/s ; at s=5 -> 1
    EXPECT_NEAR(eval_at(result, "s", 5.0), 1.0, 1e-6, "L{5}(s=5) = 5/5 = 1");
}

static void test_laplace_exp() {
    TEST_CASE("laplace_transform: e^(at) -> 1/(s-a)");

    auto t = var("t");
    auto a = var("a");
    auto at = SymbolicExpr::multiply(a, t);
    auto f = SymbolicExpr::exp(at);

    auto result = LMCAS::laplace_transform_checked(f, "t", "s").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "e^(at) Laplace is evaluated");

    std::string s = result->to_string();
    std::cout << "  Laplace(e^(at)) = " << s << std::endl;
    // ℒ{e^(at)} = 1/(s-a) ; at a=1, s=3 -> 1/2
    EXPECT_NEAR(eval_at2(result, "a", 1.0, "s", 3.0), 0.5, 1e-6,
                "L{e^t}(s=3) = 1/(3-1) = 0.5");
}

static void test_laplace_sin() {
    TEST_CASE("laplace_transform: sin(at) -> a/(s^2+a^2)");

    auto t = var("t");
    auto a = var("a");
    auto at = SymbolicExpr::multiply(a, t);
    auto f = SymbolicExpr::sin(at);

    auto result = LMCAS::laplace_transform_checked(f, "t", "s").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "sin(at) Laplace is evaluated");

    std::string s = result->to_string();
    std::cout << "  Laplace(sin(at)) = " << s << std::endl;
    // ℒ{sin(at)} = a/(s²+a²) ; at a=2, s=1 -> 2/5 = 0.4
    EXPECT_NEAR(eval_at2(result, "a", 2.0, "s", 1.0), 0.4, 1e-6,
                "L{sin(2t)}(s=1) = 2/(1+4) = 0.4");
}

static void test_laplace_t_power() {
    TEST_CASE("laplace_transform: t^n -> n!/s^(n+1)");

    auto t = var("t");
    auto t_sq = SymbolicExpr::power(t, num(2));

    auto result = LMCAS::laplace_transform_checked(t_sq, "t", "s").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "t^2 Laplace is evaluated");

    // ℒ{t²} = 2/s³
    std::string s = result->to_string();
    std::cout << "  Laplace(t^2) = " << s << std::endl;
    // at s=2 -> 2/8 = 0.25
    EXPECT_NEAR(eval_at(result, "s", 2.0), 0.25, 1e-6, "L{t^2}(s=2) = 2/8 = 0.25");
}


static void test_z_transform_constant() {
    TEST_CASE("z_transform: constant c -> c*z/(z-1)");

    auto result = LMCAS::z_transform_checked(num(3), "n", "z").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "Constant Z-transform is evaluated");

    std::string s = result->to_string();
    std::cout << "  Z(3) = " << s << std::endl;
    // Z{3} = 3z/(z-1) ; at z=2 -> 6
    EXPECT_NEAR(eval_at(result, "z", 2.0), 6.0, 1e-6, "Z{3}(z=2) = 3*2/(2-1) = 6");
}

static void test_z_transform_exp_sequence() {
    TEST_CASE("z_transform: a^n -> z/(z-a)");

    auto n = var("n");
    auto a = var("a");
    auto f = SymbolicExpr::power(a, n);

    auto result = LMCAS::z_transform_checked(f, "n", "z").value().value.expression;
    EXPECT_TRUE(result != nullptr, "Result is not null");

    auto tn = std::dynamic_pointer_cast<const TransformNode>(LMCAS::detail::node(result));
    EXPECT_TRUE(tn == nullptr, "a^n Z-transform is evaluated");

    std::string s = result->to_string();
    std::cout << "  Z(a^n) = " << s << std::endl;
    // Z{a^n} = z/(z-a) ; at a=2, z=4 -> 4/2 = 2
    EXPECT_NEAR(eval_at2(result, "a", 2.0, "z", 4.0), 2.0, 1e-6,
                "Z{2^n}(z=4) = 4/(4-2) = 2");
}

static void test_transform_checked_contracts() {
    TEST_CASE("transform_engine checked APIs: explicit errors, context, and inconclusive support domains");

    auto t = var("t");
    auto s = var("s");

    auto laplace = LMCAS::laplace_transform_checked(num(5), "t", "s");
    EXPECT_TRUE(laplace.has_value(), "checked Laplace succeeds for constants");
    if (laplace) {
        EXPECT_TRUE(laplace.value().value.expression != nullptr,
                    "checked Laplace returns an expression");
        EXPECT_TRUE(laplace.value().value.roc.size() == 1,
                    "checked Laplace reports ROC for constant input");
        if (!laplace.value().value.roc.empty()) {
            auto roc = std::dynamic_pointer_cast<const RelationalNode>(
                LMCAS::detail::node(laplace.value().value.roc[0]));
            EXPECT_TRUE(roc != nullptr && roc->op() == RelationalNode::Op::GT,
                        "checked Laplace ROC is a greater-than condition");
            auto lhs = roc ? std::dynamic_pointer_cast<const VariableNode>(roc->left()) : nullptr;
            auto rhs = roc ? std::dynamic_pointer_cast<const NumberNode>(roc->right()) : nullptr;
            EXPECT_TRUE(lhs != nullptr && lhs->name() == "s" &&
                            rhs != nullptr && rhs->is_zero(),
                        "checked Laplace constant ROC is s > 0");
        }
    }

    auto unknown = LMCAS::fourier_transform_checked(SymbolicExpr::ln(t), "t", "omega");
    EXPECT_TRUE(!unknown.has_value(),
                "checked Fourier rejects unsupported closed forms");
    EXPECT_TRUE(unknown.error().code == LMCAS::CasErrc::Inconclusive,
                "checked Fourier reports Inconclusive for unevaluated transform nodes");

    auto nested_unknown = LMCAS::fourier_transform_checked(
        SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::ln(t)),
        "t", "omega");
    EXPECT_TRUE(!nested_unknown &&
                    nested_unknown.error().code == LMCAS::CasErrc::Inconclusive,
                "checked Fourier rejects nested unevaluated transform nodes");

    auto unsupported_inverse = LMCAS::inverse_fourier_transform_checked(num(1), "omega", "t");
    EXPECT_TRUE(!unsupported_inverse.has_value(),
                "checked inverse Fourier rejects unevaluated transform nodes");
    EXPECT_TRUE(unsupported_inverse.error().code == LMCAS::CasErrc::Inconclusive,
                "checked inverse Fourier reports Inconclusive for unsupported constants");

    auto null_input = LMCAS::laplace_transform_checked(nullptr, "t", "s");
    EXPECT_TRUE(!null_input.has_value(), "checked Laplace rejects null input");
    EXPECT_TRUE(null_input.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked Laplace reports InvalidArgument for null input");

    auto empty_var = LMCAS::inverse_laplace_checked(s, "", "t");
    EXPECT_TRUE(!empty_var.has_value(), "checked inverse Laplace rejects empty variable");
    EXPECT_TRUE(empty_var.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked inverse Laplace reports InvalidArgument for empty variable");

    auto same_var = LMCAS::z_transform_checked(t, "n", "n");
    EXPECT_TRUE(!same_var.has_value(), "checked Z transform rejects same input/output variable");
    EXPECT_TRUE(same_var.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked Z transform reports InvalidArgument for same variables");

    auto z_const = LMCAS::z_transform_checked(num(5), "n", "z");
    EXPECT_TRUE(z_const.has_value(), "checked Z transform succeeds for constants");
    if (z_const) {
        EXPECT_TRUE(z_const.value().value.roc.size() == 1,
                    "checked Z constant reports one ROC condition");
        auto roc = z_const.value().value.roc.empty() ? nullptr :
            std::dynamic_pointer_cast<const RelationalNode>(
                LMCAS::detail::node(z_const.value().value.roc[0]));
        EXPECT_TRUE(roc != nullptr && roc->op() == RelationalNode::Op::GT,
                    "checked Z constant ROC is a greater-than condition");
        auto lhs_abs = roc ? std::dynamic_pointer_cast<const FunctionNode>(roc->left()) : nullptr;
        auto rhs_expr = roc ? LMCAS::detail::make_expression_ptr(roc->right()) : nullptr;
        auto abs_arg = (lhs_abs && lhs_abs->arguments().size() == 1) ?
            std::dynamic_pointer_cast<const VariableNode>(lhs_abs->arguments()[0]) : nullptr;
        EXPECT_TRUE(lhs_abs != nullptr && lhs_abs->type() == FunctionNode::FuncType::Abs &&
                        abs_arg != nullptr && abs_arg->name() == "z" &&
                        rhs_expr != nullptr && rhs_expr->is_one(),
                    "checked Z constant ROC is abs(z) > 1");
    }

    auto z_exp = LMCAS::z_transform_checked(SymbolicExpr::power(num(2), var("n")), "n", "z");
    EXPECT_TRUE(z_exp.has_value(), "checked Z transform succeeds for a^n");
    if (z_exp) {
        EXPECT_TRUE(z_exp.value().value.roc.size() == 1,
                    "checked Z exponential reports one ROC condition");
        auto roc = z_exp.value().value.roc.empty() ? nullptr :
            std::dynamic_pointer_cast<const RelationalNode>(
                LMCAS::detail::node(z_exp.value().value.roc[0]));
        auto lhs_abs = roc ? std::dynamic_pointer_cast<const FunctionNode>(roc->left()) : nullptr;
        auto rhs_abs = roc ? std::dynamic_pointer_cast<const FunctionNode>(roc->right()) : nullptr;
        auto lhs_arg = (lhs_abs && lhs_abs->arguments().size() == 1) ?
            std::dynamic_pointer_cast<const VariableNode>(lhs_abs->arguments()[0]) : nullptr;
        auto rhs_arg = (rhs_abs && rhs_abs->arguments().size() == 1) ?
            rhs_abs->arguments()[0] : nullptr;
        EXPECT_TRUE(roc != nullptr && roc->op() == RelationalNode::Op::GT &&
                        lhs_abs != nullptr && lhs_abs->type() == FunctionNode::FuncType::Abs &&
                        rhs_abs != nullptr && rhs_abs->type() == FunctionNode::FuncType::Abs &&
                        lhs_arg != nullptr && lhs_arg->name() == "z" &&
                        rhs_arg != nullptr &&
                        rhs_arg->equals(*LMCAS::detail::node(num(2))),
                    "checked Z exponential ROC is abs(z) > abs(2)");
    }

    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = LMCAS::fourier_transform_checked(t, "t", "omega", cancelled_context);
    EXPECT_TRUE(!cancelled.has_value(), "checked Fourier observes cancellation");
    EXPECT_TRUE(cancelled.error().code == LMCAS::CasErrc::Cancelled,
                "checked Fourier reports Cancelled");

    LMCAS::ResourceLimits limits;
    limits.max_steps = 0;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::convolve_checked(t, t, "t", limited_context);
    EXPECT_TRUE(!limited.has_value(), "checked convolution observes exhausted step budget");
    EXPECT_TRUE(limited.error().code == LMCAS::CasErrc::ResourceLimit,
                "checked convolution reports ResourceLimit");

    auto convolution_null = LMCAS::convolve_checked(nullptr, t, "t");
    EXPECT_TRUE(!convolution_null.has_value(), "checked convolution rejects null input");
    EXPECT_TRUE(convolution_null.error().code == LMCAS::CasErrc::InvalidArgument,
                "checked convolution reports InvalidArgument for null input");

    auto unsupported_convolution = LMCAS::convolve_checked(SymbolicExpr::ln(t), t, "t");
    EXPECT_TRUE(!unsupported_convolution.has_value(),
                "checked convolution rejects unevaluated integral results");
    EXPECT_TRUE(unsupported_convolution.error().code == LMCAS::CasErrc::Inconclusive,
                "checked convolution reports Inconclusive for unsupported integrals");
}


int main() {
    // Fourier transform tests
    test_fourier_gaussian();
    test_fourier_gaussian_with_coeff();
    test_fourier_exp_decay();
    test_fourier_causal_exp();
    test_fourier_constant_returns_unevaluated();
    test_fourier_linearity();
    test_fourier_with_constant_factor();
    test_fourier_unknown_returns_unevaluated();

    // Inverse Fourier transform tests
    test_inverse_fourier_gaussian();
    test_inverse_fourier_constant_returns_unevaluated();

    // Convolution tests
    test_convolve_returns_result();
    test_convolve_null_inputs();

    // Laplace transform tests
    test_laplace_constant();
    test_laplace_exp();
    test_laplace_sin();
    test_laplace_t_power();

    // Z transform tests
    test_z_transform_constant();
    test_z_transform_exp_sequence();
    test_transform_checked_contracts();

    std::cout << "\n===================================================\n";
    std::cout << "Results: " << g_passes << " passed, "
              << g_failures << " failed\n";
    return TEST_REPORT();
}
