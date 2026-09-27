/**
 * @file test_transform_engine.cpp
 * @brief 积分变换引擎单元测试：Fourier 变换、逆 Fourier 变换、卷积。
 */

#include "test_common.hpp"
#include "transform_engine.hpp"
#include "assumption_context.hpp"

using namespace LMCAS;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using Expr = std::shared_ptr<SymbolicExpr>;

static Expr num(int n) { return SymbolicExpr::number(n); }
static Expr var(const std::string &name) { return SymbolicExpr::variable(name); }

// Evaluate a single-variable symbolic result at a numeric point.
static double eval_at(const Expr &e, const std::string &v, double x) {
    auto sub = e->substitute(v, SymbolicExpr::number(x));
    return sub->simplify()->to_numeric();
}

// Evaluate a two-variable symbolic result at numeric points.
static double eval_at2(const Expr &e, const std::string &v1, double x1,
                       const std::string &v2, double x2) {
    auto sub = e->substitute(v1, SymbolicExpr::number(x1))
                   ->substitute(v2, SymbolicExpr::number(x2));
    return sub->simplify()->to_numeric();
}

TEST(TransformEngine, FourierGaussian) {
    auto t = var("t");
    auto t_sq = SymbolicExpr::power(t, num(2));
    auto neg_t_sq = SymbolicExpr::multiply(num(-1), t_sq);
    auto f = SymbolicExpr::exp(neg_t_sq);

    auto checked = LMCAS::fourier_transform_checked(f, "t", "omega");
    ASSERT_TRUE(checked) << "Gaussian Fourier transform succeeds";
    auto result = checked.value().value.expression;
    ASSERT_TRUE(result);

    // Verify value at omega=0 (should be sqrt(pi)) and omega=2 (sqrt(pi)*e^-1)
    {
        const double actual_value = (eval_at(result, "omega", 0.0));
        const double expected_value = (std::sqrt(M_PI));
        const double tolerance = (1e-6);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
    {
        const double actual_value = (eval_at(result, "omega", 2.0));
        const double expected_value = (std::sqrt(M_PI) * std::exp(-1.0));
        const double tolerance = (1e-6);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}

TEST(TransformEngine, FourierRequiresExactSquareExponent) {
    auto t = var("t");
    auto near_two = SymbolicExpr::number(Rational(
        BigInt("20000000000001"), BigInt("10000000000000")));
    auto input = SymbolicExpr::exp(SymbolicExpr::multiply(
        num(-1), SymbolicExpr::power(t, near_two)));

    auto result = fourier_transform_checked(input, "t", "omega");
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::Inconclusive)
        << "an exact exponent 2+10^-13 is not a Gaussian square";
}

TEST(TransformEngine, FourierExpDecay) {
    auto t = var("t");
    auto a = var("a");
    auto abs_t = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(t)}));
    auto neg_a_abs_t = SymbolicExpr::multiply(num(-1),
                                              SymbolicExpr::multiply(a, abs_t));
    auto f = SymbolicExpr::exp(neg_a_abs_t);

    auto assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE((assumptions->assume_domain_checked("a", Domain::Real).has_value())) << "decay parameter is real";
    ComputationContext context;
    EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "attach real decay assumptions";
    auto checked =
        LMCAS::fourier_transform_checked(f, "t", "omega", context);
    ASSERT_TRUE(checked) << "bilateral exponential decay transforms";
    auto result = checked.value().value.expression;
    ASSERT_TRUE(result);

    EXPECT_NEAR(eval_at2(result, "a", 1.0, "omega", 0.0), 2.0, 1e-6) << "F{e^-|t|}(0) = 2";
    {
        const double actual_value = (eval_at2(result, "a", 2.0, "omega", 2.0));
        const double expected_value = (4.0 / 8.0);
        const double tolerance = (1e-6);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}

TEST(TransformEngine, FourierBareExponentialIsNotBilateralDecay) {
    auto t = var("t");
    auto a = var("a");
    auto assumptions = std::make_shared<AssumptionContext>();
    ASSERT_TRUE(assumptions->assume_domain_checked("a", Domain::Real));
    ASSERT_TRUE(assumptions->assume_sign_checked("a", Sign::Positive));
    ComputationContext context;
    ASSERT_TRUE(context.set_assumptions(assumptions));
    auto input = SymbolicExpr::exp(SymbolicExpr::multiply(
        num(-1), SymbolicExpr::multiply(a, t)));

    auto result = fourier_transform_checked(input, "t", "omega", context);
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::Inconclusive)
        << "a bilateral Fourier transform cannot use a one-sided decay rule";
}

TEST(TransformEngine, FourierConstantReturnsUnevaluated) {
    auto result = LMCAS::fourier_transform_checked(
        num(5), "t", "omega");
    EXPECT_TRUE((!result &&
                 result.error().code == LMCAS::CasErrc::Inconclusive))
        << "constant Fourier transform is explicitly Inconclusive";
}

TEST(TransformEngine, FourierUnknownReturnsUnevaluated) {
    auto t = var("t");
    auto result = LMCAS::fourier_transform_checked(
        SymbolicExpr::ln(t), "t", "omega");
    EXPECT_TRUE((!result &&
                 result.error().code == LMCAS::CasErrc::Inconclusive))
        << "unknown Fourier transform is explicitly Inconclusive";
}

TEST(TransformEngine, InverseFourierGaussian) {
    auto w = var("omega");
    auto F = SymbolicExpr::exp(
        SymbolicExpr::multiply(
            num(-1), SymbolicExpr::power(w, num(2))));
    auto result = LMCAS::inverse_fourier_transform_checked(
        F, "omega", "t");
    EXPECT_TRUE((!result &&
                 result.error().code == LMCAS::CasErrc::Inconclusive))
        << "approximate Gaussian inverse requires a round-trip proof";
}

TEST(TransformEngine, InverseFourierConstantReturnsUnevaluated) {
    auto result = LMCAS::inverse_fourier_transform_checked(
        num(3), "omega", "t");
    EXPECT_TRUE((!result &&
                 result.error().code == LMCAS::CasErrc::Inconclusive))
        << "constant inverse Fourier transform is explicitly Inconclusive";
}

TEST(TransformEngine, ConvolveReturnsResult) {
    auto x = var("x");
    auto x_squared = SymbolicExpr::power(x, num(2));
    auto first = SymbolicExpr::multiply(
        num(2), SymbolicExpr::exp(
                    SymbolicExpr::multiply(num(-1), x_squared)));
    auto second = SymbolicExpr::multiply(
        num(3), SymbolicExpr::exp(
                    SymbolicExpr::multiply(num(-2), x_squared)));

    auto result = LMCAS::convolve_checked(first, second, "x");
    ASSERT_TRUE(result) << "positive Gaussian convolution has a closed form";
    auto expression = result.value().value.expression;
    ASSERT_TRUE(expression);
    EXPECT_TRUE(std::holds_alternative<LMCAS::ByConstructionProof>(
        result.value().certificate));

    const double scale = 6.0 * std::sqrt(M_PI / 3.0);
    EXPECT_NEAR(eval_at(expression, "x", 0.0), scale, 1e-10)
        << "Gaussian convolution has the analytic value at zero";
    EXPECT_NEAR(eval_at(expression, "x", 1.0),
                scale * std::exp(-2.0 / 3.0), 1e-10)
        << "Gaussian convolution has the analytic value at one";

    auto zero = LMCAS::convolve_checked(num(0), first, "x");
    EXPECT_TRUE((zero && zero.value().value.expression->is_zero())) << "零函数卷积应为零";

}

TEST(TransformEngine, ConvolutionRequiresExactSquareExponent) {
    auto x = var("x");
    auto near_two = SymbolicExpr::number(Rational(
        BigInt("20000000000001"), BigInt("10000000000000")));
    auto near_gaussian = SymbolicExpr::exp(SymbolicExpr::multiply(
        num(-1), SymbolicExpr::power(x, near_two)));
    auto gaussian = SymbolicExpr::exp(SymbolicExpr::multiply(
        num(-1), SymbolicExpr::power(x, num(2))));

    auto result = convolve_checked(near_gaussian, gaussian, "x");
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::Inconclusive)
        << "an exact exponent 2+10^-13 is not a Gaussian square";
}

TEST(TransformEngine, ConvolveNullInputs) {
    auto first = LMCAS::convolve_checked(nullptr, num(1), "x");
    EXPECT_TRUE((!first &&
                 first.error().code == LMCAS::CasErrc::InvalidArgument))
        << "null first convolution input is InvalidArgument";
    auto second = LMCAS::convolve_checked(num(1), nullptr, "x");
    EXPECT_TRUE((!second &&
                 second.error().code == LMCAS::CasErrc::InvalidArgument))
        << "null second convolution input is InvalidArgument";
}

TEST(TransformEngine, LaplaceConstant) {
    auto checked = LMCAS::laplace_transform_checked(num(5), "t", "s");
    ASSERT_TRUE(checked);
    auto result = checked.value().value.expression;
    ASSERT_TRUE(result);
    // ℒ{5} = 5/s ; at s=5 -> 1
    EXPECT_NEAR(eval_at(result, "s", 5.0), 1.0, 1e-6) << "L{5}(s=5) = 5/5 = 1";
}

TEST(TransformEngine, LaplaceExp) {
    auto t = var("t");
    auto a = var("a");
    auto at = SymbolicExpr::multiply(a, t);
    auto f = SymbolicExpr::exp(at);

    auto checked = LMCAS::laplace_transform_checked(f, "t", "s");
    ASSERT_TRUE(checked);
    auto result = checked.value().value.expression;
    ASSERT_TRUE(result);
    // ℒ{e^(at)} = 1/(s-a) ; at a=1, s=3 -> 1/2
    EXPECT_NEAR(eval_at2(result, "a", 1.0, "s", 3.0), 0.5, 1e-6) << "L{e^t}(s=3) = 1/(3-1) = 0.5";
}

TEST(TransformEngine, LaplaceSin) {
    auto t = var("t");
    auto a = var("a");
    auto at = SymbolicExpr::multiply(a, t);
    auto f = SymbolicExpr::sin(at);

    auto checked = LMCAS::laplace_transform_checked(f, "t", "s");
    ASSERT_TRUE(checked);
    auto result = checked.value().value.expression;
    ASSERT_TRUE(result);
    // ℒ{sin(at)} = a/(s²+a²) ; at a=2, s=1 -> 2/5 = 0.4
    EXPECT_NEAR(eval_at2(result, "a", 2.0, "s", 1.0), 0.4, 1e-6) << "L{sin(2t)}(s=1) = 2/(1+4) = 0.4";
}

TEST(TransformEngine, LaplaceTPower) {
    auto t = var("t");
    auto t_sq = SymbolicExpr::power(t, num(2));

    auto checked = LMCAS::laplace_transform_checked(t_sq, "t", "s");
    ASSERT_TRUE(checked);
    auto result = checked.value().value.expression;
    ASSERT_TRUE(result);
    // at s=2 -> 2/8 = 0.25
    EXPECT_NEAR(eval_at(result, "s", 2.0), 0.25, 1e-6) << "L{t^2}(s=2) = 2/8 = 0.25";
}

TEST(TransformEngine, ZTransformConstant) {
    auto checked = LMCAS::z_transform_checked(num(3), "n", "z");
    ASSERT_TRUE(checked);
    auto result = checked.value().value.expression;
    ASSERT_TRUE(result);
    // Z{3} = 3z/(z-1) ; at z=2 -> 6
    EXPECT_NEAR(eval_at(result, "z", 2.0), 6.0, 1e-6) << "Z{3}(z=2) = 3*2/(2-1) = 6";
}

TEST(TransformEngine, ZTransformExpSequence) {
    auto n = var("n");
    auto a = var("a");
    auto f = SymbolicExpr::power(a, n);

    auto checked = LMCAS::z_transform_checked(f, "n", "z");
    ASSERT_TRUE(checked);
    auto result = checked.value().value.expression;
    ASSERT_TRUE(result);
    // Z{a^n} = z/(z-a) ; at a=2, z=4 -> 4/2 = 2
    EXPECT_NEAR(eval_at2(result, "a", 2.0, "z", 4.0), 2.0, 1e-6) << "Z{2^n}(z=4) = 4/(4-2) = 2";
}
