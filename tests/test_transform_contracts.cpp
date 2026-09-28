#include "test_common.hpp"
#include "transform_engine.hpp"

using namespace LMCAS;

static auto num(int n) { return SymbolicExpr::number(n); }
static auto var(const std::string &name) { return SymbolicExpr::variable(name); }

TEST(TransformContracts, LaplaceRoc) {
    auto laplace = LMCAS::laplace_transform_checked(num(5), "t", "s");
    ASSERT_TRUE((laplace.has_value())) << "checked Laplace succeeds for constants";
    if (laplace) {
        EXPECT_TRUE((laplace.value().value.expression != nullptr)) << "checked Laplace returns an expression";
        EXPECT_TRUE((laplace.value().value.roc.size() == 1)) << "checked Laplace reports ROC for constant input";
        if (!laplace.value().value.roc.empty()) {
            auto roc = std::dynamic_pointer_cast<const RelationalNode>(
                LMCAS::detail::node(laplace.value().value.roc[0]));
            EXPECT_TRUE((roc != nullptr && roc->op() == RelationalNode::Op::GT)) << "checked Laplace ROC is a greater-than condition";
            auto lhs = roc ? std::dynamic_pointer_cast<const VariableNode>(roc->left()) : nullptr;
            auto rhs = roc ? std::dynamic_pointer_cast<const NumberNode>(roc->right()) : nullptr;
            EXPECT_TRUE((lhs != nullptr && lhs->name() == "s" &&
                         rhs != nullptr && rhs->is_zero()))
                << "checked Laplace constant ROC is s > 0";
        }
    }
}

TEST(TransformContracts, UnsupportedTransforms) {
    auto t = var("t");
    auto unknown = LMCAS::fourier_transform_checked(SymbolicExpr::ln(t), "t", "omega");
    EXPECT_TRUE((!unknown.has_value())) << "checked Fourier rejects unsupported closed forms";
    EXPECT_TRUE((unknown.error().code == LMCAS::CasErrc::Inconclusive)) << "checked Fourier reports Inconclusive for unevaluated transform nodes";

    auto nested_unknown = LMCAS::fourier_transform_checked(
        SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::ln(t)),
        "t", "omega");
    EXPECT_TRUE((!nested_unknown &&
                 nested_unknown.error().code == LMCAS::CasErrc::Inconclusive))
        << "checked Fourier rejects nested unevaluated transform nodes";

    auto unsupported_inverse = LMCAS::inverse_fourier_transform_checked(num(1), "omega", "t");
    EXPECT_TRUE((!unsupported_inverse.has_value())) << "checked inverse Fourier rejects unevaluated transform nodes";
    EXPECT_TRUE((unsupported_inverse.error().code == LMCAS::CasErrc::Inconclusive)) << "checked inverse Fourier reports Inconclusive for unsupported constants";
}

TEST(TransformContracts, InverseFourierExactDecayRoundTrip) {
    auto omega = var("omega");
    auto input = SymbolicExpr::divide(
        num(2), SymbolicExpr::add(num(1), SymbolicExpr::power(omega, num(2))));
    auto result = inverse_fourier_transform_checked(input, "omega", "t");
    ASSERT_TRUE(result);
    EXPECT_TRUE(std::holds_alternative<ExactRoundTripProof>(result.value().certificate));
    auto value = result.value().value.expression;
    auto abs_t = detail::make_expression_ptr(detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{detail::node(var("t"))}));
    auto expected = SymbolicExpr::exp(SymbolicExpr::multiply(num(-1), abs_t));
    EXPECT_TRUE(detail::node(value->simplify())->equals(*detail::node(expected->simplify())));
    EXPECT_TRUE(result.value().value.conditions.empty());

    auto unsupported = inverse_fourier_transform_checked(
        SymbolicExpr::divide(num(3),
            SymbolicExpr::add(num(1), SymbolicExpr::power(omega, num(2)))),
        "omega", "t");
    ASSERT_FALSE(unsupported);
    EXPECT_EQ(unsupported.error().code, CasErrc::Inconclusive);
}

TEST(TransformContracts, InvalidTransformInputs) {
    auto t = var("t");
    auto s = var("s");
    auto null_input = LMCAS::laplace_transform_checked(nullptr, "t", "s");
    EXPECT_TRUE((!null_input.has_value())) << "checked Laplace rejects null input";
    EXPECT_TRUE((null_input.error().code == LMCAS::CasErrc::InvalidArgument)) << "checked Laplace reports InvalidArgument for null input";

    auto empty_var = LMCAS::inverse_laplace_checked(s, "", "t");
    EXPECT_TRUE((!empty_var.has_value())) << "checked inverse Laplace rejects empty variable";
    EXPECT_TRUE((empty_var.error().code == LMCAS::CasErrc::InvalidArgument)) << "checked inverse Laplace reports InvalidArgument for empty variable";

    auto same_var = LMCAS::z_transform_checked(t, "n", "n");
    EXPECT_TRUE((!same_var.has_value())) << "checked Z transform rejects same input/output variable";
    EXPECT_TRUE((same_var.error().code == LMCAS::CasErrc::InvalidArgument)) << "checked Z transform reports InvalidArgument for same variables";
}

static void expect_z_roc_boundary(const EvaluatedTransform &transform, double radius) {
    EXPECT_TRUE((transform.roc.size() == 1)) << "one convergence boundary is reported";
    if (transform.roc.size() != 1)
        return;
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(detail::node(transform.roc[0]));
    EXPECT_TRUE((relation != nullptr)) << "ROC has a relational boundary";
    if (!relation)
        return;
    EXPECT_TRUE((relation->op() == RelationalNode::Op::GT)) << "ROC excludes its boundary";
    auto left = detail::make_expression_ptr(relation->left());
    auto right = detail::make_expression_ptr(relation->right());
    {
        const double actual_value = (right->simplify()->to_numeric());
        const double expected_value = (radius);
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
    EXPECT_NEAR(left->substitute("z", num(-3))->simplify()->to_numeric(), 3.0, 1e-12) << "ROC uses modulus for negative z";
    EXPECT_NEAR(left->substitute("z", num(0))->simplify()->to_numeric(), 0.0, 1e-12) << "origin remains outside the convergence region";
}

TEST(TransformContracts, ZConstantRoc) {
    auto result = z_transform_checked(num(5), "n", "z");
    ASSERT_TRUE((result.has_value())) << "constant Z transform succeeds";
    if (result)
        expect_z_roc_boundary(result.value().value, 1.0);
}

TEST(TransformContracts, ZExponentialRoc) {
    auto result = z_transform_checked(SymbolicExpr::power(num(2), var("n")), "n", "z");
    ASSERT_TRUE((result.has_value())) << "exponential Z transform succeeds";
    if (result)
        expect_z_roc_boundary(result.value().value, 2.0);
}

TEST(TransformContracts, TransformContextErrors) {
    auto t = var("t");
    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = LMCAS::fourier_transform_checked(t, "t", "omega", cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked Fourier observes cancellation";
    EXPECT_TRUE((cancelled.error().code == LMCAS::CasErrc::Cancelled)) << "checked Fourier reports Cancelled";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 0;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::convolve_checked(t, t, "t", limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked convolution observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == LMCAS::CasErrc::ResourceLimit)) << "checked convolution reports ResourceLimit";
}

TEST(TransformContracts, ConvolutionErrors) {
    auto t = var("t");
    auto convolution_null = LMCAS::convolve_checked(nullptr, t, "t");
    EXPECT_TRUE((!convolution_null.has_value())) << "checked convolution rejects null input";
    EXPECT_TRUE((convolution_null.error().code == LMCAS::CasErrc::InvalidArgument)) << "checked convolution reports InvalidArgument for null input";

    auto unsupported_convolution = LMCAS::convolve_checked(SymbolicExpr::ln(t), t, "t");
    EXPECT_TRUE((!unsupported_convolution.has_value())) << "checked convolution rejects unevaluated integral results";
    EXPECT_TRUE((unsupported_convolution.error().code == LMCAS::CasErrc::Inconclusive)) << "checked convolution reports Inconclusive for unsupported integrals";
}
