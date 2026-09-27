#include "complex_analysis.hpp"
#include "limit_result.hpp"
#include "internal/symbolic_ast.hpp"
#include "test_common.hpp"

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
static std::shared_ptr<SymbolicExpr> cmplx(int a, int b) {
    return LMCAS::detail::make_expression_ptr(
        SymbolicFactory::create_complex(
            LMCAS::detail::node(num(a)), LMCAS::detail::node(num(b))));
}

TEST(LmcasComplexAnalysis, SimplePole) {
    // ---- residue of 1/(z-a) at z=a is 1 ----
    auto z = SymbolicExpr::variable("z");
    auto a = num(2);
    auto denom = SymbolicExpr::add(z, SymbolicExpr::multiply(num(-1), a));
    auto f = SymbolicExpr::divide(num(1), denom);
    auto res = residue_checked(f, "z", a, 1).value();
    ASSERT_TRUE((res != nullptr)) << ("residue not null");
    {
        const auto actual_expr = res->simplify();
        const auto expected_expr = num(1);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Res(1/(z-2), z=2) = 1");
    }
}

TEST(LmcasComplexAnalysis, ResidueAtNonrealPole) {
    auto z = SymbolicExpr::variable("z");
    auto i = cmplx(0, 1);
    auto shifted = SymbolicExpr::add(z, SymbolicExpr::multiply(num(-1), i));
    auto result = residue_checked(SymbolicExpr::divide(num(1), shifted), "z", i, 1);
    ASSERT_TRUE(result.has_value()) << "residue at a nonreal pole succeeds";
    ASSERT_TRUE(result.value());
    auto simplified = result.value()->simplify();
    ASSERT_TRUE(simplified);
    EXPECT_EQ(simplified->to_string(), num(1)->to_string()) << "Res(1/(z-i), i) = 1";
}

TEST(LmcasComplexAnalysis, ImaginaryResidue) {
    auto z = SymbolicExpr::variable("z");
    auto i = cmplx(0, 1);
    auto result = residue_checked(SymbolicExpr::divide(i, z), "z", num(0), 1);
    ASSERT_TRUE(result.has_value()) << "nonreal residue succeeds";
    ASSERT_TRUE(result.value());
    auto simplified = result.value()->simplify();
    ASSERT_TRUE(simplified);
    EXPECT_EQ(simplified->to_string(), i->to_string()) << "Res(i/z, 0) = i";
}

TEST(LmcasComplexAnalysis, LogarithmicResidueAwayFromBranchCut) {
    auto z = SymbolicExpr::variable("z");
    auto result = residue_checked(SymbolicExpr::divide(
                                      SymbolicExpr::ln(z), SymbolicExpr::add(z, num(-2))),
                                  "z", num(2), 1);
    ASSERT_TRUE(result.has_value()) << "logarithmic residue away from branch cut succeeds";
    ASSERT_TRUE(result.value());
    auto simplified = result.value()->simplify();
    auto expected = SymbolicExpr::ln(num(2))->simplify();
    ASSERT_TRUE(simplified);
    ASSERT_TRUE(expected);
    EXPECT_EQ(simplified->to_string(), expected->to_string()) << "Res(ln(z)/(z-2), 2) = ln(2)";
}

TEST(LmcasComplexAnalysis, RadicalResidueAwayFromBranchCut) {
    auto z = SymbolicExpr::variable("z");
    auto result = residue_checked(SymbolicExpr::divide(
                                      SymbolicExpr::sqrt(z), SymbolicExpr::add(z, num(-4))),
                                  "z", num(4), 1);
    ASSERT_TRUE(result.has_value()) << "radical residue away from branch cut succeeds";
    ASSERT_TRUE(result.value());
    auto simplified = result.value()->simplify();
    ASSERT_TRUE(simplified);
    EXPECT_EQ(simplified->to_string(), num(2)->to_string()) << "Res(sqrt(z)/(z-4), 4) = 2";
}

TEST(LmcasComplexAnalysis, ComplexIdentityLimit) {
    auto z = SymbolicExpr::variable("z");
    auto i = cmplx(0, 1);
    auto result = limit_expression_checked(z, "z", i, LimitDirection::Both, Domain::Complex);
    ASSERT_TRUE(result.has_value()) << "complex identity limit succeeds at i";
    ASSERT_TRUE(result.value());
    EXPECT_EQ(result.value()->to_string(), i->to_string()) << "lim z at i = i";
}

TEST(LmcasComplexAnalysis, ComplexVaryingLimit) {
    auto z = SymbolicExpr::variable("z");
    auto i = cmplx(0, 1);
    auto product = SymbolicExpr::multiply(i, z);
    auto result = limit_expression_checked(product, "z", num(1), LimitDirection::Both, Domain::Complex);
    ASSERT_TRUE(result.has_value()) << "complex-valued varying limit succeeds";
    ASSERT_TRUE(result.value());
    auto simplified = result.value()->simplify();
    ASSERT_TRUE(simplified);
    EXPECT_EQ(simplified->to_string(), i->to_string()) << "lim i*z at 1 = i";
}

TEST(LmcasComplexAnalysis, ComplexConstantLimit) {
    auto i = cmplx(0, 1);
    auto result = limit_expression_checked(i, "z", num(0), LimitDirection::Both, Domain::Complex);
    ASSERT_TRUE(result.has_value()) << "finite complex constant is preserved";
    ASSERT_TRUE(result.value());
    EXPECT_EQ(result.value()->to_string(), i->to_string()) << "lim i = i";
}

TEST(LmcasComplexAnalysis, ComplexRemovablePoleLimit) {
    auto z = SymbolicExpr::variable("z");
    auto i = cmplx(0, 1);
    auto result = limit_expression_checked(
        SymbolicExpr::divide(SymbolicExpr::multiply(i, z), z),
        "z", num(0), LimitDirection::Both, Domain::Complex);
    ASSERT_TRUE(result.has_value()) << "complex removable pole has a finite limit";
    ASSERT_TRUE(result.value());
    auto simplified = result.value()->simplify();
    ASSERT_TRUE(simplified);
    EXPECT_EQ(simplified->to_string(), i->to_string()) << "lim (i*z)/z at 0 = i";
}

TEST(LmcasComplexAnalysis, ComplexPuncturedNormalization) {
    auto z = SymbolicExpr::variable("z");
    auto i = cmplx(0, 1);
    auto shifted = SymbolicExpr::add(z, SymbolicExpr::multiply(num(-1), i));
    auto zero_product = SymbolicExpr::multiply(num(0), SymbolicExpr::power(shifted, num(-1)));
    auto zero_limit = limit_expression_checked(zero_product, "z", i, LimitDirection::Both, Domain::Complex);
    ASSERT_TRUE((zero_limit.has_value())) << ("zero times a punctured pole has a complex limit");
    if (zero_limit) {
        {
            const auto actual_expr = zero_limit.value();
            const auto expected_expr = num(0);
            EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("punctured zero product tends to zero");
        }
    }
    auto zero_power = limit_expression_checked(SymbolicExpr::power(shifted, num(0)),
                                               "z", i, LimitDirection::Both, Domain::Complex);
    ASSERT_TRUE((zero_power.has_value())) << ("zero power is defined throughout the punctured neighborhood");
    if (zero_power) {
        {
            const auto actual_expr = zero_power.value();
            const auto expected_expr = num(1);
            EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("punctured zero power tends to one");
        }
    }
}

TEST(LmcasComplexAnalysis, ComplexLimitDomainProtection) {
    using Node = std::shared_ptr<const SymbolicNode>;
    auto z = SymbolicExpr::variable("z");
    auto zero = num(0);
    auto pole = detail::make_node<PowerNode>(detail::node(zero), detail::node(num(-1)));
    auto masked_pole = detail::make_expression_ptr(detail::make_node<MultiplyNode>(
        std::vector<Node>{detail::node(zero), pole}));
    auto invalid = limit_checked(masked_pole, "z", zero, LimitDirection::Both, Domain::Complex);
    EXPECT_TRUE((!invalid && invalid.error().code == CasErrc::Inconclusive)) << ("complex continuity cannot normalize away an undefined source");
    auto a = SymbolicExpr::variable("a");
    auto denominator = SymbolicExpr::multiply(a, z);
    auto uncertain = limit_checked(SymbolicExpr::divide(denominator, denominator),
                                   "z", zero, LimitDirection::Both, Domain::Complex);
    EXPECT_TRUE((!uncertain && uncertain.error().code == CasErrc::Inconclusive)) << ("unknown polynomial coefficient is not proof of a punctured domain");
    auto logarithm = detail::make_expression_ptr(detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Ln, std::vector<Node>{detail::node(z)}));
    auto branch_cut = limit_checked(logarithm, "z", num(-1), LimitDirection::Both, Domain::Complex);
    EXPECT_TRUE((!branch_cut && branch_cut.error().code == CasErrc::Inconclusive)) << ("complex definedness does not prove continuity across a branch cut");
    auto ordered = limit_checked(z, "z", cmplx(0, 1), LimitDirection::FromAbove, Domain::Complex);
    EXPECT_TRUE((!ordered && ordered.error().code == CasErrc::InvalidArgument)) << ("complex domain does not invent ordered one-sided approaches");
    auto essential = SymbolicExpr::exp(
        SymbolicExpr::multiply(num(-1), SymbolicExpr::power(z, num(-2))));
    auto directional_only = limit_checked(essential, "z", zero, LimitDirection::Both, Domain::Complex);
    EXPECT_TRUE((!directional_only && directional_only.error().code == CasErrc::Inconclusive)) << ("a finite limit along the real axis does not certify a complex limit");
}

TEST(LmcasComplexAnalysis, ComplexParts) {
    // ---- real_part / imag_part of (3 + 4i) ----
    auto z = cmplx(3, 4);
    {
        const auto actual_expr = real_part_checked(z).value();
        const auto expected_expr = num(3);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Re(3+4i) = 3");
    }
    {
        const auto actual_expr = imag_part_checked(z).value();
        const auto expected_expr = num(4);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Im(3+4i) = 4");
    }

    auto re_checked = real_part_checked(z);
    auto im_checked = imag_part_checked(z);
    ASSERT_TRUE(re_checked.has_value()) << "checked Re(3+4i) succeeds";
    ASSERT_TRUE(im_checked.has_value()) << "checked Im(3+4i) succeeds";
    EXPECT_TRUE(test_expression_text(re_checked.value(), num(3)))
        << "checked Re(3+4i) = 3";
    EXPECT_TRUE(test_expression_text(im_checked.value(), num(4)))
        << "checked Im(3+4i) = 4";
}

TEST(LmcasComplexAnalysis, ProductParts) {
    // ---- real/imag of (2+i)*(1+i) = 2 + 2i + i + i^2 = 1 + 3i ----
    auto prod = SymbolicExpr::multiply(cmplx(2, 1), cmplx(1, 1));
    {
        const auto actual_expr = real_part_checked(prod).value();
        const auto expected_expr = num(1);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Re((2+i)(1+i)) = 1");
    }
    {
        const auto actual_expr = imag_part_checked(prod).value();
        const auto expected_expr = num(3);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Im((2+i)(1+i)) = 3");
    }
}

TEST(LmcasComplexAnalysis, Conjugate) {
    // ---- conjugate of (3+4i) = 3-4i ----
    auto z = cmplx(3, 4);
    auto c = conjugate_checked(z).value();
    {
        const auto actual_expr = real_part_checked(c).value();
        const auto expected_expr = num(3);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Re(conj(3+4i)) = 3");
    }
    {
        const auto actual_expr = imag_part_checked(c).value()->simplify();
        const auto expected_expr = num(-4);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Im(conj(3+4i)) = -4");
    }

    auto checked = conjugate_checked(z);
    ASSERT_TRUE((checked.has_value())) << ("checked conjugate(3+4i) succeeds");
    if (checked) {
        {
            const auto actual_expr = real_part_checked(checked.value()).value();
            const auto expected_expr = num(3);
            EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("checked Re(conj(3+4i)) = 3");
        }
        {
            const auto actual_expr = imag_part_checked(checked.value()).value()->simplify();
            const auto expected_expr = num(-4);
            EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("checked Im(conj(3+4i)) = -4");
        }
    }
}

TEST(LmcasComplexAnalysis, InvalidComplexParts) {
    // ---- checked complex part APIs reject invalid expressions ----
    auto null_re = real_part_checked(nullptr);
    auto null_im = imag_part_checked(nullptr);
    auto null_conj = conjugate_checked(nullptr);
    EXPECT_TRUE((!null_re && null_re.error().code == CasErrc::InvalidArgument)) << ("checked real_part rejects null expression");
    EXPECT_TRUE((!null_im && null_im.error().code == CasErrc::InvalidArgument)) << ("checked imag_part rejects null expression");
    EXPECT_TRUE((!null_conj && null_conj.error().code == CasErrc::InvalidArgument)) << ("checked conjugate rejects null expression");

    std::shared_ptr<SymbolicExpr> empty_expr;
    auto empty_re = real_part_checked(empty_expr);
    EXPECT_TRUE((!empty_re && empty_re.error().code == CasErrc::InvalidArgument)) << ("checked real_part rejects null expression");
}

TEST(LmcasComplexAnalysis, CancelledComplexParts) {
    // ---- checked complex part APIs observe computation context cancellation ----
    CancellationToken token;
    token.cancel();
    ComputationContext context({}, token);
    auto cancelled = real_part_checked(cmplx(1, 2), context);
    EXPECT_TRUE((!cancelled && cancelled.error().code == CasErrc::Cancelled)) << ("checked real_part observes cancelled context");
}

TEST(LmcasComplexAnalysis, ComplexFunctionDomain) {
    /// 受检实部/虚部接口对复合函数参数返回明确诊断.
    auto i = cmplx(0, 1);
    auto exp_i = SymbolicExpr::exp(i);
    auto re = real_part_checked(exp_i);
    auto im = imag_part_checked(exp_i);
    auto conj = conjugate_checked(exp_i);
    EXPECT_TRUE((!re && re.error().code == CasErrc::Inconclusive)) << ("checked real_part rejects unsupported function of complex argument");
    EXPECT_TRUE((!im && im.error().code == CasErrc::Inconclusive)) << ("checked imag_part rejects unsupported function of complex argument");
    EXPECT_TRUE((!conj && conj.error().code == CasErrc::Inconclusive)) << ("checked conjugate rejects unsupported function of complex argument");
}

TEST(LmcasComplexAnalysis, UninterpretedFunctionAnalyticityUsesFreeVariables) {
    using Node = std::shared_ptr<const SymbolicNode>;
    auto z = SymbolicExpr::variable("z");
    auto w = SymbolicExpr::variable("w");
    auto f_of_z = detail::make_expression_ptr(
        detail::make_node<UninterpretedFunctionNode>(
            "f", std::vector<Node>{detail::node(z)}));
    auto f_of_w = detail::make_expression_ptr(
        detail::make_node<UninterpretedFunctionNode>(
            "f", std::vector<Node>{detail::node(w)}));

    auto dependent = is_analytic_checked(f_of_z, "z");
    ASSERT_FALSE(dependent);
    EXPECT_EQ(dependent.error().code, CasErrc::Inconclusive);
    EXPECT_EQ(dependent.error().operation, "is_analytic");

    auto constant = is_analytic_checked(f_of_w, "z");
    ASSERT_TRUE(constant);
    EXPECT_TRUE(constant.value());
}

TEST(LmcasComplexAnalysis, AnalyticityRespectsBinderShadowing) {
    using Node = std::shared_ptr<const SymbolicNode>;
    auto z = SymbolicExpr::variable("z");
    auto f_of_bound_z = detail::make_node<UninterpretedFunctionNode>(
        "f", std::vector<Node>{detail::node(z)});
    auto integral = detail::make_expression_ptr(
        detail::make_node<IntegralNode>(f_of_bound_z, "z"));

    auto result = is_analytic_checked(integral, "z");
    ASSERT_TRUE(result);
    EXPECT_TRUE(result.value());
}

TEST(LmcasComplexAnalysis, NestedContainersExposeComplexFunctionArguments) {
    using Node = std::shared_ptr<const SymbolicNode>;
    auto explicit_complex = detail::node(cmplx(2, 3));
    auto integral = detail::make_node<IntegralNode>(explicit_complex, "x");
    auto condition = detail::make_node<RelationalNode>(
        detail::node(SymbolicExpr::variable("x")), detail::node(num(0)),
        RelationOp::GT);
    auto piecewise = detail::make_node<PiecewiseNode>(
        std::vector<PiecewiseNode::Branch>{{integral, condition}},
        detail::node(num(0)));
    auto set = detail::make_node<FiniteSetNode>(std::vector<Node>{piecewise});
    auto matrix = detail::make_node<MatrixNode>(
        1, 1, MatrixNode::DenseStorage{set});
    auto expression = detail::make_expression_ptr(
        detail::make_node<UninterpretedFunctionNode>(
            "f", std::vector<Node>{matrix}));

    auto real = real_part_checked(expression);
    auto imag = imag_part_checked(expression);
    auto conjugate = conjugate_checked(expression);
    ASSERT_FALSE(real);
    ASSERT_FALSE(imag);
    ASSERT_FALSE(conjugate);
    EXPECT_EQ(real.error().code, CasErrc::Inconclusive);
    EXPECT_EQ(real.error().operation, "real_part");
    EXPECT_EQ(imag.error().code, CasErrc::Inconclusive);
    EXPECT_EQ(imag.error().operation, "imag_part");
    EXPECT_EQ(conjugate.error().code, CasErrc::Inconclusive);
    EXPECT_EQ(conjugate.error().operation, "conjugate");
}

TEST(LmcasComplexAnalysis, RawImaginaryUnitHasCanonicalPartsAndConjugate) {
    auto imaginary = SymbolicExpr::variable("I");
    auto real = real_part_checked(imaginary);
    auto imag = imag_part_checked(imaginary);
    auto conjugate = conjugate_checked(imaginary);
    ASSERT_TRUE(real);
    ASSERT_TRUE(imag);
    ASSERT_TRUE(conjugate);
    EXPECT_TRUE(test_same_expression(real.value(), num(0)));
    EXPECT_TRUE(test_same_expression(imag.value(), num(1)));

    auto conjugate_real = real_part_checked(conjugate.value());
    auto conjugate_imag = imag_part_checked(conjugate.value());
    ASSERT_TRUE(conjugate_real);
    ASSERT_TRUE(conjugate_imag);
    EXPECT_TRUE(test_same_expression(conjugate_real.value(), num(0)));
    EXPECT_TRUE(test_same_expression(conjugate_imag.value(), num(-1)));
}

TEST(LmcasComplexAnalysis, LowercaseIIsAnOrdinaryRealSymbol) {
    auto symbol = SymbolicExpr::variable("i");
    auto real = real_part_checked(symbol);
    auto imag = imag_part_checked(symbol);
    auto conjugate = conjugate_checked(symbol);
    ASSERT_TRUE(real);
    ASSERT_TRUE(imag);
    ASSERT_TRUE(conjugate);
    EXPECT_TRUE(test_same_expression(real.value(), symbol));
    EXPECT_TRUE(test_same_expression(imag.value(), num(0)));
    EXPECT_TRUE(test_same_expression(conjugate.value(), symbol));
}

TEST(LmcasComplexAnalysis, ImaginarySquare) {
    // ---- i^2 = -1 ----
    auto i = cmplx(0, 1);
    auto i2 = SymbolicExpr::multiply(i, i)->simplify();
    {
        const auto actual_expr = real_part_checked(i2).value();
        const auto expected_expr = num(-1);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Re(i^2) = -1");
    }
    {
        const auto actual_expr = imag_part_checked(i2).value()->simplify();
        const auto expected_expr = num(0);
        EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << ("Im(i^2) = 0");
    }
}

TEST(LmcasComplexAnalysis, PolynomialAnalyticity) {
    auto z = SymbolicExpr::variable("z");
    auto f = SymbolicExpr::multiply(z, z);
    auto checked = is_analytic_checked(f, "z");
    ASSERT_TRUE(checked.has_value())
        << "checked is_analytic succeeds for z^2";
    EXPECT_TRUE(checked.value()) << "checked z^2 is analytic";
}

TEST(LmcasComplexAnalysis, FunctionAnalyticity) {
    /// 受检解析性接口显式报告当前支持域之外的依赖函数.
    auto z = SymbolicExpr::variable("z");
    auto sin_z = SymbolicExpr::sin(z);
    auto checked = is_analytic_checked(sin_z, "z");
    EXPECT_TRUE((!checked && checked.error().code == CasErrc::Inconclusive)) << ("checked is_analytic reports dependent function domain as inconclusive");

    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto constant_in_z = is_analytic_checked(sin_x, "z");
    EXPECT_TRUE((constant_in_z.has_value() && constant_in_z.value())) << ("checked is_analytic treats functions independent of z as constants");
}

TEST(LmcasComplexAnalysis, ResidueAndCauchyDomains) {
    // ---- checked residue and Cauchy APIs reject invalid domains explicitly ----
    auto z = SymbolicExpr::variable("z");
    auto z0 = num(2);
    auto denom = SymbolicExpr::add(z, SymbolicExpr::multiply(num(-1), z0));
    auto f = SymbolicExpr::divide(num(1), denom);

    auto checked_residue = residue_checked(f, "z", z0, 1);
    ASSERT_TRUE(checked_residue.has_value())
        << "checked residue succeeds for simple pole";
    ASSERT_TRUE(checked_residue.value());
    EXPECT_TRUE(test_expression_text(
        checked_residue.value()->simplify(), num(1)))
        << "checked Res(1/(z-2), z=2) = 1";

    auto bad_order = residue_checked(f, "z", z0, 0);
    EXPECT_TRUE((!bad_order && bad_order.error().code == CasErrc::InvalidArgument)) << ("checked residue rejects order 0 instead of returning zero");

    auto empty_var = cauchy_integral_checked(f, "", z0, 1);
    EXPECT_TRUE((!empty_var && empty_var.error().code == CasErrc::InvalidArgument)) << ("checked Cauchy integral rejects empty variable");

    auto cauchy = cauchy_integral_checked(num(1), "z", num(0), 1);
    ASSERT_TRUE(cauchy.has_value())
        << "checked Cauchy integral constructs constant analytic formula";
    ASSERT_TRUE(cauchy.value());
}

TEST(LmcasComplexAnalysis, AnalysisContextErrors) {
    auto z = SymbolicExpr::variable("z");
    auto z0 = num(0);

    CancellationToken token;
    token.cancel();
    ComputationContext cancelled_context({}, token);
    auto cancelled = residue_checked(z, "z", z0, 1, cancelled_context);
    EXPECT_TRUE((!cancelled && cancelled.error().code == CasErrc::Cancelled)) << ("checked residue observes cancellation");

    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext limited_context(limits);
    auto limited = cauchy_integral_checked(z, "z", z0, 1, limited_context);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << ("checked Cauchy integral observes exhausted step budget");

    auto invalid_analytic = is_analytic_checked(nullptr, "z");
    EXPECT_TRUE((!invalid_analytic &&
                 invalid_analytic.error().code == CasErrc::InvalidArgument))
        << ("checked is_analytic rejects null expression");
}
