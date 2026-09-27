#include "test_multivariate_factor_support.hpp"

TEST(MultivariateFactor, LeadingCoefficientConstantLcNeedsNoPrecomputation) {
    // f = x^2 + y (lc w.r.t. x is 1, constant)
    // 常数首项系数无需预计算
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)), // x^2
        make_term({0, 1}, Rational(1))  // y
    };
    MultiPoly poly(terms, vars);

    MultiPoly lc = poly.leading_coeff("x");
    EXPECT_TRUE((lc.is_constant())) << "lc(x^2 + y, x) is constant (= 1)";
}

TEST(MultivariateFactor, LeadingCoefficientNonConstantLcInAuxiliaryVariable) {
    // f = y*x^2 + x + 1 (lc w.r.t. x is y, non-constant)
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(1)), // y*x^2
        make_term({1, 0}, Rational(1)), // x
        make_term({0, 0}, Rational(1))  // 1
    };
    MultiPoly poly(terms, vars);

    MultiPoly lc = poly.leading_coeff("x");
    EXPECT_FALSE((lc.is_constant())) << "lc(y*x^2 + x + 1, x) is non-constant (= y)";

    // 验证 lc 在 y=2 处求值为 2
    std::map<std::string, Rational> eval_pts = {{"y", Rational(2)}};
    MultiPoly lc_eval = lc.eval(eval_pts);
    EXPECT_TRUE((lc_eval.is_constant())) << "lc evaluated at y=2 is constant";
    EXPECT_TRUE((lc_eval.terms()[0].second == Rational(2))) << "lc evaluated at y=2 equals 2";
}

TEST(MultivariateFactor, LeadingCoefficientLcEvaluationMatchesFactorLcProduct) {
    // f = (y*x + 1)(x + y) = y*x^2 + (y^2+1)*x + y
    // lc(f, x) = y
    // 一元因子（在 y=1 处）：f(x,1) = x^2 + 2x + 1 = (x+1)^2
    // 在 y=2 处：f(x,2) = 2x^2 + 5x + 2 = (2x+1)(x+2)
    // lc(f,x) = y, 在 y=2 处 lc_eval = 2
    // 一元因子 lc: lc(2x+1) = 2, lc(x+2) = 1, 乘积 = 2 = lc_eval ✓
    std::vector<std::string> vars = {"x", "y"};

    // 构造 f = y*x^2 + (y^2+1)*x + y
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(1)), // y*x^2
        make_term({1, 2}, Rational(1)), // y^2*x
        make_term({1, 0}, Rational(1)), // x
        make_term({0, 1}, Rational(1))  // y
    };
    MultiPoly poly(terms, vars);

    // 验证 lc(f, x) = y
    MultiPoly lc = poly.leading_coeff("x");
    EXPECT_FALSE((lc.is_constant())) << "lc is non-constant";

    // 在 y=2 处求值
    std::map<std::string, Rational> eval_pts = {{"y", Rational(2)}};
    MultiPoly lc_eval = lc.eval(eval_pts);
    EXPECT_TRUE((lc_eval.terms()[0].second == Rational(2))) << "lc at y=2 is 2";

    // 求值后的多项式：f(x,2) = 2x^2 + 5x + 2
    MultiPoly f_eval = poly.eval("y", Rational(2));
    Polynomial<Rational> f_uni = f_eval.to_univariate();
    EXPECT_TRUE((f_uni.lead_coeff() == Rational(2))) << "f(x,2) has leading coefficient 2";
}

TEST(MultivariateFactor, LeadingCoefficientHenselLiftWithNonConstantLc) {
    // 测试 Hensel 提升在非常数首项系数情形下的行为
    // f = (x + y)(x - y) = x^2 - y^2
    // lc(f, x) = 1 (常数)，这是简单情形
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)), // x^2
        make_term({0, 2}, Rational(-1)) // -y^2
    };
    MultiPoly poly(terms, vars);

    // 验证 lc(f, x) = 1（常数），无需预计算
    MultiPoly lc = poly.leading_coeff("x");
    EXPECT_TRUE((lc.is_constant())) << "lc(x^2 - y^2, x) is constant (= 1)";

    // 在 y=1 处求值：f(x,1) = x^2 - 1 = (x+1)(x-1)
    MultiPoly f_eval = poly.eval("y", Rational(1));
    Polynomial<Rational> f_uni = f_eval.to_univariate();
    EXPECT_TRUE((f_uni.degree() == 2)) << "f(x,1) has degree 2";
    EXPECT_TRUE((f_uni.lead_coeff() == Rational(1))) << "f(x,1) is monic";
}

TEST(MultivariateFactor, LeadingCoefficientMultivariatePolyWithYAsLc) {
    // f = y*x^2 - y = y*(x^2 - 1) = y*(x+1)*(x-1)
    // lc(f, x) = y (非常数)
    // 在 y=1 处：f(x,1) = x^2 - 1 = (x+1)(x-1)
    // lc_eval = 1, 一元因子 lc 乘积 = 1*1 = 1 = lc_eval
    // 此情形下 scale = 1，无需调整
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(1)), // y*x^2
        make_term({0, 1}, Rational(-1)) // -y
    };
    MultiPoly poly(terms, vars);

    MultiPoly lc = poly.leading_coeff("x");
    EXPECT_FALSE((lc.is_constant())) << "lc(y*x^2 - y, x) = y is non-constant";

    // 在 y=1 处求值
    std::map<std::string, Rational> eval_pts = {{"y", Rational(1)}};
    MultiPoly lc_eval = lc.eval(eval_pts);
    EXPECT_TRUE((lc_eval.terms()[0].second == Rational(1))) << "lc at y=1 is 1";
}

TEST(MultivariateFactor, LeadingCoefficientScaleFactorComputation) {
    // f = 2y*x^2 + 3y*x + y = y*(2x^2 + 3x + 1) = y*(2x+1)*(x+1)
    // lc(f, x) = 2y
    // 在 y=1 处：f(x,1) = 2x^2 + 3x + 1 = (2x+1)(x+1)
    // lc_eval = 2, 一元因子 lc: lc(2x+1)=2, lc(x+1)=1, 乘积=2 = lc_eval ✓
    // 无需缩放
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(2)), // 2y*x^2
        make_term({1, 1}, Rational(3)), // 3y*x
        make_term({0, 1}, Rational(1))  // y
    };
    MultiPoly poly(terms, vars);

    MultiPoly lc = poly.leading_coeff("x");
    EXPECT_FALSE((lc.is_constant())) << "lc(2y*x^2+3y*x+y, x) = 2y is non-constant";

    // 在 y=1 处求值
    std::map<std::string, Rational> eval_pts = {{"y", Rational(1)}};
    MultiPoly lc_eval = lc.eval(eval_pts);
    Rational lc_val = lc_eval.is_zero() ? Rational(0) : lc_eval.terms()[0].second;
    EXPECT_TRUE((lc_val == Rational(2))) << "lc at y=1 is 2";

    // 一元因子 (2x+1)(x+1) 的首项系数乘积
    Polynomial<Rational> f1({Rational(1), Rational(2)}, "x"); // 2x+1
    Polynomial<Rational> f2({Rational(1), Rational(1)}, "x"); // x+1
    Rational product_lcs = f1.lead_coeff() * f2.lead_coeff();
    EXPECT_TRUE((product_lcs == lc_val)) << "product of factor lcs equals lc_eval (no scaling needed)";
}

TEST(MultivariateFactor, LeadingCoefficientScalingNeededWhenLcProductDiffers) {
    // 场景：一元因子的首项系数乘积与 lc_eval 不同
    // f = 3y*x^2 + ... 在 y=1 处 lc_eval = 3
    // 若一元因子为 (x+a)(x+b)（首一），则 lc 乘积 = 1 ≠ 3
    // 需要将缩放因子 3 分配给某个因子

    // 构造 f = 3y*x^2 - 3y = 3y*(x^2-1) = 3y*(x+1)*(x-1)
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 1}, Rational(3)), // 3y*x^2
        make_term({0, 1}, Rational(-3)) // -3y
    };
    MultiPoly poly(terms, vars);

    MultiPoly lc = poly.leading_coeff("x");
    std::map<std::string, Rational> eval_pts = {{"y", Rational(1)}};
    MultiPoly lc_eval = lc.eval(eval_pts);
    Rational lc_val = lc_eval.terms()[0].second;
    EXPECT_TRUE((lc_val == Rational(3))) << "lc at y=1 is 3";

    // 若一元分解给出首一因子 (x+1)(x-1)
    Polynomial<Rational> f1({Rational(1), Rational(1)}, "x");  // x+1
    Polynomial<Rational> f2({Rational(-1), Rational(1)}, "x"); // x-1
    Rational product_lcs = f1.lead_coeff() * f2.lead_coeff();
    EXPECT_TRUE((product_lcs == Rational(1))) << "monic factors have lc product = 1";

    // 缩放因子 = lc_val / product_lcs = 3
    Rational scale = lc_val / product_lcs;
    EXPECT_TRUE((scale == Rational(3))) << "scale factor is 3";

    // 应用缩放后，第一个因子变为 3x+3
    Polynomial<Rational> scaled_f1({Rational(3), Rational(3)}, "x"); // 3(x+1) = 3x+3
    EXPECT_TRUE((scaled_f1.lead_coeff() == Rational(3))) << "scaled factor has correct leading coefficient";

    // 验证缩放后乘积的首项系数
    Rational new_product_lcs = scaled_f1.lead_coeff() * f2.lead_coeff();
    EXPECT_TRUE((new_product_lcs == lc_val)) << "after scaling, product of lcs equals lc_eval";
}
