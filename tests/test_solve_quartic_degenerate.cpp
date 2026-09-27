#include "test_solve_quartic_support.hpp"

TEST(SolveQuarticDegenerate, SolveBiquadraticDegenerateLeadingCoefficient) {
    auto roots = LMCAS::solve_biquadratic(
        SymbolicExpr::number(0),
        SymbolicExpr::number(1),
        SymbolicExpr::number(-4),
        "x");
    EXPECT_TRUE((roots.size() == 2)) << "0*x^4+x^2-4 delegates to a quadratic in x";
    if (roots.size() == 2) {
        std::vector<double> values{
            roots[0]->to_numeric(), roots[1]->to_numeric()};
        std::sort(values.begin(), values.end());
        EXPECT_TRUE((std::abs(values[0] + 2.0) < 1e-12 &&
                     std::abs(values[1] - 2.0) < 1e-12))
            << "degenerate biquadratic roots are -2 and 2";
    }

    auto no_roots = LMCAS::solve_biquadratic(
        SymbolicExpr::number(0),
        SymbolicExpr::number(0),
        SymbolicExpr::number(1),
        "x");
    EXPECT_TRUE((no_roots.empty())) << "nonzero constant equation has no roots";

    bool rejected_indeterminate = false;
    try {
        (void)LMCAS::solve_biquadratic(
            SymbolicExpr::number(0),
            SymbolicExpr::number(0),
            SymbolicExpr::number(0),
            "x");
    } catch (const std::invalid_argument &) {
        rejected_indeterminate = true;
    }
    EXPECT_TRUE((rejected_indeterminate)) << "identically zero equation is rejected as indeterminate";
}

TEST(SolveQuarticDegenerate, SolveQuarticA0DelegatesToCubic) {
    auto a = SymbolicExpr::number(0);
    auto b = SymbolicExpr::number(1);
    auto c = SymbolicExpr::number(-6);
    auto d = SymbolicExpr::number(11);
    auto e = SymbolicExpr::number(-6);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 3)) << "Quartic with a=0 should delegate to cubic and return 3 roots";
}

TEST(SolveQuarticDegenerate, SolveQuarticBiquadraticShortcutVerifiesRootValues12) {
    auto a = SymbolicExpr::number(1);
    auto b = SymbolicExpr::number(0);
    auto c = SymbolicExpr::number(-5);
    auto d = SymbolicExpr::number(0);
    auto e = SymbolicExpr::number(4);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 4)) << "Biquadratic x^4-5x^2+4 should return 4 roots";

    std::vector<double> vals;
    for (auto &r : roots)
        vals.push_back(r->to_numeric());
    std::sort(vals.begin(), vals.end());

    double expected[] = {-2.0, -1.0, 1.0, 2.0};
    bool all_match = true;
    for (int i = 0; i < 4; ++i) {
        if (std::abs(vals[i] - expected[i]) > 1e-8) {
            all_match = false;
            break;
        }
    }
    EXPECT_TRUE((all_match)) << "Biquadratic roots should be {-2, -1, 1, 2} (got {" + std::to_string(vals[0]) + ", " + std::to_string(vals[1]) + ", " + std::to_string(vals[2]) + ", " + std::to_string(vals[3]) + "})";
}

TEST(SolveQuarticDegenerate, SolveBiquadraticExactHugeCoefficientAvoidsUnsafeNumericConversion) {
    std::string huge_digits = "1" + std::string(400, '0');
    auto roots = LMCAS::solve_biquadratic(
        SymbolicExpr::number(BigInt(huge_digits)),
        SymbolicExpr::number(0),
        SymbolicExpr::number(-1),
        "x");

    EXPECT_TRUE((roots.size() == 4)) << "huge exact biquadratic should return symbolic square-root roots";

    int zero_roots = 0;
    for (const auto &root : roots) {
        std::string text = root ? root->to_string() : "";
        EXPECT_TRUE((text.find("inf") == std::string::npos &&
                     text.find("nan") == std::string::npos))
            << "huge exact biquadratic roots should not contain fabricated inf/nan";
        if (root && root->simplify()->is_zero()) {
            ++zero_roots;
        }
    }
    EXPECT_TRUE((zero_roots < 4)) << "huge exact biquadratic must not collapse nonzero roots to zero";
}

TEST(SolveQuarticDegenerate, SolveQuarticQ0AfterDepressionNonBiquadratic) {
    auto a = SymbolicExpr::number(1);
    auto b = SymbolicExpr::number(8);
    auto c = SymbolicExpr::number(22);
    auto d = SymbolicExpr::number(24);
    auto e = SymbolicExpr::number(9);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 4)) << "q=0 quartic (x+1)^2(x+3)^2 should return 4 roots";

    for (size_t i = 0; i < roots.size(); ++i) {
        double val = roots[i]->to_numeric();
        double residual = eval_quartic(1, 8, 22, 24, 9, val);
        EXPECT_TRUE((std::abs(residual) < 1e-8)) << "q=0 quartic root " + std::to_string(i) + " residual < 1e-8 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")";
    }

    std::vector<double> vals;
    for (auto &r : roots)
        vals.push_back(r->to_numeric());
    std::sort(vals.begin(), vals.end());

    int count_neg3 = 0, count_neg1 = 0;
    for (double v : vals) {
        if (std::abs(v - (-3.0)) < 1e-8)
            count_neg3++;
        else if (std::abs(v - (-1.0)) < 1e-8)
            count_neg1++;
    }
    EXPECT_TRUE((count_neg3 == 2 && count_neg1 == 2)) << "q=0 quartic roots should be {-3, -3, -1, -1}";
}

TEST(SolveQuarticDegenerate, SolveQuarticRepeatedRootX14) {
    auto a = SymbolicExpr::number(1);
    auto b = SymbolicExpr::number(-4);
    auto c = SymbolicExpr::number(6);
    auto d = SymbolicExpr::number(-4);
    auto e = SymbolicExpr::number(1);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 4)) << "(x-1)^4 should return 4 roots";

    bool all_one = true;
    for (size_t i = 0; i < roots.size(); ++i) {
        double val = roots[i]->to_numeric();
        if (std::abs(val - 1.0) > 1e-6) {
            all_one = false;
        }
        double residual = eval_quartic(1, -4, 6, -4, 1, val);
        EXPECT_TRUE((std::abs(residual) < 1e-6)) << "(x-1)^4 root " + std::to_string(i) + " residual < 1e-6 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")";
    }
    EXPECT_TRUE((all_one)) << "(x-1)^4 all roots should equal 1";
}

TEST(SolveQuarticDegenerate, SolveQuarticRepeatedRootsX12X22) {
    auto a = SymbolicExpr::number(1);
    auto b = SymbolicExpr::number(-6);
    auto c = SymbolicExpr::number(13);
    auto d = SymbolicExpr::number(-12);
    auto e = SymbolicExpr::number(4);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 4)) << "(x-1)^2(x-2)^2 should return 4 roots";

    for (size_t i = 0; i < roots.size(); ++i) {
        double val = roots[i]->to_numeric();
        double residual = eval_quartic(1, -6, 13, -12, 4, val);
        EXPECT_TRUE((std::abs(residual) < 1e-6)) << "(x-1)^2(x-2)^2 root " + std::to_string(i) + " residual < 1e-6 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")";
    }

    std::vector<double> vals;
    for (auto &r : roots)
        vals.push_back(r->to_numeric());
    std::sort(vals.begin(), vals.end());

    int count_1 = 0, count_2 = 0;
    for (double v : vals) {
        if (std::abs(v - 1.0) < 1e-6)
            count_1++;
        else if (std::abs(v - 2.0) < 1e-6)
            count_2++;
    }
    EXPECT_TRUE((count_1 == 2 && count_2 == 2)) << "(x-1)^2(x-2)^2 roots should be {1, 1, 2, 2} (got {" + std::to_string(vals[0]) + ", " + std::to_string(vals[1]) + ", " + std::to_string(vals[2]) + ", " + std::to_string(vals[3]) + "})";
}

TEST(SolveQuarticDegenerate, SolveQuarticVietaSFormulasSumAndProductOfRoots) {
    {
        double ca = 2.0, cb = 0.0, cc = -20.0, cd = 0.0, ce = 18.0;
        auto a = SymbolicExpr::number(ca);
        auto b = SymbolicExpr::number(cb);
        auto c = SymbolicExpr::number(cc);
        auto d = SymbolicExpr::number(cd);
        auto e = SymbolicExpr::number(ce);

        auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE((roots.size() == 4)) << "Vieta quartic should return 4 roots";

        double sum = 0.0;
        double product = 1.0;
        for (auto &r : roots) {
            double val = r->to_numeric();
            sum += val;
            product *= val;
        }

        double expected_sum = -cb / ca;
        double expected_product = ce / ca;

        EXPECT_TRUE((std::abs(sum - expected_sum) < 1e-8)) << "Vieta sum: Σrᵢ = -b/a = " + std::to_string(expected_sum) + " (got " + std::to_string(sum) + ")";
        EXPECT_TRUE((std::abs(product - expected_product) < 1e-8)) << "Vieta product: ∏rᵢ = e/a = " + std::to_string(expected_product) + " (got " + std::to_string(product) + ")";
    }

    {
        double ca = 1.0, cb = -10.0, cc = 35.0, cd = -50.0, ce = 24.0;
        auto a = SymbolicExpr::number(ca);
        auto b = SymbolicExpr::number(cb);
        auto c = SymbolicExpr::number(cc);
        auto d = SymbolicExpr::number(cd);
        auto e = SymbolicExpr::number(ce);

        auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
        EXPECT_TRUE((roots.size() == 4)) << "Vieta quartic (1,2,3,4) should return 4 roots";

        double sum = 0.0;
        double product = 1.0;
        for (auto &r : roots) {
            double val = r->to_numeric();
            sum += val;
            product *= val;
        }

        double expected_sum = -cb / ca;
        double expected_product = ce / ca;

        EXPECT_TRUE((std::abs(sum - expected_sum) < 1e-8)) << "Vieta sum (1+2+3+4=10): got " + std::to_string(sum);
        EXPECT_TRUE((std::abs(product - expected_product) < 1e-8)) << "Vieta product (1*2*3*4=24): got " + std::to_string(product);
    }
}
