#include "test_solve_quartic_support.hpp"

TEST(SolveQuartic, SolveQuarticBiquadratic) {
    auto a = SymbolicExpr::number(1);
    auto b = SymbolicExpr::number(0);
    auto c = SymbolicExpr::number(-5);
    auto d = SymbolicExpr::number(0);
    auto e = SymbolicExpr::number(4);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 4)) << "Biquadratic should return 4 roots";

    for (size_t i = 0; i < roots.size(); ++i) {
        double val = roots[i]->to_numeric();
        double residual = eval_quartic(1, 0, -5, 0, 4, val);
        EXPECT_TRUE((std::abs(residual) < 1e-8)) << "Biquadratic root " + std::to_string(i) + " residual < 1e-8 (got " + std::to_string(residual) + ")";
    }
}

TEST(SolveQuartic, SolveQuarticIllConditionedBiquadratic) {
    auto roots = LMCAS::solve_quartic(
        SymbolicExpr::number(1),
        SymbolicExpr::number(0),
        SymbolicExpr::number(BigInt("-10000000000000000")),
        SymbolicExpr::number(0),
        SymbolicExpr::number(1),
        "x");
    EXPECT_TRUE((roots.size() == 4)) << "x^4-10^16*x^2+1 returns four real roots";

    std::vector<double> magnitudes;
    for (const auto &root : roots) {
        magnitudes.push_back(std::abs(root->to_numeric()));
    }
    std::sort(magnitudes.begin(), magnitudes.end());
    if (magnitudes.size() == 4) {
        EXPECT_TRUE((std::abs(magnitudes[0] - 1e-8) < 1e-20 &&
                     std::abs(magnitudes[1] - 1e-8) < 1e-20))
            << "biquadratic preserves both small roots";
        EXPECT_TRUE((std::abs(magnitudes[2] - 1e8) < 1e-6 &&
                     std::abs(magnitudes[3] - 1e8) < 1e-6))
            << "biquadratic preserves both large roots";
    }
}

TEST(SolveQuartic, SolveQuarticOverflowResistantBiquadratic) {
    auto roots = LMCAS::solve_quartic(
        SymbolicExpr::number(1.0),
        SymbolicExpr::number(0.0),
        SymbolicExpr::number(-1.0e200),
        SymbolicExpr::number(0.0),
        SymbolicExpr::number(1.0),
        "x");
    EXPECT_TRUE((roots.size() == 4)) << "x^4-10^200*x^2+1 returns four finite roots";

    std::vector<double> magnitudes;
    for (const auto &root : roots) {
        magnitudes.push_back(std::abs(root->to_numeric()));
    }
    std::sort(magnitudes.begin(), magnitudes.end());
    if (magnitudes.size() == 4) {
        EXPECT_TRUE((std::isfinite(magnitudes[0]) &&
                     std::abs(magnitudes[0] / 1e-100 - 1.0) < 1e-12 &&
                     std::abs(magnitudes[1] / 1e-100 - 1.0) < 1e-12))
            << "scaled biquadratic preserves both finite small roots";
        EXPECT_TRUE((std::isfinite(magnitudes[2]) &&
                     std::abs(magnitudes[2] / 1e100 - 1.0) < 1e-12 &&
                     std::abs(magnitudes[3] / 1e100 - 1.0) < 1e-12))
            << "scaled biquadratic preserves both finite large roots";
    }
}

TEST(SolveQuartic, SolveQuarticGeneralFerrari) {
    auto a = SymbolicExpr::number(1);
    auto b = SymbolicExpr::number(-10);
    auto c = SymbolicExpr::number(35);
    auto d = SymbolicExpr::number(-50);
    auto e = SymbolicExpr::number(24);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 4)) << "General quartic should return 4 roots";

    for (size_t i = 0; i < roots.size(); ++i) {
        double val = roots[i]->to_numeric();
        double residual = eval_quartic(1, -10, 35, -50, 24, val);
        EXPECT_TRUE((std::abs(residual) < 1e-8)) << "General quartic root " + std::to_string(i) + " residual < 1e-8 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")";
    }
}

TEST(SolveQuartic, SolveQuarticDepressedQ0Case) {
    auto a = SymbolicExpr::number(1);
    auto b = SymbolicExpr::number(2);
    auto c = SymbolicExpr::number(-7);
    auto d = SymbolicExpr::number(-8);
    auto e = SymbolicExpr::number(12);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 4)) << "Quartic (x-1)(x+2)(x-2)(x+3) should return 4 roots";

    for (size_t i = 0; i < roots.size(); ++i) {
        double val = roots[i]->to_numeric();
        double residual = eval_quartic(1, 2, -7, -8, 12, val);
        EXPECT_TRUE((std::abs(residual) < 1e-8)) << "Root " + std::to_string(i) + " residual < 1e-8 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")";
    }
}

TEST(SolveQuartic, SolveQuarticTinyNonzeroDepressedOddTerm) {
    auto roots = LMCAS::solve_quartic(
        SymbolicExpr::number(1.0),
        SymbolicExpr::number(0.0),
        SymbolicExpr::number(0.0),
        SymbolicExpr::number(1e-15),
        SymbolicExpr::number(0.0),
        "x");
    EXPECT_TRUE((roots.size() == 4)) << "x^4+10^-15*x returns four roots with multiplicity";

    bool found_nonzero_real_root = false;
    for (const auto &root : roots) {
        try {
            const double value = root->to_numeric();
            if (std::isfinite(value) && std::abs(value + 1e-5) < 1e-12) {
                found_nonzero_real_root = true;
            }
        } catch (const std::exception &) {
        }
    }
    EXPECT_TRUE((found_nonzero_real_root)) << "tiny nonzero odd term must not be discarded as zero";
}

TEST(SolveQuartic, SolveQuarticResolventScalingPreservesTinyOddTerm) {
    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    bool threw = false;
    try {
        roots = LMCAS::solve_quartic(
            SymbolicExpr::number(1.0),
            SymbolicExpr::number(0.0),
            SymbolicExpr::number(0.0),
            SymbolicExpr::number(1e-200),
            SymbolicExpr::number(0.0),
            "x");
    } catch (const std::exception &) {
        threw = true;
    }
    EXPECT_TRUE((!threw && roots.size() == 4)) << "x^4+10^-200*x must not underflow its resolvent cubic";

    bool found_nonzero_real_root = false;
    const double expected = -std::cbrt(1e-200);
    for (const auto &root : roots) {
        try {
            const double value = root->to_numeric();
            if (std::isfinite(value) &&
                std::abs((value - expected) / expected) < 1e-12) {
                found_nonzero_real_root = true;
            }
        } catch (const std::exception &) {
        }
    }
    EXPECT_TRUE((found_nonzero_real_root)) << "scaled resolvent preserves the tiny nonzero real root";
}

TEST(SolveQuartic, SolveQuarticTinyNegativeQuadraticResolventRoot) {
    auto roots = LMCAS::solve_quartic(
        SymbolicExpr::number(1.0),
        SymbolicExpr::number(4.0),
        SymbolicExpr::number(6.0 + 1e-15),
        SymbolicExpr::number(4.0 + 2e-15),
        SymbolicExpr::number(1.0 + 1e-15),
        "x");
    EXPECT_TRUE((roots.size() == 4)) << "shifted x^4+10^-15*x^2 form returns four roots with multiplicity";

    int repeated_real_roots = 0;
    int nonreal_roots = 0;
    std::ostringstream root_details;
    for (const auto &root : roots) {
        const std::string text = root->simplify()->to_string();
        root_details << " [" << text << "]";
        if (text == "-1") {
            ++repeated_real_roots;
        } else {
            ++nonreal_roots;
        }
    }
    EXPECT_TRUE((repeated_real_roots == 2 && nonreal_roots == 2)) << "small negative u root must remain a conjugate imaginary pair; roots:" +
                                                                         root_details.str();
}

TEST(SolveQuartic, SolveQuarticComplexQuadraticResolventRoots) {
    auto roots = LMCAS::solve_quartic(
        SymbolicExpr::number(1),
        SymbolicExpr::number(4),
        SymbolicExpr::number(6),
        SymbolicExpr::number(4),
        SymbolicExpr::number(2),
        "x");
    EXPECT_TRUE((roots.size() == 4)) << "(x+1)^4+1 returns four complex roots";

    bool fabricated_real_root = false;
    for (const auto &root : roots) {
        const std::string text = root->simplify()->to_string();
        if (text == "0" || text == "-2") {
            fabricated_real_root = true;
        }
    }
    EXPECT_TRUE((!fabricated_real_root)) << "negative quadratic discriminant must not fabricate real u roots";
}

TEST(SolveQuartic, SolveQuarticTinyComplexPairIsNotClampedReal) {
    const double c1 = 0.25 + 3e-14;
    const double c2 = 2.0;
    auto roots = LMCAS::solve_quartic(
        SymbolicExpr::number(1.0),
        SymbolicExpr::number(0.0),
        SymbolicExpr::number(c1 + c2 - 1.0),
        SymbolicExpr::number(c2 - c1),
        SymbolicExpr::number(c1 * c2),
        "x");
    EXPECT_TRUE((roots.size() == 4)) << "product of two complex quadratics returns four roots";

    int nonreal_roots = 0;
    std::ostringstream root_details;
    for (const auto &root : roots) {
        root_details << " [" << root->simplify()->to_string();
        try {
            root_details << " -> " << root->to_numeric();
        } catch (const std::exception &) {
            ++nonreal_roots;
            root_details << " -> nonreal";
        }
        root_details << "]";
    }
    EXPECT_TRUE((nonreal_roots == 4)) << "small negative quadratic discriminant must remain nonreal; roots:" +
                                             root_details.str();
}

TEST(SolveQuartic, SolveQuarticLeadingCoefficient1) {
    auto a = SymbolicExpr::number(2);
    auto b = SymbolicExpr::number(-20);
    auto c = SymbolicExpr::number(70);
    auto d = SymbolicExpr::number(-100);
    auto e = SymbolicExpr::number(48);

    auto roots = LMCAS::solve_quartic(a, b, c, d, e, "x");
    EXPECT_TRUE((roots.size() == 4)) << "Quartic with a=2 should return 4 roots";

    for (size_t i = 0; i < roots.size(); ++i) {
        double val = roots[i]->to_numeric();
        double residual = eval_quartic(2, -20, 70, -100, 48, val);
        EXPECT_TRUE((std::abs(residual) < 1e-8)) << "Root " + std::to_string(i) + " residual < 1e-8 (val=" + std::to_string(val) + ", res=" + std::to_string(residual) + ")";
    }
}
