#include "test_common.hpp"
#include "newton_raphson.hpp"
#include "internal/exact_sturm.hpp"

using namespace LMCAS;

using LMCAS::Polynomial;

TEST(SturmIsolation, StrictCauchyRootBound) {
    const Polynomial<Rational> monic(
        {Rational(-4), Rational(0), Rational(1)}, "x");
    EXPECT_TRUE(detail::strict_cauchy_root_bound(monic) == Rational(6));

    const Polynomial<Rational> negative_leading(
        {Rational(1), Rational(-5), Rational(-2)}, "x");
    EXPECT_TRUE(detail::strict_cauchy_root_bound(negative_leading) ==
                Rational(9, 2));
}

TEST(SturmIsolation, SturmIsolationLinearPolynomialX3) {
    Polynomial<Rational> p({Rational(-3), Rational(1)}, "x");
    auto intervals = LMCAS::isolate_real_roots_checked(p).value();
    EXPECT_TRUE((intervals.size() == 1)) << "x-3 should have 1 real root";
    if (!intervals.empty()) {
        EXPECT_TRUE((intervals[0].first <= Rational(3) && intervals[0].second >= Rational(3))) << "Root x=3 should be in the interval";
    }
}

TEST(SturmIsolation, SturmIsolationQuadraticX24RootsAt22) {
    Polynomial<Rational> p({Rational(-4), Rational(0), Rational(1)}, "x");
    auto intervals = LMCAS::isolate_real_roots_checked(p).value();
    EXPECT_TRUE((intervals.size() == 2)) << "x^2-4 should have 2 real roots";
    if (intervals.size() == 2) {

        EXPECT_TRUE((intervals[0].first <= Rational(-2) && intervals[0].second >= Rational(-2))) << "First root x=-2 should be in first interval";

        EXPECT_TRUE((intervals[1].first <= Rational(2) && intervals[1].second >= Rational(2))) << "Second root x=2 should be in second interval";
    }
}

TEST(SturmIsolation, SturmIsolationCubicX36x211x6RootsAt123) {
    Polynomial<Rational> p({Rational(-6), Rational(11), Rational(-6), Rational(1)}, "x");
    auto intervals = LMCAS::isolate_real_roots_checked(p).value();
    EXPECT_TRUE((intervals.size() == 3)) << "x^3-6x^2+11x-6 should have 3 real roots";
    if (intervals.size() == 3) {
        EXPECT_TRUE((intervals[0].first <= Rational(1) && intervals[0].second >= Rational(1))) << "Root x=1 should be in first interval";
        EXPECT_TRUE((intervals[1].first <= Rational(2) && intervals[1].second >= Rational(2))) << "Root x=2 should be in second interval";
        EXPECT_TRUE((intervals[2].first <= Rational(3) && intervals[2].second >= Rational(3))) << "Root x=3 should be in third interval";
    }
}

TEST(SturmIsolation, SturmIsolationNoRealRootsX21) {
    Polynomial<Rational> p({Rational(1), Rational(0), Rational(1)}, "x");
    auto intervals = LMCAS::isolate_real_roots_checked(p).value();
    EXPECT_TRUE((intervals.size() == 0)) << "x^2+1 should have 0 real roots";
}

TEST(SturmIsolation, SturmIsolationDoubleRootX22x1X12) {
    Polynomial<Rational> p({Rational(1), Rational(-2), Rational(1)}, "x");
    auto intervals = LMCAS::isolate_real_roots_checked(p).value();

    EXPECT_TRUE((intervals.size() == 1)) << "(x-1)^2 should have 1 distinct real root";
    if (!intervals.empty()) {
        EXPECT_TRUE((intervals[0].first <= Rational(1) && intervals[0].second >= Rational(1))) << "Root x=1 should be in the interval";
    }
}

TEST(SturmIsolation, SturmIsolationQuarticX22X23With4RealRoots) {
    Polynomial<Rational> p({Rational(6), Rational(0), Rational(-5), Rational(0), Rational(1)}, "x");
    auto intervals = LMCAS::isolate_real_roots_checked(p).value();
    EXPECT_TRUE((intervals.size() == 4)) << "(x^2-2)(x^2-3) should have 4 real roots";
}

TEST(SturmIsolation, SturmIsolationConstantPolynomial) {
    Polynomial<Rational> p({Rational(5)}, "x");
    auto intervals = LMCAS::isolate_real_roots_checked(p).value();
    EXPECT_TRUE((intervals.size() == 0)) << "Constant polynomial should have 0 roots";
}

TEST(SturmIsolation, SturmIsolationZeroPolynomial) {
    Polynomial<Rational> p("x");
    auto intervals = LMCAS::isolate_real_roots_checked(p).value();
    EXPECT_TRUE((intervals.size() == 0)) << "Zero polynomial should return empty";
}
