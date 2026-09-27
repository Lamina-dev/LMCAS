#include "test_common.hpp"
#include "rational.hpp"

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>

using namespace LMCAS;

namespace {

void expect_rejected(const std::string &literal) {
    bool rejected = false;
    try {
        (void)Rational(literal);
    } catch (const std::invalid_argument &) {
        rejected = true;
    } catch (const std::length_error &) {
        rejected = true;
    }
    std::string label = literal.size() > 80
                            ? literal.substr(0, 40) + "...(" + std::to_string(literal.size()) + " bytes)"
                            : literal;
    EXPECT_TRUE((rejected)) << "reject malformed Rational literal: '" + label + "'";
}

} // namespace

TEST(RationalSafety, RationalParser) {
    EXPECT_EQ((Rational("0").to_string()), ("0")) << "zero";
    EXPECT_EQ((Rational("-3").to_string()), ("-3")) << "negative integer";
    EXPECT_EQ((Rational("+12").to_string()), ("12")) << "explicit positive sign";
    EXPECT_EQ((Rational("1.25").to_string()), ("5/4")) << "finite decimal";
    EXPECT_EQ((Rational("1.25e-2").to_string()), ("1/80")) << "negative exponent";
    EXPECT_EQ((Rational(".5").to_string()), ("1/2")) << "leading decimal point";
    EXPECT_EQ((Rational("5.").to_string()), ("5")) << "trailing decimal point";
    EXPECT_EQ((Rational("0.(3)").to_string()), ("1/3")) << "pure repeating decimal";
    EXPECT_EQ((Rational("1.2(34)").to_string()), ("611/495")) << "mixed repeating decimal";
    EXPECT_EQ((Rational("1e+2").to_string()), ("100")) << "positive exponent with explicit sign";
    EXPECT_EQ((Rational("1e-0").to_string()), ("1")) << "negative zero exponent";
    EXPECT_EQ((Rational("0.(3)e1").to_string()), ("10/3")) << "repeating decimal with exponent";
    EXPECT_EQ((Rational("-0.0(9)").to_string()), ("-1/10")) << "negative repeating decimal";

    for (const char *literal : {
             "", "+", "-", ".", "1e", "1e-", "1a", "1..2",
             "0.()", "0.(3", "0.3)", "(3)", "1(2)", "1.(2)(3)",
             "1e1000001", "1e+1000001", "1e--1", " 1", "1 "}) {
        expect_rejected(literal);
    }

    expect_rejected("0." + std::string(1000001, '0') + "1");
    expect_rejected("0.(" + std::string(1000001, '3') + ")");
}

TEST(RationalSafety, Binary64ExactConversion) {
    EXPECT_TRUE((Rational::from_double(0.0) == Rational(0))) << "positive zero";
    EXPECT_TRUE((Rational::from_double(-0.0) == Rational(0))) << "negative zero";
    EXPECT_TRUE((Rational::from_double(std::nextafter(1.0, 2.0)) ==
                 Rational(BigInt(std::int64_t{4503599627370497}),
                          BigInt(std::int64_t{4503599627370496}))))
        << "next binary64 above one retains its low significand bit";
    EXPECT_TRUE((Rational::from_double(std::nextafter(1.0, 0.0)) ==
                 Rational(BigInt(std::int64_t{9007199254740991}),
                          BigInt(std::int64_t{9007199254740992}))))
        << "next binary64 below one uses the smaller binade spacing";
    EXPECT_TRUE((Rational::from_double(std::nextafter(-1.0, -2.0)) ==
                 Rational(BigInt(std::int64_t{-4503599627370497}),
                          BigInt(std::int64_t{4503599627370496}))))
        << "negative adjacent-to-one value retains its sign and low bit";
    EXPECT_TRUE((Rational::from_double(0.1) ==
                 Rational(BigInt(std::int64_t{3602879701896397}),
                          BigInt(std::int64_t{36028797018963968}))))
        << "binary64 one tenth is not the decimal rational one tenth";

    const BigInt subnormal_denominator = BigInt(1) << 1074;
    EXPECT_TRUE((Rational::from_double(std::numeric_limits<double>::denorm_min()) ==
                 Rational(BigInt(1), subnormal_denominator)))
        << "smallest subnormal is exactly two to the minus 1074";
    EXPECT_TRUE((Rational::from_double(-std::numeric_limits<double>::denorm_min()) ==
                 Rational(BigInt(-1), subnormal_denominator)))
        << "smallest negative subnormal retains its sign";
    EXPECT_TRUE((Rational::from_double(
                     std::nextafter(std::numeric_limits<double>::min(), 0.0)) ==
                 Rational((BigInt(1) << 52) - BigInt(1), subnormal_denominator)))
        << "largest subnormal retains every significand bit";
    EXPECT_TRUE((Rational::from_double(std::numeric_limits<double>::min()) ==
                 Rational(BigInt(1), BigInt(1) << 1022)))
        << "smallest normal is exactly two to the minus 1022";
    const BigInt max_finite = (BigInt(1) << 1024) - (BigInt(1) << 971);
    EXPECT_TRUE((Rational::from_double(std::numeric_limits<double>::max()) ==
                 Rational(max_finite)))
        << "largest finite binary64 has its exact integer value";
    EXPECT_TRUE((Rational::from_double(std::numeric_limits<double>::lowest()) ==
                 Rational(-max_finite)))
        << "lowest finite binary64 has its exact integer value";

    for (double value : {std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
        bool rejected = false;
        try {
            (void)Rational::from_double(value);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        EXPECT_TRUE((rejected)) << "nonfinite binary64 conversion is rejected";
    }
}

TEST(RationalSafety, RationalPowerAliases) {
    Rational self_power(3, 2);
    self_power.power_self(self_power, 2);
    EXPECT_TRUE((self_power == Rational(183, 100))) << "(3/2)^(3/2) truncates to 1.83 with an aliased exponent";
    Rational separate_power(3, 2);
    separate_power.power_self(Rational(3, 2), 2);
    EXPECT_TRUE((separate_power == Rational(183, 100))) << "(3/2)^(3/2) uses the same truncation without aliasing";
    Rational negative_self_power(-2, 3);
    negative_self_power.power_self(negative_self_power, 2);
    EXPECT_TRUE((negative_self_power == Rational(131, 100))) << "(-2/3)^(-2/3) truncates the cube root of 9/4 to 1.31";
}

TEST(RationalSafety, RationalRoots) {
    EXPECT_EQ((Rational(0).sqrt(8).to_string()), ("0")) << "sqrt(0)";
    EXPECT_EQ((Rational(4, 9).sqrt(8).to_string()), ("2/3")) << "exact sqrt(4/9)";
    EXPECT_EQ((Rational(2).sqrt(6).to_string()), ("1414213/1000000")) << "sqrt(2) truncated to six decimal places";

    Rational cube(8);
    cube.radicand_self(BigInt(3), 6);
    EXPECT_EQ((cube.to_string()), ("2")) << "cube root of 8";
    Rational negative_cube(-8);
    negative_cube.radicand_self(BigInt(3), 6);
    EXPECT_EQ((negative_cube.to_string()), ("-2")) << "odd root of a negative value";

    bool negative_sqrt_rejected = false;
    try {
        (void)Rational(-1).sqrt(4);
    } catch (const std::domain_error &) {
        negative_sqrt_rejected = true;
    }
    EXPECT_TRUE((negative_sqrt_rejected)) << "negative real square root is rejected";
}

TEST(RationalSafety, DecimalTruncation) {
    struct Case {
        Rational input;
        std::int64_t precision;
        Rational expected;
    };
    const Case cases[] = {
        {Rational(3, 2), -1, Rational(0)},
        {Rational(-3, 2), -1, Rational(0)},
        {Rational(3, 2), 0, Rational(1)},
        {Rational(-3, 2), 0, Rational(-1)},
        {Rational(3, 2), 1, Rational(3, 2)},
        {Rational(-3, 2), 1, Rational(-3, 2)},
        {Rational(159, 10), -1, Rational(10)},
        {Rational(-159, 10), -1, Rational(-10)},
        {Rational(-1, 2), 0, Rational(0)},
        {Rational(-7), 0, Rational(-7)},
        {Rational(0), 0, Rational(0)},
    };
    for (const auto &test : cases) {
        auto value = test.input;
        value.floor(test.precision);
        EXPECT_TRUE((value == test.expected)) << test.input.to_string() + " truncated at " + std::to_string(test.precision);
        if (test.precision <= 0) {
            EXPECT_TRUE((value.get_denominator() == BigInt(1))) << "integer truncation has a canonical unit denominator";
        }
    }
    for (const std::int64_t precision : {
             std::int64_t(-1000001), std::int64_t(1000001),
             std::numeric_limits<std::int64_t>::min(),
             std::numeric_limits<std::int64_t>::max()}) {
        Rational value(-159, 10);
        bool rejected = false;
        try {
            value.floor(precision);
        } catch (const std::length_error &) {
            rejected = true;
        }
        EXPECT_TRUE((rejected && value == Rational(-159, 10))) << "precision rejection is explicit and does not mutate the value";
    }
}
