#include "test_common.hpp"
#include "rational.hpp"

#include <functional>
#include <stdexcept>
#include <string>

using namespace LMCAS;

namespace {

void expect_rejected(const std::string& literal) {
    bool rejected = false;
    try {
        (void)Rational(literal);
    } catch (const std::invalid_argument&) {
        rejected = true;
    } catch (const std::length_error&) {
        rejected = true;
    }
    std::string label = literal.size() > 80
        ? literal.substr(0, 40) + "...(" + std::to_string(literal.size()) + " bytes)"
        : literal;
    EXPECT_TRUE(rejected, "reject malformed Rational literal: '" + label + "'");
}

} // namespace

int main() {
    TEST_CASE("Strict Rational parser");
    EXPECT_EQ_STR(Rational("0").to_string(), "0", "zero");
    EXPECT_EQ_STR(Rational("-3").to_string(), "-3", "negative integer");
    EXPECT_EQ_STR(Rational("+12").to_string(), "12", "explicit positive sign");
    EXPECT_EQ_STR(Rational("1.25").to_string(), "5/4", "finite decimal");
    EXPECT_EQ_STR(Rational("1.25e-2").to_string(), "1/80", "negative exponent");
    EXPECT_EQ_STR(Rational(".5").to_string(), "1/2", "leading decimal point");
    EXPECT_EQ_STR(Rational("5.").to_string(), "5", "trailing decimal point");
    EXPECT_EQ_STR(Rational("0.(3)").to_string(), "1/3", "pure repeating decimal");
    EXPECT_EQ_STR(Rational("1.2(34)").to_string(), "611/495", "mixed repeating decimal");
    EXPECT_EQ_STR(Rational("1e+2").to_string(), "100", "positive exponent with explicit sign");
    EXPECT_EQ_STR(Rational("1e-0").to_string(), "1", "negative zero exponent");
    EXPECT_EQ_STR(Rational("0.(3)e1").to_string(), "10/3", "repeating decimal with exponent");
    EXPECT_EQ_STR(Rational("-0.0(9)").to_string(), "-1/10", "negative repeating decimal");

    for (const std::string& literal : {
             "", "+", "-", ".", "1e", "1e-", "1a", "1..2",
             "0.()", "0.(3", "0.3)", "(3)", "1(2)", "1.(2)(3)",
             "1e1000001", "1e+1000001", "1e--1", " 1", "1 "}) {
        expect_rejected(literal);
    }

    TEST_CASE("Rational parser safety limits");
    expect_rejected("0." + std::string(1000001, '0') + "1");
    expect_rejected("0.(" + std::string(1000001, '3') + ")");

    TEST_CASE("Rational roots");
    EXPECT_EQ_STR(Rational(0).sqrt(8).to_string(), "0", "sqrt(0)");
    EXPECT_EQ_STR(Rational(4, 9).sqrt(8).to_string(), "2/3", "exact sqrt(4/9)");
    EXPECT_EQ_STR(Rational(2).sqrt(6).to_string(), "1414213/1000000",
                  "sqrt(2) truncated to six decimal places");

    Rational cube(8);
    cube.radicand_self(BigInt(3), 6);
    EXPECT_EQ_STR(cube.to_string(), "2", "cube root of 8");
    Rational negative_cube(-8);
    negative_cube.radicand_self(BigInt(3), 6);
    EXPECT_EQ_STR(negative_cube.to_string(), "-2", "odd root of a negative value");

    bool negative_sqrt_rejected = false;
    try {
        (void)Rational(-1).sqrt(4);
    } catch (const std::domain_error&) {
        negative_sqrt_rejected = true;
    }
    EXPECT_TRUE(negative_sqrt_rejected, "negative real square root is rejected");

    return TEST_REPORT();
}
