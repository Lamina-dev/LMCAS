#include "expr.hpp"
#include <gtest/gtest.h>
#include "numeric_evaluation.hpp"
#include <limits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <utility>

using namespace LMCAS;

namespace {

TEST(ExprParsing, OrdinaryAndImaginaryNames) {
    auto ordinary_i = LMCAS::sym("i");
    EXPECT_TRUE((ordinary_i && ordinary_i.value()->to_string() == "i")) << "lowercase i is an ordinary symbolic identifier";

    auto reserved_I = LMCAS::sym("I");
    EXPECT_TRUE((!reserved_I &&
                 reserved_I.error().code == LMCAS::CasErrc::InvalidArgument))
        << "uppercase I is the imaginary unit and cannot be shadowed as a symbol";
}

TEST(ExprParsing, ReservedConstantNames) {
    auto reserved_pi = LMCAS::sym("pi");
    EXPECT_TRUE((!reserved_pi &&
                 reserved_pi.error().code == LMCAS::CasErrc::InvalidArgument))
        << "pi is a constant and cannot be shadowed as a symbol";

    auto reserved_unicode_pi = LMCAS::sym("\xCF\x80");
    EXPECT_TRUE((!reserved_unicode_pi &&
                 reserved_unicode_pi.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "unicode pi is a constant alias and cannot be shadowed";

    auto reserved_e = LMCAS::sym("e");
    EXPECT_TRUE((!reserved_e &&
                 reserved_e.error().code == LMCAS::CasErrc::InvalidArgument))
        << "e is a constant and cannot be shadowed as a symbol";

    auto reserved_phi = LMCAS::sym("phi");
    EXPECT_TRUE((!reserved_phi &&
                 reserved_phi.error().code == LMCAS::CasErrc::InvalidArgument))
        << "phi is a constant and cannot be shadowed as a symbol";
}

TEST(ExprParsing, ComplexAndRelationSyntax) {
    auto parsed_complex = LMCAS::parse_expr("x + 2*I");
    EXPECT_TRUE((parsed_complex.has_value())) << "parser accepts reserved imaginary unit";
    {
        const std::string actual_text = (parsed_complex.value()->to_string());
        for (const auto &token : std::vector<std::string>{"I", "x"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "parser builds complex expression with imaginary unit" << ": missing " << token << " in " << actual_text;
        }
    }
    auto parsed_relation = LMCAS::parse_expr("x > 0 and x < 1");
    EXPECT_TRUE((parsed_relation.has_value())) << "parser accepts relational logical expressions";
    {
        const std::string actual_text = (parsed_relation.value()->to_string());
        for (const auto &token : std::vector<std::string>{"and", "x", "0", "1"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "parser prints logical expressions with ASCII operators" << ": missing " << token << " in " << actual_text;
        }
    }
}

TEST(ExprParsing, SetIntervalAndMembershipSyntax) {
    auto parsed_set_literal = LMCAS::parse_expr("{-1, 1}");
    EXPECT_TRUE((parsed_set_literal.has_value())) << "parser accepts finite set literals";
    {
        const std::string actual_text = (parsed_set_literal.value()->to_string());
        for (const auto &token : std::vector<std::string>{"{", "1", "}"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "parser prints finite set literals" << ": missing " << token << " in " << actual_text;
        }
    }

    auto parsed_interval = LMCAS::parse_expr("[0, 1)");
    ASSERT_TRUE(parsed_interval.has_value()) << "parser accepts half-open interval literals";
    ASSERT_TRUE(parsed_interval.value()) << "parsed interval has an expression";
    EXPECT_EQ(parsed_interval.value()->to_string(), "[0, 1)")
        << "parser prints half-open interval literals";

    auto parsed_membership = LMCAS::parse_expr("x in {-1, 1}");
    EXPECT_TRUE((parsed_membership.has_value())) << "parser accepts membership expressions";
    {
        const std::string actual_text = (parsed_membership.value()->to_string());
        for (const auto &token : std::vector<std::string>{"x", "in", "{"}) {
            EXPECT_NE(actual_text.find(token), std::string::npos) << "parser prints membership expressions" << ": missing " << token << " in " << actual_text;
        }
    }
}

TEST(ExprParsing, SymbolicCallsAndMalformedSyntax) {
    auto parsed_pow_star = LMCAS::parse_expr("x**2");
    EXPECT_TRUE((!parsed_pow_star &&
                 parsed_pow_star.error().code == LMCAS::CasErrc::ParseError))
        << "parser rejects ** power syntax";
    auto parsed_unknown_call = LMCAS::parse_expr("f(x, y)");
    ASSERT_TRUE(parsed_unknown_call.has_value()) << "parser accepts uninterpreted symbolic function calls";
    ASSERT_TRUE(parsed_unknown_call.value()) << "parsed call has an expression";
    EXPECT_EQ(parsed_unknown_call.value()->to_string(), "f(x, y)")
        << "parser prints uninterpreted symbolic function calls";

    auto parsed_invalid = LMCAS::parse_expr("x +");
    EXPECT_TRUE((!parsed_invalid &&
                 parsed_invalid.error().code == LMCAS::CasErrc::ParseError))
        << "parser reports malformed input as ParseError";
}

TEST(ExprParsing, PolynomialAndFunctionSyntax) {
    auto polynomial = LMCAS::parse_expr("x^2 + 1");
    ASSERT_TRUE((polynomial.has_value())) << "parser accepts powers and sums";
    if (polynomial) {
        auto value = LMCAS::evalf(*polynomial.value(), {{"x", 3.0}});
        EXPECT_TRUE((value && value.value().value == 10.0)) << "parsed polynomial preserves exponent and addition semantics";
    }
    auto trigonometric = LMCAS::parse_expr("sin(x)^2 + cos(x)^2");
    ASSERT_TRUE((trigonometric.has_value())) << "parser accepts common functions";
    if (trigonometric) {
        auto value = LMCAS::evalf(*trigonometric.value(), {{"x", 0.5}});
        ASSERT_TRUE((value.has_value())) << "parsed trigonometric expression evaluates";
        if (value) {
            EXPECT_NEAR(value.value().value, 1.0, 1e-12) << "parsed functions preserve the trigonometric identity";
        }
    }
}

TEST(ExprParsing, MultiArgumentSyntax) {
    auto maximum = LMCAS::parse_expr("max(x, y, 0)");
    ASSERT_TRUE((maximum.has_value())) << "parser accepts supported multi-argument functions";
    if (maximum) {
        auto value = LMCAS::evalf(*maximum.value(), {{"x", -2.0}, {"y", 3.0}});
        EXPECT_TRUE((value && value.value().value == 3.0)) << "parsed max compares all arguments";
    }
}

TEST(ExprParsing, LogarithmWithExplicitBase) {
    auto parsed = LMCAS::parse_expr("log(x, 10)");
    ASSERT_TRUE(parsed) << "parser accepts the two-argument logarithm";
    ASSERT_TRUE(parsed.value());

    auto value = LMCAS::evalf(*parsed.value(), {{"x", 100.0}});
    ASSERT_TRUE(value) << "parsed base-10 logarithm evaluates";
    ASSERT_TRUE(value.value().is_finite());
    EXPECT_NEAR(value.value().value, 2.0, 1e-12)
        << "log base 10 of 100 equals 2";
}

TEST(ExprParsing, NeRelationSyntax) {
    auto x = LMCAS::sym("x");
    auto rhs = LMCAS::integer(2);
    auto named = LMCAS::ne(x.value(), rhs.value());
    auto parsed = LMCAS::parse_expr("x != 2");
    EXPECT_TRUE((named && parsed &&
                 LMCAS::structurally_equal(*named.value(), *parsed.value())))
        << "named ne preserves the parsed relation operator";
}

TEST(ExprParsing, LtRelationSyntax) {
    auto x = LMCAS::sym("x");
    auto rhs = LMCAS::integer(2);
    auto named = LMCAS::lt(x.value(), rhs.value());
    auto parsed = LMCAS::parse_expr("x < 2");
    EXPECT_TRUE((named && parsed &&
                 LMCAS::structurally_equal(*named.value(), *parsed.value())))
        << "named lt preserves the parsed relation operator";
}

TEST(ExprParsing, LeRelationSyntax) {
    auto x = LMCAS::sym("x");
    auto rhs = LMCAS::integer(2);
    auto named = LMCAS::le(x.value(), rhs.value());
    auto parsed = LMCAS::parse_expr("x <= 2");
    EXPECT_TRUE((named && parsed &&
                 LMCAS::structurally_equal(*named.value(), *parsed.value())))
        << "named le preserves the parsed relation operator";
}

TEST(ExprParsing, GtRelationSyntax) {
    auto x = LMCAS::sym("x");
    auto rhs = LMCAS::integer(2);
    auto named = LMCAS::gt(x.value(), rhs.value());
    auto parsed = LMCAS::parse_expr("x > 2");
    EXPECT_TRUE((named && parsed &&
                 LMCAS::structurally_equal(*named.value(), *parsed.value())))
        << "named gt preserves the parsed relation operator";
}

TEST(ExprParsing, GeRelationSyntax) {
    auto x = LMCAS::sym("x");
    auto rhs = LMCAS::integer(2);
    auto named = LMCAS::ge(x.value(), rhs.value());
    auto parsed = LMCAS::parse_expr("x >= 2");
    EXPECT_TRUE((named && parsed &&
                 LMCAS::structurally_equal(*named.value(), *parsed.value())))
        << "named ge preserves the parsed relation operator";
}

TEST(ExprParsing, DivisionRetainsDomain) {
    for (const char *source : {"0/0", "0/x", "0/2"}) {
        auto parsed = parse_expr(source);
        ASSERT_TRUE((parsed.has_value())) << "division syntax parses independently of its value";
        if (!parsed)
            continue;
        for (double x : {0.0, 2.0}) {
            auto value = evaluate_numeric(*parsed.value(), {{"x", x}});
            const bool invalid = std::string(source) == "0/0" ||
                                 (std::string(source) == "0/x" && x == 0);
            if (invalid) {
                EXPECT_TRUE((!value && value.error().code == CasErrc::DomainError)) << "zero numerator cannot erase a denominator pole";
            } else {
                EXPECT_TRUE((value && value.value().value == 0.0)) << "legal division evaluates to zero";
            }
        }
    }
}

void expect_binary64_decimal_roundtrips() {
    const double values[] = {
        1.2345678901234567, std::nextafter(1.0, 0.0), std::nextafter(1.0, 2.0),
        1.2345678901234567e-200, std::numeric_limits<double>::denorm_min(),
        std::nextafter(std::numeric_limits<double>::min(), 0.0),
        std::numeric_limits<double>::min(), std::numeric_limits<double>::max(),
        0.0, -0.0, 0.5, 1.0};
    for (double expected : values) {
        auto original = approx_real(expected);
        ASSERT_TRUE((original.has_value())) << "finite approximate factory succeeds";
        if (!original) {
            continue;
        }
        auto parsed = parse_expr(original.value()->to_string());
        ASSERT_TRUE((parsed.has_value())) << "canonical approximate literal parses";
        if (!parsed) {
            continue;
        }
        auto actual = evalf(*parsed.value());
        ASSERT_TRUE((actual.has_value())) << "approximate literal evaluates";
        if (!actual) {
            continue;
        }
        std::uint64_t expected_bits, actual_bits;
        std::memcpy(&expected_bits, &expected, sizeof expected_bits);
        std::memcpy(&actual_bits, &actual.value().value, sizeof actual_bits);
        EXPECT_TRUE((expected_bits == actual_bits)) << "roundtrip retains every binary64 bit";
        EXPECT_TRUE((structurally_equal(*original.value(), *parsed.value()))) << "roundtrip retains approximate number identity";
        auto exact = rational(Rational::from_double(expected));
        EXPECT_TRUE((exact && !structurally_equal(*exact.value(), *parsed.value()))) << "approximate literal does not become an exact rational";
    }
}

void expect_approximate_decimal_syntax() {
    for (const char *text : {"approx( +.5 )", "approx(5.e-1)", "approx(.50)"}) {
        auto parsed = parse_expr(text);
        ASSERT_TRUE((parsed.has_value())) << "signed decimal grammar accepts valid forms";
        if (!parsed) {
            continue;
        }
        auto actual = evalf(*parsed.value());
        EXPECT_TRUE((actual && actual.value().value == 0.5)) << "valid decimals retain value";
    }
    for (const auto &test : {std::make_pair("1.25", Rational(5, 4)),
                             std::make_pair("1e-3", Rational(1, 1000))}) {
        auto parsed = parse_expr(test.first);
        auto exact = rational(test.second);
        EXPECT_TRUE((parsed && exact && structurally_equal(*parsed.value(), *exact.value()))) << "ordinary decimals remain exact";
    }
    auto identifier = parse_expr("approx");
    auto symbol = sym("approx");
    EXPECT_TRUE((identifier && symbol && structurally_equal(*identifier.value(), *symbol.value()))) << "bare approx remains an ordinary symbol";
}

TEST(ExprParsing, ApproximateDecimalRoundtrip) {
    expect_binary64_decimal_roundtrips();
    expect_approximate_decimal_syntax();
    for (const char *source : {"x + approx(1.2345678901234567)",
                               "x * approx(1.2345678901234567)",
                               "x^approx(1.2345678901234567)"}) {
        auto original = parse_expr(source);
        ASSERT_TRUE((original.has_value())) << "nested approximate expression parses";
        if (!original) {
            continue;
        }
        auto parsed = parse_expr(original.value()->to_string());
        EXPECT_TRUE((parsed && structurally_equal(*original.value(), *parsed.value()))) << "nested arithmetic retains approximate atoms and their full value";
    }
}

TEST(ExprParsing, ApproximateDecimalRejections) {
    for (const char *text : {
             "approx()", "approx(x)", "approx(1/2)", "approx(1,2)", "approx(1e)",
             "approx(nan)", "approx(inf)", "approx(-inf)", "approx(0x1p0)",
             "approx(1e400)", "approx(-1e400)", "approx(1e-4000)",
             "approx(-1e-4000)", "approx(2e-324)", "approx(1 2)",
             "approx(1e +2)", "approx(1)junk",
             "approx((1))", "approx(1+2)", "approx(--1)", "approx(1"}) {
        auto parsed = parse_expr(text);
        EXPECT_TRUE((!parsed && parsed.error().code == CasErrc::ParseError)) << "malformed approximate literal is an explicit parse error";
    }
}

} // namespace
