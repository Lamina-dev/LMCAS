#include "test_common.hpp"
#include "numeric_evaluation.hpp"
#include "symbolic_ast.hpp"

int main() {
    using namespace LMCAS;

    TEST_CASE("Explicit numeric bindings");
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(x, SymbolicExpr::number(2));
    auto bound = evaluate_numeric(*expression, {{"x", 5.0}});
    EXPECT_TRUE(bound.has_value(), "bound expression evaluates");
    if (bound) EXPECT_NEAR(bound.value().value, 7.0, 1e-12, "x + 2 at x=5");

    TEST_CASE("Unbound symbols are errors");
    auto unbound = evaluate_numeric(*expression);
    EXPECT_FALSE(unbound.has_value(), "unbound x does not become zero");
    if (!unbound) {
        EXPECT_TRUE(unbound.error().code == CasErrc::UnboundSymbol,
                    "unbound x reports UnboundSymbol");
    }

    TEST_CASE("Domain errors are explicit");
    auto invalid_log = SymbolicExpr::ln(SymbolicExpr::number(-1));
    auto log_result = evaluate_numeric(*invalid_log);
    EXPECT_FALSE(log_result.has_value(), "ln(-1) is not a real number");
    if (!log_result) {
        EXPECT_TRUE(log_result.error().code == CasErrc::DomainError,
                    "ln(-1) reports DomainError");
    }

    TEST_CASE("Constants have one numeric meaning");
    auto ascii_pi = evaluate_numeric(*SymbolicExpr::variable("pi"));
    auto unicode_pi = evaluate_numeric(*SymbolicExpr::variable("π"));
    EXPECT_TRUE(ascii_pi.has_value() && unicode_pi.has_value(), "pi constants evaluate");
    if (ascii_pi && unicode_pi) {
        EXPECT_NEAR(ascii_pi.value().value, unicode_pi.value().value, 0.0,
                    "pi and unicode pi agree");
    }

    TEST_CASE("Certified algebraic roots have real numeric values");
    auto root_variable = SymbolicExpr::variable("z");
    auto root_polynomial = SymbolicExpr::add(
        SymbolicExpr::multiply(root_variable, root_variable),
        SymbolicExpr::number(-2));
    auto root = evaluate_numeric(
        *SymbolicExpr::root_of(root_polynomial, "z", 0));
    EXPECT_TRUE(root.has_value(), "real RootOf evaluates numerically");
    if (root) {
        EXPECT_NEAR(root.value().value, -std::sqrt(2.0), 1e-12,
                    "first real root of z^2-2 is -sqrt(2)");
    }

    TEST_CASE("Real-axis complex arguments preserve the principal branch");
    auto argument_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ComplexArg,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(x)}));
    const struct {
        double input;
        double expected;
    } argument_cases[] = {
        {-4.0, std::acos(-1.0)},
        {4.0, 0.0},
        {-0.0, std::acos(-1.0)},
        {0.0, 0.0}
    };
    for (const auto& test : argument_cases) {
        auto argument = evaluate_numeric(*argument_expression, {{"x", test.input}});
        EXPECT_TRUE(argument && argument.value().is_finite(),
                    "bound real-axis argument has a finite principal value");
        if (argument) {
            EXPECT_NEAR(argument.value().value, test.expected, 1e-14,
                        "principal argument distinguishes positive and negative real axes");
        }
    }
    auto imaginary_expression = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ImagPart,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(x)}));
    auto imaginary = evaluate_numeric(*imaginary_expression, {{"x", -4.0}});
    EXPECT_TRUE(imaginary && imaginary.value().value == 0.0,
                "negative real values have zero imaginary part");

    TEST_CASE("Exact power exponents preserve parity beyond double precision");
    const BigInt odd_exponent("9007199254740993");
    for (auto exponent : {SymbolicExpr::number(odd_exponent),
                          SymbolicExpr::number(Rational(odd_exponent, BigInt(1)))}) {
        auto odd_power = SymbolicExpr::power(x, exponent);
        auto negative = evaluate_numeric(*odd_power, {{"x", -1.0}});
        EXPECT_TRUE(negative && negative.value().value == -1.0,
                    "an exact odd exponent keeps the sign of a negative base");
        auto signed_zero = evaluate_numeric(*odd_power, {{"x", -0.0}});
        EXPECT_TRUE(signed_zero && signed_zero.value().value == 0.0 &&
                        std::signbit(signed_zero.value().value),
                    "an exact odd exponent preserves negative zero");
    }
    auto even_power = SymbolicExpr::power(
        x, SymbolicExpr::number(odd_exponent + BigInt(1)));
    auto even_value = evaluate_numeric(*even_power, {{"x", -1.0}});
    EXPECT_TRUE(even_value && even_value.value().value == 1.0,
                "an adjacent exact even exponent produces a positive result");

    TEST_CASE("Negative real powers retain exact exponent domain restrictions");
    auto fractional_power = SymbolicExpr::power(
        x, SymbolicExpr::number(Rational(odd_exponent, BigInt(2))));
    auto fractional_value = evaluate_numeric(*fractional_power, {{"x", -1.0}});
    EXPECT_TRUE(!fractional_value &&
                    fractional_value.error().code == CasErrc::DomainError,
                "a fraction rounded to an integer double remains outside the real power domain");

    TEST_CASE("Resource limits are enforced");
    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext context(limits);
    auto exhausted = evaluate_numeric(*expression, {{"x", 1.0}}, context);
    EXPECT_FALSE(exhausted.has_value(), "step budget stops recursive evaluation");
    if (!exhausted) {
        EXPECT_TRUE(exhausted.error().code == CasErrc::ResourceLimit,
                    "budget exhaustion reports ResourceLimit");
    }

    return TEST_REPORT();
}
