#include "test_common.hpp"
#include "numeric_evaluation.hpp"

int main() {
    using namespace lamina;

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
