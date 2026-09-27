#include "test_common.hpp"

using namespace LMCAS;

namespace {

void expect_factor_multiset(
    const std::shared_ptr<SymbolicExpr> &actual,
    const std::vector<std::shared_ptr<SymbolicExpr>> &expected) {
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(
        detail::node(actual));
    ASSERT_NE(product, nullptr);
    ASSERT_EQ(product->operands().size(), expected.size());

    std::vector<bool> matched(expected.size(), false);
    for (const auto &operand : product->operands()) {
        auto factor = detail::make_expression_ptr(operand);
        bool found = false;
        for (std::size_t i = 0; i < expected.size(); ++i) {
            if (!matched[i] &&
                test_proved_equivalent(factor, expected[i])) {
                matched[i] = true;
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "unexpected factor " << factor->to_string();
    }
    for (bool was_matched : matched) {
        EXPECT_TRUE(was_matched);
    }
}

} // namespace

TEST(LmcasFactor, CommonTerm) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto two = SymbolicExpr::number(2);
    auto expression = SymbolicExpr::add(
        SymbolicExpr::multiply(two, x),
        SymbolicExpr::multiply(two, y));

    auto checked = expression->factor_checked();

    ASSERT_TRUE(checked) << checked.error().message;
    expect_factor_multiset(
        checked.value(),
        {two, SymbolicExpr::add(x, y)});
    EXPECT_TRUE(test_proved_equivalent(
        checked.value(), expression));
}

TEST(LmcasFactor, Quadratic) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(5), x),
            SymbolicExpr::number(6)));

    auto checked = expression->factor_checked();

    ASSERT_TRUE(checked) << checked.error().message;
    expect_factor_multiset(
        checked.value(),
        {SymbolicExpr::add(x, SymbolicExpr::number(2)),
         SymbolicExpr::add(x, SymbolicExpr::number(3))});
    EXPECT_TRUE(test_proved_equivalent(
        checked.value(), expression));
}

TEST(LmcasFactor, DifferenceOfSquares) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-4));

    auto checked = expression->factor_checked();

    ASSERT_TRUE(checked) << checked.error().message;
    expect_factor_multiset(
        checked.value(),
        {SymbolicExpr::add(x, SymbolicExpr::number(2)),
         SymbolicExpr::add(x, SymbolicExpr::number(-2))});
    EXPECT_TRUE(test_proved_equivalent(
        checked.value(), expression));
}

TEST(LmcasFactor, CheckedResourceLimit) {
    auto x = SymbolicExpr::variable("x");

    LMCAS::ResourceLimits limits;
    limits.max_steps = 0;
    LMCAS::ComputationContext context(limits);
    auto limited = x->factor_checked(context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "factor_checked preserves an exhausted budget";
}
