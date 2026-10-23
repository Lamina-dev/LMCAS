#include "expr.hpp"
#include "test_common.hpp"

#include "value.hpp"

using namespace LMCAS;

TEST(Proof, LogarithmChangeOfBase) {
    auto x = SymbolicExpr::variable("x");
    auto b = SymbolicExpr::variable("b");

    auto lhs = SymbolicExpr::log(x, b);

    auto ln_x = SymbolicExpr::ln(x);
    auto ln_b = SymbolicExpr::ln(b);
    auto rhs = SymbolicExpr::multiply(ln_x, SymbolicExpr::power(ln_b, SymbolicExpr::number(-1)));

    auto lhs_simplified = LMCAS::simplify(lhs).value();

    auto rhs_simplified = LMCAS::simplify(rhs).value();

    bool string_match = (lhs_simplified->to_string() == rhs_simplified->to_string());

    auto diff = SymbolicExpr::add(lhs_simplified, SymbolicExpr::multiply(rhs_simplified, SymbolicExpr::number(-1)));
    auto diff_simplified = LMCAS::simplify(diff).value();
    bool diff_is_zero = (diff_simplified->is_number() && diff_simplified->convert_rational() == Rational(0));

    EXPECT_TRUE(string_match || diff_is_zero) << "log_b(x) == ln(x)/ln(b)";
}
