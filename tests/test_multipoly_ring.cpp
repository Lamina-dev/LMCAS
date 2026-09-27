/**
 * @file test_multipoly_ring.cpp
 * @brief MultiPoly 环运算契约。
 */

#include "test_common.hpp"
#include "multivariate_poly.hpp"

using namespace LMCAS;

namespace {

TEST(MultipolyRing, AdditionCommutativityLinear) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(3)}, {{0, 1}, Rational(2)}}, vars);
    MultiPoly b({{{1, 0}, Rational(-1)}, {{0, 0}, Rational(5)}}, vars);
    EXPECT_TRUE((a + b == b + a)) << "commutativity: (3x+2y) + (-x+5) == (-x+5) + (3x+2y)";
}

TEST(MultipolyRing, AdditionCommutativityZero) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{2, 1}, Rational(1)}, {{0, 0}, Rational(7)}}, vars);
    MultiPoly zero;
    EXPECT_TRUE((a + zero == zero + a)) << "commutativity: a + 0 == 0 + a";
}

TEST(MultipolyRing, AdditionCommutativityDense) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{3, 0}, Rational(1)}, {{2, 1}, Rational(-2)}, {{1, 2}, Rational(3)}, {{0, 3}, Rational(-4)}}, vars);
    MultiPoly b({{{0, 0}, Rational(1)}, {{1, 0}, Rational(2)}, {{0, 1}, Rational(3)}, {{1, 1}, Rational(4)}}, vars);
    EXPECT_TRUE((a + b == b + a)) << "commutativity: dense cubic + dense linear";
}

TEST(MultipolyRing, AdditionCommutativityRational) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(1, 3)}, {{0, 1}, Rational(2, 7)}}, vars);
    MultiPoly b({{{1, 0}, Rational(5, 6)}, {{0, 0}, Rational(-1, 2)}}, vars);
    EXPECT_TRUE((a + b == b + a)) << "commutativity: rational coefficients";
}

TEST(MultipolyRing, AdditionCommutativityTrivariate) {
    std::vector<std::string> vars3 = {"x", "y", "z"};
    MultiPoly a({{{1, 1, 1}, Rational(2)}, {{2, 0, 0}, Rational(1)}}, vars3);
    MultiPoly b({{{0, 0, 2}, Rational(3)}, {{1, 1, 0}, Rational(-1)}}, vars3);
    EXPECT_TRUE((a + b == b + a)) << "commutativity: trivariate";
}

TEST(MultipolyRing, AdditionCommutativityCancellation) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{2, 0}, Rational(5)}, {{1, 1}, Rational(3)}}, vars);
    MultiPoly b({{{2, 0}, Rational(-5)}, {{0, 2}, Rational(1)}}, vars);
    EXPECT_TRUE((a + b == b + a)) << "commutativity: with cancellation";
}

TEST(MultipolyRing, AdditionCommutativityConstant) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a(Rational(42), vars);
    MultiPoly b(Rational(-17), vars);
    EXPECT_TRUE((a + b == b + a)) << "commutativity: constants";
}

TEST(MultipolyRing, AdditionCommutativityHighDegree) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{5, 0}, Rational(1)}, {{0, 5}, Rational(1)}}, vars);
    MultiPoly b({{{3, 2}, Rational(2)}, {{2, 3}, Rational(-2)}}, vars);
    EXPECT_TRUE((a + b == b + a)) << "commutativity: high degree";
}

TEST(MultipolyRing, AdditionCommutativityNegative) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(-7)}, {{0, 1}, Rational(-3)}, {{0, 0}, Rational(-1)}}, vars);
    MultiPoly b({{{1, 0}, Rational(7)}, {{0, 1}, Rational(3)}, {{0, 0}, Rational(1)}}, vars);
    EXPECT_TRUE((a + b == b + a)) << "commutativity: negatives (sum is zero)";
}

TEST(MultipolyRing, MultiplicationAssociativityLinear) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(1)}, {{0, 0}, Rational(1)}}, vars);  /**< 多项式 x + 1。 */
    MultiPoly b({{{0, 1}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); /**< 多项式 y - 1。 */
    MultiPoly c({{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}}, vars);  /**< 多项式 x + y。 */
    EXPECT_TRUE(((a * b) * c == a * (b * c))) << "associativity: (x+1)(y-1)(x+y)";
}

TEST(MultipolyRing, MultiplicationAssociativityConstant) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(2)}, {{0, 1}, Rational(3)}}, vars);  /**< 多项式 2x + 3y。 */
    MultiPoly b(Rational(5), vars);                                     /**< 多项式 5。 */
    MultiPoly c({{{1, 1}, Rational(1)}, {{0, 0}, Rational(-2)}}, vars); /**< 多项式 xy - 2。 */
    EXPECT_TRUE(((a * b) * c == a * (b * c))) << "associativity: constant factor";
}

TEST(MultipolyRing, MultiplicationAssociativityIdentity) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{2, 0}, Rational(1)}, {{0, 2}, Rational(-1)}}, vars); /**< 多项式 x^2 - y^2。 */
    MultiPoly one(Rational(1), vars);
    MultiPoly c({{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}}, vars); /**< 多项式 x + y。 */
    EXPECT_TRUE(((a * one) * c == a * (one * c))) << "associativity: identity element";
}

TEST(MultipolyRing, MultiplicationAssociativityQuadratic) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{2, 0}, Rational(1)}, {{0, 0}, Rational(1)}}, vars);  /**< 多项式 x^2 + 1。 */
    MultiPoly b({{{0, 2}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); /**< 多项式 y^2 - 1。 */
    MultiPoly c({{{1, 1}, Rational(1)}}, vars);                         /**< 多项式 xy。 */
    EXPECT_TRUE(((a * b) * c == a * (b * c))) << "associativity: quadratic factors";
}

TEST(MultipolyRing, MultiplicationAssociativityRational) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(1, 2)}, {{0, 0}, Rational(1, 3)}}, vars);
    MultiPoly b({{{0, 1}, Rational(2, 3)}, {{0, 0}, Rational(3, 4)}}, vars);
    MultiPoly c({{{1, 0}, Rational(1)}, {{0, 1}, Rational(-1)}}, vars);
    EXPECT_TRUE(((a * b) * c == a * (b * c))) << "associativity: rational coefficients";
}

TEST(MultipolyRing, MultiplicationAssociativityTrivariate) {
    std::vector<std::string> vars3 = {"x", "y", "z"};
    MultiPoly a({{{1, 0, 0}, Rational(1)}, {{0, 0, 1}, Rational(1)}}, vars3);  /**< 多项式 x + z。 */
    MultiPoly b({{{0, 1, 0}, Rational(1)}, {{0, 0, 0}, Rational(2)}}, vars3);  /**< 多项式 y + 2。 */
    MultiPoly c({{{1, 0, 0}, Rational(1)}, {{0, 1, 0}, Rational(-1)}}, vars3); /**< 多项式 x - y。 */
    EXPECT_TRUE(((a * b) * c == a * (b * c))) << "associativity: trivariate";
}

TEST(MultipolyRing, MultiplicationAssociativityZero) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(3)}, {{0, 1}, Rational(2)}}, vars);
    MultiPoly zero;
    MultiPoly c({{{1, 0}, Rational(1)}}, vars);
    EXPECT_TRUE(((a * zero) * c == a * (zero * c))) << "associativity: zero factor";
}

TEST(MultipolyRing, DistributivityLinear) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(2)}, {{0, 0}, Rational(1)}}, vars);  /**< 多项式 2x + 1。 */
    MultiPoly b({{{0, 1}, Rational(1)}, {{0, 0}, Rational(3)}}, vars);  /**< 多项式 y + 3。 */
    MultiPoly c({{{1, 0}, Rational(-1)}, {{0, 1}, Rational(2)}}, vars); /**< 多项式 -x + 2y。 */
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: (2x+1)*((y+3)+(-x+2y))";
}

TEST(MultipolyRing, DistributivityMonomial) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 1}, Rational(3)}}, vars);                         /**< 多项式 3xy。 */
    MultiPoly b({{{2, 0}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); /**< 多项式 x^2 - 1。 */
    MultiPoly c({{{0, 2}, Rational(2)}, {{1, 0}, Rational(1)}}, vars);  /**< 多项式 2y^2 + x。 */
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: monomial * sum";
}

TEST(MultipolyRing, DistributivityCancellation) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(5)}, {{0, 1}, Rational(7)}}, vars);
    MultiPoly b({{{2, 1}, Rational(3)}, {{1, 0}, Rational(1)}}, vars);
    MultiPoly c = -b;
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: b + c = 0";
}

TEST(MultipolyRing, DistributivityConstant) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a(Rational(7), vars);
    MultiPoly b({{{1, 0}, Rational(1)}, {{0, 1}, Rational(-2)}}, vars);
    MultiPoly c({{{2, 0}, Rational(3)}, {{0, 0}, Rational(4)}}, vars);
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: constant * sum";
}

TEST(MultipolyRing, DistributivityQuadratic) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{2, 0}, Rational(1)}, {{1, 1}, Rational(-1)}, {{0, 0}, Rational(2)}}, vars);
    MultiPoly b({{{0, 2}, Rational(1)}, {{1, 0}, Rational(3)}}, vars);
    MultiPoly c({{{2, 0}, Rational(-1)}, {{0, 1}, Rational(4)}}, vars);
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: quadratic * sum of quadratics";
}

TEST(MultipolyRing, DistributivityRational) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(1, 2)}, {{0, 0}, Rational(3, 4)}}, vars);
    MultiPoly b({{{0, 1}, Rational(2, 3)}, {{0, 0}, Rational(1, 5)}}, vars);
    MultiPoly c({{{1, 0}, Rational(4, 7)}, {{0, 1}, Rational(-1, 3)}}, vars);
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: rational coefficients";
}

TEST(MultipolyRing, DistributivityTrivariate) {
    std::vector<std::string> vars3 = {"x", "y", "z"};
    MultiPoly a({{{1, 0, 1}, Rational(1)}, {{0, 1, 0}, Rational(2)}}, vars3);
    MultiPoly b({{{1, 1, 0}, Rational(1)}, {{0, 0, 1}, Rational(-3)}}, vars3);
    MultiPoly c({{{0, 0, 2}, Rational(2)}, {{1, 0, 0}, Rational(1)}}, vars3);
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: trivariate";
}

TEST(MultipolyRing, DistributivityRight) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}}, vars); /**< 多项式 x + y。 */
    MultiPoly b({{{2, 0}, Rational(1)}}, vars);                        /**< 多项式 x^2。 */
    MultiPoly c({{{0, 2}, Rational(1)}}, vars);                        /**< 多项式 y^2。 */
    EXPECT_TRUE(((b + c) * a == b * a + c * a)) << "right distributive: (x^2+y^2)*(x+y)";
}

TEST(MultipolyRing, DistributivityHighDegree) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{3, 0}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); /**< 多项式 x^3 - 1。 */
    MultiPoly b({{{1, 0}, Rational(1)}, {{0, 1}, Rational(1)}}, vars);  /**< 多项式 x + y。 */
    MultiPoly c({{{0, 1}, Rational(1)}, {{0, 0}, Rational(-1)}}, vars); /**< 多项式 y - 1。 */
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: cubic * sum";
}

TEST(MultipolyRing, DistributivitySingleTerms) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly a({{{1, 0}, Rational(1)}}, vars); /**< 多项式 x。 */
    MultiPoly b({{{0, 1}, Rational(1)}}, vars); /**< 多项式 y。 */
    MultiPoly c({{{0, 0}, Rational(1)}}, vars); /**< 多项式 1。 */
    EXPECT_TRUE((a * (b + c) == a * b + a * c)) << "distributive: x*(y+1) == xy + x";
}

template <class Operation>
void expect_variable_mismatch(Operation operation, const char *message) {
    bool rejected = false;
    try {
        operation();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << message;
}

TEST(MultipolyRing, NamedVariableRingsAreDistinct) {
    const MultiPoly x({{{1}, Rational(1)}}, {"x"});
    const MultiPoly y({{{1}, Rational(1)}}, {"y"});
    const MultiPoly xy({{{1, 0}, Rational(1)}}, {"x", "y"});
    const MultiPoly yx({{{1, 0}, Rational(1)}}, {"y", "x"});
    const MultiPoly one_x(Rational(1), {"x"});
    const MultiPoly one_y(Rational(1), {"y"});
    for (const auto &operands : {std::pair{x, y}, std::pair{xy, yx},
                                 std::pair{x, xy}, std::pair{one_x, one_y}}) {
        const auto &lhs = operands.first;
        const auto &rhs = operands.second;
        EXPECT_FALSE((lhs == rhs)) << "different nonzero rings are not equal";
        EXPECT_TRUE((lhs != rhs)) << "inequality respects variable rings";
        expect_variable_mismatch([&] { (void)(lhs + rhs); }, "addition rejects different rings");
        expect_variable_mismatch([&] { (void)(lhs - rhs); }, "subtraction rejects different rings");
        expect_variable_mismatch([&] { (void)(lhs * rhs); }, "multiplication rejects different rings");
    }

    const MultiPoly common_y({{{0, 1}, Rational(1)}}, {"x", "y"});
    const std::map<std::string, Rational> point{{"x", Rational(1)}, {"y", Rational(2)}};
    EXPECT_TRUE(((xy + common_y).eval(point) == MultiPoly(Rational(3), {}))) << "common-ring x+y evaluates to 3 at x=1,y=2";
    EXPECT_TRUE(((xy * common_y).eval(point) == MultiPoly(Rational(2), {}))) << "common-ring x*y evaluates to 2 at x=1,y=2";
}

TEST(MultipolyRing, ZeroEmbedsInTheNonzeroRing) {
    const std::vector<std::string> vars{"x", "y"};
    const MultiPoly a({{{1, 0}, Rational(2)}, {{0, 2}, Rational(3)}},
                      vars, MonomialOrderType::Lex);
    const MultiPoly negative_a({{{1, 0}, Rational(-2)}, {{0, 2}, Rational(-3)}},
                               vars, MonomialOrderType::Lex);
    const MultiPoly zero;
    const MultiPoly named_zero(std::vector<MultiPoly::Term>{}, {"z"}, MonomialOrderType::Lex);
    for (const auto *z : {&zero, &named_zero}) {
        for (const auto &result : {*z + a, a + *z, a - *z}) {
            EXPECT_TRUE((result.variables() == vars && result == a)) << "zero additive identities retain the nonzero ring and value";
            EXPECT_TRUE((result.terms() == a.terms())) << "zero addition retains the nonzero ordering";
        }
        const auto difference = *z - a;
        EXPECT_TRUE((difference.variables() == vars && difference == negative_a)) << "zero-a retains the nonzero ring and negates its value";
        EXPECT_TRUE((difference.terms() == negative_a.terms())) << "zero subtraction retains the nonzero ordering";
        for (const auto &product : {*z * a, a * *z}) {
            EXPECT_TRUE((product.is_zero() && product.variables() == vars)) << "zero multiplication embeds in the nonzero ring";
        }
    }
    EXPECT_TRUE((zero == named_zero && named_zero == zero)) << "zero equality is ring independent";
    for (const auto &result : {named_zero + zero, named_zero - zero, named_zero * zero}) {
        EXPECT_TRUE((result.is_zero() && result.variables() == named_zero.variables())) << "two zero operands retain the left variable ring";
    }
}

TEST(MultipolyRing, MonomialOrderDoesNotChangeRingValues) {
    const std::vector<std::string> vars{"x", "y", "z"};
    const std::vector<MultiPoly::Term> terms{
        {{2, 0, 1}, Rational(2)}, {{1, 2, 0}, Rational(3)}, {{0, 4, 0}, Rational(5)}, {{0, 0, 0}, Rational(7)}};
    const MultiPoly lex(terms, vars, MonomialOrderType::Lex);
    for (const auto order : {MonomialOrderType::GrevLex, MonomialOrderType::DegLex,
                             MonomialOrderType::DegRevLex}) {
        const MultiPoly other(terms, vars, order);
        EXPECT_TRUE((lex == other && other == lex)) << "equal coefficients compare across orderings";
        auto changed_terms = terms;
        changed_terms[1].second = Rational(4);
        EXPECT_FALSE((lex == MultiPoly(changed_terms, vars, order))) << "different coefficients remain unequal across orderings";
        changed_terms = terms;
        changed_terms[1].first = Monomial{1, 1, 1};
        EXPECT_FALSE((lex == MultiPoly(changed_terms, vars, order))) << "different monomials remain unequal across orderings";
    }

    const std::vector<std::string> xy{"x", "y"};
    const std::vector<MultiPoly::Term> a_terms{{{1, 0}, Rational(1)}, {{0, 2}, Rational(1)}};
    const std::vector<MultiPoly::Term> square_terms{
        {{2, 0}, Rational(1)}, {{1, 2}, Rational(2)}, {{0, 4}, Rational(1)}};
    const MultiPoly a(a_terms, xy, MonomialOrderType::Lex);
    const MultiPoly b(a_terms, xy, MonomialOrderType::GrevLex);
    EXPECT_TRUE((a - b == MultiPoly())) << "subtraction cancels equal differently ordered values";
    for (const auto order : {MonomialOrderType::Lex, MonomialOrderType::GrevLex}) {
        const auto sum = order == MonomialOrderType::Lex ? a + b : b + a;
        const auto product = order == MonomialOrderType::Lex ? a * b : b * a;
        const MultiPoly expected_sum({{{1, 0}, Rational(2)}, {{0, 2}, Rational(2)}}, xy, order);
        const MultiPoly expected_product(square_terms, xy, order);
        EXPECT_TRUE((sum.terms() == expected_sum.terms() && sum.variables() == xy)) << "mixed-order addition returns the correct sum in the left ordering";
        EXPECT_TRUE((product.terms() == expected_product.terms() && product.variables() == xy)) << "mixed-order multiplication returns the correct product in the left ordering";
    }
}

} // namespace
