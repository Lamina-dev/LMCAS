#include "test_common.hpp"

int main() {
    TEST_CASE("Matrix Determinant");

    {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m1_data = {
            {SymbolicExpr::number(1), SymbolicExpr::number(0)},
            {SymbolicExpr::number(0), SymbolicExpr::number(1)}
        };
        auto m1 = SymbolicExpr::matrix(m1_data);
        auto det1 = SymbolicExpr::determinant(m1);
        EXPECT_EQ_EXPR(det1, SymbolicExpr::number(1), "Det(I_2)");
    }

    {
        auto a = SymbolicExpr::variable("a");
        auto b = SymbolicExpr::variable("b");
        auto c = SymbolicExpr::variable("c");
        auto d = SymbolicExpr::variable("d");

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> m2_data = {
            {a, b},
            {c, d}
        };
        auto m2 = SymbolicExpr::matrix(m2_data);
        auto det2 = SymbolicExpr::determinant(m2);

        auto ad = SymbolicExpr::multiply(a, d);
        auto bc = SymbolicExpr::multiply(b, c);
        auto neg_bc = SymbolicExpr::multiply(bc, SymbolicExpr::number(-1));
        auto expected2 = SymbolicExpr::add(ad, neg_bc)->simplify();

        EXPECT_EQ_EXPR(det2, expected2, "Det(Symbolic 2x2)");
    }

    return TEST_REPORT();
}
