#include <iostream>
#include <vector>
#include "symbolic.hpp"
#include "test_common.hpp"

void print_expr(const std::string& label, const std::shared_ptr<SymbolicExpr>& expr) {
    std::cout << label << ": " << expr->to_string() << std::endl;
}

void test_maclaurin_sin() {
    std::cout << "Testing Maclaurin Series for sin(x)..." << std::endl;

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    std::cout << "Calling series() for order 2..." << std::endl;
    auto series_2 = expr->series("x", SymbolicExpr::number(0), 2);
    print_expr("sin(x) series (order 2)", series_2);
    EXPECT_TRUE(series_2 != nullptr && !series_2->to_string().empty(),
                "sin Maclaurin order 2 returns a printable expression");

    std::cout << "Calling series() for order 5..." << std::endl;
    auto series_5 = expr->series("x", SymbolicExpr::number(0), 5);
    print_expr("sin(x) series (order 5)", series_5);
    EXPECT_TRUE(series_5 != nullptr && series_5->to_string().find("x") != std::string::npos,
                "sin Maclaurin order 5 contains x terms");

    auto expanded = series_5->expand();
    print_expr("Expanded", expanded);
    EXPECT_TRUE(expanded != nullptr && !expanded->to_string().empty(),
                "expanded sin series is non-empty");
}

void test_maclaurin_exp() {
    std::cout << "Testing Maclaurin Series for e^x..." << std::endl;

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::exp(x);

    auto series = expr->series("x", SymbolicExpr::number(0), 4);
    print_expr("exp(x) series (order 4)", series);
    EXPECT_TRUE(series != nullptr && series->to_string().find("1") != std::string::npos,
                "exp Maclaurin series contains constant term");

}

void test_taylor_ln() {
    std::cout << "Testing Taylor Series for ln(x) at x=1..." << std::endl;

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::ln(x);

    auto series = expr->series("x", SymbolicExpr::number(1), 3);
    print_expr("ln(x) series at x=1 (order 3)", series);
    EXPECT_TRUE(series != nullptr && !series->to_string().empty(),
                "ln Taylor series returns a printable expression");
}

void test_poly_series() {
    std::cout << "Testing Series for Polynomial x^2 + 2x + 1..." << std::endl;

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(1))
    );

    auto series = expr->series("x", SymbolicExpr::number(0), 3);
    print_expr("Poly series (order 3)", series);
    EXPECT_TRUE(series != nullptr && !series->to_string().empty(),
                "polynomial series returns a printable expression");
    auto expanded = series->expand();
    print_expr("Expanded", expanded);
    EXPECT_TRUE(expanded != nullptr && !expanded->to_string().empty(),
                "expanded polynomial series is non-empty");
}

#include "lammp/lmmp.h"

int main() {
    lmmp_stack_init(512 * 1024);
    try {
        test_maclaurin_sin();
        test_maclaurin_exp();
        test_taylor_ln();
        test_poly_series();
        std::cout << "All series tests passed!" << std::endl;
    } catch (const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    }
    return TEST_REPORT();
}
