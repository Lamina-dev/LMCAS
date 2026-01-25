#include "cas.hpp"
#include "symbolic.hpp"
#include "value.hpp"
#include <iostream>

int main() {
    std::cout << "Prove: log_b(x) = ln(x)/ln(b)" << std::endl;

    // 1. 定义变量 x 和底数 b
    auto x = SymbolicExpr::variable("x");
    auto b = SymbolicExpr::variable("b");

    // 2. 构造左边 (LHS): log_b(x)
    // 使用我们新添加的 SymbolicExpr::log(value, base)
    auto lhs = SymbolicExpr::log(x, b);
    std::cout << "LHS: " << lhs->to_string() << std::endl;

    // 3. 构造右边 (RHS): ln(x) / ln(b)
    // 在系统内部表示为 ln(x) * (ln(b))^-1
    auto ln_x = SymbolicExpr::ln(x);
    auto ln_b = SymbolicExpr::ln(b);
    auto rhs = SymbolicExpr::multiply(ln_x, SymbolicExpr::power(ln_b, SymbolicExpr::number(-1)));
    std::cout << "RHS: " << rhs->to_string() << std::endl;

    // 4. 化简两者
    // 根据我们在 symbolic.cpp 中实现的 simplify()，Type::Log 会被自动转换为 ln(x)/ln(b)
    std::cout << "" << std::endl;
    auto lhs_simplified = lhs->simplify();
    std::cout << "LHS Simplified: " << lhs_simplified->to_string() << std::endl;

    std::cout << "" << std::endl;
    auto rhs_simplified = rhs->simplify();
    std::cout << "RHS Simplified: " << rhs_simplified->to_string() << std::endl;

    // 5. 验证相等性
    // 方法 1: 字符串比较
    bool string_match = (lhs_simplified->to_string() == rhs_simplified->to_string());
    
    // 方法 2: 作差判零 (LHS - RHS = 0?)
    auto diff = SymbolicExpr::add(lhs_simplified, SymbolicExpr::multiply(rhs_simplified, SymbolicExpr::number(-1)));
    auto diff_simplified = diff->simplify();
    std::cout << "Difference: " << diff_simplified->to_string() << std::endl;

    bool diff_is_zero = (diff_simplified->is_number() && diff_simplified->convert_rational() == ::Rational(0));

    if (string_match || diff_is_zero) {
        std::cout << "\nQ.E.D." << std::endl;
    } else {
        std::cout << "\nCAS could not verify the identity." << std::endl;
    }

    return 0;
}
