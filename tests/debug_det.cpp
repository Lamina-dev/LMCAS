#include <iostream>
#include <memory>
#include "../include/symbolic.hpp"

int main() {

    auto one = SymbolicExpr::number(1);
    auto zero = SymbolicExpr::number(0);

    std::cout << "one: " << one->to_string() << " type=" << one->get_type_name() << " is_zero=" << one->is_zero() << std::endl;
    std::cout << "zero: " << zero->to_string() << " type=" << zero->get_type_name() << " is_zero=" << zero->is_zero() << std::endl;

    auto prod = SymbolicExpr::multiply(one, one)->simplify();
    std::cout << "1*1 = " << prod->to_string() << " is_zero=" << prod->is_zero() << std::endl;

    auto prod0 = SymbolicExpr::multiply(zero, one)->simplify();
    std::cout << "0*1 = " << prod0->to_string() << " is_zero=" << prod0->is_zero() << std::endl;

    auto b = SymbolicExpr::number(0);
    auto c = SymbolicExpr::number(0);
    auto bc = SymbolicExpr::multiply(b, c)->simplify();
    std::cout << "0*0 = " << bc->to_string() << " is_zero=" << bc->is_zero() << std::endl;

    auto neg_bc = SymbolicExpr::multiply(bc, SymbolicExpr::number(-1))->simplify();
    std::cout << "(0*0)*(-1) = " << neg_bc->to_string() << " is_zero=" << neg_bc->is_zero() << std::endl;

    auto sum = SymbolicExpr::add(one, zero)->simplify();
    std::cout << "1 + 0 = " << sum->to_string() << " is_zero=" << sum->is_zero() << std::endl;

    auto m1_data = std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>{
        {SymbolicExpr::number(1), SymbolicExpr::number(0)},
        {SymbolicExpr::number(0), SymbolicExpr::number(1)}
    };
    auto m1 = SymbolicExpr::matrix(m1_data);

    std::cout << "\nMatrix: " << m1->to_string() << std::endl;

    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(m1->root);
    if (mat_node) {
        std::cout << "MatrixNode rows=" << mat_node->rows << " cols=" << mat_node->cols << std::endl;
        std::cout << "is_sparse=" << mat_node->is_sparse() << std::endl;
        for (size_t i = 0; i < mat_node->rows; i++) {
            for (size_t j = 0; j < mat_node->cols; j++) {
                auto elem = std::make_shared<SymbolicExpr>(mat_node->get(i,j));
                std::cout << "  get(" << i << "," << j << "): " << elem->to_string()
                          << " type=" << elem->get_type_name()
                          << " is_zero=" << elem->is_zero() << std::endl;
            }
        }
    }

    auto det1 = SymbolicExpr::determinant(m1);
    std::cout << "det: " << det1->to_string() << " is_zero=" << det1->is_zero() << std::endl;

    return 0;
}
