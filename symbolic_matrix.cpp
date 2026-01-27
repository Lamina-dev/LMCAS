#include "symbolic.hpp"
#include "symbolic_internal.hpp"
#include <vector>

// Helper for matrix addition
std::shared_ptr<SymbolicExpr> add_matrices(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    if (a->operands.size() != b->operands.size()) return nullptr; // Dimension mismatch
    if (a->operands.empty()) return a;
    if (a->operands[0]->operands.size() != b->operands[0]->operands.size()) return nullptr;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
    for (size_t i = 0; i < a->operands.size(); ++i) {
        std::vector<std::shared_ptr<SymbolicExpr>> row;
        auto vec_a = a->operands[i];
        auto vec_b = b->operands[i];
        for (size_t j = 0; j < vec_a->operands.size(); ++j) {
            row.push_back(SymbolicExpr::add(vec_a->operands[j], vec_b->operands[j])->simplify());
        }
        res_mat.push_back(row);
    }
    return SymbolicExpr::matrix(res_mat);
}

std::shared_ptr<SymbolicExpr> multiply_matrices(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    if (a->type == SymbolicExpr::Type::Matrix && b->type == SymbolicExpr::Type::Matrix) {
        if (a->operands.empty() || b->operands.empty()) return nullptr;
        
        size_t rowsA = a->operands.size();
        size_t colsA = 0;
        if (!a->operands[0]->operands.empty()) colsA = a->operands[0]->operands.size();
        
        size_t rowsB = b->operands.size();
        size_t colsB = 0;
        if(!b->operands.empty() && !b->operands[0]->operands.empty()) colsB = b->operands[0]->operands.size();
        
        if (colsA != rowsB) return nullptr;

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
        res_mat.resize(rowsA);

        for(size_t i=0; i<rowsA; ++i) {
            for(size_t j=0; j<colsB; ++j) {
                std::shared_ptr<SymbolicExpr> sum = SymbolicExpr::number(0);
                for(size_t k=0; k<colsA; ++k) {
                    auto term = SymbolicExpr::multiply(a->operands[i]->operands[k], b->operands[k]->operands[j]);
                    sum = SymbolicExpr::add(sum, term);
                }
                auto sim = sum->simplify();
                res_mat[i].push_back(sim);
            }
        }
        return SymbolicExpr::matrix(res_mat);
    }
    auto scalar_mul = [](const std::shared_ptr<SymbolicExpr>& scalar, const std::shared_ptr<SymbolicExpr>& mat) {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> res_mat;
        for(const auto& row_vec : mat->operands) {
             std::vector<std::shared_ptr<SymbolicExpr>> new_row;
             for(const auto& elem : row_vec->operands) {
                 new_row.push_back(SymbolicExpr::multiply(scalar, elem)->simplify());
             }
             res_mat.push_back(new_row);
        }
        return SymbolicExpr::matrix(res_mat);
    };
    if (a->type != SymbolicExpr::Type::Matrix && b->type == SymbolicExpr::Type::Matrix) return scalar_mul(a, b);
    if (a->type == SymbolicExpr::Type::Matrix && b->type != SymbolicExpr::Type::Matrix) return scalar_mul(b, a);
    return nullptr;
}
