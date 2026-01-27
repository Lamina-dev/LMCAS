#pragma once
#include "symbolic.hpp"
#include <vector>
#include <memory>

// Shared internal helper functions

// Check if two expressions are structurally equal
bool symbolic_equal(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

// Helper for matrix addition
std::shared_ptr<SymbolicExpr> add_matrices(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

// Helper for matrix multiplication
std::shared_ptr<SymbolicExpr> multiply_matrices(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

// Helper for distribution
std::shared_ptr<SymbolicExpr> distribute_multiply(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);
