#pragma once

#include <memory>
#include <string>

#include "bigint.hpp"
#include "conditional_result.hpp"
#include "numeric_evaluation.hpp"
#include "rational.hpp"
#include "solve_strategies.hpp"
#include "symbolic.hpp"

namespace lamina::lsr {

using Expr = SymbolicExpr;
using ExprPtr = std::shared_ptr<SymbolicExpr>;
using ExprResult = Result<ExprPtr>;

LAMINA_API ExprResult sym(const std::string& name);
LAMINA_API ExprResult integer(long long value);
LAMINA_API ExprResult integer(const BigInt& value);
LAMINA_API ExprResult rational(const Rational& value);
LAMINA_API ExprResult approx_real(double value);

LAMINA_API ExprResult imaginary_unit();
LAMINA_API ExprResult complex(ExprPtr real, ExprPtr imag);

LAMINA_API Result<ApproxReal> evalf(const SymbolicExpr& expression,
                                    const NumericBindings& bindings,
                                    ComputationContext& context);

LAMINA_API Result<ApproxReal> evalf(const SymbolicExpr& expression,
                                    const NumericBindings& bindings = {});

LAMINA_API SolveResult solve_set(const ExprPtr& equation,
                                 const std::string& variable,
                                 ComputationContext& context,
                                 const SolveOptions& options = {});

LAMINA_API SolveResult solve_set(const ExprPtr& equation,
                                 const std::string& variable,
                                 const SolveOptions& options = {});

LAMINA_API bool structurally_equal(const SymbolicExpr& lhs,
                                   const SymbolicExpr& rhs);

LAMINA_API Result<bool> equivalent_core(const SymbolicExpr& lhs,
                                        const SymbolicExpr& rhs,
                                        ComputationContext& context);

} // namespace lamina::lsr
