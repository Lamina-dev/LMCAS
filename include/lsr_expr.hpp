#pragma once

#include <memory>
#include <string>
#include <vector>

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

struct ApproxComplex {
    ApproxReal real;
    ApproxReal imag;

    bool is_finite() const noexcept {
        return real.is_finite() && imag.is_finite();
    }
};

class LAMINA_API ExprSet {
public:
    ExprSet() = default;

    static Result<ExprSet> make(std::vector<ExprPtr> elements);

    bool empty() const noexcept { return elements_.empty(); }
    std::size_t size() const noexcept { return elements_.size(); }
    const std::vector<ExprPtr>& elements() const noexcept { return elements_; }

    bool contains(const SymbolicExpr& expression) const;
    bool subset_of(const ExprSet& other) const;

    ExprSet set_union(const ExprSet& other) const;
    ExprSet intersection(const ExprSet& other) const;
    ExprSet difference(const ExprSet& other) const;
    ExprSet symmetric_difference(const ExprSet& other) const;

private:
    explicit ExprSet(std::vector<ExprPtr> elements);

    std::vector<ExprPtr> elements_;
};

using ExprSetResult = Result<ExprSet>;

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

LAMINA_API Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                              const NumericBindings& bindings,
                                              ComputationContext& context);

LAMINA_API Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                              const NumericBindings& bindings = {});

LAMINA_API SolveResult solve_set(const ExprPtr& equation,
                                 const std::string& variable,
                                 ComputationContext& context,
                                 const SolveOptions& options = {});

LAMINA_API SolveResult solve_set(const ExprPtr& equation,
                                 const std::string& variable,
                                 const SolveOptions& options = {});

LAMINA_API ExprSetResult expr_set(std::vector<ExprPtr> elements);

LAMINA_API ExprSetResult solve_expr_set(const ExprPtr& equation,
                                        const std::string& variable,
                                        ComputationContext& context,
                                        const SolveOptions& options = {});

LAMINA_API ExprSetResult solve_expr_set(const ExprPtr& equation,
                                        const std::string& variable,
                                        const SolveOptions& options = {});

LAMINA_API bool structurally_equal(const SymbolicExpr& lhs,
                                   const SymbolicExpr& rhs);

LAMINA_API Result<bool> equivalent_core(const SymbolicExpr& lhs,
                                        const SymbolicExpr& rhs,
                                        ComputationContext& context);

} // namespace lamina::lsr
