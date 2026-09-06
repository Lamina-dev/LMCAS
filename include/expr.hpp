#pragma once

#include <memory>
#include <string>
#include <vector>

#include "bigint.hpp"
#include "conditional_result.hpp"
#include "numeric_evaluation.hpp"
#include "rational.hpp"
#include "solve_strategies.hpp"
#include "quantity.hpp"
#include "symbolic_set.hpp"
#include "symbolic.hpp"

namespace LMCAS {
class ComputationContext;
}

namespace LMCAS {

using Expr = SymbolicExpr;
using ExprResult = Result<ExprPtr>;

struct Binding {
    ExprPtr symbol;
    ExprPtr value;
};

using BindingResult = Result<Binding>;

struct ExprMatchBinding {
    std::string name;
    ExprPtr value;
};

struct ExprMatch {
    bool matched = false;
    std::vector<ExprMatchBinding> bindings;
};

using ExprMatchResult = Result<ExprMatch>;


enum class EqvProfile {
    Core,
    TrigBasic,
    ExpLogBasic
};

struct EqvBudget {
    std::size_t max_rewrite_steps = 256;
    std::size_t max_rewrite_depth = 64;
    std::size_t max_node_growth_factor = 4;
};

struct EqvOptions {
    EqvProfile profile = EqvProfile::Core;
    EqvBudget budget = {};
};

LMCAS_API Result<EqvProfile> eqv_profile_from_name(const std::string& name);
LMCAS_API Result<void> set_eqv_profile(EqvOptions& options,
                                        const std::string& name);
LMCAS_API Result<void> set_eqv_budget(EqvOptions& options,
                                       std::size_t steps,
                                       std::size_t depth,
                                       std::size_t growth);

class LMCAS_API ExprSet {
public:
    ExprSet() = default;

    static Result<ExprSet> make(std::vector<ExprPtr> elements);

    bool empty() const noexcept { return elements_.empty(); }
    std::size_t size() const noexcept { return elements_.size(); }
    const std::vector<ExprPtr>& elements() const noexcept { return elements_; }
    const ExprPtr& expression() const noexcept { return expression_; }

    bool contains(const SymbolicExpr& expression) const;
    bool subset_of(const ExprSet& other) const;

    ExprSet set_union(const ExprSet& other) const;
    ExprSet intersection(const ExprSet& other) const;
    ExprSet difference(const ExprSet& other) const;
    ExprSet symmetric_difference(const ExprSet& other) const;

private:
    ExprSet(std::vector<ExprPtr> elements, ExprPtr expression);

    std::vector<ExprPtr> elements_;
    ExprPtr expression_;
};

using ExprSetResult = Result<ExprSet>;

enum class NumberDomain {
    Integers,
    Rationals,
    Reals,
    Complexes,
    Expressions
};

class LMCAS_API NumberDomainSet {
public:
    explicit NumberDomainSet(NumberDomain domain) : domain_(domain) {}

    NumberDomain domain() const noexcept { return domain_; }
    const char* name() const noexcept;
    bool subset_of(const NumberDomainSet& other) const noexcept;

    Result<bool> contains(const ExprPtr& element) const;

private:
    NumberDomain domain_;
};

LMCAS_API ExprResult sym(const std::string& name);
LMCAS_API ExprResult parse_expr(const std::string& source);
LMCAS_API ExprResult parse_expr(const std::string& source,
                                ComputationContext& context);
LMCAS_API ExprResult integer(long long value);
LMCAS_API ExprResult integer(const BigInt& value);
LMCAS_API ExprResult rational(const Rational& value);
LMCAS_API ExprResult approx_real(double value);

LMCAS_API ExprResult pi();
LMCAS_API ExprResult e();
LMCAS_API ExprResult phi();

LMCAS_API ExprResult I();
LMCAS_API ExprResult imaginary_unit();
LMCAS_API ExprResult complex(ExprPtr real, ExprPtr imag);
LMCAS_API ExprResult function(const std::string& name,
                               std::vector<ExprPtr> arguments);
LMCAS_API ExprResult finite_set(std::vector<ExprPtr> elements);
LMCAS_API ExprResult interval(ExprPtr lower, ExprPtr upper,
                               bool lower_closed, bool upper_closed);
LMCAS_API ExprResult relation(const ExprPtr& lhs, const ExprPtr& rhs,
                               RelationOp op);
LMCAS_API ExprResult ne(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult lt(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult le(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult gt(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult ge(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult logical_and(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult logical_or(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult logical_not(const ExprPtr& expression);
LMCAS_API ExprResult membership(const ExprPtr& element,
                                 const ExprPtr& set, bool negated = false);
LMCAS_API ExprResult with_unit(const ExprPtr& value,
                                const std::string& unit,
                                ComputationContext& context);
LMCAS_API ExprResult with_unit_definition(const ExprPtr& value,
                                           std::string display_unit,
                                           UnitDefinition definition,
                                           ComputationContext& context);
LMCAS_API ExprResult convert_to_unit(const ExprPtr& quantity,
                                      const std::string& unit,
                                      ComputationContext& context);
LMCAS_API ExprResult convert_to_unit_definition(const ExprPtr& quantity,
                                                 std::string display_unit,
                                                 UnitDefinition definition,
                                                 ComputationContext& context);
LMCAS_API ExprResult strip_to_base_value(const ExprPtr& quantity,
                                          ComputationContext& context);
LMCAS_API ExprResult strip_to_display_value(const ExprPtr& quantity,
                                             ComputationContext& context);
LMCAS_API ExprResult finite_set(std::vector<ExprPtr> elements,
                                 ComputationContext& context);
LMCAS_API ExprResult interval(const ExprPtr& lower, const ExprPtr& upper,
                               bool lower_closed, bool upper_closed,
                               ComputationContext& context);
LMCAS_API ExprResult member(const ExprPtr& element, const ExprPtr& set,
                             ComputationContext& context);
LMCAS_API ExprResult add(const ExprPtr& lhs,
                          const ExprPtr& rhs,
                          ComputationContext& context);
LMCAS_API ExprResult add(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult sub(const ExprPtr& lhs,
                          const ExprPtr& rhs,
                          ComputationContext& context);
LMCAS_API ExprResult sub(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult mul(const ExprPtr& lhs,
                          const ExprPtr& rhs,
                          ComputationContext& context);
LMCAS_API ExprResult mul(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult div(const ExprPtr& numerator,
                          const ExprPtr& denominator,
                          ComputationContext& context);
LMCAS_API ExprResult div(const ExprPtr& numerator,
                          const ExprPtr& denominator);
LMCAS_API ExprResult neg(const ExprPtr& expression,
                          ComputationContext& context);
LMCAS_API ExprResult neg(const ExprPtr& expression);
LMCAS_API ExprResult eq(const ExprPtr& lhs,
                         const ExprPtr& rhs,
                         ComputationContext& context);
LMCAS_API ExprResult eq(const ExprPtr& lhs, const ExprPtr& rhs);
LMCAS_API ExprResult sqrt(const ExprPtr& expression,
                           ComputationContext& context);
LMCAS_API ExprResult sqrt(const ExprPtr& expression);
LMCAS_API ExprResult pow(const ExprPtr& base,
                          const ExprPtr& exponent,
                          ComputationContext& context);
LMCAS_API ExprResult pow(const ExprPtr& base, const ExprPtr& exponent);
LMCAS_API ExprResult sin(const ExprPtr& expression,
                          ComputationContext& context);
LMCAS_API ExprResult sin(const ExprPtr& expression);
LMCAS_API ExprResult cos(const ExprPtr& expression,
                          ComputationContext& context);
LMCAS_API ExprResult cos(const ExprPtr& expression);
LMCAS_API ExprResult tan(const ExprPtr& expression,
                          ComputationContext& context);
LMCAS_API ExprResult tan(const ExprPtr& expression);
LMCAS_API ExprResult asin(const ExprPtr& expression,
                           ComputationContext& context);
LMCAS_API ExprResult asin(const ExprPtr& expression);
LMCAS_API ExprResult acos(const ExprPtr& expression,
                           ComputationContext& context);
LMCAS_API ExprResult acos(const ExprPtr& expression);
LMCAS_API ExprResult atan(const ExprPtr& expression,
                           ComputationContext& context);
LMCAS_API ExprResult atan(const ExprPtr& expression);
LMCAS_API ExprResult exp(const ExprPtr& expression,
                          ComputationContext& context);
LMCAS_API ExprResult exp(const ExprPtr& expression);
LMCAS_API ExprResult log(const ExprPtr& expression,
                          ComputationContext& context);
LMCAS_API ExprResult log(const ExprPtr& expression);
LMCAS_API ExprResult log10(const ExprPtr& expression,
                            ComputationContext& context);
LMCAS_API ExprResult log10(const ExprPtr& expression);
LMCAS_API ExprResult floor(const ExprPtr& expression,
                            ComputationContext& context);
LMCAS_API ExprResult floor(const ExprPtr& expression);
LMCAS_API ExprResult ceil(const ExprPtr& expression,
                           ComputationContext& context);
LMCAS_API ExprResult ceil(const ExprPtr& expression);
LMCAS_API ExprResult round(const ExprPtr& expression,
                            ComputationContext& context);
LMCAS_API ExprResult round(const ExprPtr& expression);
LMCAS_API ExprResult clamp(const ExprPtr& expression,
                            const ExprPtr& lower,
                            const ExprPtr& upper,
                            ComputationContext& context);
LMCAS_API ExprResult clamp(const ExprPtr& expression,
                            const ExprPtr& lower,
                            const ExprPtr& upper);
LMCAS_API ExprResult real(const ExprPtr& expression,
                           ComputationContext& context);
LMCAS_API ExprResult real(const ExprPtr& expression);
LMCAS_API ExprResult imag(const ExprPtr& expression,
                           ComputationContext& context);
LMCAS_API ExprResult imag(const ExprPtr& expression);
LMCAS_API ExprResult conj(const ExprPtr& expression,
                           ComputationContext& context);
LMCAS_API ExprResult conj(const ExprPtr& expression);
/**
 * @brief Construct the complex modulus in square-sum form.
 *
 * Exact numeric components retain exact simplification. Approximate components
 * and bound real components use a scaled norm without forming their squares.
 */
LMCAS_API ExprResult abs(const ExprPtr& expression,
                          ComputationContext& context);
LMCAS_API ExprResult abs(const ExprPtr& expression);
LMCAS_API ExprResult simplify(const ExprPtr& expression,
                               ComputationContext& context);
LMCAS_API ExprResult simplify(const ExprPtr& expression);
LMCAS_API ExprResult expand(const ExprPtr& expression,
                             ComputationContext& context);
LMCAS_API ExprResult expand(const ExprPtr& expression);
LMCAS_API ExprResult differentiate(const ExprPtr& expression,
                                    const std::string& variable,
                                    ComputationContext& context);
LMCAS_API ExprResult differentiate(const ExprPtr& expression,
                                    const std::string& variable);

LMCAS_API ExprResult substitute(const ExprPtr& expression,
                                 const std::string& variable,
                                 const ExprPtr& value,
                                 ComputationContext& context);

LMCAS_API ExprResult substitute(const ExprPtr& expression,
                                 const std::string& variable,
                                 const ExprPtr& value);

LMCAS_API BindingResult binding(ExprPtr symbol, ExprPtr value);

LMCAS_API ExprResult substitute(const ExprPtr& expression,
                                 const Binding& binding,
                                 ComputationContext& context);

LMCAS_API ExprResult substitute(const ExprPtr& expression,
                                 const Binding& binding);

LMCAS_API ExprResult substitute(const ExprPtr& expression,
                                 const std::vector<Binding>& bindings,
                                 ComputationContext& context);

LMCAS_API ExprResult substitute(const ExprPtr& expression,
                                 const std::vector<Binding>& bindings);

LMCAS_API ExprMatchResult expr_match(const ExprPtr& pattern,
                                      const ExprPtr& target,
                                      const std::vector<std::string>& wildcards,
                                      ComputationContext& context);

LMCAS_API ExprMatchResult expr_match(const ExprPtr& pattern,
                                      const ExprPtr& target,
                                      const std::vector<std::string>& wildcards);

LMCAS_API Result<ApproxReal> evalf(const SymbolicExpr& expression,
                                    const NumericBindings& bindings,
                                    ComputationContext& context);

LMCAS_API Result<ApproxReal> evalf(const SymbolicExpr& expression,
                                    const NumericBindings& bindings = {});

/**
 * @brief Explicitly evaluate an expression as a finite complex value.
 *
 * Powers with integer exponents in [-64, 64] use exponentiation by squaring.
 * Exact exponent expressions are normalized before numeric conversion, with
 * rational exponents classified in their exact number domain.
 * Negative powers first form the reciprocal base, keeping the computation
 * on the magnitude scale of the requested result.
 * Complex traversal and real-component evaluation share the context's active
 * recursion budget. Each frame releases its depth on return or stack unwinding.
 * Complex products use LMMC's component-wise binary scaling at exponent
 * boundaries and retain the propagated component error estimates.
 */
LMCAS_API Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                              const NumericBindings& bindings,
                                              ComputationContext& context);

LMCAS_API Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                              const NumericBindings& bindings = {});

LMCAS_API SolveResult solve_set(const ExprPtr& equation,
                                 const std::string& variable,
                                 ComputationContext& context,
                                 const SolveOptions& options = {});

LMCAS_API SolveResult solve_set(const ExprPtr& equation,
                                 const std::string& variable,
                                 const SolveOptions& options = {});

LMCAS_API ExprSetResult expr_set(std::vector<ExprPtr> elements);
LMCAS_API NumberDomainSet integers();
LMCAS_API NumberDomainSet rationals();
LMCAS_API NumberDomainSet reals();
LMCAS_API NumberDomainSet complexes();
LMCAS_API NumberDomainSet expressions();
LMCAS_API Result<bool> domain_contains(const NumberDomainSet& domain,
                                        const ExprPtr& element);
LMCAS_API Result<bool> domain_subset(const NumberDomainSet& lhs,
                                      const NumberDomainSet& rhs);
LMCAS_API Result<bool> expr_set_contains(const ExprSet& set,
                                          const ExprPtr& element);
LMCAS_API Result<bool> expr_set_not_contains(const ExprSet& set,
                                              const ExprPtr& element);
LMCAS_API Result<bool> expr_set_subset(const ExprSet& lhs,
                                        const ExprSet& rhs);
LMCAS_API Result<bool> expr_set_subset_domain(const ExprSet& set,
                                               const NumberDomainSet& domain);
LMCAS_API ExprSetResult expr_set_union(const ExprSet& lhs,
                                        const ExprSet& rhs);
LMCAS_API ExprSetResult expr_set_intersection(const ExprSet& lhs,
                                               const ExprSet& rhs);
LMCAS_API ExprSetResult expr_set_difference(const ExprSet& lhs,
                                             const ExprSet& rhs);
LMCAS_API ExprSetResult expr_set_symmetric_difference(const ExprSet& lhs,
                                                       const ExprSet& rhs);

LMCAS_API ExprSetResult solve_expr_set(const ExprPtr& equation,
                                        const std::string& variable,
                                        ComputationContext& context,
                                        const SolveOptions& options = {});

LMCAS_API ExprSetResult solve_expr_set(const ExprPtr& equation,
                                        const std::string& variable,
                                        const SolveOptions& options = {});

LMCAS_API ExprSetResult roots(const ExprPtr& expression,
                               const std::string& variable,
                               ComputationContext& context,
                               const SolveOptions& options = {});

LMCAS_API ExprSetResult roots(const ExprPtr& expression,
                               const std::string& variable,
                               const SolveOptions& options = {});

LMCAS_API ExprSetResult solve(const ExprPtr& equation,
                               const std::string& variable,
                               ComputationContext& context,
                               const SolveOptions& options = {});

LMCAS_API ExprSetResult solve(const ExprPtr& equation,
                               const std::string& variable,
                               const SolveOptions& options = {});

LMCAS_API const char* error_name(CasErrc code) noexcept;
LMCAS_API const char* error_name(const CasError& error) noexcept;

LMCAS_API bool structurally_equal(const SymbolicExpr& lhs,
                                   const SymbolicExpr& rhs);

LMCAS_API Result<bool> equivalent_core(const SymbolicExpr& lhs,
                                        const SymbolicExpr& rhs,
                                        ComputationContext& context);

LMCAS_API Result<bool> equivalent_core(const SymbolicExpr& lhs,
                                        const SymbolicExpr& rhs,
                                        ComputationContext& context,
                                        const EqvOptions& options);

LMCAS_API Result<bool> equivalent(const SymbolicExpr& lhs,
                                   const SymbolicExpr& rhs,
                                   ComputationContext& context);

LMCAS_API Result<bool> equivalent(const SymbolicExpr& lhs,
                                   const SymbolicExpr& rhs,
                                   ComputationContext& context,
                                   const EqvOptions& options);

} // namespace LMCAS
