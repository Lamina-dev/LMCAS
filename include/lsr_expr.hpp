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

namespace lamina::lsr {

using Expr = SymbolicExpr;
using ExprPtr = std::shared_ptr<SymbolicExpr>;
using ExprResult = Result<ExprPtr>;

struct ExprMatchBinding {
    std::string name;
    ExprPtr value;
};

struct ExprMatch {
    bool matched = false;
    std::vector<ExprMatchBinding> bindings;
};

using ExprMatchResult = Result<ExprMatch>;

struct ApproxComplex {
    ApproxReal real;
    ApproxReal imag;

    bool is_finite() const noexcept {
        return real.is_finite() && imag.is_finite();
    }
};

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

LAMINA_API Result<EqvProfile> eqv_profile_from_name(const std::string& name);
LAMINA_API Result<void> set_eqv_profile(EqvOptions& options,
                                        const std::string& name);
LAMINA_API Result<void> set_eqv_budget(EqvOptions& options,
                                       std::size_t steps,
                                       std::size_t depth,
                                       std::size_t growth);

class LAMINA_API ExprSet {
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

class LAMINA_API NumberDomainSet {
public:
    explicit NumberDomainSet(NumberDomain domain) : domain_(domain) {}

    NumberDomain domain() const noexcept { return domain_; }
    const char* name() const noexcept;
    bool subset_of(const NumberDomainSet& other) const noexcept;

    Result<bool> contains(const ExprPtr& element) const;

private:
    NumberDomain domain_;
};

LAMINA_API ExprResult sym(const std::string& name);
LAMINA_API ExprResult parse_expr(const std::string& source);
LAMINA_API ExprResult integer(long long value);
LAMINA_API ExprResult integer(const BigInt& value);
LAMINA_API ExprResult rational(const Rational& value);
LAMINA_API ExprResult approx_real(double value);

LAMINA_API ExprResult pi();
LAMINA_API ExprResult e();
LAMINA_API ExprResult phi();

LAMINA_API ExprResult i();
LAMINA_API ExprResult I();
LAMINA_API ExprResult imaginary_unit();
LAMINA_API ExprResult complex(ExprPtr real, ExprPtr imag);
LAMINA_API ExprResult with_unit(const ExprPtr& value,
                                const std::string& unit,
                                ComputationContext& context);
LAMINA_API ExprResult convert_to_unit(const ExprPtr& quantity,
                                      const std::string& unit,
                                      ComputationContext& context);
LAMINA_API ExprResult strip_to_base_value(const ExprPtr& quantity,
                                          ComputationContext& context);
LAMINA_API ExprResult strip_to_display_value(const ExprPtr& quantity,
                                             ComputationContext& context);
LAMINA_API ExprResult finite_set(std::vector<ExprPtr> elements,
                                 ComputationContext& context);
LAMINA_API ExprResult interval(const ExprPtr& lower, const ExprPtr& upper,
                               bool lower_closed, bool upper_closed,
                               ComputationContext& context);
LAMINA_API ExprResult member(const ExprPtr& element, const ExprPtr& set,
                             ComputationContext& context);
LAMINA_API ExprResult add(const ExprPtr& lhs,
                          const ExprPtr& rhs,
                          ComputationContext& context);
LAMINA_API ExprResult add(const ExprPtr& lhs, const ExprPtr& rhs);
LAMINA_API ExprResult sub(const ExprPtr& lhs,
                          const ExprPtr& rhs,
                          ComputationContext& context);
LAMINA_API ExprResult sub(const ExprPtr& lhs, const ExprPtr& rhs);
LAMINA_API ExprResult mul(const ExprPtr& lhs,
                          const ExprPtr& rhs,
                          ComputationContext& context);
LAMINA_API ExprResult mul(const ExprPtr& lhs, const ExprPtr& rhs);
LAMINA_API ExprResult div(const ExprPtr& numerator,
                          const ExprPtr& denominator,
                          ComputationContext& context);
LAMINA_API ExprResult div(const ExprPtr& numerator,
                          const ExprPtr& denominator);
LAMINA_API ExprResult neg(const ExprPtr& expression,
                          ComputationContext& context);
LAMINA_API ExprResult neg(const ExprPtr& expression);
LAMINA_API ExprResult eq(const ExprPtr& lhs,
                         const ExprPtr& rhs,
                         ComputationContext& context);
LAMINA_API ExprResult eq(const ExprPtr& lhs, const ExprPtr& rhs);
LAMINA_API ExprResult sqrt(const ExprPtr& expression,
                           ComputationContext& context);
LAMINA_API ExprResult sqrt(const ExprPtr& expression);
LAMINA_API ExprResult pow(const ExprPtr& base,
                          const ExprPtr& exponent,
                          ComputationContext& context);
LAMINA_API ExprResult pow(const ExprPtr& base, const ExprPtr& exponent);
LAMINA_API ExprResult sin(const ExprPtr& expression,
                          ComputationContext& context);
LAMINA_API ExprResult sin(const ExprPtr& expression);
LAMINA_API ExprResult cos(const ExprPtr& expression,
                          ComputationContext& context);
LAMINA_API ExprResult cos(const ExprPtr& expression);
LAMINA_API ExprResult tan(const ExprPtr& expression,
                          ComputationContext& context);
LAMINA_API ExprResult tan(const ExprPtr& expression);
LAMINA_API ExprResult asin(const ExprPtr& expression,
                           ComputationContext& context);
LAMINA_API ExprResult asin(const ExprPtr& expression);
LAMINA_API ExprResult acos(const ExprPtr& expression,
                           ComputationContext& context);
LAMINA_API ExprResult acos(const ExprPtr& expression);
LAMINA_API ExprResult atan(const ExprPtr& expression,
                           ComputationContext& context);
LAMINA_API ExprResult atan(const ExprPtr& expression);
LAMINA_API ExprResult exp(const ExprPtr& expression,
                          ComputationContext& context);
LAMINA_API ExprResult exp(const ExprPtr& expression);
LAMINA_API ExprResult log(const ExprPtr& expression,
                          ComputationContext& context);
LAMINA_API ExprResult log(const ExprPtr& expression);
LAMINA_API ExprResult log10(const ExprPtr& expression,
                            ComputationContext& context);
LAMINA_API ExprResult log10(const ExprPtr& expression);
LAMINA_API ExprResult floor(const ExprPtr& expression,
                            ComputationContext& context);
LAMINA_API ExprResult floor(const ExprPtr& expression);
LAMINA_API ExprResult ceil(const ExprPtr& expression,
                           ComputationContext& context);
LAMINA_API ExprResult ceil(const ExprPtr& expression);
LAMINA_API ExprResult round(const ExprPtr& expression,
                            ComputationContext& context);
LAMINA_API ExprResult round(const ExprPtr& expression);
LAMINA_API ExprResult clamp(const ExprPtr& expression,
                            const ExprPtr& lower,
                            const ExprPtr& upper,
                            ComputationContext& context);
LAMINA_API ExprResult clamp(const ExprPtr& expression,
                            const ExprPtr& lower,
                            const ExprPtr& upper);
LAMINA_API ExprResult real(const ExprPtr& expression,
                           ComputationContext& context);
LAMINA_API ExprResult real(const ExprPtr& expression);
LAMINA_API ExprResult imag(const ExprPtr& expression,
                           ComputationContext& context);
LAMINA_API ExprResult imag(const ExprPtr& expression);
LAMINA_API ExprResult conj(const ExprPtr& expression,
                           ComputationContext& context);
LAMINA_API ExprResult conj(const ExprPtr& expression);
LAMINA_API ExprResult abs(const ExprPtr& expression,
                          ComputationContext& context);
LAMINA_API ExprResult abs(const ExprPtr& expression);
LAMINA_API ExprResult simplify(const ExprPtr& expression,
                               ComputationContext& context);
LAMINA_API ExprResult simplify(const ExprPtr& expression);
LAMINA_API ExprResult expand(const ExprPtr& expression,
                             ComputationContext& context);
LAMINA_API ExprResult expand(const ExprPtr& expression);
LAMINA_API ExprResult differentiate(const ExprPtr& expression,
                                    const std::string& variable,
                                    ComputationContext& context);
LAMINA_API ExprResult differentiate(const ExprPtr& expression,
                                    const std::string& variable);

LAMINA_API ExprResult substitute(const ExprPtr& expression,
                                 const std::string& variable,
                                 const ExprPtr& value,
                                 ComputationContext& context);

LAMINA_API ExprResult substitute(const ExprPtr& expression,
                                 const std::string& variable,
                                 const ExprPtr& value);

LAMINA_API ExprMatchResult expr_match(const ExprPtr& pattern,
                                      const ExprPtr& target,
                                      const std::vector<std::string>& wildcards,
                                      ComputationContext& context);

LAMINA_API ExprMatchResult expr_match(const ExprPtr& pattern,
                                      const ExprPtr& target,
                                      const std::vector<std::string>& wildcards);

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
LAMINA_API NumberDomainSet integers();
LAMINA_API NumberDomainSet rationals();
LAMINA_API NumberDomainSet reals();
LAMINA_API NumberDomainSet complexes();
LAMINA_API NumberDomainSet expressions();
LAMINA_API Result<bool> domain_contains(const NumberDomainSet& domain,
                                        const ExprPtr& element);
LAMINA_API Result<bool> domain_subset(const NumberDomainSet& lhs,
                                      const NumberDomainSet& rhs);
LAMINA_API Result<bool> expr_set_contains(const ExprSet& set,
                                          const ExprPtr& element);
LAMINA_API Result<bool> expr_set_not_contains(const ExprSet& set,
                                              const ExprPtr& element);
LAMINA_API Result<bool> expr_set_subset(const ExprSet& lhs,
                                        const ExprSet& rhs);
LAMINA_API Result<bool> expr_set_subset_domain(const ExprSet& set,
                                               const NumberDomainSet& domain);
LAMINA_API ExprSetResult expr_set_union(const ExprSet& lhs,
                                        const ExprSet& rhs);
LAMINA_API ExprSetResult expr_set_intersection(const ExprSet& lhs,
                                               const ExprSet& rhs);
LAMINA_API ExprSetResult expr_set_difference(const ExprSet& lhs,
                                             const ExprSet& rhs);
LAMINA_API ExprSetResult expr_set_symmetric_difference(const ExprSet& lhs,
                                                       const ExprSet& rhs);

LAMINA_API ExprSetResult solve_expr_set(const ExprPtr& equation,
                                        const std::string& variable,
                                        ComputationContext& context,
                                        const SolveOptions& options = {});

LAMINA_API ExprSetResult solve_expr_set(const ExprPtr& equation,
                                        const std::string& variable,
                                        const SolveOptions& options = {});

LAMINA_API ExprSetResult roots(const ExprPtr& expression,
                               const std::string& variable,
                               ComputationContext& context,
                               const SolveOptions& options = {});

LAMINA_API ExprSetResult roots(const ExprPtr& expression,
                               const std::string& variable,
                               const SolveOptions& options = {});

LAMINA_API ExprSetResult solve(const ExprPtr& equation,
                               const std::string& variable,
                               ComputationContext& context,
                               const SolveOptions& options = {});

LAMINA_API ExprSetResult solve(const ExprPtr& equation,
                               const std::string& variable,
                               const SolveOptions& options = {});

LAMINA_API const char* error_name(CasErrc code) noexcept;
LAMINA_API const char* error_name(const CasError& error) noexcept;

LAMINA_API bool structurally_equal(const SymbolicExpr& lhs,
                                   const SymbolicExpr& rhs);

LAMINA_API Result<bool> equivalent_core(const SymbolicExpr& lhs,
                                        const SymbolicExpr& rhs,
                                        ComputationContext& context);

LAMINA_API Result<bool> equivalent_core(const SymbolicExpr& lhs,
                                        const SymbolicExpr& rhs,
                                        ComputationContext& context,
                                        const EqvOptions& options);

LAMINA_API Result<bool> equivalent(const SymbolicExpr& lhs,
                                   const SymbolicExpr& rhs,
                                   ComputationContext& context);

LAMINA_API Result<bool> equivalent(const SymbolicExpr& lhs,
                                   const SymbolicExpr& rhs,
                                   ComputationContext& context,
                                   const EqvOptions& options);

} // namespace lamina::lsr
