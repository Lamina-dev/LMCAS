#pragma once

#include "internal/integration_support.hpp"
#include "internal/exact_algebraic.hpp"

namespace LMCAS::detail::definite {
using Node = std::shared_ptr<const SymbolicNode>;
using Poly = Polynomial<Rational>;
using Algebraic = ExactRealAlgebraic;
inline constexpr const char* operation = "integrate.definite";

struct Boundary {
    int infinity = 0;
    Algebraic finite;
};

struct RationalFunction { Poly numerator; Poly denominator; };
using RationalResult = Result<std::optional<RationalFunction>>;

template<class T>
Result<T> undecided() {
    return Result<T>::failure(CasErrc::Inconclusive,
        "interval certificate is unavailable", operation);
}

template<class T>
Result<T> divergent() {
    return Result<T>::failure(CasErrc::DomainError,
        "ordinary improper integral does not exist",
        "integrate.definite.convergence");
}

inline bool is_undecided(const CasError& error) {
    return error.code == CasErrc::Inconclusive ||
        error.code == CasErrc::UnsupportedExpression;
}

std::optional<Rational> exact_number(const Node& node, bool allow_binary = false);
Node absolute_argument(const SymbolicNode& node);
int infinity_sign(const Node& node);
Result<Boundary> boundary_checked(const SymbolicExpr& expression, ComputationContext& context);
Result<int> compare_boundaries(const Boundary& a, const Boundary& b, ComputationContext& context);
Result<Poly> product_checked(const Poly& a, const Poly& b, ComputationContext& context);
RationalResult rational_checked(const Node& node, const std::string& var, ComputationContext& context);
Result<std::vector<Algebraic>> roots_checked(const Poly& polynomial, ComputationContext& context);
Result<std::optional<int>> polynomial_sign_checked(const Poly& polynomial,
    const Boundary& a, const Boundary& b, ComputationContext& context);
Result<int> rational_divergence_checked(RationalFunction& rational,
    const Boundary& a, const Boundary& b, ComputationContext& context);

Result<AssumptionContext> real_assumptions(const std::string& var, ComputationContext& context);
bool needs_segment_expression(const Node& node);
Result<SymbolicExpr> segment_expression_checked(const SymbolicExpr& expression,
    const std::string& var, const SymbolicExpr& lower, const SymbolicExpr& upper,
    ComputationContext& context);
Result<std::vector<SymbolicExpr>> definite_partition_checked(const SymbolicExpr& original,
    const std::string& var, const SymbolicExpr& lower, const SymbolicExpr& upper,
    ComputationContext& context);
Result<Tribool> certify_continuity_checked(const SymbolicExpr& expression,
    const std::string& var, const SymbolicExpr& lower, const SymbolicExpr& upper,
    ComputationContext& context);
Result<SymbolicExpr> definite_segment_checked(const SymbolicExpr& original,
    const SymbolicExpr& primitive, const std::string& var,
    const SymbolicExpr& a, const SymbolicExpr& b, ComputationContext& context);
}
