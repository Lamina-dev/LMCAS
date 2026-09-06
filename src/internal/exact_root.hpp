#pragma once

#include "internal/exact_root_id.hpp"
#include "exact_algebraic.hpp"
#include "numeric_evaluation.hpp"
#include "result.hpp"

#include <variant>
#include <vector>

namespace LMCAS::detail {

struct RealIsolation {
    ExactRealAlgebraic value;
};

struct ComplexIsolation {
    Rational real_lower;
    Rational real_upper;
    Rational imaginary_lower;
    Rational imaginary_upper;
    ExactRealAlgebraic real_projection;
    ExactRealAlgebraic imaginary_projection;
};

using RootIsolation = std::variant<RealIsolation, ComplexIsolation>;

struct NumericEvaluationOptions {
    double absolute_tolerance = 1e-12;
    double relative_tolerance = 1e-12;
};

Result<ExactRootId> make_exact_root_id(
    Polynomial<Rational> polynomial,
    std::size_t index,
    ComputationContext& context,
    const std::string& operation = "rootof.construct");

Result<std::vector<RootIsolation>> isolate_exact_roots(
    const Polynomial<Rational>& polynomial,
    ComputationContext& context,
    const std::string& operation = "rootof.isolate");

Result<RootIsolation> isolate_exact_root(
    const ExactRootId& root,
    ComputationContext& context,
    const std::string& operation = "rootof.isolate");

Result<void> refine_complex_isolation(
    const Polynomial<Rational>& polynomial,
    ComplexIsolation& isolation,
    const NumericEvaluationOptions& options,
    ComputationContext& context,
    const std::string& operation = "rootof.refine_complex");

Result<ApproxReal> evaluate_root_real(
    const ExactRootId& root,
    const NumericEvaluationOptions& options,
    ComputationContext& context);

Result<ApproxComplex> evaluate_root_complex(
    const ExactRootId& root,
    const NumericEvaluationOptions& options,
    ComputationContext& context);

} // namespace LMCAS::detail
