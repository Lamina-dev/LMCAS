#include "internal/integration_definite_support.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {
namespace {
using namespace detail::definite;

struct DefiniteInterval {
    SymbolicExpr lo;
    SymbolicExpr hi;
    Result<Boundary> a;
    Result<Boundary> b;
    bool reversed = false;

    DefiniteInterval(const SymbolicExpr& lower, const SymbolicExpr& upper,
        ComputationContext& context)
        : lo(lower), hi(upper), a(boundary_checked(lo, context)), b(boundary_checked(hi, context)) {}
};

Result<SymbolicExpr> unresolved_integral(const SymbolicExpr& expr,
    const std::string& var_name, const SymbolicExpr& lower, const SymbolicExpr& upper) {
    return detail::expression_from_node(detail::make_node<IntegralNode>(detail::node(expr),
        var_name, detail::node(lower), detail::node(upper)));
}

Result<SymbolicExpr> equal_endpoint_integral(const SymbolicExpr& expr,
    const std::string& var_name, const SymbolicExpr& lower, const SymbolicExpr& upper,
    const DefiniteInterval& interval, ComputationContext& context) {
    if (interval.a.value().infinity) { return unresolved_integral(expr, var_name, lower, upper); }
    auto assumptions = real_assumptions(var_name, context);
    if (!assumptions) { return Result<SymbolicExpr>::failure(assumptions.error()); }
    detail::AssumptionFacts facts(assumptions.value());
    auto value = expr.substitute(var_name, detail::make_expression_ptr(interval.lo));
    auto defined = detail::query_definedness(detail::node(value), facts, Domain::Real, context);
    if (!defined) { return Result<SymbolicExpr>::failure(defined.error()); }
    if (defined.value() == Tribool::False) { return divergent<SymbolicExpr>(); }
    return defined.value() == Tribool::True ? Result<SymbolicExpr>(*SymbolicExpr::number(0)) : unresolved_integral(expr, var_name, lower, upper);
}

Result<Tribool> certify_partition(const SymbolicExpr& expr, const std::string& var_name,
    const std::vector<SymbolicExpr>& partition, ComputationContext& context) {
    for (std::size_t i = 1; i < partition.size(); ++i) {
        auto continuous = certify_continuity_checked(expr, var_name,
            partition[i - 1], partition[i], context);
        if (!continuous || continuous.value() != Tribool::True) { return continuous; }
    }
    return Tribool::True;
}

struct RationalIntegrand {
    SymbolicExpr expression;
    int divergence = 0;
};

Result<RationalIntegrand> reduced_integrand_checked(const SymbolicExpr& expr,
    const std::string& var_name, const DefiniteInterval& interval, ComputationContext& context) {
    auto rational = rational_checked(detail::node(expr), var_name, context);
    if (!rational) { return Result<RationalIntegrand>::failure(rational.error()); }
    RationalIntegrand result{expr};
    if (rational.value() && interval.a && interval.b) {
        auto divergence = rational_divergence_checked(
            *rational.value(), interval.a.value(), interval.b.value(), context);
        if (!divergence) { return Result<RationalIntegrand>::failure(divergence.error()); }
        result.divergence = divergence.value();
        if (result.divergence) { return result; }
        result.expression = *SymbolicExpr::divide(poly_to_symbolic(rational.value()->numerator),
            poly_to_symbolic(rational.value()->denominator));
    }
    return result;
}

Result<SymbolicExpr> integrate_segment_primitive(Integrator& integrator,
    const SymbolicExpr& value_expression, const std::string& var_name,
    const SymbolicExpr& lower, const SymbolicExpr& upper, ComputationContext& context) {
    auto selected = segment_expression_checked(value_expression, var_name, lower, upper, context);
    if (!selected) { return selected; }
    auto integrated = integrator.integrate_checked(selected.value(), var_name, context);
    if (!integrated) { return integrated; }
    if (contains_unevaluated_integral(detail::node(integrated.value()))) { return undecided<SymbolicExpr>(); }
    return integrated;
}

Result<SymbolicExpr> accumulate_segments(Integrator& integrator, const SymbolicExpr& expr,
    const SymbolicExpr& value_expression, const std::string& var_name,
    const std::vector<SymbolicExpr>& partition, const std::optional<SymbolicExpr>& primitive,
    bool reversed, ComputationContext& context) {
    SymbolicExpr sum = *SymbolicExpr::number(0);
    int infinity = 0;
    for (std::size_t i = 1; i < partition.size(); ++i) {
        std::optional<SymbolicExpr> local_primitive;
        if (!primitive) {
            auto integrated = integrate_segment_primitive(integrator, value_expression, var_name,
                partition[i - 1], partition[i], context);
            if (!integrated) { return integrated; }
            local_primitive = std::move(integrated.value());
        }
        auto segment = definite_segment_checked(expr, primitive ? *primitive : *local_primitive, var_name,
            partition[i - 1], partition[i], context);
        if (!segment) { return segment; }
        int contribution = infinity_sign(detail::node(segment.value()));
        if (contribution) {
            if (infinity && infinity != contribution) { return divergent<SymbolicExpr>(); }
            infinity = contribution;
        } else { sum = *SymbolicExpr::add(detail::make_expression_ptr(sum), detail::make_expression_ptr(segment.value())); }
    }
    if (infinity) { return *SymbolicExpr::infinity(infinity * (reversed ? -1 : 1)); }
    if (reversed) { sum = *SymbolicExpr::multiply(SymbolicExpr::number(-1), detail::make_expression_ptr(sum)); }
    return *sum.simplify();
}

Result<SymbolicExpr> integrate_ordered_interval(Integrator& integrator, const SymbolicExpr& expr,
    const std::string& var_name, const SymbolicExpr& lower, const SymbolicExpr& upper,
    const DefiniteInterval& interval, ComputationContext& context) {
    const auto unresolved = [&]() { return unresolved_integral(expr, var_name, lower, upper); };
    auto partition = definite_partition_checked(expr, var_name, interval.lo, interval.hi, context);
    if (!partition) { return is_undecided(partition.error()) ? unresolved() : Result<SymbolicExpr>::failure(partition.error()); }
    auto continuous = certify_partition(expr, var_name, partition.value(), context);
    if (!continuous) { return Result<SymbolicExpr>::failure(continuous.error()); }
    if (continuous.value() == Tribool::False) { return divergent<SymbolicExpr>(); }
    if (continuous.value() != Tribool::True) { return unresolved(); }
    auto reduced = reduced_integrand_checked(expr, var_name, interval, context);
    if (!reduced) { return Result<SymbolicExpr>::failure(reduced.error()); }
    if (reduced.value().divergence) {
        return *SymbolicExpr::infinity(reduced.value().divergence * (interval.reversed ? -1 : 1));
    }
    const auto& value_expression = reduced.value().expression;
    std::optional<SymbolicExpr> primitive;
    if (!needs_segment_expression(detail::node(value_expression))) {
        auto integrated = integrator.integrate_checked(value_expression, var_name, context);
        if (!integrated) { return is_undecided(integrated.error()) ? unresolved() : integrated; }
        if (contains_unevaluated_integral(detail::node(integrated.value()))) { return unresolved(); }
        primitive = std::move(integrated.value());
    }
    auto sum = accumulate_segments(integrator, expr, value_expression, var_name,
        partition.value(), primitive, interval.reversed, context);
    if (!sum && is_undecided(sum.error())) { return unresolved(); }
    return sum;
}

Result<SymbolicExpr> integrate_interval_checked(Integrator& integrator, const SymbolicExpr& expr,
    const std::string& var_name, const SymbolicExpr& lower, const SymbolicExpr& upper,
    ComputationContext& context) {
    DefiniteInterval interval(lower, upper, context);
    if (!interval.a && !is_undecided(interval.a.error())) { return Result<SymbolicExpr>::failure(interval.a.error()); }
    if (!interval.b && !is_undecided(interval.b.error())) { return Result<SymbolicExpr>::failure(interval.b.error()); }
    if (interval.a && interval.b) {
        auto order = compare_boundaries(interval.a.value(), interval.b.value(), context);
        if (!order) { return Result<SymbolicExpr>::failure(order.error()); }
        if (!order.value()) {
            return equal_endpoint_integral(expr, var_name, lower, upper, interval, context);
        }
        if (order.value() > 0) {
            std::swap(interval.lo, interval.hi);
            std::swap(interval.a.value(), interval.b.value());
            interval.reversed = true;
        }
    }
    return integrate_ordered_interval(integrator, expr, var_name, lower, upper, interval, context);
}
}

Result<SymbolicExpr> Integrator::integrate_def_checked(const SymbolicExpr& expr,
    const std::string& var_name, const SymbolicExpr& lower, const SymbolicExpr& upper,
    ComputationContext& context) {
    auto access = context.consume_steps(1, operation);
    if (!access) { return Result<SymbolicExpr>::failure(access.error()); }
    if (var_name.empty() || depends_on_integration_variable(lower, var_name) ||
        depends_on_integration_variable(upper, var_name))
        { return Result<SymbolicExpr>::failure(CasErrc::InvalidArgument,
            "integration endpoints must be independent of the integration variable", operation); }
    try {
        return integrate_interval_checked(*this, expr, var_name, lower, upper, context);
    } catch (const std::bad_alloc&) {
        return Result<SymbolicExpr>::failure(CasErrc::ResourceLimit,
            "definite integration allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<SymbolicExpr>::failure(CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<SymbolicExpr> Integrator::integrate_def(const SymbolicExpr& expr,
    const std::string& var_name, const SymbolicExpr& lower, const SymbolicExpr& upper) {
    ComputationContext context;
    return integrate_def_checked(expr, var_name, lower, upper, context);
}
}
