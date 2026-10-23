#include "expr.hpp"
#include "test_common.hpp"
#include "integration.hpp"
#include "matcher.hpp"
#include "assumption_context.hpp"
#include "residual_verification.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>

using namespace LMCAS;

using LMCAS::IntegrationEntry;
using LMCAS::IntegrationTable;
using LMCAS::Matcher;
using LMCAS::MatchMap;

namespace {

constexpr const char *kVarName = "x";
constexpr double kTolerance = 1e-10;

// Concrete values to bind to the well-known wildcard names.
//   _u  -> integration variable x
//   _a  -> 2          (non-zero positive constant, avoids sqrt(0) issues)
//   _n  -> 3          (integer != -1, x^n rule rejects n = -1)
//   any other wildcard -> 1 as a defensive fallback
std::shared_ptr<SymbolicExpr> concrete_for_wildcard(const std::string &wname) {
    if (wname == "_u") {
        return SymbolicExpr::variable(kVarName);
    }
    if (wname == "_a") {
        return SymbolicExpr::number(2);
    }
    if (wname == "_n") {
        return SymbolicExpr::number(3);
    }
    return SymbolicExpr::number(1);
}

MatchMap make_bindings(const IntegrationEntry &entry) {
    MatchMap m;
    for (const auto &w : entry.wildcards) {
        auto val = concrete_for_wildcard(w);
        m.emplace(w, *val);
    }
    return m;
}

const std::vector<double> &sample_points() {
    static const std::vector<double> pts = {0.5, 1.0, 1.5, 2.0, 2.5};
    return pts;
}

struct EntryReport {
    bool checked = false;
    bool failed = false;
    std::string failure_detail;
};

bool evaluate_roundtrip_point(
    const std::shared_ptr<SymbolicExpr> &pattern,
    const std::shared_ptr<SymbolicExpr> &derivative,
    double point, double &pattern_value, double &derivative_value) {
    auto x_value = SymbolicExpr::number(point);
    auto pattern_at = LMCAS::substitute(pattern, kVarName, x_value).value();
    auto derivative_at = LMCAS::substitute(derivative, kVarName, x_value).value();
    if (!pattern_at || !derivative_at) {
        return false;
    }
    pattern_at = LMCAS::simplify(pattern_at).value();
    derivative_at = LMCAS::simplify(derivative_at).value();
    auto pv = test_numeric_eval(pattern_at);
    auto dv = test_numeric_eval(derivative_at);
    if (!pv || !dv || !std::isfinite(*pv) || !std::isfinite(*dv)) {
        return false;
    }
    pattern_value = *pv;
    derivative_value = *dv;
    return true;
}

EntryReport verify_entry(const IntegrationEntry &entry) {
    EntryReport rep;

    MatchMap bindings = make_bindings(entry);

    SymbolicExpr pat_inst = Matcher::replace(entry.pattern, bindings, false);
    auto pattern = LMCAS::detail::make_expression_ptr(pat_inst);
    Integrator integrator;
    ComputationContext context;
    auto assumptions = std::make_shared<AssumptionContext>();
    auto real_variable = assumptions->assume_domain_checked(kVarName, Domain::Real);
    if (!real_variable) {
        rep.failed = true;
        rep.failure_detail = real_variable.error().message;
        return rep;
    }
    auto attached = context.set_assumptions(assumptions);
    if (!attached) {
        rep.failed = true;
        rep.failure_detail = attached.error().message;
        return rep;
    }
    auto entry_result = detail::make_expression_ptr(
        Matcher::replace(entry.result, bindings, false));
    auto entry_derivative = LMCAS::differentiate(entry_result, kVarName, context);
    if (!entry_derivative) {
        rep.failed = true;
        rep.failure_detail = entry_derivative.error().message;
        return rep;
    }
    auto entry_proof = check_equivalent(entry_derivative.value(), pattern, context);
    if (!entry_proof || !std::holds_alternative<ProvedZeroResidual>(entry_proof.value())) {
        rep.failed = true;
        rep.failure_detail = entry_proof ? "table entry residual was not proved zero" : entry_proof.error().message;
        return rep;
    }
    auto integrated = integrator.integrate_checked(*pattern, kVarName, context);
    if (!integrated) {
        rep.failed = true;
        rep.failure_detail = integrated.error().message;
        return rep;
    }
    auto result_simp = detail::make_expression_ptr(integrated.value());
    if (detail::contains_node_type<IntegralNode>(detail::node(result_simp))) {
        rep.failed = true;
        rep.failure_detail = "integrator left an unevaluated integral: " + result_simp->to_string();
        return rep;
    }
    auto differentiated = LMCAS::differentiate(result_simp, kVarName, context);
    if (!differentiated) {
        rep.failed = true;
        rep.failure_detail = differentiated.error().message;
        return rep;
    }
    auto deriv_simp = differentiated.value();
    auto proof = check_equivalent(deriv_simp, pattern, context);
    if (!proof || !std::holds_alternative<ProvedZeroResidual>(proof.value())) {
        rep.failed = true;
        rep.failure_detail = proof ? "derivative residual was not proved zero" : proof.error().message;
        return rep;
    }
    rep.checked = true;
    auto pattern_simp = pattern;

    for (double xv : sample_points()) {
        double pattern_value = 0.0;
        double derivative_value = 0.0;
        if (!evaluate_roundtrip_point(pattern_simp, deriv_simp, xv,
                                      pattern_value, derivative_value)) {
            continue;
        }

        double delta = std::abs(pattern_value - derivative_value);
        if (delta > kTolerance) {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv
                << ": d/dx(result)=" << derivative_value
                << " vs pattern=" << pattern_value
                << " |delta|=" << delta
                << " | pattern_expr=" << pattern_simp->to_string()
                << " | result_expr=" << result_simp->to_string()
                << " | derivative_expr=" << deriv_simp->to_string();
            rep.failure_detail = oss.str();
            break;
        }
    }
    return rep;
}

struct CategorySpec {
    IntegrationTable::Category cat;
    const char *name;
};

const std::vector<CategorySpec> &all_categories() {
    static const std::vector<CategorySpec> cats = {
        {IntegrationTable::Category::Polynomial, "Polynomial"},
        {IntegrationTable::Category::Exponential, "Exponential"},
        {IntegrationTable::Category::Logarithmic, "Logarithmic"},
        {IntegrationTable::Category::Trigonometric, "Trigonometric"},
        {IntegrationTable::Category::InverseTrig, "InverseTrig"},
        {IntegrationTable::Category::Hyperbolic, "Hyperbolic"},
        {IntegrationTable::Category::Algebraic, "Algebraic"},
        {IntegrationTable::Category::Special, "Special"},
        {IntegrationTable::Category::UserDefined, "UserDefined"},
    };
    return cats;
}

} // anonymous namespace

TEST(IntegrationTableRoundtrip, EntriesDifferentiateToTheirPatterns) {
    IntegrationTable table;
    int verified_entries = 0;
    for (const auto &spec : all_categories()) {
        for (const auto &entry : table.get_entries(spec.cat)) {
            SCOPED_TRACE(std::string(spec.name) + "/" + entry.name);
            const EntryReport rep = verify_entry(entry);
            EXPECT_FALSE(rep.failed) << rep.failure_detail;
            if (!rep.failed && rep.checked) {
                ++verified_entries;
            }
        }
    }
    EXPECT_GT(verified_entries, 0)
        << "at least one entry verified by exact derivative residual";
}
