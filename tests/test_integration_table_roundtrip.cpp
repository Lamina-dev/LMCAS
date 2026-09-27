
#include "test_common.hpp"
#include "integration.hpp"
#include "matcher.hpp"

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
    bool checked = false; // at least one (pattern, derivative) pair was evaluated
    bool failed = false;  // a numeric mismatch was observed
    int matches = 0;
    std::string failure_detail;
};

bool evaluate_roundtrip_point(
    const std::shared_ptr<SymbolicExpr> &pattern,
    const std::shared_ptr<SymbolicExpr> &derivative,
    double point, double &pattern_value, double &derivative_value) {
    auto x_value = SymbolicExpr::number(point);
    auto pattern_at = pattern->substitute(kVarName, x_value);
    auto derivative_at = derivative->substitute(kVarName, x_value);
    if (!pattern_at || !derivative_at) {
        return false;
    }
    pattern_at = pattern_at->simplify();
    derivative_at = derivative_at->simplify();
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
    SymbolicExpr res_inst = Matcher::replace(entry.result, bindings, false);

    auto pattern = LMCAS::detail::make_expression_ptr(pat_inst);
    auto result = LMCAS::detail::make_expression_ptr(res_inst);

    auto pattern_simp = pattern->simplify();
    auto result_simp = result->simplify();
    if (!pattern_simp || !result_simp) {
        rep.failure_detail = "instantiation/simplification returned null";
        return rep;
    }

    auto deriv = result_simp->differentiate(kVarName);
    if (!deriv) {
        rep.failure_detail = "differentiation returned null";
        return rep;
    }
    auto deriv_simp = deriv->simplify();
    if (!deriv_simp) {
        deriv_simp = deriv;
    }

    for (double xv : sample_points()) {
        double pattern_value = 0.0;
        double derivative_value = 0.0;
        if (!evaluate_roundtrip_point(pattern_simp, deriv_simp, xv,
                                      pattern_value, derivative_value)) {
            continue;
        }

        rep.checked = true;
        double delta = std::abs(pattern_value - derivative_value);
        if (delta <= kTolerance) {
            ++rep.matches;
        } else {
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
                EXPECT_GT(rep.matches, 0);
                ++verified_entries;
            }
        }
    }
    EXPECT_GT(verified_entries, 0)
        << "at least one entry verified by numeric round-trip";
}
