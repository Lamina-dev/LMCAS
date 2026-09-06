
#include "test_common.hpp"
#include "integration.hpp"
#include "matcher.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>

using namespace LMCAS;

using LMCAS::IntegrationTable;
using LMCAS::IntegrationEntry;
using LMCAS::Matcher;
using LMCAS::MatchMap;

namespace {

constexpr const char* kVarName = "x";
constexpr double kTolerance = 1e-10;

// Concrete values to bind to the well-known wildcard names.
//   _u  -> integration variable x
//   _a  -> 2          (non-zero positive constant, avoids sqrt(0) issues)
//   _n  -> 3          (integer != -1, x^n rule rejects n = -1)
//   any other wildcard -> 1 as a defensive fallback
std::shared_ptr<SymbolicExpr> concrete_for_wildcard(const std::string& wname) {
    if (wname == "_u") return SymbolicExpr::variable(kVarName);
    if (wname == "_a") return SymbolicExpr::number(2);
    if (wname == "_n") return SymbolicExpr::number(3);
    return SymbolicExpr::number(1);
}

MatchMap make_bindings(const IntegrationEntry& entry) {
    MatchMap m;
    for (const auto& w : entry.wildcards) {
        auto val = concrete_for_wildcard(w);
        m.emplace(w, *val);
    }
    return m;
}

const std::vector<double>& sample_points() {
    static const std::vector<double> pts = {0.5, 1.0, 1.5, 2.0, 2.5};
    return pts;
}

struct EntryReport {
    bool checked = false;       // at least one (pattern, derivative) pair was evaluated
    bool failed  = false;       // a numeric mismatch was observed
    int  matches = 0;
    int  skipped_points = 0;
    std::string failure_detail;
};

EntryReport verify_entry(const IntegrationEntry& entry) {
    EntryReport rep;

    MatchMap bindings = make_bindings(entry);

    SymbolicExpr pat_inst = Matcher::replace(entry.pattern, bindings, false);
    SymbolicExpr res_inst = Matcher::replace(entry.result,  bindings, false);

    auto pattern = LMCAS::detail::make_expression_ptr(pat_inst);
    auto result  = LMCAS::detail::make_expression_ptr(res_inst);

    auto pattern_simp = pattern->simplify();
    auto result_simp  = result->simplify();
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
    if (!deriv_simp) deriv_simp = deriv;

    for (double xv : sample_points()) {
        auto x_val = SymbolicExpr::number(xv);

        auto pat_at = pattern_simp->substitute(kVarName, x_val);
        auto der_at = deriv_simp->substitute(kVarName, x_val);
        if (!pat_at || !der_at) {
            ++rep.skipped_points;
            continue;
        }
        pat_at = pat_at->simplify();
        der_at = der_at->simplify();

        auto pv = test_numeric_eval(pat_at);
        auto dv = test_numeric_eval(der_at);

        if (!pv || !dv || !std::isfinite(*pv) || !std::isfinite(*dv)) {
            ++rep.skipped_points;
            continue;
        }

        rep.checked = true;
        double delta = std::abs(*pv - *dv);
        if (delta <= kTolerance) {
            ++rep.matches;
        } else {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv
                << ": d/dx(result)=" << *dv
                << " vs pattern=" << *pv
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
    const char* name;
};

const std::vector<CategorySpec>& all_categories() {
    static const std::vector<CategorySpec> cats = {
        {IntegrationTable::Category::Polynomial,    "Polynomial"},
        {IntegrationTable::Category::Exponential,   "Exponential"},
        {IntegrationTable::Category::Logarithmic,   "Logarithmic"},
        {IntegrationTable::Category::Trigonometric, "Trigonometric"},
        {IntegrationTable::Category::InverseTrig,   "InverseTrig"},
        {IntegrationTable::Category::Hyperbolic,    "Hyperbolic"},
        {IntegrationTable::Category::Algebraic,     "Algebraic"},
        {IntegrationTable::Category::Special,       "Special"},
        {IntegrationTable::Category::UserDefined,   "UserDefined"},
    };
    return cats;
}

}// anonymous namespace

int main() {
    TEST_CASE("Integration table entry round-trip");

    IntegrationTable table;

    int total_entries = 0;
    int verified_entries = 0;
    int skipped_entries = 0;
    int failed_entries = 0;

    for (const auto& spec : all_categories()) {
        const auto& entries = table.get_entries(spec.cat);
        for (const auto& entry : entries) {
            ++total_entries;
            EntryReport rep = verify_entry(entry);

            std::string prefix = std::string(spec.name) + "/" + entry.name;

            if (rep.failed) {
                ++failed_entries;
                std::cerr << "[FAIL] " << prefix << " : " << rep.failure_detail << std::endl;
                EXPECT_TRUE(false, prefix + ": numeric round-trip mismatch");
            } else if (!rep.checked) {
                ++skipped_entries;
                std::cout << "[SKIP] " << prefix
                          << " (no closed-form numeric eval at any sample point)"
                          << std::endl;
            } else {
                ++verified_entries;
                std::ostringstream oss;
                oss << prefix << ": " << rep.matches << " match(es), "
                    << rep.skipped_points << " skipped point(s)";
                EXPECT_TRUE(rep.checked && !rep.failed && rep.matches > 0, oss.str());
            }
        }
    }

    std::cout << "\nSummary: total=" << total_entries
              << ", verified=" << verified_entries
              << ", skipped=" << skipped_entries
              << ", failed=" << failed_entries << std::endl;

    EXPECT_TRUE(verified_entries > 0,
                "at least one entry verified by numeric round-trip");

    return TEST_REPORT();
}
