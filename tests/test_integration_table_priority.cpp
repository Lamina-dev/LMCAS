
#include "test_common.hpp"
#include "integration.hpp"
#include "expr.hpp"
#include "residual_verification.hpp"

#include <string>
#include <vector>

using namespace LMCAS;

using LMCAS::IntegrationEntry;
using LMCAS::IntegrationTable;

namespace {

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

// Returns true iff the entries are in non-decreasing order of priority.
bool is_non_decreasing_priority(const std::vector<IntegrationEntry> &entries) {
    for (size_t i = 1; i < entries.size(); ++i) {
        if (entries[i].priority < entries[i - 1].priority) {
            return false;
        }
    }
    return true;
}

void check_priority_order(const IntegrationTable &table,
                          const CategorySpec &spec,
                          const std::string &phase) {
    const auto &entries = table.get_entries(spec.cat);
    bool ok = is_non_decreasing_priority(entries);
    std::string msg = phase + ": category " + spec.name + " entries (n=" + std::to_string(entries.size()) + ") sorted by priority";
    if (!ok) {
        std::cerr << "Priority sequence for " << spec.name << ":";
        for (const auto &e : entries) {
            std::cerr << " [" << e.name << ":" << e.priority << "]";
        }
        std::cerr << std::endl;
    }
    EXPECT_TRUE((ok)) << msg;
}

} // anonymous namespace

TEST(LmcasIntegrationTablePriority, DefaultPriorityOrdering) {
    IntegrationTable table;
    for (const auto &spec : all_categories()) {
        check_priority_order(table, spec, "after load_defaults");
    }
}

TEST(LmcasIntegrationTablePriority, InsertedPriorityOrdering) {
    IntegrationTable table;
    const std::vector<int> priorities = {
        100, 5, 250, 1, 75, 30, 999, 12, 60, 7, 200, 42, 8, 33, 150, 0};

    auto x_var = SymbolicExpr::variable("x");
    auto pat = *SymbolicExpr::power(x_var, SymbolicExpr::number(2));
    auto res = *SymbolicExpr::power(x_var, SymbolicExpr::number(3));

    for (const auto &spec : all_categories()) {
        for (size_t i = 0; i < priorities.size(); ++i) {
            IntegrationEntry e(
                "synthetic_" + std::string(spec.name) + "_" + std::to_string(i),
                pat, res, {"_u"}, nullptr, priorities[i]);
            table.add_entry(spec.cat, e);
            check_priority_order(table, spec, "after add_entry #" + std::to_string(i));
        }
    }
}

TEST(LmcasIntegrationTablePriority, IncorrectSamePatternFallsBack) {
    for (bool affine : {false, true}) {
        SCOPED_TRACE(affine ? "affine substitution" : "direct lookup");
        Integrator integrator;
        for (const auto &spec : all_categories()) {
            integrator.table().clear_category(spec.cat);
        }
        auto u = SymbolicExpr::variable("_u");
        auto pattern = SymbolicExpr::sin(u);
        auto correct = SymbolicExpr::multiply(
            SymbolicExpr::number(-1), SymbolicExpr::cos(u));
        integrator.table().add_entry(IntegrationTable::Category::UserDefined,
            IntegrationEntry("wrong sine", *pattern, *u, {"_u"}, nullptr, -100));
        integrator.table().add_entry(IntegrationTable::Category::UserDefined,
            IntegrationEntry("correct sine", *pattern, *correct, {"_u"}, nullptr, -99));
        auto x = SymbolicExpr::variable("x");
        auto argument = affine
            ? SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x),
                                SymbolicExpr::number(1))
            : x;
        auto integrand = SymbolicExpr::sin(argument);
        ComputationContext context;
        TableLookupStrategy direct;
        LinearSubstitutionStrategy substitution;
        auto attempted = affine
            ? substitution.try_integrate(*integrand, "x", integrator, context)
            : direct.try_integrate(*integrand, "x", integrator, context);
        ASSERT_TRUE(attempted);
        auto *candidate = std::get_if<IntegrationCandidate>(&attempted.value());
        ASSERT_NE(candidate, nullptr);
        auto candidate_derivative = LMCAS::differentiate(candidate->expression, "x", context);
        ASSERT_TRUE(candidate_derivative);
        auto candidate_proof = check_equivalent(candidate_derivative.value(), integrand, context);
        ASSERT_TRUE(candidate_proof);
        EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(candidate_proof.value()));
        auto result = integrator.integrate_checked(*integrand, "x", context);
        ASSERT_TRUE(result) << result.error().message;
        ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
        auto derivative = LMCAS::differentiate(
            detail::make_expression_ptr(result.value()), "x", context);
        ASSERT_TRUE(derivative);
        auto proof = check_equivalent(derivative.value(), integrand, context);
        ASSERT_TRUE(proof);
        EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()));
    }
}
