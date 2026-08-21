
#include "test_common.hpp"
#include "integration.hpp"

#include <string>
#include <vector>

using lamina::IntegrationTable;
using lamina::IntegrationEntry;

namespace {

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

// Returns true iff the entries are in non-decreasing order of priority.
bool is_non_decreasing_priority(const std::vector<IntegrationEntry>& entries) {
    for (size_t i = 1; i < entries.size(); ++i) {
        if (entries[i].priority < entries[i - 1].priority) {
            return false;
        }
    }
    return true;
}

void check_priority_order(const IntegrationTable& table,
                          const CategorySpec& spec,
                          const std::string& phase) {
    const auto& entries = table.get_entries(spec.cat);
    bool ok = is_non_decreasing_priority(entries);
    std::string msg = phase + ": category " + spec.name
                    + " entries (n=" + std::to_string(entries.size())
                    + ") sorted by priority";
    if (!ok) {
        std::cerr << "Priority sequence for " << spec.name << ":";
        for (const auto& e : entries) {
            std::cerr << " [" << e.name << ":" << e.priority << "]";
        }
        std::cerr << std::endl;
    }
    EXPECT_TRUE(ok, msg);
}

}// anonymous namespace

int main() {
    TEST_CASE("Integration table priority ordering invariant");

    // Phase 1: defaults loaded by the IntegrationTable ctor.
    {
        IntegrationTable table;
        for (const auto& spec : all_categories()) {
            check_priority_order(table, spec, "after load_defaults");
        }
    }

    // Phase 2: stress the invariant by inserting entries with arbitrary
    // priorities (interleaved high/low values) into every category and
    // re-checking. This validates that add_entry preserves the ordering even
    // when called many times with non-monotonic priorities, which is the
    {
        IntegrationTable table;
        const std::vector<int> priorities = {
            100, 5, 250, 1, 75, 30, 999, 12, 60, 7, 200, 42, 8, 33, 150, 0
        };

        auto x_var = SymbolicExpr::variable("x");
        auto pat = *SymbolicExpr::power(x_var, SymbolicExpr::number(2));
        auto res = *SymbolicExpr::power(x_var, SymbolicExpr::number(3));

        for (const auto& spec : all_categories()) {
            for (size_t i = 0; i < priorities.size(); ++i) {
                IntegrationEntry e(
                    "synthetic_" + std::string(spec.name) + "_" + std::to_string(i),
                    pat, res, {"_u"}, nullptr, priorities[i]);
                table.add_entry(spec.cat, e);
                check_priority_order(table, spec, "after add_entry #" + std::to_string(i));
            }
        }
    }

    return TEST_REPORT();
}
