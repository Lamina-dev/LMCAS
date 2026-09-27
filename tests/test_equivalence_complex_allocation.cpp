#include "test_common.hpp"

#include "expr.hpp"
#include "internal/equivalence_engine.hpp"
#include "internal/equivalence_support.hpp"
#include "internal/rewrite_budget.hpp"
#include "internal/symbolic_ast.hpp"

#include <cstddef>
#include <cstdlib>
#include <limits>
#include <new>
#include <optional>

namespace {

struct AllocationProbeState {
    std::size_t ordinal;
    std::size_t fail_at;
    bool enabled;
};

thread_local AllocationProbeState allocation_probe{};

void* allocate_with_probe(std::size_t size) {
    if (allocation_probe.enabled) {
        ++allocation_probe.ordinal;
        if (allocation_probe.ordinal == allocation_probe.fail_at) {
            allocation_probe.enabled = false;
            throw std::bad_alloc{};
        }
    }
    if (void* memory = std::malloc(size == 0 ? 1 : size)) return memory;
    throw std::bad_alloc{};
}

class AllocationProbeScope final {
public:
    explicit AllocationProbeScope(std::size_t fail_at) noexcept {
        allocation_probe = AllocationProbeState{0, fail_at, true};
    }

    ~AllocationProbeScope() { allocation_probe.enabled = false; }

    std::size_t count() const noexcept { return allocation_probe.ordinal; }

    AllocationProbeScope(const AllocationProbeScope&) = delete;
    AllocationProbeScope& operator=(const AllocationProbeScope&) = delete;
};

void expect_allocation_failures_surface(
    const LMCAS::ExprPtr& lhs, const LMCAS::ExprPtr& rhs,
    const LMCAS::EqvOptions& options = {}) {
    {
        LMCAS::ComputationContext context;
        const auto warmup = LMCAS::equivalent_core(
            *lhs, *rhs, context, options);
        ASSERT_TRUE(warmup);
    }

    std::optional<LMCAS::Result<bool>> counted;
    std::size_t allocation_count = 0;
    {
        LMCAS::ComputationContext context;
        AllocationProbeScope probe(std::numeric_limits<std::size_t>::max());
        counted.emplace(LMCAS::equivalent_core(
            *lhs, *rhs, context, options));
        allocation_count = probe.count();
    }
    ASSERT_TRUE(counted.has_value());
    ASSERT_TRUE(*counted);
    ASSERT_GT(allocation_count, 0U);

    for (std::size_t ordinal = 1; ordinal <= allocation_count; ++ordinal) {
        LMCAS::ComputationContext context;
        std::optional<LMCAS::Result<bool>> outcome;
        {
            AllocationProbeScope probe(ordinal);
            try {
                outcome.emplace(LMCAS::equivalent_core(
                    *lhs, *rhs, context, options));
            } catch (const std::bad_alloc&) {
                FAIL() << "allocation ordinal " << ordinal
                       << " escaped equivalent_core";
            }
        }

        ASSERT_TRUE(outcome.has_value()) << "allocation ordinal " << ordinal;
        ASSERT_FALSE(*outcome)
            << "allocation ordinal " << ordinal << " produced boolean success";
        EXPECT_EQ(outcome->error().code, LMCAS::CasErrc::ResourceLimit)
            << "allocation ordinal " << ordinal;
        EXPECT_EQ(outcome->error().operation, "LMCAS.equivalent_core")
            << "allocation ordinal " << ordinal;
    }
}

} // namespace

void* operator new(std::size_t size) {
    return allocate_with_probe(size);
}

void* operator new[](std::size_t size) {
    return allocate_with_probe(size);
}

void operator delete(void* memory) noexcept {
    std::free(memory);
}

void operator delete[](void* memory) noexcept {
    std::free(memory);
}

void operator delete(void* memory, std::size_t) noexcept {
    std::free(memory);
}

void operator delete[](void* memory, std::size_t) noexcept {
    std::free(memory);
}

TEST(EquivalenceComplexAllocation, ImaginaryUnitFailuresSurface) {
    const auto imaginary = LMCAS::SymbolicExpr::variable("I");
    expect_allocation_failures_surface(imaginary, imaginary);
}

TEST(EquivalenceComplexAllocation, ExplicitComplexFailuresSurface) {
    const auto expression = LMCAS::detail::make_expression_ptr(
        LMCAS::SymbolicFactory::create_complex(
            LMCAS::detail::node(LMCAS::SymbolicExpr::number(2)),
            LMCAS::detail::node(LMCAS::SymbolicExpr::number(3))));
    expect_allocation_failures_surface(expression, expression);
}

TEST(EquivalenceComplexAllocation, InvalidBudgetReturnsResourceLimit) {
    const auto x = LMCAS::SymbolicExpr::variable("x");
    LMCAS::EqvOptions invalid_options;
    invalid_options.budget.max_rewrite_steps = 0;
    LMCAS::ComputationContext context;
    const auto outcome = LMCAS::equivalent_core(
        *x, *x, context, invalid_options);

    ASSERT_FALSE(outcome);
    EXPECT_EQ(outcome.error().code, LMCAS::CasErrc::ResourceLimit);
    EXPECT_EQ(outcome.error().operation, "LMCAS.equivalent_core");
}

TEST(EquivalenceComplexAllocation, ProofProfileFailuresSurface) {
    const auto x = LMCAS::SymbolicExpr::variable("x");
    const auto two = LMCAS::SymbolicExpr::number(2);
    const auto one = LMCAS::SymbolicExpr::number(1);
    const auto trig_identity = LMCAS::SymbolicExpr::add(
        LMCAS::SymbolicExpr::power(LMCAS::SymbolicExpr::sin(x), two),
        LMCAS::SymbolicExpr::power(LMCAS::SymbolicExpr::cos(x), two));
    LMCAS::EqvOptions trig_options;
    trig_options.profile = LMCAS::EqvProfile::TrigBasic;
    expect_allocation_failures_surface(
        trig_identity, one, trig_options);

    const auto exp_identity =
        LMCAS::SymbolicExpr::exp(LMCAS::SymbolicExpr::number(0));
    LMCAS::EqvOptions exp_log_options;
    exp_log_options.profile = LMCAS::EqvProfile::ExpLogBasic;
    expect_allocation_failures_surface(
        exp_identity, one, exp_log_options);
}

TEST(EquivalenceComplexAllocation, InvalidComponentPreservesOriginal) {
    const auto original = LMCAS::SymbolicExpr::variable("x");
    const auto real = LMCAS::SymbolicExpr::number(1);
    LMCAS::ComputationContext context;
    LMCAS::detail::RewriteBudget budget(
        context, 64, 1024, "test.rewritten_complex");

    const auto rewritten = LMCAS::equivalence_detail::rewritten_complex(
        real, {}, LMCAS::detail::node(original), budget);

    ASSERT_TRUE(rewritten);
    EXPECT_TRUE(LMCAS::structurally_equal(*rewritten, *original));
}
