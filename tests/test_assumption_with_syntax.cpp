
#include "test_common.hpp"
#include "assumption_context.hpp"
#include <stdexcept>

using namespace LMCAS;

static void declare_parent_fact(AssumptionContext &ctx) {
    ASSERT_TRUE(ctx.assume_sign("parent", Sign::Negative));
}

static void expect_parent_restored(const AssumptionContext &ctx) {
    EXPECT_EQ(ctx.depth(), 1);
    EXPECT_TRUE(ctx.has_sign("parent", Sign::Negative));
    EXPECT_FALSE(ctx.has_sign("temporary", Sign::Positive));
    EXPECT_FALSE(ctx.has_sign("nested", Sign::NonZero));
}

static void expect_payload(const CasError &error) {
    EXPECT_EQ(error.code, CasErrc::DomainError);
    EXPECT_EQ(error.message, "callable domain payload");
    EXPECT_EQ(error.operation, "test.callable");
}

TEST(LmcasAssumptionWithSyntax, CasErrorPayloadAndScopeCleanup) {
    AssumptionContext value_context;
    declare_parent_fact(value_context);
    auto value_result = with_assumptions(
        value_context, {AssumptionDecl::make_sign("temporary", Sign::Positive)},
        [&]() -> int {
            EXPECT_EQ(value_context.depth(), 2);
            EXPECT_TRUE(value_context.has_sign("temporary", Sign::Positive));
            throw CasError{CasErrc::DomainError, "callable domain payload",
                           "test.callable"};
        });
    ASSERT_FALSE(value_result);
    expect_payload(value_result.error());
    expect_parent_restored(value_context);

    AssumptionContext void_context;
    declare_parent_fact(void_context);
    auto void_result = with_assumptions(
        void_context, {AssumptionDecl::make_sign("temporary", Sign::Positive)},
        [&]() {
            EXPECT_EQ(void_context.depth(), 2);
            EXPECT_TRUE(void_context.has_sign("temporary", Sign::Positive));
            throw CasError{CasErrc::DomainError, "callable domain payload",
                           "test.callable"};
        });
    ASSERT_FALSE(void_result);
    expect_payload(void_result.error());
    expect_parent_restored(void_context);
}

TEST(LmcasAssumptionWithSyntax, StandardExceptionPayloadAndScopeCleanup) {
    AssumptionContext value_context;
    declare_parent_fact(value_context);
    auto value_result = with_assumptions(
        value_context, {AssumptionDecl::make_sign("temporary", Sign::Positive)},
        []() -> int { throw std::runtime_error("value callable failure"); });
    ASSERT_FALSE(value_result);
    EXPECT_EQ(value_result.error().code, CasErrc::InternalInvariant);
    EXPECT_EQ(value_result.error().message, "value callable failure");
    EXPECT_EQ(value_result.error().operation, "with_assumptions");
    expect_parent_restored(value_context);

    AssumptionContext void_context;
    declare_parent_fact(void_context);
    auto void_result = with_assumptions(
        void_context, {AssumptionDecl::make_sign("temporary", Sign::Positive)},
        []() { throw std::logic_error("void callable failure"); });
    ASSERT_FALSE(void_result);
    EXPECT_EQ(void_result.error().code, CasErrc::InternalInvariant);
    EXPECT_EQ(void_result.error().message, "void callable failure");
    EXPECT_EQ(void_result.error().operation, "with_assumptions");
    expect_parent_restored(void_context);
}

TEST(LmcasAssumptionWithSyntax, NonstandardExceptionRethrowsAfterScopeCleanup) {
    AssumptionContext value_context;
    declare_parent_fact(value_context);
    try {
        (void)with_assumptions(
            value_context,
            {AssumptionDecl::make_sign("temporary", Sign::Positive)},
            []() -> int { throw 17; });
        FAIL() << "value callable must rethrow its nonstandard exception";
    } catch (int payload) {
        EXPECT_EQ(payload, 17);
    } catch (...) {
        FAIL() << "value callable changed the nonstandard exception type";
    }
    expect_parent_restored(value_context);

    AssumptionContext void_context;
    declare_parent_fact(void_context);
    try {
        (void)with_assumptions(
            void_context,
            {AssumptionDecl::make_sign("temporary", Sign::Positive)},
            []() { throw 23; });
        FAIL() << "void callable must rethrow its nonstandard exception";
    } catch (int payload) {
        EXPECT_EQ(payload, 23);
    } catch (...) {
        FAIL() << "void callable changed the nonstandard exception type";
    }
    expect_parent_restored(void_context);
}

TEST(LmcasAssumptionWithSyntax, NormalReturnRestoresParentScope) {
    AssumptionContext value_context;
    declare_parent_fact(value_context);
    auto value_result = with_assumptions(
        value_context, {AssumptionDecl::make_sign("temporary", Sign::Positive)},
        [&]() -> int {
            EXPECT_TRUE(value_context.has_sign("parent", Sign::Negative));
            EXPECT_TRUE(value_context.has_sign("temporary", Sign::Positive));
            return 42;
        });
    ASSERT_TRUE(value_result);
    EXPECT_EQ(value_result.value(), 42);
    expect_parent_restored(value_context);

    AssumptionContext void_context;
    declare_parent_fact(void_context);
    bool called = false;
    auto void_result = with_assumptions(
        void_context, {AssumptionDecl::make_sign("temporary", Sign::Positive)},
        [&]() {
            called = true;
            EXPECT_TRUE(void_context.has_sign("parent", Sign::Negative));
            EXPECT_TRUE(void_context.has_sign("temporary", Sign::Positive));
        });
    ASSERT_TRUE(void_result);
    EXPECT_TRUE(called);
    expect_parent_restored(void_context);
}

TEST(LmcasAssumptionWithSyntax, NestedScopesRestoreEachParent) {
    AssumptionContext value_context;
    declare_parent_fact(value_context);
    auto outer_value = with_assumptions(
        value_context, {AssumptionDecl::make_sign("temporary", Sign::Positive)},
        [&]() -> int {
            auto inner = with_assumptions(
                value_context,
                {AssumptionDecl::make_sign("nested", Sign::NonZero)},
                [&]() -> int {
                    EXPECT_EQ(value_context.depth(), 3);
                    EXPECT_TRUE(value_context.has_sign("temporary", Sign::Positive));
                    throw CasError{CasErrc::DomainError,
                                   "callable domain payload", "test.callable"};
                });
            if (inner) {
                ADD_FAILURE() << "nested value callable must preserve its CasError";
                return -1;
            }
            expect_payload(inner.error());
            EXPECT_EQ(value_context.depth(), 2);
            EXPECT_TRUE(value_context.has_sign("temporary", Sign::Positive));
            EXPECT_FALSE(value_context.has_sign("nested", Sign::NonZero));
            return 7;
        });
    ASSERT_TRUE(outer_value);
    EXPECT_EQ(outer_value.value(), 7);
    expect_parent_restored(value_context);

    AssumptionContext void_context;
    declare_parent_fact(void_context);
    auto outer_void = with_assumptions(
        void_context, {AssumptionDecl::make_sign("temporary", Sign::Positive)},
        [&]() {
            auto inner = with_assumptions(
                void_context,
                {AssumptionDecl::make_sign("nested", Sign::NonZero)},
                [&]() {
                    EXPECT_EQ(void_context.depth(), 3);
                    throw CasError{CasErrc::DomainError,
                                   "callable domain payload", "test.callable"};
                });
            ASSERT_FALSE(inner);
            expect_payload(inner.error());
            EXPECT_EQ(void_context.depth(), 2);
            EXPECT_TRUE(void_context.has_sign("temporary", Sign::Positive));
            EXPECT_FALSE(void_context.has_sign("nested", Sign::NonZero));
        });
    ASSERT_TRUE(outer_void);
    expect_parent_restored(void_context);
}
