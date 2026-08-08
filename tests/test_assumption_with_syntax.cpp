
#include "test_common.hpp"
#include "rapidcheck/rapidcheck.h"
#include "assumption_context.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <stdexcept>
#include <vector>
#include <string>

using namespace lamina;


/// Generate a random Domain from the valid set.
static Domain random_domain() {
    std::vector<Domain> domains = {
        Domain::Real, Domain::Integer, Domain::Rational,
        Domain::Natural, Domain::PositiveInt
    };
    return rc::gen::elementOf(domains);
}

/// Generate a random Sign from the valid set.
static Sign random_sign() {
    std::vector<Sign> signs = {
        Sign::Positive, Sign::Negative, Sign::NonNegative,
        Sign::NonPositive, Sign::NonZero
    };
    return rc::gen::elementOf(signs);
}

/// Generate a random variable name.
static std::string random_var_name() {
    return "var_" + std::to_string(rc::gen::inRange(0, 999));
}

/// Generate a random set of AssumptionDecl (domain and sign only, to avoid
/// complex expression construction).
static std::vector<AssumptionDecl> random_decls() {
    int count = rc::gen::inRange(0, 4);
    std::vector<AssumptionDecl> decls;
    for (int i = 0; i < count; ++i) {
        std::string var = "d_" + std::to_string(i) + "_" + std::to_string(rc::gen::inRange(0, 99));
        if (rc::gen::boolean()) {
            decls.push_back(AssumptionDecl::make_domain(var, random_domain()));
        } else {
            decls.push_back(AssumptionDecl::make_sign(var, random_sign()));
        }
    }
    return decls;
}


static void test_with_assumptions_preserves_depth_normal() {
    TEST_CASE("with_assumptions preserves depth on normal return");

    rc::check("with_assumptions preserves context depth on normal callable return", []() {
        AssumptionContext ctx;

        // Randomly push some scopes to start at a non-trivial depth
        int initial_pushes = rc::gen::inRange(0, 5);
        for (int i = 0; i < initial_pushes; ++i) {
            ctx.push();
        }
        int depth_before = ctx.depth();

        // Generate random declarations
        auto decls = random_decls();

        // Call with_assumptions with a normal callable
        int result = with_assumptions(ctx, decls, [&]() -> int {
            // Inside the scope, depth should be one more
            RC_ASSERT(ctx.depth() == depth_before + 1);
            return 42;
        });

        // After with_assumptions, depth must be restored
        RC_ASSERT(ctx.depth() == depth_before);
        RC_ASSERT(result == 42);
    });
}

static void test_with_assumptions_preserves_depth_exception() {
    TEST_CASE("with_assumptions preserves depth on exception");

    rc::check("with_assumptions preserves context depth even when callable throws", []() {
        AssumptionContext ctx;

        // Randomly push some scopes
        int initial_pushes = rc::gen::inRange(0, 5);
        for (int i = 0; i < initial_pushes; ++i) {
            ctx.push();
        }
        int depth_before = ctx.depth();

        // Generate random declarations
        auto decls = random_decls();

        // Call with_assumptions with a throwing callable
        bool caught = false;
        try {
            with_assumptions(ctx, decls, [&]() -> int {
                RC_ASSERT(ctx.depth() == depth_before + 1);
                throw std::runtime_error("test exception");
                return 0; // unreachable
            });
        } catch (const std::runtime_error& e) {
            caught = true;
            RC_ASSERT(std::string(e.what()) == "test exception");
        }

        // Exception must have been caught
        RC_ASSERT(caught);
        // Depth must be restored despite the exception
        RC_ASSERT(ctx.depth() == depth_before);
    });
}

/// Test: with_assumptions preserves depth with void callable on normal return
static void test_with_assumptions_void_preserves_depth_normal() {
    TEST_CASE("void with_assumptions preserves depth on normal return");

    rc::check("void with_assumptions preserves context depth on normal return", []() {
        AssumptionContext ctx;

        int initial_pushes = rc::gen::inRange(0, 5);
        for (int i = 0; i < initial_pushes; ++i) {
            ctx.push();
        }
        int depth_before = ctx.depth();

        auto decls = random_decls();

        // Call with void callable
        with_assumptions(ctx, decls, [&]() {
            RC_ASSERT(ctx.depth() == depth_before + 1);
        });

        RC_ASSERT(ctx.depth() == depth_before);
    });
}

/// Test: with_assumptions preserves depth with void callable on exception
static void test_with_assumptions_void_preserves_depth_exception() {
    TEST_CASE("void with_assumptions preserves depth on exception");

    rc::check("void with_assumptions preserves context depth when void callable throws", []() {
        AssumptionContext ctx;

        int initial_pushes = rc::gen::inRange(0, 5);
        for (int i = 0; i < initial_pushes; ++i) {
            ctx.push();
        }
        int depth_before = ctx.depth();

        auto decls = random_decls();

        bool caught = false;
        try {
            with_assumptions(ctx, decls, [&]() {
                RC_ASSERT(ctx.depth() == depth_before + 1);
                throw std::logic_error("void exception");
            });
        } catch (const std::logic_error& e) {
            caught = true;
            RC_ASSERT(std::string(e.what()) == "void exception");
        }

        RC_ASSERT(caught);
        RC_ASSERT(ctx.depth() == depth_before);
    });
}

/// Test: with_assumptions with vector overload preserves depth
static void test_with_assumptions_vector_preserves_depth() {
    TEST_CASE("vector overload preserves depth");

    rc::check("with_assumptions (vector overload) preserves context depth on normal and exception paths", []() {
        AssumptionContext ctx;

        int initial_pushes = rc::gen::inRange(0, 3);
        for (int i = 0; i < initial_pushes; ++i) {
            ctx.push();
        }
        int depth_before = ctx.depth();

        // Build a vector of declarations
        std::vector<AssumptionDecl> decls = random_decls();

        // Normal path
        int result = with_assumptions(ctx, decls, [&]() -> int {
            RC_ASSERT(ctx.depth() == depth_before + 1);
            return 99;
        });
        RC_ASSERT(ctx.depth() == depth_before);
        RC_ASSERT(result == 99);

        // Exception path
        bool caught = false;
        try {
            with_assumptions(ctx, decls, [&]() -> int {
                throw std::runtime_error("vec exception");
                return 0;
            });
        } catch (const std::runtime_error&) {
            caught = true;
        }
        RC_ASSERT(caught);
        RC_ASSERT(ctx.depth() == depth_before);
    });
}

/// Test: assumptions applied inside with_assumptions are visible during callable
static void test_with_assumptions_applies_declarations() {
    TEST_CASE("declarations are applied inside scope");

    rc::check("with_assumptions applies declarations that are visible inside the callable", []() {
        AssumptionContext ctx;
        std::string var = random_var_name();
        Domain dom = random_domain();

        int depth_before = ctx.depth();

        with_assumptions(ctx, {AssumptionDecl::make_domain(var, dom)}, [&]() {
            // The domain should be visible inside
            RC_ASSERT(ctx.has_domain(var, dom));
        });

        // After with_assumptions, the declaration should NOT be visible
        // (unless it was already declared in the parent scope)
        RC_ASSERT(ctx.depth() == depth_before);
    });
}


int main() {
    test_with_assumptions_preserves_depth_normal();
    test_with_assumptions_preserves_depth_exception();
    test_with_assumptions_void_preserves_depth_normal();
    test_with_assumptions_void_preserves_depth_exception();
    test_with_assumptions_vector_preserves_depth();
    test_with_assumptions_applies_declarations();

    return TEST_REPORT();
}
