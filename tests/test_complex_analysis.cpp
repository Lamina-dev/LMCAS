#include "complex_analysis.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"

using namespace lamina;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
static std::shared_ptr<SymbolicExpr> cmplx(int a, int b) {
    return lamina::detail::make_expression_ptr(
        SymbolicFactory::create_complex(
            lamina::detail::node(num(a)), lamina::detail::node(num(b))));
}

int main() {
    // ---- residue of 1/(z-a) at z=a is 1 ----
    {
        auto z = SymbolicExpr::variable("z");
        auto a = num(2);
        auto denom = SymbolicExpr::add(z, SymbolicExpr::multiply(num(-1), a));
        auto f = SymbolicExpr::divide(num(1), denom);
        auto res = residue(f, "z", a, 1);
        EXPECT_TRUE(res != nullptr, "residue not null");
        EXPECT_EQ_EXPR(res->simplify(), num(1), "Res(1/(z-2), z=2) = 1");
    }

    // ---- real_part / imag_part of (3 + 4i) ----
    {
        auto z = cmplx(3, 4);
        EXPECT_EQ_EXPR(real_part(z), num(3), "Re(3+4i) = 3");
        EXPECT_EQ_EXPR(imag_part(z), num(4), "Im(3+4i) = 4");

        auto re_checked = real_part_checked(z);
        auto im_checked = imag_part_checked(z);
        EXPECT_TRUE(re_checked.has_value(), "checked Re(3+4i) succeeds");
        EXPECT_TRUE(im_checked.has_value(), "checked Im(3+4i) succeeds");
        if (re_checked && im_checked) {
            EXPECT_EQ_EXPR(re_checked.value(), num(3), "checked Re(3+4i) = 3");
            EXPECT_EQ_EXPR(im_checked.value(), num(4), "checked Im(3+4i) = 4");
        }
    }

    // ---- real/imag of (2+i)*(1+i) = 2 + 2i + i + i^2 = 1 + 3i ----
    {
        auto prod = SymbolicExpr::multiply(cmplx(2,1), cmplx(1,1));
        EXPECT_EQ_EXPR(real_part(prod), num(1), "Re((2+i)(1+i)) = 1");
        EXPECT_EQ_EXPR(imag_part(prod), num(3), "Im((2+i)(1+i)) = 3");
    }

    // ---- conjugate of (3+4i) = 3-4i ----
    {
        auto z = cmplx(3, 4);
        auto c = conjugate(z);
        EXPECT_EQ_EXPR(real_part(c), num(3), "Re(conj(3+4i)) = 3");
        EXPECT_EQ_EXPR(imag_part(c)->simplify(), num(-4), "Im(conj(3+4i)) = -4");

        auto checked = conjugate_checked(z);
        EXPECT_TRUE(checked.has_value(), "checked conjugate(3+4i) succeeds");
        if (checked) {
            EXPECT_EQ_EXPR(real_part(checked.value()), num(3),
                           "checked Re(conj(3+4i)) = 3");
            EXPECT_EQ_EXPR(imag_part(checked.value())->simplify(), num(-4),
                           "checked Im(conj(3+4i)) = -4");
        }
    }

    // ---- checked complex part APIs reject invalid expressions ----
    {
        auto null_re = real_part_checked(nullptr);
        auto null_im = imag_part_checked(nullptr);
        auto null_conj = conjugate_checked(nullptr);
        EXPECT_TRUE(!null_re && null_re.error().code == CasErrc::InvalidArgument,
                    "checked real_part rejects null expression");
        EXPECT_TRUE(!null_im && null_im.error().code == CasErrc::InvalidArgument,
                    "checked imag_part rejects null expression");
        EXPECT_TRUE(!null_conj && null_conj.error().code == CasErrc::InvalidArgument,
                    "checked conjugate rejects null expression");

        std::shared_ptr<SymbolicExpr> empty_expr;
        auto empty_re = real_part_checked(empty_expr);
        EXPECT_TRUE(!empty_re && empty_re.error().code == CasErrc::InvalidArgument,
                    "checked real_part rejects null expression");

        EXPECT_TRUE(real_part(nullptr) == nullptr,
                    "legacy real_part unwraps invalid input to nullptr");
        EXPECT_TRUE(imag_part(nullptr) == nullptr,
                    "legacy imag_part unwraps invalid input to nullptr");
        EXPECT_TRUE(conjugate(nullptr) == nullptr,
                    "legacy conjugate unwraps invalid input to nullptr");
    }

    // ---- checked complex part APIs observe computation context cancellation ----
    {
        CancellationToken token;
        token.cancel();
        ComputationContext context({}, token);
        auto cancelled = real_part_checked(cmplx(1, 2), context);
        EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                    "checked real_part observes cancelled context");
    }

    // ---- checked complex part APIs do not fabricate splits for complex function arguments ----
    {
        auto i = cmplx(0, 1);
        auto exp_i = SymbolicExpr::exp(i);
        auto re = real_part_checked(exp_i);
        auto im = imag_part_checked(exp_i);
        auto conj = conjugate_checked(exp_i);
        EXPECT_TRUE(!re && re.error().code == CasErrc::Inconclusive,
                    "checked real_part rejects unsupported function of complex argument");
        EXPECT_TRUE(!im && im.error().code == CasErrc::Inconclusive,
                    "checked imag_part rejects unsupported function of complex argument");
        EXPECT_TRUE(!conj && conj.error().code == CasErrc::Inconclusive,
                    "checked conjugate rejects unsupported function of complex argument");
    }

    // ---- i^2 = -1 ----
    {
        auto i = cmplx(0, 1);
        auto i2 = SymbolicExpr::multiply(i, i)->simplify();
        EXPECT_EQ_EXPR(real_part(i2), num(-1), "Re(i^2) = -1");
        EXPECT_EQ_EXPR(imag_part(i2)->simplify(), num(0), "Im(i^2) = 0");
    }

    // ---- is_analytic: f(z) = z^2 is analytic ----
    {
        auto z = SymbolicExpr::variable("z");
        auto f = SymbolicExpr::multiply(z, z);
        EXPECT_TRUE(is_analytic(f, "z"), "z^2 is analytic");

        auto checked = is_analytic_checked(f, "z");
        EXPECT_TRUE(checked.has_value(), "checked is_analytic succeeds for z^2");
        if (checked) {
            EXPECT_TRUE(checked.value(), "checked z^2 is analytic");
        }
    }

    // ---- checked analyticity reports unsupported dependent functions explicitly ----
    {
        auto z = SymbolicExpr::variable("z");
        auto sin_z = SymbolicExpr::sin(z);
        auto checked = is_analytic_checked(sin_z, "z");
        EXPECT_TRUE(!checked && checked.error().code == CasErrc::Inconclusive,
                    "checked is_analytic reports dependent function domain as inconclusive");

        auto x = SymbolicExpr::variable("x");
        auto sin_x = SymbolicExpr::sin(x);
        auto constant_in_z = is_analytic_checked(sin_x, "z");
        EXPECT_TRUE(constant_in_z.has_value() && constant_in_z.value(),
                    "checked is_analytic treats functions independent of z as constants");
    }

    // ---- checked residue and Cauchy APIs reject invalid domains explicitly ----
    {
        auto z = SymbolicExpr::variable("z");
        auto z0 = num(2);
        auto denom = SymbolicExpr::add(z, SymbolicExpr::multiply(num(-1), z0));
        auto f = SymbolicExpr::divide(num(1), denom);

        auto checked_residue = residue_checked(f, "z", z0, 1);
        EXPECT_TRUE(checked_residue.has_value(),
                    "checked residue succeeds for simple pole");
        if (checked_residue) {
            EXPECT_EQ_EXPR(checked_residue.value()->simplify(), num(1),
                           "checked Res(1/(z-2), z=2) = 1");
        }

        auto bad_order = residue_checked(f, "z", z0, 0);
        EXPECT_TRUE(!bad_order && bad_order.error().code == CasErrc::InvalidArgument,
                    "checked residue rejects order 0 instead of returning zero");

        auto empty_var = cauchy_integral_checked(f, "", z0, 1);
        EXPECT_TRUE(!empty_var && empty_var.error().code == CasErrc::InvalidArgument,
                    "checked Cauchy integral rejects empty variable");

        auto cauchy = cauchy_integral_checked(num(1), "z", num(0), 1);
        EXPECT_TRUE(cauchy.has_value(),
                    "checked Cauchy integral constructs formula for constant analytic part");
        if (cauchy) {
            EXPECT_TRUE(cauchy.value() != nullptr,
                        "checked Cauchy integral result is non-null");
        }
    }

    // ---- checked complex-analysis APIs observe context errors ----
    {
        auto z = SymbolicExpr::variable("z");
        auto z0 = num(0);

        CancellationToken token;
        token.cancel();
        ComputationContext cancelled_context({}, token);
        auto cancelled = residue_checked(z, "z", z0, 1, cancelled_context);
        EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                    "checked residue observes cancellation");

        ResourceLimits limits;
        limits.max_steps = 1;
        ComputationContext limited_context(limits);
        auto limited = cauchy_integral_checked(z, "z", z0, 1, limited_context);
        EXPECT_TRUE(!limited && limited.error().code == CasErrc::ResourceLimit,
                    "checked Cauchy integral observes exhausted step budget");

        auto invalid_analytic = is_analytic_checked(nullptr, "z");
        EXPECT_TRUE(!invalid_analytic &&
                    invalid_analytic.error().code == CasErrc::InvalidArgument,
                    "checked is_analytic rejects null expression");
    }

    return TEST_REPORT();
}
