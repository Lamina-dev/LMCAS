#include "differential_geometry.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"

using namespace lamina;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

int main() {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    std::vector<std::string> vars = {"x", "y", "z"};

    // ---- Lie derivative of f = x^2 along X = (1,0,0) is 2x ----
    {
        auto f = SymbolicExpr::multiply(x, x);
        std::vector<std::shared_ptr<SymbolicExpr>> X = {num(1), num(0), num(0)};
        auto L = lie_derivative(f, X, vars, 1);
        auto expected = SymbolicExpr::multiply(num(2), x);
        EXPECT_EQ_EXPR(L->simplify(), expected->simplify(), "L_(1,0,0)(x^2) = 2x");
    }

    // ---- Lie derivative of f = x*y along X = (y, x, 0) = y*y + x*x ----
    {
        auto f = SymbolicExpr::multiply(x, y);
        std::vector<std::shared_ptr<SymbolicExpr>> X = {y, x, num(0)};
        auto L = lie_derivative(f, X, vars, 1);
        // d/dx(xy)=y, d/dy(xy)=x => y*y + x*x
        auto expected = SymbolicExpr::add(SymbolicExpr::multiply(y, y), SymbolicExpr::multiply(x, x));
        EXPECT_EQ_EXPR(L->simplify(), expected->simplify(), "L_(y,x,0)(xy) = y^2 + x^2");
    }

    // ---- exterior derivative of 0-form f = x^2 + y^2 + z^2 equals gradient (2x,2y,2z) ----
    {
        auto f = SymbolicExpr::add(SymbolicExpr::multiply(x, x),
                 SymbolicExpr::add(SymbolicExpr::multiply(y, y), SymbolicExpr::multiply(z, z)));
        auto d = exterior_derivative({f}, 0, vars);
        EXPECT_TRUE(d.size() == 3, "d(0-form) has 3 components");
        EXPECT_EQ_EXPR(d[0]->simplify(), SymbolicExpr::multiply(num(2), x)->simplify(), "d f / dx = 2x");
        EXPECT_EQ_EXPR(d[1]->simplify(), SymbolicExpr::multiply(num(2), y)->simplify(), "d f / dy = 2y");
        EXPECT_EQ_EXPR(d[2]->simplify(), SymbolicExpr::multiply(num(2), z)->simplify(), "d f / dz = 2z");
    }

    // ---- d^2 = 0: exterior derivative of (d of a 0-form) is the zero 2-form ----
    {
        auto f = SymbolicExpr::multiply(x, y); // smooth scalar
        auto d1 = exterior_derivative({f}, 0, vars); // 1-form coeffs
        auto d2 = exterior_derivative(d1, 1, vars);   // 2-form coeffs, should all be 0
        EXPECT_TRUE(d2.size() == 3, "2-form has 3 components in 3D");
        bool all_zero = true;
        for (auto& c : d2) {
            auto simplified = c->simplify();
            if (!(lamina::detail::node(simplified) &&
                  lamina::detail::node(simplified)->is_zero())) {
                all_zero = false;
            }
        }
        EXPECT_TRUE(all_zero, "d^2(f) = 0");
    }

    // ---- arbitrary-degree exterior derivative uses lexicographic basis order ----
    {
        // omega = z dx^dy + 0 dx^dz + x dy^dz.
        // d omega = (d(x)/dx - d(0)/dy + d(z)/dz) dx^dy^dz = 2 volume.
        auto d2 = exterior_derivative_checked({z, num(0), x}, 2, vars);
        EXPECT_TRUE(d2.has_value(), "checked exterior derivative supports 2-forms");
        if (d2) {
            EXPECT_TRUE(d2.value().size() == 1,
                        "2-form in three dimensions differentiates to one component");
            if (d2.value().size() == 1) {
                EXPECT_EQ_EXPR(d2.value()[0], num(2),
                               "2-form exterior derivative has the alternating-sign coefficient");
            }
        }

        // alpha = x dy + xy dz has d alpha = dx^dy + y dx^dz + x dy^dz.
        auto d_alpha = exterior_derivative_checked({num(0), x, SymbolicExpr::multiply(x, y)},
                                                   1,
                                                   vars);
        EXPECT_TRUE(d_alpha.has_value(), "nontrivial 1-form exterior derivative succeeds");
        if (d_alpha) {
            auto d_squared = exterior_derivative_checked(d_alpha.value(), 2, vars);
            EXPECT_TRUE(d_squared.has_value(), "d of a 2-form succeeds");
            if (d_squared && d_squared.value().size() == 1) {
                EXPECT_TRUE(lamina::detail::node(d_squared.value()[0]->simplify())->is_zero(),
                            "d squared is zero for a nontrivial 1-form");
            }
        }

        auto top = exterior_derivative_checked({x}, 3, vars);
        EXPECT_TRUE(top.has_value() && top.value().empty(),
                    "exterior derivative of a top-degree form is the zero higher form");

        std::vector<std::string> vars4 = {"x", "y", "z", "w"};
        // Input basis: 01, 02, 03, 12, 13, 23. Output basis: 012, 013, 023, 123.
        auto d4 = exterior_derivative_checked(
            {z, num(0), num(0), x, num(0), x}, 2, vars4);
        EXPECT_TRUE(d4.has_value() && d4.value().size() == 4,
                    "2-form in four dimensions has four derivative components");
        if (d4 && d4.value().size() == 4) {
            EXPECT_EQ_EXPR(d4.value()[0], num(2), "4D component 012 is ordered correctly");
            EXPECT_EQ_EXPR(d4.value()[1], num(0), "4D component 013 is ordered correctly");
            EXPECT_EQ_EXPR(d4.value()[2], num(1), "4D component 023 is ordered correctly");
            EXPECT_EQ_EXPR(d4.value()[3], num(0), "4D component 123 is ordered correctly");
        }
    }

    // ---- existing Christoffel smoke check still works ----
    {
        // metric_inverse expects a matrix; use a 2x2 identity metric.
        auto g = SymbolicExpr::matrix({{num(1), num(0)}, {num(0), num(1)}});
        auto mi = metric_inverse(g);
        EXPECT_TRUE(mi != nullptr, "metric_inverse returns non-null");
    }

    // ---- checked differential-geometry APIs preserve success paths ----
    {
        auto g = SymbolicExpr::matrix({{num(1), num(0)}, {num(0), num(1)}});
        std::vector<std::string> coords2 = {"x", "y"};
        auto inv = metric_inverse_checked(g);
        EXPECT_TRUE(inv.has_value(), "checked metric_inverse succeeds for square metric");
        if (inv) {
            EXPECT_TRUE(inv.value() != nullptr, "checked metric_inverse returns non-null");
        }

        auto gamma1 = christoffel_first_kind_checked(g, coords2, 0, 0, 0);
        EXPECT_TRUE(gamma1.has_value(), "checked first-kind Christoffel succeeds");
        if (gamma1) {
            auto simplified_gamma = gamma1.value()->simplify();
            EXPECT_TRUE(lamina::detail::node(simplified_gamma) &&
                        lamina::detail::node(simplified_gamma)->is_zero(),
                        "flat first-kind Christoffel is zero");
        }

        auto gamma2 = christoffel_second_kind_checked(g, g, coords2, 0, 0, 0);
        EXPECT_TRUE(gamma2.has_value(), "checked second-kind Christoffel succeeds");

        auto riemann = riemann_curvature_tensor_checked(g, coords2, 0, 0, 0, 1);
        EXPECT_TRUE(riemann.has_value(), "checked Riemann tensor succeeds");

        std::vector<std::shared_ptr<SymbolicExpr>> X = {num(1), num(0), num(0)};
        auto lie = lie_derivative_checked(SymbolicExpr::multiply(x, x), X, vars, 1);
        EXPECT_TRUE(lie.has_value(), "checked Lie derivative succeeds");
        if (lie) {
            EXPECT_EQ_EXPR(lie.value()->simplify(),
                           SymbolicExpr::multiply(num(2), x)->simplify(),
                           "checked Lie derivative computes 2x");
        }

        auto d = exterior_derivative_checked({x}, 0, vars);
        EXPECT_TRUE(d.has_value(), "checked exterior derivative succeeds");
        if (d) {
            EXPECT_TRUE(d.value().size() == 3,
                        "checked exterior derivative of 0-form has 3 components");
        }
    }

    // ---- checked differential-geometry APIs reject invalid inputs explicitly ----
    {
        auto g2 = SymbolicExpr::matrix({{num(1), num(0)}}); // 1x2, not square
        auto bad_metric = metric_inverse_checked(g2);
        EXPECT_TRUE(!bad_metric && bad_metric.error().code == CasErrc::InvalidArgument,
                    "checked metric_inverse rejects non-square metric");

        auto scalar_metric = metric_inverse_checked(num(1));
        EXPECT_TRUE(!scalar_metric && scalar_metric.error().code == CasErrc::InvalidArgument,
                    "checked metric_inverse rejects non-matrix input");

        auto g = SymbolicExpr::matrix({{num(1), num(0)}, {num(0), num(1)}});
        std::vector<std::string> coords2 = {"x", "y"};
        auto singular_metric = SymbolicExpr::matrix({{num(1), num(2)}, {num(2), num(4)}});
        auto singular_inverse = metric_inverse_checked(singular_metric);
        EXPECT_TRUE(!singular_inverse &&
                        singular_inverse.error().code == CasErrc::DomainError,
                    "checked metric_inverse rejects singular metric as domain error");
        auto singular_riemann = riemann_curvature_tensor_checked(
            singular_metric, coords2, 0, 0, 0, 1);
        EXPECT_TRUE(!singular_riemann &&
                        singular_riemann.error().code == CasErrc::DomainError,
                    "checked Riemann rejects singular metric instead of fabricating zero");

        auto bad_index = christoffel_first_kind_checked(g, coords2, 2, 0, 0);
        EXPECT_TRUE(!bad_index && bad_index.error().code == CasErrc::InvalidArgument,
                    "checked Christoffel rejects out-of-range index");

        auto bad_coords = riemann_curvature_tensor_checked(g, {"x"}, 0, 0, 0, 0);
        EXPECT_TRUE(!bad_coords && bad_coords.error().code == CasErrc::InvalidArgument,
                    "checked Riemann rejects coordinate/metric dimension mismatch");

        auto bad_lie = lie_derivative_checked(x, {num(1)}, vars, 1);
        EXPECT_TRUE(!bad_lie && bad_lie.error().code == CasErrc::InvalidArgument,
                    "checked Lie derivative rejects vector dimension mismatch");

        auto bad_degree = exterior_derivative_checked({}, 4, vars);
        EXPECT_TRUE(!bad_degree && bad_degree.error().code == CasErrc::InvalidArgument,
                    "checked exterior derivative rejects degree above the dimension");

        auto negative_degree = exterior_derivative_checked({}, -1, vars);
        EXPECT_TRUE(!negative_degree && negative_degree.error().code == CasErrc::InvalidArgument,
                    "checked exterior derivative rejects negative degree");

        auto bad_count = exterior_derivative_checked({x}, 2, vars);
        EXPECT_TRUE(!bad_count && bad_count.error().code == CasErrc::InvalidArgument,
                    "checked exterior derivative validates high-degree coefficient count");

        auto duplicate_coords = exterior_derivative_checked({x}, 0, {"x", "x"});
        EXPECT_TRUE(!duplicate_coords &&
                        duplicate_coords.error().code == CasErrc::InvalidArgument,
                    "checked exterior derivative rejects duplicate coordinates");

        auto bad_form = exterior_derivative_checked({}, 0, vars);
        EXPECT_TRUE(!bad_form && bad_form.error().code == CasErrc::InvalidArgument,
                    "checked exterior derivative rejects missing 0-form coefficient");
    }

    // ---- checked differential-geometry APIs observe context cancellation and budgets ----
    {
        CancellationToken token;
        token.cancel();
        ComputationContext cancelled_context({}, token);
        auto cancelled = metric_inverse_checked(
            SymbolicExpr::matrix({{num(1), num(0)}, {num(0), num(1)}}),
            cancelled_context);
        EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                    "checked metric_inverse observes cancellation");

        ResourceLimits limits;
        limits.max_steps = 1;
        ComputationContext limited_context(limits);
        auto limited = exterior_derivative_checked({x}, 0, vars, limited_context);
        EXPECT_TRUE(!limited && limited.error().code == CasErrc::ResourceLimit,
                    "checked exterior derivative observes exhausted step budget");

        ResourceLimits expansion_limits;
        expansion_limits.max_expansion_terms = 2;
        ComputationContext expansion_context(expansion_limits);
        auto expansion_limited = exterior_derivative_checked({x, y, z},
                                                             1,
                                                             vars,
                                                             expansion_context);
        EXPECT_TRUE(!expansion_limited &&
                        expansion_limited.error().code == CasErrc::ResourceLimit,
                    "checked exterior derivative bounds combinatorial component growth");

        ResourceLimits zero_expansion_limits;
        zero_expansion_limits.max_expansion_terms = 0;
        ComputationContext zero_expansion_context(zero_expansion_limits);
        auto zero_expansion = exterior_derivative_checked({x}, 0, vars,
                                                          zero_expansion_context);
        EXPECT_TRUE(!zero_expansion && zero_expansion.error().code == CasErrc::ResourceLimit,
                    "checked exterior derivative enforces a zero expansion budget");

        ResourceLimits node_limits;
        node_limits.max_ast_nodes = 0;
        ComputationContext node_context(node_limits);
        auto node_limited = exterior_derivative_checked({x}, 0, vars, node_context);
        EXPECT_TRUE(!node_limited && node_limited.error().code == CasErrc::ResourceLimit,
                    "checked exterior derivative observes AST node budget");

        CancellationToken exterior_token;
        exterior_token.cancel();
        ComputationContext exterior_cancelled({}, exterior_token);
        auto cancelled_exterior = exterior_derivative_checked({x}, 0, vars,
                                                              exterior_cancelled);
        EXPECT_TRUE(!cancelled_exterior &&
                        cancelled_exterior.error().code == CasErrc::Cancelled,
                    "checked exterior derivative observes cancellation");
    }

    return TEST_REPORT();
}
