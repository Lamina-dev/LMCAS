#include "test_common.hpp"
#include "expr.hpp"
#include "vector_calculus.hpp"
#include "internal/symbolic_ast.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> get_mat_entry(
    const std::shared_ptr<SymbolicExpr> &mat, size_t r, size_t c) {
    if (!mat || !LMCAS::detail::node(mat))
        return nullptr;
    auto mn = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(mat));
    if (!mn)
        return nullptr;
    auto node = mn->get(r, c);
    if (!node)
        return SymbolicExpr::number(0);
    return LMCAS::detail::make_expression_ptr(node);
}

static size_t get_mat_rows(const std::shared_ptr<SymbolicExpr> &mat) {
    if (!mat || !LMCAS::detail::node(mat))
        return 0;
    auto mn = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(mat));
    return mn ? mn->rows() : 0;
}

static size_t get_mat_cols(const std::shared_ptr<SymbolicExpr> &mat) {
    if (!mat || !LMCAS::detail::node(mat))
        return 0;
    auto mn = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(mat));
    return mn ? mn->cols() : 0;
}

TEST(VectorCalculusJacobianHessian, JacobianSquare) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f1 = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto f2 = SymbolicExpr::multiply(x, y);

    auto J = jacobian({f1, f2}, {"x", "y"});

    EXPECT_TRUE((J != nullptr)) << "Jacobian is not null";
    EXPECT_TRUE((get_mat_rows(J) == 2)) << "Jacobian has 2 rows";
    EXPECT_TRUE((get_mat_cols(J) == 2)) << "Jacobian has 2 columns";

    auto j00 = get_mat_entry(J, 0, 0);
    auto j00_at = j00->substitute("x", SymbolicExpr::number(3));
    j00_at = j00_at->simplify();
    auto v00 = test_numeric_eval(j00_at);
    EXPECT_TRUE((v00.has_value() && std::abs(*v00 - 6.0) < 1e-10)) << "J[0][0] = 2x => 6 at x=3";

    auto j01 = get_mat_entry(J, 0, 1);
    auto j01_at = j01->substitute("y", SymbolicExpr::number(5));
    j01_at = j01_at->simplify();
    auto v01 = test_numeric_eval(j01_at);
    EXPECT_TRUE((v01.has_value() && std::abs(*v01 - 10.0) < 1e-10)) << "J[0][1] = 2y => 10 at y=5";

    EXPECT_TRUE(test_expression_text(get_mat_entry(J, 1, 0), y)) << "J[1][0] = y";
    EXPECT_TRUE(test_expression_text(get_mat_entry(J, 1, 1), x)) << "J[1][1] = x";
}

TEST(VectorCalculusJacobianHessian, JacobianNonSquare) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f1 = SymbolicExpr::add(x, y);
    auto f2 = SymbolicExpr::multiply(x, y);
    auto f3 = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto J = jacobian({f1, f2, f3}, {"x", "y"});

    EXPECT_TRUE((get_mat_rows(J) == 3)) << "3x2 Jacobian has 3 rows";
    EXPECT_TRUE((get_mat_cols(J) == 2)) << "3x2 Jacobian has 2 columns";

    auto j00 = get_mat_entry(J, 0, 0);
    auto v00 = test_numeric_eval(j00);
    EXPECT_TRUE((v00.has_value() && std::abs(*v00 - 1.0) < 1e-10)) << "J[0][0] = 1";

    auto j01 = get_mat_entry(J, 0, 1);
    auto v01 = test_numeric_eval(j01);
    EXPECT_TRUE((v01.has_value() && std::abs(*v01 - 1.0) < 1e-10)) << "J[0][1] = 1";

    EXPECT_TRUE(test_expression_text(get_mat_entry(J, 1, 0), y)) << "J[1][0] = y";
    EXPECT_TRUE(test_expression_text(get_mat_entry(J, 1, 1), x)) << "J[1][1] = x";

    auto j20 = get_mat_entry(J, 2, 0);
    auto j20_at = j20->substitute("x", SymbolicExpr::number(4));
    j20_at = j20_at->simplify();
    auto v20 = test_numeric_eval(j20_at);
    EXPECT_TRUE((v20.has_value() && std::abs(*v20 - 8.0) < 1e-10)) << "J[2][0] = 2x => 8 at x=4";

    auto j21 = get_mat_entry(J, 2, 1);
    EXPECT_TRUE((j21->is_zero())) << "J[2][1] = 0";
}

TEST(VectorCalculusJacobianHessian, JacobianWide) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto f = SymbolicExpr::add(x,
                               SymbolicExpr::add(
                                   SymbolicExpr::multiply(SymbolicExpr::number(2), y),
                                   SymbolicExpr::multiply(SymbolicExpr::number(3), z)));

    auto J = jacobian({f}, {"x", "y", "z"});

    EXPECT_TRUE((get_mat_rows(J) == 1)) << "1x3 Jacobian has 1 row";
    EXPECT_TRUE((get_mat_cols(J) == 3)) << "1x3 Jacobian has 3 columns";

    auto v0 = test_numeric_eval(get_mat_entry(J, 0, 0));
    auto v1 = test_numeric_eval(get_mat_entry(J, 0, 1));
    auto v2 = test_numeric_eval(get_mat_entry(J, 0, 2));
    EXPECT_TRUE((v0.has_value() && std::abs(*v0 - 1.0) < 1e-10)) << "J[0][0] = 1";
    EXPECT_TRUE((v1.has_value() && std::abs(*v1 - 2.0) < 1e-10)) << "J[0][1] = 2";
    EXPECT_TRUE((v2.has_value() && std::abs(*v2 - 3.0) < 1e-10)) << "J[0][2] = 3";
}

TEST(VectorCalculusJacobianHessian, HessianQuadratic) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2),
                                   SymbolicExpr::multiply(x, y)),
            SymbolicExpr::multiply(SymbolicExpr::number(3),
                                   SymbolicExpr::power(y, SymbolicExpr::number(2)))));

    auto H = hessian(f, {"x", "y"});

    EXPECT_TRUE((H != nullptr)) << "Hessian is not null";
    EXPECT_TRUE((get_mat_rows(H) == 2)) << "Hessian has 2 rows";
    EXPECT_TRUE((get_mat_cols(H) == 2)) << "Hessian has 2 columns";

    auto v00 = test_numeric_eval(get_mat_entry(H, 0, 0));
    auto v01 = test_numeric_eval(get_mat_entry(H, 0, 1));
    auto v10 = test_numeric_eval(get_mat_entry(H, 1, 0));
    auto v11 = test_numeric_eval(get_mat_entry(H, 1, 1));

    EXPECT_TRUE((v00.has_value() && std::abs(*v00 - 2.0) < 1e-10)) << "H[0][0] = 2";
    EXPECT_TRUE((v01.has_value() && std::abs(*v01 - 2.0) < 1e-10)) << "H[0][1] = 2";
    EXPECT_TRUE((v10.has_value() && std::abs(*v10 - 2.0) < 1e-10)) << "H[1][0] = 2 (symmetric)";
    EXPECT_TRUE((v11.has_value() && std::abs(*v11 - 6.0) < 1e-10)) << "H[1][1] = 6";
}

TEST(VectorCalculusJacobianHessian, HessianSingleVar) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::power(x, SymbolicExpr::number(3));

    auto H = hessian(f, {"x"});

    EXPECT_TRUE((get_mat_rows(H) == 1)) << "1x1 Hessian has 1 row";
    EXPECT_TRUE((get_mat_cols(H) == 1)) << "1x1 Hessian has 1 column";

    auto h00 = get_mat_entry(H, 0, 0);
    auto h00_at = h00->substitute("x", SymbolicExpr::number(2));
    h00_at = h00_at->simplify();
    auto v = test_numeric_eval(h00_at);
    EXPECT_TRUE((v.has_value() && std::abs(*v - 12.0) < 1e-10)) << "H[0][0] = 6x => 12 at x=2";
}

TEST(VectorCalculusJacobianHessian, Hessian3dDiagonal) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto f = SymbolicExpr::multiply(
        SymbolicExpr::multiply(x, y), z);

    auto H = hessian(f, {"x", "y", "z"});

    EXPECT_TRUE((get_mat_rows(H) == 3)) << "3D Hessian has 3 rows";
    EXPECT_TRUE((get_mat_cols(H) == 3)) << "3D Hessian has 3 columns";

    auto d00 = test_numeric_eval(get_mat_entry(H, 0, 0));
    auto d11 = test_numeric_eval(get_mat_entry(H, 1, 1));
    auto d22 = test_numeric_eval(get_mat_entry(H, 2, 2));
    EXPECT_TRUE((d00.has_value() && std::abs(*d00) < 1e-10)) << "H[0][0] = 0";
    EXPECT_TRUE((d11.has_value() && std::abs(*d11) < 1e-10)) << "H[1][1] = 0";
    EXPECT_TRUE((d22.has_value() && std::abs(*d22) < 1e-10)) << "H[2][2] = 0";
}

TEST(VectorCalculusJacobianHessian, Hessian3dMixed) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");

    auto f = SymbolicExpr::multiply(
        SymbolicExpr::multiply(x, y), z);

    auto H = hessian(f, {"x", "y", "z"});

    auto h01 = get_mat_entry(H, 0, 1);
    auto h01_at = h01->substitute("x", SymbolicExpr::number(2));
    h01_at = h01_at->substitute("y", SymbolicExpr::number(3));
    h01_at = h01_at->substitute("z", SymbolicExpr::number(5));
    h01_at = h01_at->simplify();
    auto v01 = test_numeric_eval(h01_at);
    EXPECT_TRUE((v01.has_value() && std::abs(*v01 - 5.0) < 1e-10)) << "H[0][1] = z => 5 at z=5";

    auto h02 = get_mat_entry(H, 0, 2);
    auto h02_at = h02->substitute("x", SymbolicExpr::number(2));
    h02_at = h02_at->substitute("y", SymbolicExpr::number(3));
    h02_at = h02_at->substitute("z", SymbolicExpr::number(5));
    h02_at = h02_at->simplify();
    auto v02 = test_numeric_eval(h02_at);
    EXPECT_TRUE((v02.has_value() && std::abs(*v02 - 3.0) < 1e-10)) << "H[0][2] = y => 3 at y=3";

    auto h12 = get_mat_entry(H, 1, 2);
    auto h12_at = h12->substitute("x", SymbolicExpr::number(2));
    h12_at = h12_at->substitute("y", SymbolicExpr::number(3));
    h12_at = h12_at->substitute("z", SymbolicExpr::number(5));
    h12_at = h12_at->simplify();
    auto v12 = test_numeric_eval(h12_at);
    EXPECT_TRUE((v12.has_value() && std::abs(*v12 - 2.0) < 1e-10)) << "H[1][2] = x => 2 at x=2";

    auto h10 = get_mat_entry(H, 1, 0);
    auto h10_at = h10->substitute("z", SymbolicExpr::number(5));
    h10_at = h10_at->simplify();
    auto v10 = test_numeric_eval(h10_at);
    EXPECT_TRUE((v10.has_value() && std::abs(*v10 - 5.0) < 1e-10)) << "H[1][0] = z => 5 (symmetric)";
}

TEST(VectorCalculusJacobianHessian, JacobianCheckedShape) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    auto checked_J = jacobian_checked({f}, {"x", "y"});
    ASSERT_TRUE((checked_J.has_value())) << "checked Jacobian succeeds";
    if (checked_J) {
        EXPECT_TRUE((get_mat_rows(checked_J.value()) == 1)) << "checked Jacobian has one row";
        EXPECT_TRUE((get_mat_cols(checked_J.value()) == 2)) << "checked Jacobian has two columns";
    }
}

TEST(VectorCalculusJacobianHessian, HessianCheckedShape) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    auto checked_H = hessian_checked(f, {"x", "y"});
    ASSERT_TRUE((checked_H.has_value())) << "checked Hessian succeeds";
    if (checked_H) {
        EXPECT_TRUE((get_mat_rows(checked_H.value()) == 2)) << "checked Hessian has two rows";
        EXPECT_TRUE((get_mat_cols(checked_H.value()) == 2)) << "checked Hessian has two columns";
    }
}

TEST(VectorCalculusJacobianHessian, JacobianCheckedInvalidFunctions) {
    auto empty_functions = jacobian_checked({}, {"x"});
    EXPECT_TRUE((!empty_functions.has_value())) << "checked Jacobian rejects empty function list";
    EXPECT_TRUE((empty_functions.error().code == CasErrc::InvalidArgument)) << "checked Jacobian reports InvalidArgument for empty functions";

    std::shared_ptr<SymbolicExpr> null_root;
    auto null_function = jacobian_checked({null_root}, {"x"});
    EXPECT_TRUE((!null_function.has_value())) << "checked Jacobian rejects null function";
    EXPECT_TRUE((null_function.error().code == CasErrc::InvalidArgument)) << "checked Jacobian reports InvalidArgument for null function";
}

TEST(VectorCalculusJacobianHessian, HessianCheckedEmptyVariables) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    auto empty_vars = hessian_checked(f, {});
    EXPECT_TRUE((!empty_vars.has_value())) << "checked Hessian rejects empty variables";
    EXPECT_TRUE((empty_vars.error().code == CasErrc::InvalidArgument)) << "checked Hessian reports InvalidArgument for empty variables";
}

TEST(VectorCalculusJacobianHessian, JacobianCheckedCancellation) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = jacobian_checked({f}, {"x"}, cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked Jacobian observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked Jacobian reports Cancelled";
}

TEST(VectorCalculusJacobianHessian, HessianCheckedBudget) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = hessian_checked(f, {"x"}, limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked Hessian observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked Hessian reports ResourceLimit";
}

TEST(VectorCalculusJacobianHessian, UnsupportedDerivativeErrors) {
    const auto floor_result = LMCAS::floor(SymbolicExpr::variable("x"));
    ASSERT_TRUE(floor_result.has_value());
    const auto& floor_x = floor_result.value();
    auto jac = jacobian_checked({floor_x}, {"x"});
    ASSERT_FALSE(jac.has_value());
    EXPECT_EQ(jac.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(jac.error().operation, "jacobian");

    auto hess = hessian_checked(floor_x, {"x"});
    ASSERT_FALSE(hess.has_value());
    EXPECT_EQ(hess.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(hess.error().operation, "hessian");

    const auto integral = SymbolicExpr::make_integral(floor_x, "x");
    auto second = hessian_checked(integral, {"x"});
    ASSERT_FALSE(second.has_value());
    EXPECT_EQ(second.error().code, CasErrc::UnsupportedExpression);
    EXPECT_EQ(second.error().operation, "hessian");
}
