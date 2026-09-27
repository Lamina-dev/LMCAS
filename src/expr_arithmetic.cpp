#include "expr.hpp"
#include "complex_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expr_common.hpp"
#include <cmath>
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace LMCAS {

using namespace expr_detail::expr_common;

ExprResult add(const ExprPtr& lhs,
               const ExprPtr& rhs,
               ComputationContext& context) {
    return quantity_add(lhs, rhs, context);
}

ExprResult add(const ExprPtr& lhs, const ExprPtr& rhs) {
    ComputationContext context;
    return add(lhs, rhs, context);
}

ExprResult mul(const ExprPtr& lhs,
               const ExprPtr& rhs,
               ComputationContext& context) {
    return quantity_multiply(lhs, rhs, context);
}

ExprResult mul(const ExprPtr& lhs, const ExprPtr& rhs) {
    ComputationContext context;
    return mul(lhs, rhs, context);
}

ExprResult div(const ExprPtr& numerator,
               const ExprPtr& denominator,
               ComputationContext& context) {
    return quantity_divide(numerator, denominator, context);
}

ExprResult div(const ExprPtr& numerator, const ExprPtr& denominator) {
    ComputationContext context;
    return div(numerator, denominator, context);
}

ExprResult neg(const ExprPtr& expression, ComputationContext& context) {
    auto minus_one = integer(-1);
    if (!minus_one) { return ExprResult::failure(minus_one.error()); }
    return mul(minus_one.value(), expression, context);
}

ExprResult neg(const ExprPtr& expression) {
    ComputationContext context;
    return neg(expression, context);
}

ExprResult sub(const ExprPtr& lhs,
               const ExprPtr& rhs,
               ComputationContext& context) {
    return quantity_subtract(lhs, rhs, context);
}

ExprResult sub(const ExprPtr& lhs, const ExprPtr& rhs) {
    ComputationContext context;
    return sub(lhs, rhs, context);
}

ExprResult eq(const ExprPtr& lhs,
              const ExprPtr& rhs,
              ComputationContext& context) {
    if (!lhs || !rhs) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "eq operands cannot be null", kExprOperation);
    }
    auto left_dimension = dimension_of(*lhs);
    if (!left_dimension) { return ExprResult::failure(left_dimension.error()); }
    auto right_dimension = dimension_of(*rhs);
    if (!right_dimension) { return ExprResult::failure(right_dimension.error()); }
    if (left_dimension.value() != right_dimension.value()) { return integer(0); }
    auto left = comparison_value(lhs, context);
    if (!left) { return left; }
    auto right = comparison_value(rhs, context);
    if (!right) { return right; }
    return make_binary_expr(left.value(), right.value(), context, "eq", SymbolicExpr::eq);
}

ExprResult eq(const ExprPtr& lhs, const ExprPtr& rhs) {
    ComputationContext context;
    return eq(lhs, rhs, context);
}


ExprResult real(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(real_part_checked(expression, context),
                                    kRealOperation);
}

ExprResult real(const ExprPtr& expression) {
    ComputationContext context;
    return real(expression, context);
}

ExprResult imag(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(imag_part_checked(expression, context),
                                    kImagOperation);
}

ExprResult imag(const ExprPtr& expression) {
    ComputationContext context;
    return imag(expression, context);
}

ExprResult conj(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(conjugate_checked(expression, context),
                                    kConjOperation);
}

ExprResult conj(const ExprPtr& expression) {
    ComputationContext context;
    return conj(expression, context);
}

ExprResult abs(const ExprPtr& expression, ComputationContext& context) {
    auto step = context.consume_steps(1, kAbsOperation);
    if (!step) {
        return ExprResult::failure(step.error());
    }
    auto re = real(expression, context);
    if (!re) {
        return re;
    }
    auto im = imag(expression, context);
    if (!im) {
        return im;
    }
    try {
        auto re_squared = SymbolicExpr::power(re.value(), SymbolicExpr::number(2));
        auto im_squared = SymbolicExpr::power(im.value(), SymbolicExpr::number(2));
        auto sum = SymbolicExpr::add(re_squared, im_squared);
        auto result = SymbolicExpr::sqrt(sum)->simplify();
        if (!result || !LMCAS::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "complex absolute value construction failed",
                                      kAbsOperation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "complex absolute value allocation failed",
                                  kAbsOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InternalInvariant, error.what(),
                                  kAbsOperation);
    }
}

ExprResult abs(const ExprPtr& expression) {
    ComputationContext context;
    return abs(expression, context);
}
}
