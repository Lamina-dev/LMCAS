#include "expr.hpp"
#include "symbolic_ast.hpp"
#include "internal/expr_common.hpp"
#include <cmath>
#include <exception>
#include <memory>
#include <string>
#include <utility>

namespace LMCAS {

using namespace expr_detail::expr_common;

ExprResult sym(const std::string& name) {
    if (name.empty()) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "symbol name cannot be empty", kSymOperation);
    }
    if (is_reserved_symbol_name(name)) {
        if (is_imaginary_unit_name(name)) {
            return expression_failure(CasErrc::InvalidArgument,
                                      "imaginary unit symbol is reserved",
                                      kSymOperation);
        }
        return expression_failure(CasErrc::InvalidArgument,
                                  "reserved mathematical constants cannot be shadowed",
                                  kSymOperation);
    }
    try {
        auto expression = SymbolicExpr::variable(name);
        if (!expression) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "symbol factory returned null", kSymOperation);
        }
        return ExprResult::success(std::move(expression));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "symbol allocation failed", kSymOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kSymOperation);
    }
}

ExprResult integer(long long value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "integer allocation failed", kIntegerOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kIntegerOperation);
    }
}

ExprResult integer(const BigInt& value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "integer allocation failed", kIntegerOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kIntegerOperation);
    }
}

ExprResult rational(const Rational& value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "rational allocation failed", kRationalOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kRationalOperation);
    }
}

ExprResult approx_real(double value) {
    if (!std::isfinite(value)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "approximate real must be finite", kApproxOperation);
    }
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "approximate real allocation failed", kApproxOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kApproxOperation);
    }
}

ExprResult constant_symbol(const char* name) {
    try {
        auto expression = SymbolicExpr::variable(name);
        if (!expression || !LMCAS::detail::node(expression)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "constant factory returned null",
                                      kConstantOperation);
        }
        return ExprResult::success(std::move(expression));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "constant allocation failed",
                                  kConstantOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kConstantOperation);
    }
}

ExprResult pi() {
    return constant_symbol("pi");
}

ExprResult e() {
    return constant_symbol("e");
}

ExprResult phi() {
    return constant_symbol("phi");
}

ExprResult I() {
    return imaginary_unit();
}

ExprResult imaginary_unit() {
    auto zero = integer(0);
    if (!zero) return ExprResult::failure(zero.error());
    auto one = integer(1);
    if (!one) return ExprResult::failure(one.error());
    return complex(zero.value(), one.value());
}

ExprResult complex(ExprPtr real, ExprPtr imag) {
    if (!real || !imag) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "complex expression parts cannot be null",
                                  kComplexOperation);
    }
    try {
        auto node = SymbolicFactory::create_complex(LMCAS::detail::node(real),
                                                    LMCAS::detail::node(imag));
        return ExprResult::success(LMCAS::detail::make_expression_ptr(std::move(node)));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "complex expression allocation failed",
                                  kComplexOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kComplexOperation);
    }
}

} // namespace LMCAS
