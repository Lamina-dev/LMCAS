#include "limit_result.hpp"
/**
 * @file complex_analysis.cpp
 * @brief 复变函数分析实现。
 */
#include "complex_analysis.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include <cmath>
#include <optional>

namespace LMCAS {

namespace {

Result<void> validate_complex_expr_input(const std::shared_ptr<SymbolicExpr>& expr,
                                         ComputationContext& context,
                                         const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "expression cannot be null", operation);
    }
    return Result<void>::success();
}

Result<void> validate_complex_expr_point_input(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!expr || !LMCAS::detail::node(expr) || !z0 || !LMCAS::detail::node(z0)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "complex-analysis expressions cannot be null",
                                     operation);
    }
    if (z.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "complex variable name cannot be empty",
                                     operation);
    }
    if (order < 1 || order > 13) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "complex-analysis order must be between 1 and 13",
                                     operation);
    }
    return Result<void>::success();
}

class FunctionDependencyFinder final : public detail::RecursiveSymbolicVisitor {
public:
    explicit FunctionDependencyFinder(const std::string& variable)
        : variable_(variable) {}

    bool found() const noexcept { return found_; }

    void visit(const FunctionNode& node) override {
        visit_function(node.arguments());
    }

    void visit(const UninterpretedFunctionNode& node) override {
        visit_function(node.arguments());
    }

    void visit(const SummationNode& node) override { visit_binder(node); }
    void visit(const ProductNode& node) override { visit_binder(node); }
    void visit(const IntegralNode& node) override { visit_binder(node); }
    void visit(const TransformNode& node) override { visit_binder(node); }
    void visit(const QuantifierNode& node) override { visit_binder(node); }
    void visit(const SetBuilderNode& node) override { visit_binder(node); }
    void visit(const LimitNode& node) override { visit_binder(node); }

protected:
    void visit_child(const detail::SymbolicNodePtr& child) override {
        if (!found_) RecursiveSymbolicVisitor::visit_child(child);
    }

private:
    void visit_function(const std::vector<detail::SymbolicNodePtr>& arguments) {
        if (!is_bound_) {
            for (const auto& argument : arguments) {
                if (expression_depends_on_variable(argument, variable_)) {
                    found_ = true;
                    return;
                }
            }
        }
        visit_children(arguments);
    }

    void visit_binder(const SymbolicNode& node) {
        const auto binder = detail::binder_view(node);
        visit_children(binder->outside_scope);
        const bool previous = is_bound_;
        is_bound_ = is_bound_ || binder->bound_name == variable_;
        visit_child(binder->scoped_body);
        is_bound_ = previous;
    }

    const std::string& variable_;
    bool is_bound_ = false;
    bool found_ = false;
};

class ExplicitComplexFinder final : public detail::RecursiveSymbolicVisitor {
public:
    bool found() const noexcept { return found_; }

    void visit(const VariableNode& node) override {
        found_ = found_ || detail::is_imaginary_unit_name(node.name());
    }

    void visit(const ComplexNode&) override { found_ = true; }

protected:
    void visit_child(const detail::SymbolicNodePtr& child) override {
        if (!found_) RecursiveSymbolicVisitor::visit_child(child);
    }

private:
    bool found_ = false;
};

bool contains_explicit_complex(const detail::SymbolicNodePtr& node) {
    if (!node) return false;
    ExplicitComplexFinder finder;
    node->accept(finder);
    return finder.found();
}

class FunctionOfExplicitComplexFinder final : public detail::RecursiveSymbolicVisitor {
public:
    bool found() const noexcept { return found_; }

    void visit(const FunctionNode& node) override {
        visit_function(node.arguments());
    }

    void visit(const UninterpretedFunctionNode& node) override {
        visit_function(node.arguments());
    }

protected:
    void visit_child(const detail::SymbolicNodePtr& child) override {
        if (!found_) RecursiveSymbolicVisitor::visit_child(child);
    }

private:
    void visit_function(const std::vector<detail::SymbolicNodePtr>& arguments) {
        for (const auto& argument : arguments) {
            if (contains_explicit_complex(argument)) {
                found_ = true;
                return;
            }
        }
        visit_children(arguments);
    }

    bool found_ = false;
};

bool has_z_dependent_function(const detail::SymbolicNodePtr& node,
                              const std::string& z) {
    if (!node) return false;
    FunctionDependencyFinder finder(z);
    node->accept(finder);
    return finder.found();
}

bool has_function_of_explicit_complex(const detail::SymbolicNodePtr& node) {
    if (!node) return false;
    FunctionOfExplicitComplexFinder finder;
    node->accept(finder);
    return finder.found();
}

} // namespace

static ExpressionResult calculate_residue_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&,
    const std::shared_ptr<SymbolicExpr>&,
    int,
    ComputationContext&);
static ExpressionResult cauchy_integral_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&,
    const std::shared_ptr<SymbolicExpr>&,
    int, ComputationContext&);
static bool is_analytic_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&);

ExpressionResult calculate_residue_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order,
    ComputationContext& context)
{
    const std::string operation = "calculate_residue";
    auto input = validate_complex_expr_point_input(f, z, z0, order, context, operation);
    if (!input) return ExpressionResult::failure(input.error());
    auto budget = context.consume_steps(static_cast<std::size_t>(order) * 12 + 12,
                                        operation);
    if (!budget) return ExpressionResult::failure(budget.error());

    try {
        auto calculated =
            calculate_residue_impl(f, z, z0, order, context);
        if (!calculated) return calculated;
        auto result = std::move(calculated.value());
        if (!result || !LMCAS::detail::node(result)) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "residue could not be constructed in the supported symbolic domain",
                operation);
        }
        return result;
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                          "allocation failed while calculating residue",
                                          operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                          ex.what(),
                                          operation);
    }
}

ExpressionResult calculate_residue_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order)
{
    ComputationContext context;
    return calculate_residue_checked(f, z, z0, order, context);
}

static ExpressionResult calculate_residue_impl(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order,
    ComputationContext& context) {
    if (order < 1) return SymbolicExpr::number(0);
    
    auto var_z = SymbolicExpr::variable(z);
    auto z_minus_z0 = SymbolicExpr::add(var_z, SymbolicExpr::multiply(z0, SymbolicExpr::number(-1)));
    auto pow_term = SymbolicExpr::power(z_minus_z0, SymbolicExpr::number(order));
    
    auto F = SymbolicExpr::multiply(pow_term, f);
    
    for (int i = 0; i < order - 1; i++) {
        F = F->differentiate(z);
    }
    
    auto fact = BigInt::factorial_checked(BigInt(order - 1), context);
    if (!fact) return ExpressionResult::failure(fact.error());

    F = SymbolicExpr::divide(F, SymbolicExpr::number(fact.value()));
    
    return limit_expression_checked(
        F, z, z0, LimitDirection::Both, context, Domain::Complex);
}

ExpressionResult cauchy_integral_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int n,
    ComputationContext& context)
{
    const std::string operation = "cauchy_integral";
    auto input = validate_complex_expr_point_input(f, z, z0, n, context, operation);
    if (!input) return ExpressionResult::failure(input.error());
    auto budget = context.consume_steps(static_cast<std::size_t>(n) * 12 + 12,
                                        operation);
    if (!budget) return ExpressionResult::failure(budget.error());

    try {
        auto calculated = cauchy_integral_impl(f, z, z0, n, context);
        if (!calculated) return calculated;
        auto result = std::move(calculated.value());
        if (!result || !LMCAS::detail::node(result)) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "Cauchy integral formula could not be constructed in the supported symbolic domain",
                operation);
        }
        return ExpressionResult::success(result);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                          "allocation failed while applying Cauchy integral formula",
                                          operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                          ex.what(),
                                          operation);
    }
}

ExpressionResult cauchy_integral_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int n)
{
    ComputationContext context;
    return cauchy_integral_checked(f, z, z0, n, context);
}

static ExpressionResult cauchy_integral_impl(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int n, ComputationContext& context) {
    
    /// n is the power in the denominator: \oint f(z)/(z-z0)^n dz
    if (n < 1) return SymbolicExpr::number(0);
    
    auto deriv = f;
    for (int i = 0; i < n - 1; i++) {
        deriv = deriv->differentiate(z);
    }
    
    auto f_n_minus_1_z0 = deriv->substitute(z, z0);
    
    auto fact = BigInt::factorial_checked(BigInt(n - 1), context);
    if (!fact) return ExpressionResult::failure(fact.error());

    auto term = SymbolicExpr::divide(f_n_minus_1_z0, SymbolicExpr::number(fact.value()));
    
    auto pi_node = LMCAS::detail::make_node<VariableNode>("pi");
    auto i_node = SymbolicFactory::create_complex(
        LMCAS::detail::node(SymbolicExpr::number(0)),
        LMCAS::detail::node(SymbolicExpr::number(1)));
    
    auto pi_expr = LMCAS::detail::make_expression_ptr(pi_node);
    auto i_expr = LMCAS::detail::make_expression_ptr(i_node);
    
    auto two_pi_i = SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::multiply(pi_expr, i_expr));
    
    return SymbolicExpr::multiply(two_pi_i, term)->simplify();
}

std::shared_ptr<SymbolicExpr> analytic_continuation(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string&) {
    /// Basic analytic continuation via symbolic simplification
    /// Simplification often reduces locally defined series (if represented)
    /// to their global analytic forms (e.g. geometric series).
    if (!f) return f;
    return f->simplify();
}

namespace {

std::optional<int> bounded_complex_power(
    const std::shared_ptr<const SymbolicNode>& node) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;
    BigInt exponent;
    if (const auto* integer = std::get_if<BigInt>(&number->value())) {
        exponent = *integer;
    } else if (const auto* rational = std::get_if<Rational>(&number->value())) {
        if (!rational->is_integer()) return std::nullopt;
        exponent = rational->get_numerator();
    } else {
        const double value = std::get<lmmc_real_t>(number->value());
        if (!std::isfinite(value) || value != std::floor(value) ||
            value < 0 || value > 16) return std::nullopt;
        return static_cast<int>(value);
    }
    if (exponent < BigInt(0) || exponent > BigInt(16)) return std::nullopt;
    return static_cast<int>(*exponent.try_to_int64());
}

/// 递归地将表达式分解为 (实部, 虚部)，把 ComplexNode 视为 a+bi。
/// 仅处理加法、乘法、数值与 ComplexNode 组合；其余子表达式视为实值。
void split_real_imag(const std::shared_ptr<const SymbolicNode>& node,
                     std::shared_ptr<SymbolicExpr>& re,
                     std::shared_ptr<SymbolicExpr>& im) {
    re = SymbolicExpr::number(0);
    im = SymbolicExpr::number(0);
    if (!node) return;
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node);
        variable && detail::is_imaginary_unit_name(variable->name())) {
        im = SymbolicExpr::number(1);
        return;
    }


    if (auto cn = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        re = LMCAS::detail::make_expression_ptr(cn->real());
        im = LMCAS::detail::make_expression_ptr(cn->imag());
        return;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            std::shared_ptr<SymbolicExpr> r2, i2;
            split_real_imag(op, r2, i2);
            re = SymbolicExpr::add(re, r2);
            im = SymbolicExpr::add(im, i2);
        }
        re = re->simplify();
        im = im->simplify();
        return;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        /// 累乘：(a+bi)(c+di) = (ac-bd) + (ad+bc)i
        std::shared_ptr<SymbolicExpr> accR = SymbolicExpr::number(1);
        std::shared_ptr<SymbolicExpr> accI = SymbolicExpr::number(0);
        for (const auto& op : mul->operands()) {
            std::shared_ptr<SymbolicExpr> r2, i2;
            split_real_imag(op, r2, i2);
            auto ac = SymbolicExpr::multiply(accR, r2);
            auto bd = SymbolicExpr::multiply(accI, i2);
            auto ad = SymbolicExpr::multiply(accR, i2);
            auto bc = SymbolicExpr::multiply(accI, r2);
            accR = SymbolicExpr::add(ac, SymbolicExpr::multiply(SymbolicExpr::number(-1), bd))->simplify();
            accI = SymbolicExpr::add(ad, bc)->simplify();
        }
        re = accR;
        im = accI;
        return;
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        /// 对整数次幂，展开为重复乘法以分离实/虚部。
        auto exponent = bounded_complex_power(pw->exponent());
        std::shared_ptr<SymbolicExpr> baseR, baseI;
        split_real_imag(pw->base(), baseR, baseI);
        bool base_real = LMCAS::detail::node(baseI) && LMCAS::detail::node(baseI)->is_zero();
        if (exponent && !base_real) {
            /// (a+bi)^e via repeated complex multiplication
            std::shared_ptr<SymbolicExpr> accR = SymbolicExpr::number(1);
            std::shared_ptr<SymbolicExpr> accI = SymbolicExpr::number(0);
            for (int k = 0; k < *exponent; ++k) {
                auto ac = SymbolicExpr::multiply(accR, baseR);
                auto bd = SymbolicExpr::multiply(accI, baseI);
                auto ad = SymbolicExpr::multiply(accR, baseI);
                auto bc = SymbolicExpr::multiply(accI, baseR);
                accR = SymbolicExpr::add(ac, SymbolicExpr::multiply(SymbolicExpr::number(-1), bd))->simplify();
                accI = SymbolicExpr::add(ad, bc)->simplify();
            }
            re = accR; im = accI;
            return;
        }
        /// 实底数或非整数指数：视为实值
        if (base_real) {
            re = LMCAS::detail::make_expression_ptr(node);
            im = SymbolicExpr::number(0);
            return;
        }
        /// 退化情形：原样返回为实部
        re = LMCAS::detail::make_expression_ptr(node);
        im = SymbolicExpr::number(0);
        return;
    }
    /// 默认：视为实值表达式
    re = LMCAS::detail::make_expression_ptr(node);
    im = SymbolicExpr::number(0);
}

} // namespace

ExpressionResult real_part_checked(const std::shared_ptr<SymbolicExpr>& expr,
                                    ComputationContext& context) {
    const std::string operation = "real_part";
    auto input = validate_complex_expr_input(expr, context, operation);
    if (!input) return ExpressionResult::failure(input.error());
    if (has_function_of_explicit_complex(LMCAS::detail::node(expr))) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "real part of functions with complex arguments is outside the current support domain",
            operation);
    }

    std::shared_ptr<SymbolicExpr> re, im;
    split_real_imag(LMCAS::detail::node(expr), re, im);
    auto simplified = re ? re->simplify() : nullptr;
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                          "real part construction failed",
                                          operation);
    }
    return ExpressionResult::success(simplified);
}

ExpressionResult real_part_checked(const std::shared_ptr<SymbolicExpr>& expr) {
    ComputationContext context;
    return real_part_checked(expr, context);
}


ExpressionResult imag_part_checked(const std::shared_ptr<SymbolicExpr>& expr,
                                    ComputationContext& context) {
    const std::string operation = "imag_part";
    auto input = validate_complex_expr_input(expr, context, operation);
    if (!input) return ExpressionResult::failure(input.error());
    if (has_function_of_explicit_complex(LMCAS::detail::node(expr))) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "imaginary part of functions with complex arguments is outside the current support domain",
            operation);
    }

    std::shared_ptr<SymbolicExpr> re, im;
    split_real_imag(LMCAS::detail::node(expr), re, im);
    auto simplified = im ? im->simplify() : nullptr;
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                          "imaginary part construction failed",
                                          operation);
    }
    return ExpressionResult::success(simplified);
}

ExpressionResult imag_part_checked(const std::shared_ptr<SymbolicExpr>& expr) {
    ComputationContext context;
    return imag_part_checked(expr, context);
}


ExpressionResult conjugate_checked(const std::shared_ptr<SymbolicExpr>& expr,
                                    ComputationContext& context) {
    const std::string operation = "conjugate";
    auto input = validate_complex_expr_input(expr, context, operation);
    if (!input) return ExpressionResult::failure(input.error());
    if (has_function_of_explicit_complex(LMCAS::detail::node(expr))) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "conjugate of functions with complex arguments is outside the current support domain",
            operation);
    }

    std::shared_ptr<SymbolicExpr> re, im;
    split_real_imag(LMCAS::detail::node(expr), re, im);
    /// conj(a+bi) = a - bi
    auto neg_im = SymbolicExpr::multiply(SymbolicExpr::number(-1), im)->simplify();
    if (!re || !LMCAS::detail::node(re) || !neg_im || !LMCAS::detail::node(neg_im)) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                          "conjugate construction failed",
                                          operation);
    }
    if (LMCAS::detail::node(neg_im) && LMCAS::detail::node(neg_im)->is_zero()) {
        auto simplified_re = re->simplify();
        if (!simplified_re || !LMCAS::detail::node(simplified_re)) {
            return ExpressionResult::failure(CasErrc::InternalInvariant,
                                              "conjugate construction failed",
                                              operation);
        }
        return ExpressionResult::success(simplified_re);
    }
    auto cn = SymbolicFactory::create_complex(
        LMCAS::detail::node(re->simplify()), LMCAS::detail::node(neg_im));
    if (!cn) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                          "conjugate construction failed",
                                          operation);
    }
    return ExpressionResult::success(LMCAS::detail::make_expression_ptr(cn));
}

ExpressionResult conjugate_checked(const std::shared_ptr<SymbolicExpr>& expr) {
    ComputationContext context;
    return conjugate_checked(expr, context);
}


ComplexBoolResult is_analytic_checked(const std::shared_ptr<SymbolicExpr>& f,
                                      const std::string& z,
                                      ComputationContext& context) {
    const std::string operation = "is_analytic";
    auto input = validate_complex_expr_input(f, context, operation);
    if (!input) return ComplexBoolResult::failure(input.error());
    if (z.empty()) {
        return ComplexBoolResult::failure(CasErrc::InvalidArgument,
                                          "complex variable name cannot be empty",
                                          operation);
    }
    try {
        if (has_z_dependent_function(LMCAS::detail::node(f), z)) {
            return ComplexBoolResult::failure(
                CasErrc::Inconclusive,
                "analyticity of functions depending on the complex variable is outside the current support domain",
                operation);
        }
        auto budget = context.consume_steps(24, operation);
        if (!budget) return ComplexBoolResult::failure(budget.error());
        return ComplexBoolResult::success(is_analytic_impl(f, z));
    } catch (const detail::UnsupportedDifferentiation& ex) {
        return ComplexBoolResult::failure(CasErrc::Inconclusive,
                                          ex.what(),
                                          operation);
    } catch (const std::bad_alloc&) {
        return ComplexBoolResult::failure(CasErrc::ResourceLimit,
                                          "allocation failed while checking analyticity",
                                          operation);
    } catch (const std::exception& ex) {
        return ComplexBoolResult::failure(CasErrc::InternalInvariant,
                                          ex.what(),
                                          operation);
    }
}

ComplexBoolResult is_analytic_checked(const std::shared_ptr<SymbolicExpr>& f,
                                      const std::string& z) {
    ComputationContext context;
    return is_analytic_checked(f, z, context);
}

static std::shared_ptr<SymbolicExpr> differentiate_if_dependent(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable) {
    if (!expression_depends_on_variable(detail::node(expression), variable)) {
        return SymbolicExpr::number(0);
    }
    return expression->differentiate(variable);
}

static bool is_analytic_impl(const std::shared_ptr<SymbolicExpr>& f, const std::string& z) {
    if (!f) return false;
    /// 将 z 替换为 (z_re + i·z_im)，分离 u、v，检验 Cauchy-Riemann。
    std::string xr = z + "_re";
    std::string xi = z + "_im";
    auto zr = SymbolicExpr::variable(xr);
    auto zi = SymbolicExpr::variable(xi);
    auto i_unit = LMCAS::detail::make_expression_ptr(
        SymbolicFactory::create_complex(
            LMCAS::detail::node(SymbolicExpr::number(0)),
            LMCAS::detail::node(SymbolicExpr::number(1))));
    auto z_sub = SymbolicExpr::add(zr, SymbolicExpr::multiply(i_unit, zi));

    auto fz = f->substitute(z, z_sub);
    if (!fz) return false;
    fz = fz->simplify();

    std::shared_ptr<SymbolicExpr> u, v;
    split_real_imag(LMCAS::detail::node(fz), u, v);

    auto ux = differentiate_if_dependent(u, xr);
    auto uy = differentiate_if_dependent(u, xi);
    auto vx = differentiate_if_dependent(v, xr);
    auto vy = differentiate_if_dependent(v, xi);

    /// CR1: ux - vy == 0 ; CR2: uy + vx == 0
    auto cr1 = SymbolicExpr::add(ux, SymbolicExpr::multiply(SymbolicExpr::number(-1), vy))->simplify();
    auto cr2 = SymbolicExpr::add(uy, vx)->simplify();

    return LMCAS::detail::node(cr1) && LMCAS::detail::node(cr1)->is_zero() && LMCAS::detail::node(cr2) && LMCAS::detail::node(cr2)->is_zero();
}

ExpressionResult residue_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order,
    ComputationContext& context) {
    return calculate_residue_checked(f, z, z0, order, context);
}

ExpressionResult residue_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order) {
    return calculate_residue_checked(f, z, z0, order);
}


} // namespace LMCAS
