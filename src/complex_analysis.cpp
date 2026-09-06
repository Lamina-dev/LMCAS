/**
 * @file complex_analysis.cpp
 * @brief 复变函数分析实现。
 */
#include "complex_analysis.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

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

bool has_z_dependent_function(const std::shared_ptr<const SymbolicNode>& node,
                              const std::string& z)
{
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return expression_depends_on_variable(node, z);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (has_z_dependent_function(op, z)) return true;
        }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (has_z_dependent_function(op, z)) return true;
        }
        return false;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return has_z_dependent_function(pow->base(), z) ||
               has_z_dependent_function(pow->exponent(), z);
    }
    if (auto matrix = std::dynamic_pointer_cast<const MatrixNode>(node)) {
        if (std::holds_alternative<MatrixNode::DenseStorage>(matrix->storage())) {
            for (const auto& item : std::get<MatrixNode::DenseStorage>(matrix->storage())) {
                if (has_z_dependent_function(item, z)) return true;
            }
        } else {
            for (const auto& [idx, item] :
                 std::get<MatrixNode::SparseStorage>(matrix->storage())) {
                (void)idx;
                if (has_z_dependent_function(item, z)) return true;
            }
        }
        return false;
    }
    if (auto rel = std::dynamic_pointer_cast<const RelationalNode>(node)) {
        return has_z_dependent_function(rel->left(), z) ||
               has_z_dependent_function(rel->right(), z);
    }
    if (auto logic = std::dynamic_pointer_cast<const LogicalNode>(node)) {
        return has_z_dependent_function(logic->left(), z) ||
               has_z_dependent_function(logic->right(), z);
    }
    if (auto piecewise = std::dynamic_pointer_cast<const PiecewiseNode>(node)) {
        for (const auto& branch : piecewise->branches()) {
            if (has_z_dependent_function(branch.expression, z) ||
                has_z_dependent_function(branch.condition, z)) {
                return true;
            }
        }
        return has_z_dependent_function(piecewise->default_expr(), z);
    }
    if (auto sum = std::dynamic_pointer_cast<const SummationNode>(node)) {
        return has_z_dependent_function(sum->lower_bound(), z) ||
               has_z_dependent_function(sum->upper_bound(), z) ||
               (sum->index_var() != z && has_z_dependent_function(sum->body(), z));
    }
    if (auto product = std::dynamic_pointer_cast<const ProductNode>(node)) {
        return has_z_dependent_function(product->lower_bound(), z) ||
               has_z_dependent_function(product->upper_bound(), z) ||
               (product->index_var() != z && has_z_dependent_function(product->body(), z));
    }
    if (auto transform = std::dynamic_pointer_cast<const TransformNode>(node)) {
        return transform->source_var() != z &&
               has_z_dependent_function(transform->body(), z);
    }
    if (auto quantifier = std::dynamic_pointer_cast<const QuantifierNode>(node)) {
        return has_z_dependent_function(quantifier->domain(), z) ||
               (quantifier->bound_var() != z &&
                has_z_dependent_function(quantifier->predicate(), z));
    }
    if (auto set_builder = std::dynamic_pointer_cast<const SetBuilderNode>(node)) {
        return has_z_dependent_function(set_builder->domain(), z) ||
               (set_builder->element_var() != z &&
                has_z_dependent_function(set_builder->predicate(), z));
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return has_z_dependent_function(complex->real(), z) ||
               has_z_dependent_function(complex->imag(), z);
    }
    return false;
}

bool contains_explicit_complex(const std::shared_ptr<const SymbolicNode>& node)
{
    if (!node) return false;
    if (std::dynamic_pointer_cast<const ComplexNode>(node)) return true;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        for (const auto& arg : fn->arguments()) {
            if (contains_explicit_complex(arg)) return true;
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (contains_explicit_complex(op)) return true;
        }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (contains_explicit_complex(op)) return true;
        }
        return false;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return contains_explicit_complex(pow->base()) ||
               contains_explicit_complex(pow->exponent());
    }
    return false;
}

bool has_function_of_explicit_complex(const std::shared_ptr<const SymbolicNode>& node)
{
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        for (const auto& arg : fn->arguments()) {
            if (contains_explicit_complex(arg) ||
                has_function_of_explicit_complex(arg)) {
                return true;
            }
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (has_function_of_explicit_complex(op)) return true;
        }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (has_function_of_explicit_complex(op)) return true;
        }
        return false;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return has_function_of_explicit_complex(pow->base()) ||
               has_function_of_explicit_complex(pow->exponent());
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return has_function_of_explicit_complex(complex->real()) ||
               has_function_of_explicit_complex(complex->imag());
    }
    return false;
}

} // namespace

static ExpressionResult calculate_residue_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&,
    const std::shared_ptr<SymbolicExpr>&,
    int,
    ComputationContext&);
static std::shared_ptr<SymbolicExpr> cauchy_integral_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&,
    const std::shared_ptr<SymbolicExpr>&,
    int);
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
    
    int fact = 1;
    for (int i = 2; i <= order - 1; i++) fact *= i;
    
    F = SymbolicExpr::divide(F, SymbolicExpr::number(fact));
    
    return limit_expression_checked(
        F, z, z0, LimitDirection::Both, context);
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
        auto result = cauchy_integral_impl(f, z, z0, n);
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

static std::shared_ptr<SymbolicExpr> cauchy_integral_impl(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int n) {
    
    /// n is the power in the denominator: \oint f(z)/(z-z0)^n dz
    if (n < 1) return SymbolicExpr::number(0);
    
    auto deriv = f;
    for (int i = 0; i < n - 1; i++) {
        deriv = deriv->differentiate(z);
    }
    
    auto f_n_minus_1_z0 = deriv->substitute(z, z0);
    
    int fact = 1;
    for (int i = 2; i <= n - 1; i++) fact *= i;
    
    auto term = SymbolicExpr::divide(f_n_minus_1_z0, SymbolicExpr::number(fact));
    
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

/// 递归地将表达式分解为 (实部, 虚部)，把 ComplexNode 视为 a+bi。
/// 仅处理加法、乘法、数值与 ComplexNode 组合；其余子表达式视为实值。
void split_real_imag(const std::shared_ptr<const SymbolicNode>& node,
                     std::shared_ptr<SymbolicExpr>& re,
                     std::shared_ptr<SymbolicExpr>& im) {
    re = SymbolicExpr::number(0);
    im = SymbolicExpr::number(0);
    if (!node) return;

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
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
        long long e = 0;
        bool int_exp = false;
        if (exp_num) {
            if (std::holds_alternative<BigInt>(exp_num->value())) {
                e = (long long)std::get<BigInt>(exp_num->value()).to_int();
                int_exp = true;
            } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
                double dv = std::get<lmmc_real_t>(exp_num->value());
                if (dv == (long long)dv) { e = (long long)dv; int_exp = true; }
            }
        }
        std::shared_ptr<SymbolicExpr> baseR, baseI;
        split_real_imag(pw->base(), baseR, baseI);
        bool base_real = LMCAS::detail::node(baseI) && LMCAS::detail::node(baseI)->is_zero();
        if (int_exp && e >= 0 && e <= 16 && !base_real) {
            /// (a+bi)^e via repeated complex multiplication
            std::shared_ptr<SymbolicExpr> accR = SymbolicExpr::number(1);
            std::shared_ptr<SymbolicExpr> accI = SymbolicExpr::number(0);
            for (long long k = 0; k < e; ++k) {
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
    if (has_z_dependent_function(LMCAS::detail::node(f), z)) {
        return ComplexBoolResult::failure(
            CasErrc::Inconclusive,
            "analyticity of functions depending on the complex variable is outside the current support domain",
            operation);
    }
    auto budget = context.consume_steps(24, operation);
    if (!budget) return ComplexBoolResult::failure(budget.error());

    try {
        return ComplexBoolResult::success(is_analytic_impl(f, z));
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

    auto ux = u->differentiate(xr);
    auto uy = u->differentiate(xi);
    auto vx = v->differentiate(xr);
    auto vy = v->differentiate(xi);

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
