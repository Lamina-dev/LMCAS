/**
 * @file complex_analysis.cpp
 * @brief 复变函数分析实现。
 */
#include "complex_analysis.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

namespace lamina {

std::shared_ptr<SymbolicExpr> calculate_residue(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order) {
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
    
    return F->limit(z, z0);
}

std::shared_ptr<SymbolicExpr> cauchy_integral(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int n) {
    
    // n is the power in the denominator: \oint f(z)/(z-z0)^n dz
    if (n < 1) return SymbolicExpr::number(0);
    
    auto deriv = f;
    for (int i = 0; i < n - 1; i++) {
        deriv = deriv->differentiate(z);
    }
    
    auto f_n_minus_1_z0 = deriv->substitute(z, z0);
    
    int fact = 1;
    for (int i = 2; i <= n - 1; i++) fact *= i;
    
    auto term = SymbolicExpr::divide(f_n_minus_1_z0, SymbolicExpr::number(fact));
    
    auto pi_node = std::make_shared<VariableNode>("pi");
    auto i_node = SymbolicFactory::create_complex(SymbolicExpr::number(0)->root, SymbolicExpr::number(1)->root);
    
    auto pi_expr = std::make_shared<SymbolicExpr>(pi_node);
    auto i_expr = std::make_shared<SymbolicExpr>(i_node);
    
    auto two_pi_i = SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::multiply(pi_expr, i_expr));
    
    return SymbolicExpr::multiply(two_pi_i, term)->simplify();
}

std::shared_ptr<SymbolicExpr> analytic_continuation(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z) {
    // Basic analytic continuation via symbolic simplification
    // Simplification often reduces locally defined series (if represented)
    // to their global analytic forms (e.g. geometric series).
    if (!f) return f;
    return f->simplify();
}

namespace {

// 递归地将表达式分解为 (实部, 虚部)，把 ComplexNode 视为 a+bi。
// 仅处理加法、乘法、数值与 ComplexNode 组合；其余子表达式视为实值。
void split_real_imag(const std::shared_ptr<SymbolicNode>& node,
                     std::shared_ptr<SymbolicExpr>& re,
                     std::shared_ptr<SymbolicExpr>& im) {
    re = SymbolicExpr::number(0);
    im = SymbolicExpr::number(0);
    if (!node) return;

    if (auto cn = std::dynamic_pointer_cast<ComplexNode>(node)) {
        re = std::make_shared<SymbolicExpr>(cn->real);
        im = std::make_shared<SymbolicExpr>(cn->imag);
        return;
    }
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (const auto& op : add->operands) {
            std::shared_ptr<SymbolicExpr> r2, i2;
            split_real_imag(op, r2, i2);
            re = SymbolicExpr::add(re, r2);
            im = SymbolicExpr::add(im, i2);
        }
        re = re->simplify();
        im = im->simplify();
        return;
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        // 累乘：(a+bi)(c+di) = (ac-bd) + (ad+bc)i
        std::shared_ptr<SymbolicExpr> accR = SymbolicExpr::number(1);
        std::shared_ptr<SymbolicExpr> accI = SymbolicExpr::number(0);
        for (const auto& op : mul->operands) {
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
    if (auto pw = std::dynamic_pointer_cast<PowerNode>(node)) {
        // 对整数次幂，展开为重复乘法以分离实/虚部。
        auto exp_num = std::dynamic_pointer_cast<NumberNode>(pw->exponent);
        long long e = 0;
        bool int_exp = false;
        if (exp_num) {
            if (std::holds_alternative<BigInt>(exp_num->value)) {
                e = (long long)std::get<BigInt>(exp_num->value).to_int();
                int_exp = true;
            } else if (std::holds_alternative<lmmc_real_t>(exp_num->value)) {
                double dv = std::get<lmmc_real_t>(exp_num->value);
                if (dv == (long long)dv) { e = (long long)dv; int_exp = true; }
            }
        }
        std::shared_ptr<SymbolicExpr> baseR, baseI;
        split_real_imag(pw->base, baseR, baseI);
        bool base_real = baseI->root && baseI->root->is_zero();
        if (int_exp && e >= 0 && e <= 16 && !base_real) {
            // (a+bi)^e via repeated complex multiplication
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
        // 实底数或非整数指数：视为实值
        if (base_real) {
            re = std::make_shared<SymbolicExpr>(node);
            im = SymbolicExpr::number(0);
            return;
        }
        // 退化情形：原样返回为实部
        re = std::make_shared<SymbolicExpr>(node);
        im = SymbolicExpr::number(0);
        return;
    }
    // 默认：视为实值表达式
    re = std::make_shared<SymbolicExpr>(node);
    im = SymbolicExpr::number(0);
}

} // anonymous namespace

std::shared_ptr<SymbolicExpr> real_part(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return nullptr;
    std::shared_ptr<SymbolicExpr> re, im;
    split_real_imag(expr->root, re, im);
    return re->simplify();
}

std::shared_ptr<SymbolicExpr> imag_part(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return nullptr;
    std::shared_ptr<SymbolicExpr> re, im;
    split_real_imag(expr->root, re, im);
    return im->simplify();
}

std::shared_ptr<SymbolicExpr> conjugate(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return nullptr;
    std::shared_ptr<SymbolicExpr> re, im;
    split_real_imag(expr->root, re, im);
    // conj(a+bi) = a - bi
    auto neg_im = SymbolicExpr::multiply(SymbolicExpr::number(-1), im)->simplify();
    if (neg_im->root && neg_im->root->is_zero()) {
        return re->simplify();
    }
    auto cn = SymbolicFactory::create_complex(re->simplify()->root, neg_im->root);
    return std::make_shared<SymbolicExpr>(cn);
}

bool is_analytic(const std::shared_ptr<SymbolicExpr>& f, const std::string& z) {
    if (!f) return false;
    // 将 z 替换为 (z_re + i·z_im)，分离 u、v，检验 Cauchy-Riemann。
    std::string xr = z + "_re";
    std::string xi = z + "_im";
    auto zr = SymbolicExpr::variable(xr);
    auto zi = SymbolicExpr::variable(xi);
    auto i_unit = std::make_shared<SymbolicExpr>(
        SymbolicFactory::create_complex(SymbolicExpr::number(0)->root, SymbolicExpr::number(1)->root));
    auto z_sub = SymbolicExpr::add(zr, SymbolicExpr::multiply(i_unit, zi));

    auto fz = f->substitute(z, z_sub);
    if (!fz) return false;
    fz = fz->simplify();

    std::shared_ptr<SymbolicExpr> u, v;
    split_real_imag(fz->root, u, v);

    auto ux = u->differentiate(xr);
    auto uy = u->differentiate(xi);
    auto vx = v->differentiate(xr);
    auto vy = v->differentiate(xi);

    // CR1: ux - vy == 0 ; CR2: uy + vx == 0
    auto cr1 = SymbolicExpr::add(ux, SymbolicExpr::multiply(SymbolicExpr::number(-1), vy))->simplify();
    auto cr2 = SymbolicExpr::add(uy, vx)->simplify();

    return cr1->root && cr1->root->is_zero() && cr2->root && cr2->root->is_zero();
}

std::shared_ptr<SymbolicExpr> residue(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& z,
    const std::shared_ptr<SymbolicExpr>& z0,
    int order) {
    return calculate_residue(f, z, z0, order);
}

} // namespace lamina
