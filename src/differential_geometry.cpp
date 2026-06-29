/**
 * @file differential_geometry.cpp
 * @brief 微分几何实现。
 */
#include "differential_geometry.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

namespace lamina {

std::shared_ptr<SymbolicExpr> metric_inverse(
    const std::shared_ptr<SymbolicExpr>& g_ij) {
    return SymbolicExpr::inverse(g_ij);
}

std::shared_ptr<SymbolicExpr> christoffel_first_kind(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j) {
    
    auto mat = std::dynamic_pointer_cast<MatrixNode>(g_ij->root);
    if (!mat) return SymbolicExpr::number(0);
    
    auto g_jk = std::make_shared<SymbolicExpr>(mat->get(j, k));
    auto g_ik = std::make_shared<SymbolicExpr>(mat->get(i, k));
    auto g_ij_comp = std::make_shared<SymbolicExpr>(mat->get(i, j));
    
    auto term1 = g_jk->differentiate(coords[i]);
    auto term2 = g_ik->differentiate(coords[j]);
    auto term3 = g_ij_comp->differentiate(coords[k]);
    
    auto sum12 = SymbolicExpr::add(term1, term2);
    auto diff = SymbolicExpr::add(sum12, SymbolicExpr::multiply(term3, SymbolicExpr::number(-1)));
    
    return SymbolicExpr::multiply(SymbolicExpr::number(0.5), diff)->simplify();
}

std::shared_ptr<SymbolicExpr> christoffel_second_kind(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::shared_ptr<SymbolicExpr>& g_up_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j) {
    
    auto mat_up = std::dynamic_pointer_cast<MatrixNode>(g_up_ij->root);
    if (!mat_up) return SymbolicExpr::number(0);
    
    size_t dim = coords.size();
    auto result = SymbolicExpr::number(0);
    
    for (size_t m = 0; m < dim; m++) {
        auto g_up_km = std::make_shared<SymbolicExpr>(mat_up->get(k, m));
        auto gamma_first = christoffel_first_kind(g_ij, coords, m, i, j);
        auto term = SymbolicExpr::multiply(g_up_km, gamma_first);
        result = SymbolicExpr::add(result, term);
    }
    
    return result->simplify();
}

std::shared_ptr<SymbolicExpr> riemann_curvature_tensor(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int rho, int sigma, int mu, int nu) {
    
    auto g_up_ij = metric_inverse(g_ij);
    if (!g_up_ij) return SymbolicExpr::number(0);
    
    auto gamma_rho_nu_sigma = christoffel_second_kind(g_ij, g_up_ij, coords, rho, nu, sigma);
    auto gamma_rho_mu_sigma = christoffel_second_kind(g_ij, g_up_ij, coords, rho, mu, sigma);
    
    auto term1 = gamma_rho_nu_sigma->differentiate(coords[mu]);
    auto term2 = gamma_rho_mu_sigma->differentiate(coords[nu]);
    
    auto diff = SymbolicExpr::add(term1, SymbolicExpr::multiply(term2, SymbolicExpr::number(-1)));
    
    size_t dim = coords.size();
    auto sum_term = SymbolicExpr::number(0);
    
    for (size_t lambda = 0; lambda < dim; lambda++) {
        auto gamma_rho_mu_lambda = christoffel_second_kind(g_ij, g_up_ij, coords, rho, mu, lambda);
        auto gamma_lambda_nu_sigma = christoffel_second_kind(g_ij, g_up_ij, coords, lambda, nu, sigma);
        
        auto gamma_rho_nu_lambda = christoffel_second_kind(g_ij, g_up_ij, coords, rho, nu, lambda);
        auto gamma_lambda_mu_sigma = christoffel_second_kind(g_ij, g_up_ij, coords, lambda, mu, sigma);
        
        auto prod1 = SymbolicExpr::multiply(gamma_rho_mu_lambda, gamma_lambda_nu_sigma);
        auto prod2 = SymbolicExpr::multiply(gamma_rho_nu_lambda, gamma_lambda_mu_sigma);
        
        auto sub = SymbolicExpr::add(prod1, SymbolicExpr::multiply(prod2, SymbolicExpr::number(-1)));
        sum_term = SymbolicExpr::add(sum_term, sub);
    }
    
    return SymbolicExpr::add(diff, sum_term)->simplify();
}

std::shared_ptr<SymbolicExpr> lie_derivative(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& X,
    const std::vector<std::string>& vars,
    int order) {
    if (!f || X.size() != vars.size()) return SymbolicExpr::number(0);
    auto result = f;
    for (int it = 0; it < order; ++it) {
        auto acc = SymbolicExpr::number(0);
        for (size_t i = 0; i < vars.size(); ++i) {
            auto partial = result->differentiate(vars[i]);
            acc = SymbolicExpr::add(acc, SymbolicExpr::multiply(X[i], partial));
        }
        result = acc->simplify();
    }
    return result;
}

std::vector<std::shared_ptr<SymbolicExpr>> exterior_derivative(
    const std::vector<std::shared_ptr<SymbolicExpr>>& form_coeffs,
    int degree,
    const std::vector<std::string>& vars) {
    std::vector<std::shared_ptr<SymbolicExpr>> result;
    size_t n = vars.size();

    if (degree == 0) {
        /// d(f) = ∑ (∂f/∂xᵢ) dxᵢ  → 返回梯度分量
        if (form_coeffs.empty()) return result;
        auto f = form_coeffs[0];
        for (size_t i = 0; i < n; ++i) {
            result.push_back(f->differentiate(vars[i])->simplify());
        }
        return result;
    }

    if (degree == 1) {
        /// d(∑ωᵢ dxᵢ) = ∑_{i<j} (∂ωⱼ/∂xᵢ - ∂ωᵢ/∂xⱼ) dxᵢ∧dxⱼ
        if (form_coeffs.size() != n) return result;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                auto dwj_dxi = form_coeffs[j]->differentiate(vars[i]);
                auto dwi_dxj = form_coeffs[i]->differentiate(vars[j]);
                auto comp = SymbolicExpr::add(dwj_dxi,
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), dwi_dxj))->simplify();
                result.push_back(comp);
            }
        }
        return result;
    }

    return result;
}

} // namespace lamina
