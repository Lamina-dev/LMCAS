/**
 * @file differential_geometry.cpp
 * @brief 微分几何实现。
 */
#include "differential_geometry.hpp"
#include "symbolic_matrix.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

#include <algorithm>
#include <map>
#include <numeric>

namespace lamina {

namespace {

Result<std::shared_ptr<const MatrixNode>> validate_metric_matrix(
    const std::shared_ptr<SymbolicExpr>& metric,
    const std::vector<std::string>* coords,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return Result<std::shared_ptr<const MatrixNode>>::failure(step.error());
    if (!metric || !lamina::detail::node(metric)) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument,
            "metric tensor cannot be null",
            operation);
    }
    auto matrix = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(metric));
    if (!matrix) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument,
            "metric tensor must be a matrix",
            operation);
    }
    if (matrix->rows() != matrix->cols()) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument,
            "metric tensor must be square",
            operation);
    }
    if (coords) {
        if (coords->empty()) {
            return Result<std::shared_ptr<const MatrixNode>>::failure(
                CasErrc::InvalidArgument,
                "coordinate list cannot be empty",
                operation);
        }
        if (coords->size() != matrix->rows()) {
            return Result<std::shared_ptr<const MatrixNode>>::failure(
                CasErrc::InvalidArgument,
                "coordinate count must match metric dimension",
                operation);
        }
        for (const auto& coord : *coords) {
            if (coord.empty()) {
                return Result<std::shared_ptr<const MatrixNode>>::failure(
                    CasErrc::InvalidArgument,
                    "coordinate names cannot be empty",
                    operation);
            }
        }
    }
    return Result<std::shared_ptr<const MatrixNode>>::success(matrix);
}

Result<void> validate_index(int idx,
                            std::size_t dim,
                            const std::string& name,
                            const std::string& operation)
{
    if (idx < 0 || static_cast<std::size_t>(idx) >= dim) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     name + " index is out of bounds",
                                     operation);
    }
    return Result<void>::success();
}

Result<void> validate_expr_vector(const std::vector<std::shared_ptr<SymbolicExpr>>& values,
                                  const std::vector<std::string>& vars,
                                  const std::string& value_name,
                                  const std::string& operation)
{
    if (vars.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "variable list cannot be empty",
                                     operation);
    }
    if (values.size() != vars.size()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     value_name + " dimension must match variable count",
                                     operation);
    }
    for (const auto& var : vars) {
        if (var.empty()) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "variable names cannot be empty",
                                         operation);
        }
    }
    for (const auto& value : values) {
        if (!value || !lamina::detail::node(value)) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         value_name + " components cannot be null",
                                         operation);
        }
    }
    return Result<void>::success();
}

Result<std::size_t> bounded_binomial(std::size_t n,
                                     std::size_t k,
                                     std::size_t limit,
                                     const std::string& operation)
{
    if (k > n) return Result<std::size_t>::success(0);
    k = std::min(k, n - k);
    std::size_t value = 1;
    for (std::size_t i = 1; i <= k; ++i) {
        std::size_t numerator = n - k + i;
        std::size_t denominator = i;
        const auto fraction_gcd = std::gcd(numerator, denominator);
        numerator /= fraction_gcd;
        denominator /= fraction_gcd;
        const auto value_gcd = std::gcd(value, denominator);
        value /= value_gcd;
        denominator /= value_gcd;
        if (denominator != 1) {
            return Result<std::size_t>::failure(
                CasErrc::InternalInvariant,
                "binomial coefficient reduction was not exact",
                operation);
        }
        if (numerator != 0 && value > limit / numerator) {
            return Result<std::size_t>::failure(
                CasErrc::ResourceLimit,
                "differential-form component count exceeds the expansion budget",
                operation);
        }
        value *= numerator;
    }
    if (value > limit) {
        return Result<std::size_t>::failure(
            CasErrc::ResourceLimit,
            "differential-form component count exceeds the expansion budget",
            operation);
    }
    return Result<std::size_t>::success(value);
}

std::vector<std::vector<std::size_t>> index_combinations(std::size_t n,
                                                          std::size_t k,
                                                          std::size_t count)
{
    std::vector<std::vector<std::size_t>> combinations;
    combinations.reserve(count);
    if (k > n) return combinations;
    if (k == 0) {
        combinations.emplace_back();
        return combinations;
    }

    std::vector<std::size_t> current(k);
    std::iota(current.begin(), current.end(), 0);
    while (true) {
        combinations.push_back(current);
        std::size_t position = k;
        while (position > 0) {
            --position;
            if (current[position] < n - k + position) break;
        }
        if (position == 0 && current[0] == n - k) break;
        ++current[position];
        for (std::size_t i = position + 1; i < k; ++i) {
            current[i] = current[i - 1] + 1;
        }
    }
    return combinations;
}

} // namespace

static std::shared_ptr<SymbolicExpr> christoffel_first_kind_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::vector<std::string>&, int, int, int);
static std::shared_ptr<SymbolicExpr> christoffel_second_kind_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::vector<std::string>&, int, int, int);
static std::shared_ptr<SymbolicExpr> lie_derivative_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::vector<std::shared_ptr<SymbolicExpr>>&,
    const std::vector<std::string>&, int);

DifferentialGeometryExprResult metric_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    ComputationContext& context)
{
    const std::string operation = "metric_inverse";
    auto metric = validate_metric_matrix(g_ij, nullptr, context, operation);
    if (!metric) return DifferentialGeometryExprResult::failure(metric.error());
    auto budget = context.consume_steps(12, operation);
    if (!budget) return DifferentialGeometryExprResult::failure(budget.error());
    try {
        auto result = matrix_inverse_checked(g_ij, context);
        if (!result) {
            return DifferentialGeometryExprResult::failure(result.error());
        }
        return DifferentialGeometryExprResult::success(result.value());
    } catch (const std::bad_alloc&) {
        return DifferentialGeometryExprResult::failure(CasErrc::ResourceLimit,
                                                       "allocation failed while inverting metric",
                                                       operation);
    } catch (const std::exception& ex) {
        return DifferentialGeometryExprResult::failure(CasErrc::InternalInvariant,
                                                       ex.what(),
                                                       operation);
    }
}

DifferentialGeometryExprResult metric_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij)
{
    ComputationContext context;
    return metric_inverse_checked(g_ij, context);
}


DifferentialGeometryExprResult christoffel_first_kind_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j,
    ComputationContext& context)
{
    const std::string operation = "christoffel_first_kind";
    auto metric = validate_metric_matrix(g_ij, &coords, context, operation);
    if (!metric) return DifferentialGeometryExprResult::failure(metric.error());
    const std::size_t dim = metric.value()->rows();
    auto k_check = validate_index(k, dim, "k", operation);
    if (!k_check) return DifferentialGeometryExprResult::failure(k_check.error());
    auto i_check = validate_index(i, dim, "i", operation);
    if (!i_check) return DifferentialGeometryExprResult::failure(i_check.error());
    auto j_check = validate_index(j, dim, "j", operation);
    if (!j_check) return DifferentialGeometryExprResult::failure(j_check.error());
    auto budget = context.consume_steps(16, operation);
    if (!budget) return DifferentialGeometryExprResult::failure(budget.error());
    try {
        auto result = christoffel_first_kind_impl(g_ij, coords, k, i, j);
        if (!result || !lamina::detail::node(result)) {
            return DifferentialGeometryExprResult::failure(
                CasErrc::Inconclusive,
                "Christoffel symbol could not be constructed",
                operation);
        }
        return DifferentialGeometryExprResult::success(result);
    } catch (const std::bad_alloc&) {
        return DifferentialGeometryExprResult::failure(CasErrc::ResourceLimit,
                                                       "allocation failed while calculating Christoffel symbol",
                                                       operation);
    } catch (const std::exception& ex) {
        return DifferentialGeometryExprResult::failure(CasErrc::InternalInvariant,
                                                       ex.what(),
                                                       operation);
    }
}

DifferentialGeometryExprResult christoffel_first_kind_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j)
{
    ComputationContext context;
    return christoffel_first_kind_checked(g_ij, coords, k, i, j, context);
}

static std::shared_ptr<SymbolicExpr> christoffel_first_kind_impl(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j) {
    
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(g_ij));
    if (!mat) return SymbolicExpr::number(0);
    
    auto g_jk = lamina::detail::make_expression_ptr(mat->get(j, k));
    auto g_ik = lamina::detail::make_expression_ptr(mat->get(i, k));
    auto g_ij_comp = lamina::detail::make_expression_ptr(mat->get(i, j));
    
    auto term1 = g_jk->differentiate(coords[i]);
    auto term2 = g_ik->differentiate(coords[j]);
    auto term3 = g_ij_comp->differentiate(coords[k]);
    
    auto sum12 = SymbolicExpr::add(term1, term2);
    auto diff = SymbolicExpr::add(sum12, SymbolicExpr::multiply(term3, SymbolicExpr::number(-1)));
    
    return SymbolicExpr::multiply(
        SymbolicExpr::number(Rational(BigInt(1), BigInt(2))),
        diff)->simplify();
}

DifferentialGeometryExprResult christoffel_second_kind_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::shared_ptr<SymbolicExpr>& g_up_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j,
    ComputationContext& context)
{
    const std::string operation = "christoffel_second_kind";
    auto metric = validate_metric_matrix(g_ij, &coords, context, operation);
    if (!metric) return DifferentialGeometryExprResult::failure(metric.error());
    auto inverse_metric = validate_metric_matrix(g_up_ij, &coords, context, operation);
    if (!inverse_metric) return DifferentialGeometryExprResult::failure(inverse_metric.error());
    const std::size_t dim = metric.value()->rows();
    auto k_check = validate_index(k, dim, "k", operation);
    if (!k_check) return DifferentialGeometryExprResult::failure(k_check.error());
    auto i_check = validate_index(i, dim, "i", operation);
    if (!i_check) return DifferentialGeometryExprResult::failure(i_check.error());
    auto j_check = validate_index(j, dim, "j", operation);
    if (!j_check) return DifferentialGeometryExprResult::failure(j_check.error());
    auto budget = context.consume_steps(dim * 20 + 12, operation);
    if (!budget) return DifferentialGeometryExprResult::failure(budget.error());
    try {
        auto result = christoffel_second_kind_impl(g_ij, g_up_ij, coords, k, i, j);
        if (!result || !lamina::detail::node(result)) {
            return DifferentialGeometryExprResult::failure(
                CasErrc::Inconclusive,
                "Christoffel symbol could not be constructed",
                operation);
        }
        return DifferentialGeometryExprResult::success(result);
    } catch (const std::bad_alloc&) {
        return DifferentialGeometryExprResult::failure(CasErrc::ResourceLimit,
                                                       "allocation failed while calculating Christoffel symbol",
                                                       operation);
    } catch (const std::exception& ex) {
        return DifferentialGeometryExprResult::failure(CasErrc::InternalInvariant,
                                                       ex.what(),
                                                       operation);
    }
}

DifferentialGeometryExprResult christoffel_second_kind_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::shared_ptr<SymbolicExpr>& g_up_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j)
{
    ComputationContext context;
    return christoffel_second_kind_checked(g_ij, g_up_ij, coords, k, i, j, context);
}

static std::shared_ptr<SymbolicExpr> christoffel_second_kind_impl(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::shared_ptr<SymbolicExpr>& g_up_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j) {
    
    auto mat_up = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(g_up_ij));
    if (!mat_up) return SymbolicExpr::number(0);
    
    size_t dim = coords.size();
    auto result = SymbolicExpr::number(0);
    
    for (size_t m = 0; m < dim; m++) {
        auto g_up_km = lamina::detail::make_expression_ptr(mat_up->get(k, m));
        auto gamma_first = christoffel_first_kind_impl(g_ij, coords, m, i, j);
        auto term = SymbolicExpr::multiply(g_up_km, gamma_first);
        result = SymbolicExpr::add(result, term);
    }
    
    return result->simplify();
}

static std::shared_ptr<SymbolicExpr> riemann_curvature_tensor_with_inverse(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::shared_ptr<SymbolicExpr>& g_up_ij,
    const std::vector<std::string>& coords,
    int rho, int sigma, int mu, int nu) {
    auto gamma_rho_nu_sigma =
        christoffel_second_kind_impl(g_ij, g_up_ij, coords, rho, nu, sigma);
    auto gamma_rho_mu_sigma =
        christoffel_second_kind_impl(g_ij, g_up_ij, coords, rho, mu, sigma);

    auto term1 = gamma_rho_nu_sigma->differentiate(coords[mu]);
    auto term2 = gamma_rho_mu_sigma->differentiate(coords[nu]);
    auto result = SymbolicExpr::add(
        term1,
        SymbolicExpr::multiply(term2, SymbolicExpr::number(-1)));

    for (std::size_t lambda = 0; lambda < coords.size(); ++lambda) {
        auto gamma_rho_mu_lambda =
            christoffel_second_kind_impl(g_ij, g_up_ij, coords, rho, mu, lambda);
        auto gamma_lambda_nu_sigma =
            christoffel_second_kind_impl(g_ij, g_up_ij, coords, lambda, nu, sigma);
        auto gamma_rho_nu_lambda =
            christoffel_second_kind_impl(g_ij, g_up_ij, coords, rho, nu, lambda);
        auto gamma_lambda_mu_sigma =
            christoffel_second_kind_impl(g_ij, g_up_ij, coords, lambda, mu, sigma);
        auto product_difference = SymbolicExpr::add(
            SymbolicExpr::multiply(
                gamma_rho_mu_lambda, gamma_lambda_nu_sigma),
            SymbolicExpr::multiply(
                SymbolicExpr::multiply(
                    gamma_rho_nu_lambda, gamma_lambda_mu_sigma),
                SymbolicExpr::number(-1)));
        result = SymbolicExpr::add(result, product_difference);
    }
    return result->simplify();
}

DifferentialGeometryExprResult riemann_curvature_tensor_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int rho, int sigma, int mu, int nu,
    ComputationContext& context)
{
    const std::string operation = "riemann_curvature_tensor";
    auto metric = validate_metric_matrix(g_ij, &coords, context, operation);
    if (!metric) return DifferentialGeometryExprResult::failure(metric.error());
    const std::size_t dim = metric.value()->rows();
    auto rho_check = validate_index(rho, dim, "rho", operation);
    if (!rho_check) return DifferentialGeometryExprResult::failure(rho_check.error());
    auto sigma_check = validate_index(sigma, dim, "sigma", operation);
    if (!sigma_check) return DifferentialGeometryExprResult::failure(sigma_check.error());
    auto mu_check = validate_index(mu, dim, "mu", operation);
    if (!mu_check) return DifferentialGeometryExprResult::failure(mu_check.error());
    auto nu_check = validate_index(nu, dim, "nu", operation);
    if (!nu_check) return DifferentialGeometryExprResult::failure(nu_check.error());
    auto budget = context.consume_steps(dim * 80 + 40, operation);
    if (!budget) return DifferentialGeometryExprResult::failure(budget.error());
    try {
        auto inverse_metric = metric_inverse_checked(g_ij, context);
        if (!inverse_metric) {
            return DifferentialGeometryExprResult::failure(inverse_metric.error());
        }
        auto result = riemann_curvature_tensor_with_inverse(
            g_ij, inverse_metric.value(), coords, rho, sigma, mu, nu);
        if (!result || !lamina::detail::node(result)) {
            return DifferentialGeometryExprResult::failure(
                CasErrc::Inconclusive,
                "Riemann curvature component could not be constructed",
                operation);
        }
        return DifferentialGeometryExprResult::success(result);
    } catch (const std::bad_alloc&) {
        return DifferentialGeometryExprResult::failure(CasErrc::ResourceLimit,
                                                       "allocation failed while calculating Riemann curvature",
                                                       operation);
    } catch (const std::exception& ex) {
        return DifferentialGeometryExprResult::failure(CasErrc::InternalInvariant,
                                                       ex.what(),
                                                       operation);
    }
}

DifferentialGeometryExprResult riemann_curvature_tensor_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int rho, int sigma, int mu, int nu)
{
    ComputationContext context;
    return riemann_curvature_tensor_checked(g_ij, coords, rho, sigma, mu, nu, context);
}


DifferentialGeometryExprResult lie_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& X,
    const std::vector<std::string>& vars,
    int order,
    ComputationContext& context)
{
    const std::string operation = "lie_derivative";
    auto step = context.consume_steps(1, operation);
    if (!step) return DifferentialGeometryExprResult::failure(step.error());
    if (!f || !lamina::detail::node(f)) {
        return DifferentialGeometryExprResult::failure(CasErrc::InvalidArgument,
                                                       "scalar function cannot be null",
                                                       operation);
    }
    if (order < 0 || order > 64) {
        return DifferentialGeometryExprResult::failure(CasErrc::InvalidArgument,
                                                       "Lie derivative order must be between 0 and 64",
                                                       operation);
    }
    auto vector_check = validate_expr_vector(X, vars, "vector field", operation);
    if (!vector_check) return DifferentialGeometryExprResult::failure(vector_check.error());
    auto budget = context.consume_steps(static_cast<std::size_t>(order + 1) * vars.size() * 8 + 8,
                                        operation);
    if (!budget) return DifferentialGeometryExprResult::failure(budget.error());
    try {
        auto result = lie_derivative_impl(f, X, vars, order);
        if (!result || !lamina::detail::node(result)) {
            return DifferentialGeometryExprResult::failure(
                CasErrc::Inconclusive,
                "Lie derivative could not be constructed",
                operation);
        }
        return DifferentialGeometryExprResult::success(result);
    } catch (const std::bad_alloc&) {
        return DifferentialGeometryExprResult::failure(CasErrc::ResourceLimit,
                                                       "allocation failed while calculating Lie derivative",
                                                       operation);
    } catch (const std::exception& ex) {
        return DifferentialGeometryExprResult::failure(CasErrc::InternalInvariant,
                                                       ex.what(),
                                                       operation);
    }
}

DifferentialGeometryExprResult lie_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& X,
    const std::vector<std::string>& vars,
    int order)
{
    ComputationContext context;
    return lie_derivative_checked(f, X, vars, order, context);
}

static std::shared_ptr<SymbolicExpr> lie_derivative_impl(
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

DifferentialGeometryVectorResult exterior_derivative_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& form_coeffs,
    int degree,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "exterior_derivative";
    auto step = context.consume_steps(1, operation);
    if (!step) return DifferentialGeometryVectorResult::failure(step.error());
    if (vars.empty()) {
        return DifferentialGeometryVectorResult::failure(CasErrc::InvalidArgument,
                                                         "variable list cannot be empty",
                                                         operation);
    }
    for (std::size_t i = 0; i < vars.size(); ++i) {
        if (vars[i].empty()) {
            return DifferentialGeometryVectorResult::failure(CasErrc::InvalidArgument,
                                                             "variable names cannot be empty",
                                                             operation);
        }
        for (std::size_t j = 0; j < i; ++j) {
            if (vars[i] == vars[j]) {
                return DifferentialGeometryVectorResult::failure(
                    CasErrc::InvalidArgument,
                    "coordinate variables must be unique",
                    operation);
            }
        }
    }
    if (degree < 0 || static_cast<std::size_t>(degree) > vars.size()) {
        return DifferentialGeometryVectorResult::failure(CasErrc::InvalidArgument,
                                                         "form degree must be between zero and the coordinate dimension",
                                                         operation);
    }
    const auto input_count = bounded_binomial(vars.size(),
                                              static_cast<std::size_t>(degree),
                                              context.limits().max_expansion_terms,
                                              operation);
    if (!input_count) return DifferentialGeometryVectorResult::failure(input_count.error());
    if (form_coeffs.size() != input_count.value()) {
        return DifferentialGeometryVectorResult::failure(CasErrc::InvalidArgument,
                                                         "form coefficient count does not match degree and dimension",
                                                         operation);
    }
    for (const auto& coeff : form_coeffs) {
        if (!coeff || !lamina::detail::node(coeff)) {
            return DifferentialGeometryVectorResult::failure(CasErrc::InvalidArgument,
                                                             "form coefficients cannot be null",
                                                             operation);
        }
    }
    const auto output_count = bounded_binomial(vars.size(),
                                               static_cast<std::size_t>(degree) + 1,
                                               context.limits().max_expansion_terms,
                                               operation);
    if (!output_count) return DifferentialGeometryVectorResult::failure(output_count.error());
    try {
        const auto input_indices = index_combinations(vars.size(),
                                                      static_cast<std::size_t>(degree),
                                                      input_count.value());
        const auto output_indices = index_combinations(vars.size(),
                                                       static_cast<std::size_t>(degree) + 1,
                                                       output_count.value());
        std::map<std::vector<std::size_t>, std::size_t> coefficient_index;
        for (std::size_t i = 0; i < input_indices.size(); ++i) {
            coefficient_index.emplace(input_indices[i], i);
        }

        std::vector<std::shared_ptr<SymbolicExpr>> result;
        result.reserve(output_count.value());
        for (const auto& output_index : output_indices) {
            auto component = SymbolicExpr::number(0);
            for (std::size_t removed = 0; removed < output_index.size(); ++removed) {
                auto term_budget = context.consume_steps(4, operation);
                if (!term_budget) {
                    return DifferentialGeometryVectorResult::failure(term_budget.error());
                }
                auto node_budget = context.reserve_nodes(4, operation);
                if (!node_budget) {
                    return DifferentialGeometryVectorResult::failure(node_budget.error());
                }

                auto source_index = output_index;
                const std::size_t derivative_variable = source_index[removed];
                source_index.erase(source_index.begin() + static_cast<std::ptrdiff_t>(removed));
                const auto source = coefficient_index.find(source_index);
                if (source == coefficient_index.end()) {
                    return DifferentialGeometryVectorResult::failure(
                        CasErrc::InternalInvariant,
                        "exterior derivative could not locate an input coefficient",
                        operation);
                }
                auto partial = form_coeffs[source->second]->differentiate(
                    vars[derivative_variable]);
                if (!partial || !lamina::detail::node(partial)) {
                    return DifferentialGeometryVectorResult::failure(
                        CasErrc::InternalInvariant,
                        "exterior derivative produced a null partial derivative",
                        operation);
                }
                if (removed % 2 != 0) {
                    partial = SymbolicExpr::multiply(SymbolicExpr::number(-1), partial);
                }
                component = SymbolicExpr::add(component, partial);
            }
            component = component->simplify();
            if (!component || !lamina::detail::node(component)) {
                return DifferentialGeometryVectorResult::failure(
                    CasErrc::InternalInvariant,
                    "exterior derivative produced a null component",
                    operation);
            }
            result.push_back(std::move(component));
        }
        return DifferentialGeometryVectorResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return DifferentialGeometryVectorResult::failure(CasErrc::ResourceLimit,
                                                         "allocation failed while calculating exterior derivative",
                                                         operation);
    } catch (const std::exception& ex) {
        return DifferentialGeometryVectorResult::failure(CasErrc::InternalInvariant,
                                                         ex.what(),
                                                         operation);
    }
}

DifferentialGeometryVectorResult exterior_derivative_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& form_coeffs,
    int degree,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return exterior_derivative_checked(form_coeffs, degree, vars, context);
}


} // namespace lamina
