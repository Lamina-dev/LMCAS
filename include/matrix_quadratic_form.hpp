/**
 * @file matrix_quadratic_form.hpp
 * @brief Extraction of symmetric matrices from homogeneous quadratic forms.
 */
#pragma once

#include "symbolic.hpp"

#include <memory>
#include <string>
#include <vector>

namespace LMCAS {

/**
 * @brief Extract the symmetric matrix A from the quadratic form x^T A x.
 *
 * Coefficients may depend on parameters outside @p vars. The expression must
 * reconstruct exactly from the extracted matrix; nonhomogeneous expressions,
 * higher-degree terms, unsupported derivatives, and invalid variable lists
 * are rejected.
 *
 * @return The extracted matrix, or nullptr when validation fails.
 */
LMCAS_API std::shared_ptr<SymbolicExpr> quadratic_form_matrix(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& vars);

}
