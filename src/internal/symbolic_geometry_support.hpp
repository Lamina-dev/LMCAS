#pragma once

#include "symbolic_geometry_line_plane.hpp"
#include "symbolic_geometry_quadric.hpp"

namespace LMCAS::geometry_detail {

Result<void> validate_symbolic_vector(
    const std::vector<std::shared_ptr<SymbolicExpr>>& vector,
    const std::string& name,
    const std::string& operation);

Result<void> validate_same_dimension_vectors(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context,
    const std::string& operation);

Result<void> validate_surface_point(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    ComputationContext& context,
    const std::string& operation);

ExpressionResult simplify_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& operation,
    const std::string& message);

Result<double> numeric_vector_scale_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& expressions,
    ComputationContext& context,
    const std::string& operation,
    std::vector<double>* values = nullptr);

Result<void> checked_nonzero_numeric_vector(
    const std::vector<std::shared_ptr<SymbolicExpr>>& expressions,
    ComputationContext& context,
    const std::string& operation,
    const std::string& domain_message,
    const std::string& inconclusive_message);

Result<void> checked_nonzero_numeric_or_exact(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation,
    const std::string& domain_message,
    const std::string& inconclusive_message);

Result<void> validate_line_checked(const LineSymbolic& line,
                                   ComputationContext& context,
                                   const std::string& operation);

Result<void> validate_plane_checked(const PlaneSymbolic& plane,
                                    ComputationContext& context,
                                    const std::string& operation);

}
