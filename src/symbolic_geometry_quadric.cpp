#include "symbolic_geometry_quadric.hpp"
#include "internal/symbolic_geometry_support.hpp"
#include "internal/symbolic_quadric_classification.hpp"
#include <stdexcept>
#include <utility>

namespace LMCAS {
using namespace geometry_detail;

VectorStringResult classify_quadric_checked(
    const SurfaceSymbolic& surf,
    ComputationContext& context) {
    const std::string operation = "classify_quadric";
    try {
        auto valid = validate_surface_point(
            surf,
            std::vector<std::shared_ptr<SymbolicExpr>>(
                surf.vars.size(), SymbolicExpr::number(0)),
            context,
            operation);
        if (!valid) return VectorStringResult::failure(valid.error());
        if (surf.vars.size() != 3) {
            return VectorStringResult::failure(
                CasErrc::InvalidArgument,
                "quadric classification requires exactly three variables",
                operation);
        }
        auto step = context.consume_steps(12, operation);
        if (!step) return VectorStringResult::failure(step.error());
        return detail::classify_quadric_impl(surf, context, operation);
    } catch (const std::bad_alloc&) {
        return VectorStringResult::failure(CasErrc::ResourceLimit,
                                           "vector geometry allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return VectorStringResult::failure(CasErrc::Inconclusive,
                                           e.what(), operation);
    }
}

VectorStringResult classify_quadric_checked(const SurfaceSymbolic& surf) {
    ComputationContext context;
    return classify_quadric_checked(surf, context);
}

}
