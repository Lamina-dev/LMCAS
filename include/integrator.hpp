#pragma once

#include "symbolic.hpp"
#include "computation_context.hpp"
#include <cstddef>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace LMCAS {
namespace detail {
struct IntegrationEntry;
struct IntegralParts;
}

class LMCAS_API Integrator {
public:
    Integrator();
    ~Integrator();
    Integrator(Integrator&&) noexcept;
    Integrator& operator=(Integrator&&) noexcept;
    Integrator(const Integrator&) = delete;
    Integrator& operator=(const Integrator&) = delete;

    Result<SymbolicExpr> integrate_checked(
        const SymbolicExpr& expr, const std::string& var_name,
        ComputationContext& context);

    Result<SymbolicExpr> integrate_def_checked(
        const SymbolicExpr& expr, const std::string& var_name,
        const SymbolicExpr& lower, const SymbolicExpr& upper,
        ComputationContext& context);

private:
    Result<detail::IntegralParts> integrate_recursive(
        const SymbolicExpr&, const std::string&, ComputationContext&, int,
        const SymbolicExpr* original = nullptr);
    Result<std::optional<detail::IntegralParts>> apply_linearity(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> dispatch_builtins(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_power_rule(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_table_lookup(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_linear_substitution(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_substitution(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_trig_substitution(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_trig_combination(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_weierstrass(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_rational_decomposition(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_special_function(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);
    Result<std::optional<detail::IntegralParts>> try_integration_by_parts(
        const SymbolicExpr&, const std::string&, ComputationContext&, int);

    std::vector<detail::IntegrationEntry> entries_;
    SymbolicExpr zero_;
    struct CycleState {
        struct Entry {
            SymbolicExpr integrand;
            SymbolicExpr original;
            std::string placeholder;
        };
        std::vector<Entry> history;
        std::set<std::string> occupied_names;
        std::size_t next_placeholder = 0;
    };
    CycleState cycle_state_;
    std::size_t query_depth_ = 0;

    std::shared_ptr<SymbolicExpr> check_cycle(
        const SymbolicExpr&, const std::string&, ComputationContext&);
    bool resolve_cycle(detail::IntegralParts&, const std::string&,
        const std::string&, ComputationContext&);
    bool has_cycle_placeholder(const SymbolicExpr&) const;
};

} // namespace LMCAS
