#pragma once

#include "internal/facts_query.hpp"
#include "inference_engine.hpp"
#include <set>

namespace LMCAS {

class AssumptionContext;

namespace detail {
class RewriteBudget;
class AssumptionFacts final : public FactsQuery {
public:
    explicit AssumptionFacts(const AssumptionContext& assumptions)
        : owned_engine_(std::in_place, assumptions), engine_(&*owned_engine_) {}
    explicit AssumptionFacts(const InferenceEngine& engine) noexcept : engine_(&engine) {}
    AssumptionFacts(const AssumptionFacts&) = delete;
    AssumptionFacts& operator=(const AssumptionFacts&) = delete;
    static void collect_names(const AssumptionContext& assumptions,
        std::set<std::string>& names, RewriteBudget& budget);

    Result<Tribool> is_positive(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;
    Result<Tribool> is_negative(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;
    Result<Tribool> is_nonnegative(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;
    Result<Tribool> is_nonzero(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;
    Result<Tribool> is_real(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;

private:
    std::optional<InferenceEngine> owned_engine_;
    const InferenceEngine* engine_;
};

}
}
