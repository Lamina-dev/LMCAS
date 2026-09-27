#include "internal/facts_query.hpp"
#include "computation_context.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/rewrite_budget.hpp"

#include <limits>
#include <new>

namespace LMCAS::detail {

Result<bool> ScopedFacts::shadows(const std::shared_ptr<const SymbolicNode>& node,
                                 ComputationContext& context) const {
    auto access = context.consume_steps(0, "facts.scoped");
    if (!access) { return Result<bool>::failure(access.error()); }
    try {
        RewriteBudget budget(context, context.limits().max_recursion_depth,
                             std::numeric_limits<std::size_t>::max(), "facts.scoped");
        return expression_depends_on_variable(node, bound_name_, &budget);
    } catch (const CasError& error) {
        return Result<bool>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<bool>::failure(CasErrc::ResourceLimit,
                                    "scoped facts allocation failed", "facts.scoped");
    }
}

Result<Tribool> ScopedFacts::is_positive(const std::shared_ptr<const SymbolicNode>& node,
                                       ComputationContext& context) const {
    auto hidden = shadows(node, context);
    if (!hidden) { return Result<Tribool>::failure(hidden.error()); }
    return hidden.value() ? Result<Tribool>(Tribool::Unknown) : outer_.is_positive(node, context);
}

Result<Tribool> ScopedFacts::is_negative(const std::shared_ptr<const SymbolicNode>& node,
                                       ComputationContext& context) const {
    auto hidden = shadows(node, context);
    if (!hidden) { return Result<Tribool>::failure(hidden.error()); }
    return hidden.value() ? Result<Tribool>(Tribool::Unknown) : outer_.is_negative(node, context);
}

Result<Tribool> ScopedFacts::is_nonnegative(const std::shared_ptr<const SymbolicNode>& node,
                                          ComputationContext& context) const {
    auto hidden = shadows(node, context);
    if (!hidden) { return Result<Tribool>::failure(hidden.error()); }
    return hidden.value() ? Result<Tribool>(Tribool::Unknown) : outer_.is_nonnegative(node, context);
}

Result<Tribool> ScopedFacts::is_nonzero(const std::shared_ptr<const SymbolicNode>& node,
                                      ComputationContext& context) const {
    auto hidden = shadows(node, context);
    if (!hidden) { return Result<Tribool>::failure(hidden.error()); }
    return hidden.value() ? Result<Tribool>(Tribool::Unknown) : outer_.is_nonzero(node, context);
}

Result<Tribool> ScopedFacts::is_real(const std::shared_ptr<const SymbolicNode>& node,
                                   ComputationContext& context) const {
    auto hidden = shadows(node, context);
    if (!hidden) { return Result<Tribool>::failure(hidden.error()); }
    return hidden.value() ? Result<Tribool>(Tribool::Unknown) : outer_.is_real(node, context);
}

namespace {

class NoFacts final : public FactsQuery {
public:
    Result<Tribool> is_positive(const std::shared_ptr<const SymbolicNode>&, ComputationContext& context) const override {
        return unknown(context);
    }

    Result<Tribool> is_negative(const std::shared_ptr<const SymbolicNode>&, ComputationContext& context) const override {
        return unknown(context);
    }

    Result<Tribool> is_nonnegative(const std::shared_ptr<const SymbolicNode>&, ComputationContext& context) const override {
        return unknown(context);
    }

    Result<Tribool> is_nonzero(const std::shared_ptr<const SymbolicNode>&, ComputationContext& context) const override {
        return unknown(context);
    }

    Result<Tribool> is_real(const std::shared_ptr<const SymbolicNode>&, ComputationContext& context) const override {
        return unknown(context);
    }

private:
    static Result<Tribool> unknown(ComputationContext& context) {
        auto access = context.consume_steps(0, "facts.empty");
        if (!access) { return Result<Tribool>::failure(access.error()); }
        return Tribool::Unknown;
    }
};

}

const FactsQuery& no_facts() noexcept {
    static const NoFacts facts;
    return facts;
}

}
