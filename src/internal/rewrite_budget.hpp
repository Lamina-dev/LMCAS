#pragma once

#include "computation_context.hpp"

#include <cstddef>
#include <memory>
#include <string>

namespace LMCAS {
class SymbolicNode;
}

namespace LMCAS::detail {
class RewriteBudget {
public:
    RewriteBudget(ComputationContext& context, std::size_t max_depth,
                  std::size_t max_nodes, const char* operation);
    RewriteBudget(ComputationContext& context, std::size_t max_steps,
                  std::size_t max_depth, std::size_t max_nodes,
                  const char* operation);
    RewriteBudget(const RewriteBudget&) = delete;
    RewriteBudget& operator=(const RewriteBudget&) = delete;
    ComputationContext& context() const noexcept { return context_; }
    std::size_t max_nodes() const noexcept { return max_nodes_; }
    void consume();

    void enter();
    void leave() noexcept;
    void require_nodes(std::size_t nodes) const;
    std::size_t measure(const std::shared_ptr<const SymbolicNode>& node) const;
    std::size_t append_size(std::size_t current, std::size_t child);
    static std::size_t input_nodes(
        const std::shared_ptr<const SymbolicNode>& lhs,
        const std::shared_ptr<const SymbolicNode>& rhs,
        ComputationContext& context, std::size_t max_depth, const char* operation);

private:
    void check_step_access() const;
    void commit_step();

    ComputationContext& context_;
    std::size_t remaining_steps_;
    std::size_t max_depth_;
    std::size_t max_nodes_;
    std::size_t depth_ = 0;
    bool limit_steps_;
    const std::string operation_;
};
class RewriteScope {
public:
    explicit RewriteScope(RewriteBudget* budget) : budget_(budget) {
        if (budget_) { budget_->enter(); }
    }
    ~RewriteScope() {
        if (budget_) { budget_->leave(); }
    }
    RewriteScope(const RewriteScope&) = delete;
    RewriteScope& operator=(const RewriteScope&) = delete;

private:
    RewriteBudget* budget_;
};

}
