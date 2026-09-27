#include "internal/rewrite_budget.hpp"
#include "internal/symbolic_ast.hpp"

#include <limits>

namespace LMCAS::detail {
namespace {

constexpr const char* kDepthMessage = "equivalence rewrite depth budget exhausted";
constexpr const char* kGrowthMessage = "equivalence node growth budget exhausted";
constexpr const char* kStepMessage = "equivalence rewrite step budget exhausted";
class MeasurementRecursion {
public:
    MeasurementRecursion(ComputationContext& context, const std::string& operation)
        : context_(context) {
        auto entered = context_.enter_recursion(operation);
        if (!entered) { throw entered.error(); }
    }
    ~MeasurementRecursion() { context_.leave_recursion(); }
    MeasurementRecursion(const MeasurementRecursion&) = delete;
    MeasurementRecursion& operator=(const MeasurementRecursion&) = delete;

private:
    ComputationContext& context_;
};
class NodeCounter final : public RecursiveSymbolicVisitor {
public:
    NodeCounter(ComputationContext& context, std::size_t max_depth,
                std::size_t max_nodes, const std::string& operation)
        : context_(context), max_depth_(max_depth), max_nodes_(max_nodes),
          operation_(operation) {}

    std::size_t count(const SymbolicNodePtr& lhs, const SymbolicNodePtr& rhs = {}) {
        visit_child(lhs);
        visit_child(rhs);
        return nodes_;
    }

protected:
    void visit_child(const SymbolicNodePtr& child) override {
        if (!child) { return; }
        if (static_cast<std::size_t>(current_depth) >= max_depth_) {
            throw CasError{CasErrc::ResourceLimit, kDepthMessage, operation_};
        }
        if (nodes_ >= max_nodes_) {
            throw CasError{CasErrc::ResourceLimit, kGrowthMessage, operation_};
        }
        MeasurementRecursion recursion(context_, operation_);
        ++nodes_;
        RecursiveSymbolicVisitor::visit_child(child);
    }

private:
    ComputationContext& context_;
    std::size_t max_depth_;
    std::size_t max_nodes_;
    const std::string& operation_;
    std::size_t nodes_ = 0;
};

}

RewriteBudget::RewriteBudget(ComputationContext& context, std::size_t max_depth,
                             std::size_t max_nodes, const char* operation)
    : context_(context), remaining_steps_(std::numeric_limits<std::size_t>::max()),
      max_depth_(max_depth), max_nodes_(max_nodes), limit_steps_(false),
      operation_(operation) {}

RewriteBudget::RewriteBudget(ComputationContext& context, std::size_t max_steps,
                             std::size_t max_depth, std::size_t max_nodes,
                             const char* operation)
    : context_(context), remaining_steps_(max_steps), max_depth_(max_depth),
      max_nodes_(max_nodes), limit_steps_(true), operation_(operation) {}

void RewriteBudget::check_step_access() const {
    auto access = context_.consume_steps(0, operation_);
    if (!access) { throw access.error(); }
    if (limit_steps_ && remaining_steps_ == 0) {
        throw CasError{CasErrc::ResourceLimit, kStepMessage, operation_};
    }
}

void RewriteBudget::commit_step() {
    auto step = context_.consume_steps(1, operation_);
    if (!step) { throw step.error(); }
    if (limit_steps_) { --remaining_steps_; }
}

void RewriteBudget::consume() {
    check_step_access();
    commit_step();
}

void RewriteBudget::enter() {
    check_step_access();
    if (depth_ >= max_depth_) {
        throw CasError{CasErrc::ResourceLimit, kDepthMessage, operation_};
    }
    auto entered = context_.enter_recursion(operation_);
    if (!entered) { throw entered.error(); }
    if (limit_steps_) { --remaining_steps_; }
    ++depth_;
}

void RewriteBudget::leave() noexcept {
    --depth_;
    context_.leave_recursion();
}

void RewriteBudget::require_nodes(std::size_t nodes) const {
    if (nodes > max_nodes_) {
        throw CasError{CasErrc::ResourceLimit, kGrowthMessage, operation_};
    }
}

std::size_t RewriteBudget::measure(const std::shared_ptr<const SymbolicNode>& node) const {
    NodeCounter counter(context_, max_depth_, max_nodes_, operation_);
    return counter.count(node);
}

std::size_t RewriteBudget::append_size(std::size_t current, std::size_t child) {
    if (child != 0 && current == 0) { consume(); }
    if (current > max_nodes_ || child > max_nodes_ - current) {
        throw CasError{CasErrc::ResourceLimit, kGrowthMessage, operation_};
    }
    return current + child;
}

std::size_t RewriteBudget::input_nodes(
    const std::shared_ptr<const SymbolicNode>& lhs,
    const std::shared_ptr<const SymbolicNode>& rhs,
    ComputationContext& context, std::size_t max_depth, const char* operation) {
    const std::string operation_name(operation);
    NodeCounter counter(context, max_depth, std::numeric_limits<std::size_t>::max(),
                        operation_name);
    return counter.count(lhs, rhs);
}

}
