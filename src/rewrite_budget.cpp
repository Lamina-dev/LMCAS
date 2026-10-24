#include "internal/rewrite_budget.hpp"
#include "internal/symbolic_ast.hpp"

#include <algorithm>
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
}

class RewriteBudget::Measurement final : public RecursiveSymbolicVisitor {
public:
    // A fresh visitor handles each node; MeasurementRecursion checks the full path.
    Measurement(RewriteBudget& budget, std::size_t path_depth, std::size_t active_depth)
        : budget_(budget), path_depth_(path_depth), active_depth_(active_depth) {}

    Shape shape{1, 1};

    void visit(const NumberNode& number) override {
        if (const auto* integer = std::get_if<BigInt>(&number.value())) {
            require_bits(*integer);
        } else {
            const auto& rational = std::get<Rational>(number.value());
            require_bits(rational.get_numerator());
            require_bits(rational.get_denominator());
        }
    }

protected:
    void visit_child(const SymbolicNodePtr& child) override {
        const auto measured = budget_.measure_shape(child, path_depth_ + 1, active_depth_);
        shape.nodes = budget_.append_size(shape.nodes, measured.nodes);
        if (measured.height == std::numeric_limits<std::size_t>::max()) {
            throw CasError{CasErrc::ResourceLimit, kDepthMessage, budget_.operation_};
        }
        shape.height = std::max(shape.height, measured.height + 1);
    }

private:
    void require_bits(const BigInt& integer) {
        auto checked = budget_.context_.require_integer_bits(
            integer.bit_length(), budget_.operation_);
        if (!checked) { throw checked.error(); }
    }

    RewriteBudget& budget_;
    const std::size_t path_depth_;
    const std::size_t active_depth_;
};


RewriteBudget::RewriteBudget(ComputationContext& context, std::size_t max_depth,
                             std::size_t max_nodes, const char* operation)
    : context_(context), remaining_steps_(std::numeric_limits<std::size_t>::max()),
      max_depth_(max_depth), max_nodes_(std::min(max_nodes, context.limits().max_ast_nodes)),
      limit_steps_(false),
      operation_(operation) {}

RewriteBudget::RewriteBudget(ComputationContext& context, std::size_t max_steps,
                             std::size_t max_depth, std::size_t max_nodes,
                             const char* operation)
    : context_(context), remaining_steps_(max_steps), max_depth_(max_depth),
      max_nodes_(std::min(max_nodes, context.limits().max_ast_nodes)),
      limit_steps_(true), operation_(operation) {}

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
    auto access = context_.require_ast_nodes(nodes, operation_);
    if (!access) { throw access.error(); }
    if (nodes > max_nodes_) {
        throw CasError{CasErrc::ResourceLimit, kGrowthMessage, operation_};
    }
}

void RewriteBudget::require_shape(const Shape& shape, std::size_t path_depth,
                                  std::size_t active_depth) const {
    require_nodes(shape.nodes);
    const auto context_depth = context_.recursion_depth();
    const auto context_limit = context_.limits().max_recursion_depth;
    if (active_depth > max_depth_ || shape.height > max_depth_ - active_depth ||
        path_depth > max_depth_ || shape.height > max_depth_ - path_depth ||
        context_depth > context_limit || shape.height > context_limit - context_depth) {
        throw CasError{CasErrc::ResourceLimit, kDepthMessage, operation_};
    }
}

RewriteBudget::Shape RewriteBudget::measure_shape(
    const SymbolicNodePtr& node, std::size_t path_depth, std::size_t active_depth) {
    check_step_access();
    auto step = context_.consume_steps(1, operation_);
    if (!step) { throw step.error(); }
    if (!node) { return {}; }
    {
        const auto found = shapes_.find(node.get());
        if (found != shapes_.end()) {
            const auto owner = found->second.owner.lock();
            if (owner && owner.get() == node.get()) {
                const auto shape = found->second.shape;
                require_shape(shape, path_depth, active_depth);
                return shape;
            }
            shapes_.erase(found);
        }
    }
    if (path_depth >= max_depth_) {
        throw CasError{CasErrc::ResourceLimit, kDepthMessage, operation_};
    }
    Shape shape;
    {
        MeasurementRecursion recursion(context_, operation_);
        Measurement measurement(*this, path_depth, active_depth);
        node->accept(measurement);
        shape = measurement.shape;
    }
    require_shape(shape, path_depth, active_depth);
    if (max_nodes_ != 0) {
        if (shapes_.size() >= max_nodes_) {
            auto cleared = context_.consume_steps(shapes_.size(), operation_);
            if (!cleared) { throw cleared.error(); }
            shapes_.clear();
        }
        shapes_.emplace(node.get(), Entry{node, shape});
    }
    return shape;
}

std::size_t RewriteBudget::measure(const SymbolicNodePtr& node) {
    return measure_shape(node, 0, depth_).nodes;
}

std::size_t RewriteBudget::measure_result(const SymbolicNodePtr& node) {
    // The result root replaces the current rewrite frame.
    return measure_shape(node, 0, depth_ == 0 ? 0 : depth_ - 1).nodes;
}

std::size_t RewriteBudget::append_size(std::size_t current, std::size_t child) {
    auto access = context_.consume_steps(0, operation_);
    if (!access) { throw access.error(); }
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
    RewriteBudget budget(context, max_depth, context.limits().max_ast_nodes, operation);
    const auto left = budget.measure(lhs);
    const auto right = budget.measure(rhs);
    if (right > std::numeric_limits<std::size_t>::max() - left) {
        throw CasError{CasErrc::ResourceLimit, kGrowthMessage, operation};
    }
    return left + right;
}

}
