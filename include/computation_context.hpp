#pragma once

#include <atomic>
#include <cstddef>
#include <functional>
#include <memory>
#include <new>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "result.hpp"
#include "unit.hpp"

namespace LMCAS {

class AssumptionContext;

/**
 * Per-operation limits for APIs that accept a ComputationContext.
 *
 * Context-free BigInt and symbolic constructors are intentionally unbounded;
 * use checked, context-aware entry points when processing untrusted input.
 */
struct ResourceLimits {
    std::size_t max_ast_nodes = 1'000'000;       ///< Cumulative reserved nodes.
    std::size_t max_steps = 10'000'000;          ///< Cumulative checked steps.
    std::size_t max_recursion_depth = 256;       ///< Simultaneous checked depth.
    std::size_t max_integer_bits = 1'000'000;    ///< One checked integer result.
    std::size_t max_expansion_terms = 100'000;   ///< One checked expansion.
    std::size_t max_input_bytes = 1'048'576;     ///< One checked input buffer.
    std::size_t max_diagnostics = 10'000;        ///< Stored diagnostics.
};

enum class DiagnosticSeverity { Info, Warning, Error };

struct Diagnostic {
    DiagnosticSeverity severity = DiagnosticSeverity::Info;
    std::string operation;
    std::string message;
};

class DiagnosticEngine {
public:
    using Consumer = std::function<void(const Diagnostic&)>;

    explicit DiagnosticEngine(std::size_t limit, Consumer consumer = {})
        : limit_(limit), consumer_(std::move(consumer)) {}

    Result<void> emit(Diagnostic diagnostic) {
        if (diagnostics_.size() >= limit_) {
            return Result<void>::failure(CasErrc::ResourceLimit,
                                         "diagnostic budget exhausted",
                                         diagnostic.operation);
        }
        bool inserted = false;
        try {
            diagnostics_.push_back(std::move(diagnostic));
            inserted = true;
            if (consumer_) consumer_(diagnostics_.back());
        } catch (const std::bad_alloc&) {
            if (inserted) diagnostics_.pop_back();
            return Result<void>::failure(CasErrc::ResourceLimit,
                                         "diagnostic allocation failed");
        } catch (...) {
            if (inserted) diagnostics_.pop_back();
            return Result<void>::failure(CasErrc::InternalInvariant,
                                         "diagnostic consumer failed");
        }
        return Result<void>::success();
    }

    void set_consumer(Consumer consumer) { consumer_ = std::move(consumer); }

    const std::vector<Diagnostic>& diagnostics() const noexcept {
        return diagnostics_;
    }

private:
    std::size_t limit_;
    Consumer consumer_;
    std::vector<Diagnostic> diagnostics_;
};

class CancellationToken {
public:
    CancellationToken() : cancelled_(std::make_shared<std::atomic<bool>>(false)) {}

    void cancel() const noexcept { cancelled_->store(true, std::memory_order_release); }
    bool is_cancelled() const noexcept {
        return cancelled_->load(std::memory_order_acquire);
    }

private:
    std::shared_ptr<std::atomic<bool>> cancelled_;
};

class ComputationContext {
public:
    explicit ComputationContext(ResourceLimits limits = {},
                                CancellationToken cancellation = {})
        : limits_(limits), cancellation_(std::move(cancellation)),
          owner_thread_(std::this_thread::get_id()),
          diagnostic_engine_(limits.max_diagnostics) {}

    ComputationContext(const ComputationContext&) = delete;
    ComputationContext& operator=(const ComputationContext&) = delete;
    ComputationContext(ComputationContext&&) = delete;
    ComputationContext& operator=(ComputationContext&&) = delete;

    Result<void> consume_steps(std::size_t amount = 1,
                               const std::string& operation = {}) {
        auto access = check_access(operation);
        if (!access) return access;
        if (steps_used_ > limits_.max_steps ||
            amount > limits_.max_steps - steps_used_) {
            return Result<void>::failure(CasErrc::ResourceLimit,
                                         "computation step budget exhausted", operation);
        }
        steps_used_ += amount;
        return Result<void>::success();
    }

    Result<void> enter_recursion(const std::string& operation = {}) {
        auto step = consume_steps(1, operation);
        if (!step) return step;
        if (recursion_depth_ >= limits_.max_recursion_depth) {
            return Result<void>::failure(CasErrc::ResourceLimit,
                                         "recursion depth limit exceeded", operation);
        }
        ++recursion_depth_;
        return Result<void>::success();
    }

    void leave_recursion() noexcept {
        if (recursion_depth_ > 0) --recursion_depth_;
    }

    Result<void> reserve_nodes(std::size_t amount,
                               const std::string& operation = {}) {
        auto access = check_access(operation);
        if (!access) return access;
        if (nodes_created_ > limits_.max_ast_nodes ||
            amount > limits_.max_ast_nodes - nodes_created_) {
            return Result<void>::failure(CasErrc::ResourceLimit,
                                         "AST node budget exhausted", operation);
        }
        nodes_created_ += amount;
        return Result<void>::success();
    }

    Result<void> require_integer_bits(
        std::size_t amount, const std::string& operation = {}) const {
        auto access = check_access(operation);
        if (!access) return access;
        if (amount > limits_.max_integer_bits) {
            return Result<void>::failure(
                CasErrc::ResourceLimit,
                "integer bit budget exceeded", operation);
        }
        return Result<void>::success();
    }

    Result<void> require_expansion_terms(
        std::size_t amount, const std::string& operation = {}) const {
        auto access = check_access(operation);
        if (!access) return access;
        if (amount > limits_.max_expansion_terms) {
            return Result<void>::failure(
                CasErrc::ResourceLimit,
                "expansion term budget exceeded", operation);
        }
        return Result<void>::success();
    }

    Result<void> require_input_bytes(
        std::size_t amount, const std::string& operation = {}) const {
        auto access = check_access(operation);
        if (!access) return access;
        if (amount > limits_.max_input_bytes) {
            return Result<void>::failure(
                CasErrc::ResourceLimit,
                "input byte budget exceeded", operation);
        }
        return Result<void>::success();
    }

    Result<void> add_diagnostic(Diagnostic diagnostic) {
        auto access = check_access(diagnostic.operation);
        if (!access) return access;
        return diagnostic_engine_.emit(std::move(diagnostic));
    }

    const ResourceLimits& limits() const noexcept { return limits_; }
    std::size_t steps_used() const noexcept { return steps_used_; }
    std::size_t nodes_created() const noexcept { return nodes_created_; }
    std::size_t recursion_depth() const noexcept { return recursion_depth_; }
    const std::vector<Diagnostic>& diagnostics() const noexcept {
        return diagnostic_engine_.diagnostics();
    }
    const DiagnosticEngine& diagnostic_engine() const noexcept {
        return diagnostic_engine_;
    }
    Result<void> set_diagnostic_consumer(
        DiagnosticEngine::Consumer consumer,
        const std::string& operation = {}) {
        auto access = check_access(operation);
        if (!access) return access;
        diagnostic_engine_.set_consumer(std::move(consumer));
        return Result<void>::success();
    }
    const CancellationToken& cancellation() const noexcept { return cancellation_; }
    UnitSystem& units() noexcept { return units_; }
    const UnitSystem& units() const noexcept { return units_; }

    Result<void> set_assumptions(
        std::shared_ptr<const AssumptionContext> assumptions,
        const std::string& operation = {}) {
        auto access = check_access(operation);
        if (!access) return access;
        assumptions_ = std::move(assumptions);
        return Result<void>::success();
    }
    const std::shared_ptr<const AssumptionContext>& assumptions() const noexcept {
        return assumptions_;
    }

private:
    Result<void> check_access(const std::string& operation) const {
        if (std::this_thread::get_id() != owner_thread_) {
            return Result<void>::failure(
                CasErrc::InternalInvariant,
                "ComputationContext cannot be shared across threads",
                operation);
        }
        if (cancellation_.is_cancelled()) {
            return Result<void>::failure(CasErrc::Cancelled,
                                         "computation was cancelled", operation);
        }
        return Result<void>::success();
    }

    ResourceLimits limits_;
    CancellationToken cancellation_;
    const std::thread::id owner_thread_;
    std::shared_ptr<const AssumptionContext> assumptions_;
    std::size_t steps_used_ = 0;
    std::size_t nodes_created_ = 0;
    std::size_t recursion_depth_ = 0;
    DiagnosticEngine diagnostic_engine_;
    UnitSystem units_;
};

} // namespace LMCAS
