#pragma once
#include "symbolic.hpp"
#include "matcher.hpp"
#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <unordered_map>

namespace lamina {

// ============================================================================
// Integration Table Entry
// ============================================================================

// Represents a single integration rule with pattern matching and optional conditions.
// This is the fundamental building block of the table-driven integration system.
struct LAMINA_API IntegrationEntry {
    // Human-readable name for debugging
    std::string name;
    
    // The pattern to match against the integrand (uses wildcards)
    SymbolicExpr pattern;
    
    // The antiderivative template (uses same wildcards as pattern)
    SymbolicExpr result;
    
    // Wildcard names used in pattern/result
    std::unordered_set<std::string> wildcards;
    
    // Optional condition: return true if this rule applies given the bindings
    std::function<bool(const MatchMap&, const std::string& var)> condition;
    
    // Priority: lower = tried first (allows overriding)
    int priority = 100;
    
    IntegrationEntry() = default;
    IntegrationEntry(std::string name, SymbolicExpr pat, SymbolicExpr res,
                     std::unordered_set<std::string> wc,
                     std::function<bool(const MatchMap&, const std::string& var)> cond = nullptr,
                     int prio = 100)
        : name(std::move(name)), pattern(std::move(pat)), result(std::move(res)),
          wildcards(std::move(wc)), condition(std::move(cond)), priority(prio) {}
};

// ============================================================================
// Integration Table
// ============================================================================

// A categorized collection of integration rules.
// Rules are organized by category for efficient lookup and easy extension.
class LAMINA_API IntegrationTable {
public:
    enum class Category {
        Polynomial,       // x^n, constants
        Exponential,      // e^x, a^x
        Logarithmic,      // ln(x), log_a(x)
        Trigonometric,    // sin, cos, tan, etc.
        InverseTrig,      // arcsin, arccos, arctan
        Hyperbolic,       // sinh, cosh, tanh
        Algebraic,        // 1/(ax+b), 1/sqrt(...)
        Special,          // LambertW, etc.
        UserDefined       // User-added rules
    };
    
    IntegrationTable();
    
    // Add a rule to a specific category
    void add_entry(Category cat, const IntegrationEntry& entry);
    
    // Remove all rules in a category (for replacement)
    void clear_category(Category cat);
    
    // Get all entries (sorted by priority)
    const std::vector<IntegrationEntry>& get_entries(Category cat) const;
    
    // Get all entries across all categories, sorted by priority
    std::vector<const IntegrationEntry*> get_all_sorted() const;
    
    // Load the default built-in rules
    void load_defaults();
    
private:
    std::unordered_map<int, std::vector<IntegrationEntry>> entries_;
    static const std::vector<IntegrationEntry> empty_entries_;
};

// ============================================================================
// Integration Strategy Interface
// ============================================================================

// Forward declaration
class Integrator;

// Abstract base for integration strategies.
// Each strategy encapsulates one approach to integration.
class LAMINA_API IntegrationStrategy {
public:
    virtual ~IntegrationStrategy() = default;
    
    // Attempt to integrate expr w.r.t. var.
    // Returns nullptr if this strategy cannot handle the expression.
    // @param depth Current recursion depth (strategies should pass depth+1
    //              when calling ctx.integrate_recursive to ensure the global
    //              depth limit is respected across nested strategy invocations).
    virtual std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr,
        const std::string& var,
        Integrator& ctx,
        int depth = 0) = 0;
    
    // Human-readable name for debugging
    virtual std::string name() const = 0;
};

// ============================================================================
// Concrete Strategies
// ============================================================================

// Strategy: Table lookup via pattern matching
class LAMINA_API TableLookupStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "TableLookup"; }
};

// Strategy: Power rule (x^n -> x^(n+1)/(n+1))
class LAMINA_API PowerRuleStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "PowerRule"; }
};

// Strategy: U-substitution
class LAMINA_API SubstitutionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "Substitution"; }
};

// Strategy: Partial fractions for rational functions
class LAMINA_API PartialFractionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "PartialFraction"; }
};

// Strategy: Integration by parts (LIATE heuristic)
class LAMINA_API IBPStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "IntegrationByParts"; }
};

// ============================================================================
// Integrator (Main Class - Redesigned)
// ============================================================================

class LAMINA_API Integrator {
public:
    Integrator();
    
    // Main entry point: integrate expr w.r.t. var_name
    SymbolicExpr integrate(const SymbolicExpr& expr, const std::string& var_name);
    
    // Definite integration
    SymbolicExpr integrate_def(const SymbolicExpr& expr, const std::string& var_name,
                               const SymbolicExpr& lower, const SymbolicExpr& upper);
    
    // --- Extension API ---
    
    // Add a custom integration strategy (inserted at given position in pipeline)
    void add_strategy(std::unique_ptr<IntegrationStrategy> strategy, int position = -1);
    
    // Access the integration table for adding custom rules
    IntegrationTable& table() { return table_; }
    const IntegrationTable& table() const { return table_; }
    
    // --- Internal API (used by strategies) ---
    
    // Recursive integration with cycle detection and depth limiting
    std::shared_ptr<SymbolicExpr> integrate_recursive(
        const SymbolicExpr& expr, const std::string& var, int depth = 0);
    
    // Check if expression depends on variable
    static bool depends_on(const SymbolicExpr& expr, const std::string& var);
    
    // Get max recursion depth
    int max_depth() const { return max_depth_; }
    void set_max_depth(int d) { max_depth_ = d; }
    
private:
    IntegrationTable table_;
    std::vector<std::unique_ptr<IntegrationStrategy>> strategies_;
    
    // Cycle detection state (per-integration-call)
    struct CycleState {
        std::vector<SymbolicExpr> history;
    };
    CycleState cycle_state_;
    
    int max_depth_ = 8;
    
    // Apply linearity: split sums and factor out constants
    std::shared_ptr<SymbolicExpr> apply_linearity(
        const SymbolicExpr& expr, const std::string& var);
    
    // Core dispatch: try each strategy in order
    std::shared_ptr<SymbolicExpr> dispatch_strategies(
        const SymbolicExpr& expr, const std::string& var, int depth);
    
    // Create a symbolic "unevaluated integral" node
    static std::shared_ptr<SymbolicExpr> make_integral_node(
        const SymbolicExpr& expr, const std::string& var);
    
    // Cycle detection and resolution
    std::shared_ptr<SymbolicExpr> check_cycle(
        const SymbolicExpr& expr, const std::string& var);
    void resolve_cycle(std::shared_ptr<SymbolicExpr>& result, size_t cycle_idx);
};

} // namespace lamina
