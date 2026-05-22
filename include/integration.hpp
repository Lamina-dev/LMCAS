#pragma once
#include "symbolic.hpp"
#include "matcher.hpp"
#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <unordered_map>

namespace lamina {

struct LAMINA_API IntegrationEntry {

    std::string name;

    SymbolicExpr pattern;

    SymbolicExpr result;

    std::unordered_set<std::string> wildcards;

    std::function<bool(const MatchMap&, const std::string& var)> condition;

    int priority = 100;

    IntegrationEntry() = default;
    IntegrationEntry(std::string name, SymbolicExpr pat, SymbolicExpr res,
                     std::unordered_set<std::string> wc,
                     std::function<bool(const MatchMap&, const std::string& var)> cond = nullptr,
                     int prio = 100)
        : name(std::move(name)), pattern(std::move(pat)), result(std::move(res)),
          wildcards(std::move(wc)), condition(std::move(cond)), priority(prio) {}
};

class LAMINA_API IntegrationTable {
public:
    enum class Category {
        Polynomial,
        Exponential,
        Logarithmic,
        Trigonometric,
        InverseTrig,
        Hyperbolic,
        Algebraic,
        Special,
        UserDefined
    };

    IntegrationTable();

    void add_entry(Category cat, const IntegrationEntry& entry);

    void clear_category(Category cat);

    const std::vector<IntegrationEntry>& get_entries(Category cat) const;

    std::vector<const IntegrationEntry*> get_all_sorted() const;

    void load_defaults();

private:
    std::unordered_map<int, std::vector<IntegrationEntry>> entries_;
    static const std::vector<IntegrationEntry> empty_entries_;
};

class Integrator;

class LAMINA_API IntegrationStrategy {
public:
    virtual ~IntegrationStrategy() = default;

    virtual std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr,
        const std::string& var,
        Integrator& ctx,
        int depth = 0) = 0;

    virtual std::string name() const = 0;
};

class LAMINA_API TableLookupStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "TableLookup"; }
};

class LAMINA_API PowerRuleStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "PowerRule"; }
};

class LAMINA_API SubstitutionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "Substitution"; }
};

class LAMINA_API PartialFractionStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "PartialFraction"; }
};

class LAMINA_API IBPStrategy : public IntegrationStrategy {
public:
    std::shared_ptr<SymbolicExpr> try_integrate(
        const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth = 0) override;
    std::string name() const override { return "IntegrationByParts"; }
};

class LAMINA_API Integrator {
public:
    Integrator();

    SymbolicExpr integrate(const SymbolicExpr& expr, const std::string& var_name);

    SymbolicExpr integrate_def(const SymbolicExpr& expr, const std::string& var_name,
                               const SymbolicExpr& lower, const SymbolicExpr& upper);

    void add_strategy(std::unique_ptr<IntegrationStrategy> strategy, int position = -1);

    IntegrationTable& table() { return table_; }
    const IntegrationTable& table() const { return table_; }

    std::shared_ptr<SymbolicExpr> integrate_recursive(
        const SymbolicExpr& expr, const std::string& var, int depth = 0);

    static bool depends_on(const SymbolicExpr& expr, const std::string& var);

    int max_depth() const { return max_depth_; }
    void set_max_depth(int d) { max_depth_ = d; }

private:
    IntegrationTable table_;
    std::vector<std::unique_ptr<IntegrationStrategy>> strategies_;

    struct CycleState {
        std::vector<SymbolicExpr> history;
    };
    CycleState cycle_state_;

    int max_depth_ = 8;

    std::shared_ptr<SymbolicExpr> apply_linearity(
        const SymbolicExpr& expr, const std::string& var);

    std::shared_ptr<SymbolicExpr> dispatch_strategies(
        const SymbolicExpr& expr, const std::string& var, int depth);

    static std::shared_ptr<SymbolicExpr> make_integral_node(
        const SymbolicExpr& expr, const std::string& var);

    std::shared_ptr<SymbolicExpr> check_cycle(
        const SymbolicExpr& expr, const std::string& var);
    void resolve_cycle(std::shared_ptr<SymbolicExpr>& result, size_t cycle_idx);
};

}
