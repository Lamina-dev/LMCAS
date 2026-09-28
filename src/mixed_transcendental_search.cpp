#include "solve_mixed_transcendental.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include <cmath>
#include <limits>

namespace LMCAS {
namespace {
static lmmc_real_t extract_real_value(const std::shared_ptr<const SymbolicNode>& node) {
    auto num = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!num) {
        return std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }

    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        return std::get<lmmc_real_t>(num->value());
    }
    if (std::holds_alternative<Rational>(num->value())) {
        return static_cast<lmmc_real_t>(std::get<Rational>(num->value()).to_double());
    }
    if (std::holds_alternative<BigInt>(num->value())) {
        return static_cast<lmmc_real_t>(std::get<BigInt>(num->value()).to_double());
    }
    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

// Accept k*x in either operand order. The caller also handles x and k*x+c.
static lmmc_real_t product_linear_coefficient(
    const MultiplyNode& product, const std::string& var) {
    if (product.operands().size() != 2) {
        return std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }
    std::shared_ptr<const SymbolicNode> number;
    bool has_variable = false;
    for (const auto& operand : product.operands()) {
        if (operand->is_number()) {
            number = operand;
        } else if (auto variable = std::dynamic_pointer_cast<const VariableNode>(operand)) {
            if (!variable->is_constant() && variable->name() == var) {
                has_variable = true;
            }
        }
    }
    if (number && has_variable) {
        return extract_real_value(number);
    }
    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

static lmmc_real_t extract_linear_coefficient(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var)
{
    if (!node || !expression_depends_on_variable(node, var)) {
        return 0.0;
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (!v->is_constant() && v->name() == var) {
            return 1.0;
        }
        return std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return product_linear_coefficient(*mul, var);
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        lmmc_real_t coeff = std::numeric_limits<lmmc_real_t>::quiet_NaN();
        bool found_var_term = false;

        for (const auto& op : add->operands()) {
            if (!expression_depends_on_variable(op, var)) {
                continue;
            }
            if (found_var_term) {
                return std::numeric_limits<lmmc_real_t>::quiet_NaN();
            }
            found_var_term = true;
            coeff = extract_linear_coefficient(op, var);
        }

        return found_var_term ? coeff : std::numeric_limits<lmmc_real_t>::quiet_NaN();
    }

    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

static bool is_periodic_func(FunctionNode::FuncType t) {
    switch (t) {
    case FunctionNode::FuncType::Sin:
    case FunctionNode::FuncType::Cos:
    case FunctionNode::FuncType::Tan:
        return true;
    default:
        return false;
    }
}

// Record the periods of trigonometric nodes with affine arguments in target_var.
// Sine and cosine use 2pi/|k|; tangent uses pi/|k|.
struct PeriodicCollector : public LMCAS::detail::RecursiveSymbolicVisitor {
    const std::string& target_var;
    lmmc_real_t max_period = 0.0;
    bool found_periodic = false;
    bool has_nonlinear_periodic = false;

    explicit PeriodicCollector(const std::string& v) : target_var(v) {}

    void visit(const MatrixNode&) override {}

    void visit(const TransformNode& node) override {
        visit_child(node.body());
    }

    void collect_period(const FunctionNode& n) {
        if (!is_periodic_func(n.type()) || n.arguments().empty()) {
            return;
        }
        if (!expression_depends_on_variable(n.arguments()[0], target_var)) {
            return;
        }
        found_periodic = true;
        const auto k = extract_linear_coefficient(n.arguments()[0], target_var);
        if (std::isnan(k) || k == 0.0) {
            has_nonlinear_periodic = true;
            return;
        }
        const auto period = n.type() == FunctionNode::FuncType::Tan
            ? LMMC_CONST_PI / std::fabs(k)
            : 2.0 * LMMC_CONST_PI / std::fabs(k);
        if (period > max_period) {
            max_period = period;
        }
    }

    void visit(const FunctionNode& node) override {
        collect_period(node);
        LMCAS::detail::RecursiveSymbolicVisitor::visit(node);
    }
};

}

std::optional<SearchInterval> determine_search_interval(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts)
{
    if (opts.has_search_interval) {
        if (opts.search_lo >= opts.search_hi) {
            return std::nullopt;
        }
        if ((opts.search_hi - opts.search_lo) <= opts.tolerance) {
            return std::nullopt;
        }
        return SearchInterval{opts.search_lo, opts.search_hi};
    }

    lmmc_real_t lo = -10.0;
    lmmc_real_t hi = 10.0;

    if (expr && LMCAS::detail::node(expr)) {
        PeriodicCollector collector(var);
        LMCAS::detail::node(expr)->accept(collector);

        if (collector.found_periodic && !collector.has_nonlinear_periodic && collector.max_period > 0.0) {
            // [-period, period] covers two complete periods around zero.
            lmmc_real_t half_span = collector.max_period;
            if (half_span > lo * -1.0 || half_span > hi) {
                lo = -half_span;
                hi = half_span;
            }
            // Cap automatic interval growth to keep numeric scanning bounded.
            if (lo < -100.0) {
                lo = -100.0;
            }
            if (hi > 100.0) {
                hi = 100.0;
            }
        }
    }

    if (lo >= hi) {
        return std::nullopt;
    }
    if ((hi - lo) <= opts.tolerance) {
        return std::nullopt;
    }

    return SearchInterval{lo, hi};
}

}
