#include "solver.hpp"
#include "solve_strategies.hpp"
#include "symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "assumption_context.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <queue>
#include <unordered_map>
#include <optional>
#include <limits>

namespace lamina {
namespace {

/**
 * @brief Visitor that checks if an expression tree contains an imaginary unit.
 *
 * The imaginary unit in LMCAS is represented as sqrt(-1), i.e., a FunctionNode
 * with type Sqrt and a single argument that is a NumberNode with value -1.
 */
class ContainsImaginaryVisitor : public lamina::detail::RecursiveSymbolicVisitor {
public:
    bool found = false;

    void visit(const NumberNode&) override {}
    void visit(const VariableNode&) override {}

    void visit(const AddNode& node) override {
        for (auto& op : node.operands()) {
            if (found) return;
            op->accept(*this);
        }
    }

    void visit(const MultiplyNode& node) override {
        for (auto& op : node.operands()) {
            if (found) return;
            op->accept(*this);
        }
    }

    void visit(const PowerNode& node) override {
        if (found) return;

        // Check if this is a negative base raised to a fractional exponent
        // (e.g., (-4)^0.5 which produces an imaginary result)
        auto base_num = std::dynamic_pointer_cast<const NumberNode>(node.base());
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(node.exponent());
        if (base_num && exp_num) {
            double base_val = 0.0;
            bool base_is_numeric = false;
            if (std::holds_alternative<BigInt>(base_num->value())) {
                base_val = std::get<BigInt>(base_num->value()).to_double();
                base_is_numeric = true;
            } else if (std::holds_alternative<Rational>(base_num->value())) {
                base_val = std::get<Rational>(base_num->value()).to_double();
                base_is_numeric = true;
            } else if (std::holds_alternative<lmmc_real_t>(base_num->value())) {
                base_val = std::get<lmmc_real_t>(base_num->value());
                base_is_numeric = true;
            }

            double exp_val = 0.0;
            bool exp_is_numeric = false;
            if (std::holds_alternative<BigInt>(exp_num->value())) {
                exp_val = std::get<BigInt>(exp_num->value()).to_double();
                exp_is_numeric = true;
            } else if (std::holds_alternative<Rational>(exp_num->value())) {
                exp_val = std::get<Rational>(exp_num->value()).to_double();
                exp_is_numeric = true;
            } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
                exp_val = std::get<lmmc_real_t>(exp_num->value());
                exp_is_numeric = true;
            }

            if (base_is_numeric && exp_is_numeric) {
                // Negative base with non-integer exponent → imaginary
                if (base_val < 0.0) {
                    double rounded_exp = std::round(exp_val);
                    if (std::abs(exp_val - rounded_exp) > 1e-12) {
                        found = true;
                        return;
                    }
                }
            }
        }

        node.base()->accept(*this);
        if (found) return;
        node.exponent()->accept(*this);
    }

    void visit(const FunctionNode& node) override {
        // Check if this is sqrt(-1) — the imaginary unit
        if (node.type() == FunctionNode::FuncType::Sqrt && node.arguments().size() == 1) {
            auto num = std::dynamic_pointer_cast<const NumberNode>(node.arguments()[0]);
            if (num) {
                if (std::holds_alternative<BigInt>(num->value())) {
                    if (std::get<BigInt>(num->value()) == BigInt(-1)) {
                        found = true;
                        return;
                    }
                } else if (std::holds_alternative<Rational>(num->value())) {
                    if (std::get<Rational>(num->value()) == Rational(-1)) {
                        found = true;
                        return;
                    }
                } else if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    int eq;
                    lmmc_double_nearly_equal(v, -1.0, &eq);
                    if (eq) {
                        found = true;
                        return;
                    }
                }
            }
        }
        // Also check for sqrt of any negative number (produces imaginary result)
        if (node.type() == FunctionNode::FuncType::Sqrt && node.arguments().size() == 1) {
            auto num = std::dynamic_pointer_cast<const NumberNode>(node.arguments()[0]);
            if (num) {
                if (std::holds_alternative<BigInt>(num->value())) {
                    if (std::get<BigInt>(num->value()).IsNegative()) {
                        found = true;
                        return;
                    }
                } else if (std::holds_alternative<Rational>(num->value())) {
                    if (std::get<Rational>(num->value()) < Rational(0)) {
                        found = true;
                        return;
                    }
                } else if (std::holds_alternative<lmmc_real_t>(num->value())) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                    if (std::isfinite(v) && v < 0.0) {
                        found = true;
                        return;
                    }
                }
            }
        }
        // Recurse into arguments
        for (auto& arg : node.arguments()) {
            if (found) return;
            arg->accept(*this);
        }
    }

    void visit(const MatrixNode& node) override {
        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
            for (auto& item : std::get<MatrixNode::DenseStorage>(node.storage())) {
                if (found) return;
                if (item) item->accept(*this);
            }
        } else {
            for (auto& [idx, item] : std::get<MatrixNode::SparseStorage>(node.storage())) {
                (void)idx;
                if (found) return;
                if (item) item->accept(*this);
            }
        }
    }

    void visit(const RelationalNode& node) override {
        node.left()->accept(*this);
        if (!found) node.right()->accept(*this);
    }

    void visit(const LogicalNode& node) override {
        node.left()->accept(*this);
        if (!found && node.right()) node.right()->accept(*this);
    }

    void visit(const PiecewiseNode& node) override {
        for (const auto& branch : node.branches()) {
            if (found) return;
            branch.expression->accept(*this);
            if (!found) branch.condition->accept(*this);
        }
        if (!found && node.default_expr()) node.default_expr()->accept(*this);
    }

    void visit(const SummationNode& node) override {
        node.body()->accept(*this);
        if (!found) node.lower_bound()->accept(*this);
        if (!found) node.upper_bound()->accept(*this);
    }

    void visit(const ProductNode& node) override {
        node.body()->accept(*this);
        if (!found) node.lower_bound()->accept(*this);
        if (!found) node.upper_bound()->accept(*this);
    }

    void visit(const TransformNode& node) override {
        node.body()->accept(*this);
    }

    void visit(const QuantifierNode& node) override {
        node.domain()->accept(*this);
        if (!found) node.predicate()->accept(*this);
    }

    void visit(const SetBuilderNode& node) override {
        node.domain()->accept(*this);
        if (!found) node.predicate()->accept(*this);
    }

    void visit(const ComplexNode&) override {
        found = true;
    }
    void visit(const FiniteSetNode& node) override {
        for (const auto& element : node.elements()) {
            if (found) return;
            element->accept(*this);
        }
    }
    void visit(const IntervalNode& node) override {
        node.lower()->accept(*this);
        if (!found) node.upper()->accept(*this);
    }
    void visit(const MembershipNode& node) override {
        node.element()->accept(*this);
        if (!found) node.set()->accept(*this);
    }
    void visit(const QuantityNode& node) override { node.value()->accept(*this); }
};

/// Check if a solution expression contains imaginary components (sqrt of negative).
static bool contains_imaginary(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return false;
    ContainsImaginaryVisitor visitor;
    lamina::detail::node(expr)->accept(visitor);
    return visitor.found;
}

/// Try to extract a numeric double value from a solution expression.
/// Returns true if the expression is a pure numeric value, and sets out_value.
static bool try_get_numeric_value(const std::shared_ptr<SymbolicExpr>& expr, double& out_value) {
    if (!expr || !lamina::detail::node(expr)) return false;
    auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr));
    if (!num) return false;

    if (std::holds_alternative<BigInt>(num->value())) {
        out_value = std::get<BigInt>(num->value()).to_double();
        return true;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        out_value = std::get<Rational>(num->value()).to_double();
        return true;
    }
    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        out_value = std::get<lmmc_real_t>(num->value());
        return true;
    }
    return false;
}

/// Check if a numeric value is an integer (within tolerance).
static bool is_integer_value(double v) {
    if (!std::isfinite(v)) return false;
    int eq;
    lmmc_double_nearly_equal_tol(v, std::round(v), 1e-12, 1e-12, &eq);
    return eq != 0;
}

} // anonymous namespace

static std::vector<std::shared_ptr<SymbolicExpr>> solve_with_assumptions_impl(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* ctx,
    ComputationContext& context)
{
    // Step 1: Compute all solutions using the existing solver
    auto all_solutions = detail::propagate_result(
        solve_finite_checked(
            std::const_pointer_cast<SymbolicExpr>(equation), variable, context));

    // Step 2: If no context provided, return all solutions unfiltered
    if (!ctx) {
        return all_solutions;
    }

    // Step 3: Get the domain for the variable from the context
    Domain domain = ctx->get_domain(variable);

    // Step 4: Apply domain filtering as post-processing
    std::vector<std::shared_ptr<SymbolicExpr>> filtered;
    filtered.reserve(all_solutions.size());

    // Check if we have any sign constraints to apply
    bool has_sign_constraint = ctx->has_sign(variable, Sign::NonNegative) ||
                               ctx->has_sign(variable, Sign::Positive) ||
                               ctx->has_sign(variable, Sign::Negative) ||
                               ctx->has_sign(variable, Sign::NonPositive);

    // If domain is Complex (default/no restriction) and no sign constraints, return all
    if (domain == Domain::Complex && !has_sign_constraint) {
        return all_solutions;
    }

    for (const auto& sol : all_solutions) {
        if (!sol) continue;

        bool exclude = false;

        // --- Real domain filter ---
        // Exclude solutions containing imaginary components
        if (domain == Domain::Real || domain == Domain::Algebraic ||
            domain == Domain::Rational ||
            domain == Domain::Integer || domain == Domain::Natural ||
            domain == Domain::PositiveInt) {
            if (contains_imaginary(sol)) {
                exclude = true;
            }
        }

        if (!exclude) {
            double numeric_val = 0.0;
            bool is_numeric = try_get_numeric_value(sol, numeric_val);

            // --- Integer domain filter ---
            // Exclude non-integer numeric solutions
            if (domain == Domain::Integer || domain == Domain::Natural ||
                domain == Domain::PositiveInt) {
                if (is_numeric && !is_integer_value(numeric_val)) {
                    exclude = true;
                }
            }

            // --- Natural domain filter (non-negative integers) ---
            // Exclude negative solutions
            if (!exclude && (domain == Domain::Natural || domain == Domain::PositiveInt)) {
                if (is_numeric && numeric_val < 0.0) {
                    exclude = true;
                }
            }

            // --- PositiveInt domain filter ---
            // Exclude zero and negative solutions
            if (!exclude && domain == Domain::PositiveInt) {
                if (is_numeric) {
                    int eq;
                    lmmc_double_nearly_equal(numeric_val, 0.0, &eq);
                    if (eq || numeric_val <= 0.0) {
                        exclude = true;
                    }
                }
            }

            // --- NonNegative sign check (for Real domain with NonNegative sign) ---
            // Check if the variable has NonNegative sign assumption
            if (!exclude && is_numeric) {
                if (ctx->has_sign(variable, Sign::NonNegative) && numeric_val < 0.0) {
                    exclude = true;
                }
                if (ctx->has_sign(variable, Sign::Positive)) {
                    int eq;
                    lmmc_double_nearly_equal(numeric_val, 0.0, &eq);
                    if (eq || numeric_val <= 0.0) {
                        exclude = true;
                    }
                }
                if (ctx->has_sign(variable, Sign::Negative) && numeric_val >= 0.0) {
                    exclude = true;
                }
                if (ctx->has_sign(variable, Sign::NonPositive) && numeric_val > 0.0) {
                    exclude = true;
                }
            }
        }

        if (!exclude) {
            filtered.push_back(sol);
        }
    }

    return filtered;
}

AssumptionSolveResult solve_with_assumptions_checked(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* assumptions,
    ComputationContext& context)
{
    constexpr const char* operation = "solve_with_assumptions";
    if (!equation || variable.empty()) {
        return AssumptionSolveResult::failure(
            CasErrc::InvalidArgument,
            "assumption-aware solve requires an equation and variable",
            operation);
    }
    auto budget = context.consume_steps(1, operation);
    if (!budget) return AssumptionSolveResult::failure(budget.error());
    try {
        return AssumptionSolveResult::success(
            solve_with_assumptions_impl(
                equation, variable, assumptions, context));
    } catch (const detail::ResultPropagation& propagation) {
        return AssumptionSolveResult::failure(propagation.error());
    } catch (const std::bad_alloc&) {
        return AssumptionSolveResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving with assumptions",
            operation);
    } catch (const std::exception& ex) {
        return AssumptionSolveResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

AssumptionSolveResult solve_with_assumptions_checked(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* assumptions)
{
    ComputationContext context;
    return solve_with_assumptions_checked(
        equation, variable, assumptions, context);
}


} // namespace lamina
