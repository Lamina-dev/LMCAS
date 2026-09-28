#include "internal/inference_engine_impl.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/assumption_simplification.hpp"
#include "internal/expression_construction.hpp"
#include "internal/rewrite_budget.hpp"
#include <limits>

namespace LMCAS {
namespace {
constexpr const char* periodic_operation = "inference.periodic";
using Node = std::shared_ptr<const SymbolicNode>;

struct PeriodEvidence {
    Tribool periodic = Tribool::Unknown;
    std::optional<SymbolicExpr> period;
    bool real_valued = false;
};

class PeriodInference {
public:
    PeriodInference(const AssumptionContext& assumptions, const InferenceEngine& engine,
                    const std::string& variable, ComputationContext& context)
        : assumptions_(assumptions), facts_(engine), variable_(variable), context_(context),
          max_depth_(engine.get_max_depth()) {}

    Result<PeriodEvidence> infer(const Node& node) {
        auto entered = context_.enter_recursion(periodic_operation);
        if (!entered) { return Result<PeriodEvidence>::failure(entered.error()); }
        struct Leave {
            ComputationContext& context;
            ~Leave() { context.leave_recursion(); }
        } leave{context_};
        if (depth_ >= max_depth_) return PeriodEvidence{};
        ++depth_;
        struct LeaveProof {
            int& depth;
            ~LeaveProof() { --depth; }
        } leave_proof{depth_};

        if (const auto* symbol = dynamic_cast<const VariableNode*>(node.get())) {
            if (!symbol->is_constant() &&
                assumptions_.has_period_declarations(symbol->name())) {
                return declared_period(node, *symbol);
            }
        }
        auto constant_real = real_defined(node);
        if (!constant_real) { return Result<PeriodEvidence>::failure(constant_real.error()); }
        if (constant_real.value() == Tribool::True && independent(node))
            { return PeriodEvidence{Tribool::True, std::nullopt, true}; }

        const auto* function = dynamic_cast<const FunctionNode*>(node.get());
        if (!function) { return PeriodEvidence{}; }
        if (function->type() == FunctionNode::FuncType::Infinity)
            { return PeriodEvidence{Tribool::False, std::nullopt, false}; }
        if (function->arguments().size() != 1) { return PeriodEvidence{}; }
        const Node& argument = function->arguments().front();
        if (function->type() == FunctionNode::FuncType::ArcTan) {
            return arctangent_period(argument);
        }
        const auto type = function->type();
        if (type != FunctionNode::FuncType::Sin && type != FunctionNode::FuncType::Cos &&
            type != FunctionNode::FuncType::Tan) { return PeriodEvidence{}; }
        return affine_period(node, argument, type);
    }

private:
    Result<PeriodEvidence> declared_period(const Node& node, const VariableNode& symbol) {
        auto period = assumptions_.get_period(symbol.name(), variable_);
        if (!period) { return PeriodEvidence{}; }
        auto defined = detail::query_definedness(
            detail::node(*period), facts_, Domain::Real, context_);
        if (!defined) { return Result<PeriodEvidence>::failure(defined.error()); }
        auto positive = detail::query_positive_value(
            detail::node(*period), facts_, context_);
        if (!positive) { return Result<PeriodEvidence>::failure(positive.error()); }
        if (defined.value() != Tribool::True || positive.value() != Tribool::True ||
            !independent(detail::node(*period)))
            { return PeriodEvidence{}; }
        auto real = detail::query_real_value(node, facts_, context_);
        if (!real) { return Result<PeriodEvidence>::failure(real.error()); }
        return PeriodEvidence{
            Tribool::True, *period, real.value() == Tribool::True};
    }

    Result<PeriodEvidence> arctangent_period(const Node& argument) {
        const auto* symbol = dynamic_cast<const VariableNode*>(argument.get());
        if (symbol && !symbol->is_constant() && symbol->name() == variable_ &&
            !assumptions_.has_period_declarations(symbol->name())) {
            return PeriodEvidence{Tribool::False, std::nullopt, true};
        }
        auto inner = infer(argument);
        if (!inner) { return inner; }
        if (inner.value().periodic == Tribool::True && inner.value().real_valued)
            { return inner; }
        return PeriodEvidence{};
    }

    Result<PeriodEvidence> affine_period(const Node& node, const Node& argument, FunctionNode::FuncType type) {
        auto defined = detail::query_definedness(argument, facts_, Domain::Real, context_);
        if (!defined) { return Result<PeriodEvidence>::failure(defined.error()); }
        if (defined.value() != Tribool::True) { return PeriodEvidence{}; }
        auto affine = detail::recognize_affine(detail::expression_from_node(argument), variable_, context_);
        if (!affine) { return Result<PeriodEvidence>::failure(affine.error()); }
        if (!affine.value()) { return PeriodEvidence{}; }
        for (auto* coefficient : {&affine.value()->slope, &affine.value()->offset}) {
            auto normalized = detail::simplify_expression(*coefficient, context_);
            if (!normalized) { return Result<PeriodEvidence>::failure(normalized.error()); }
            *coefficient = std::move(normalized.value());
        }
        const auto& slope = affine.value()->slope;
        const auto& offset = affine.value()->offset;
        for (const auto& coefficient : {slope, offset}) {
            auto real = real_defined(detail::node(*coefficient));
            if (!real) { return Result<PeriodEvidence>::failure(real.error()); }
            if (real.value() != Tribool::True || !independent(detail::node(*coefficient)))
                { return PeriodEvidence{}; }
        }
        auto nonzero = detail::query_nonzero_value(detail::node(*slope), facts_, Domain::Real, context_);
        if (!nonzero) { return Result<PeriodEvidence>::failure(nonzero.error()); }
        return period_from_slope(node, type, slope, nonzero.value());
    }

    Result<PeriodEvidence> period_from_slope(const Node& node, FunctionNode::FuncType type, const std::shared_ptr<SymbolicExpr>& slope, Tribool nonzero) {
        if (nonzero == Tribool::False) {
            auto constant_domain = detail::query_definedness(node, facts_, Domain::Real, context_);
            if (!constant_domain) { return Result<PeriodEvidence>::failure(constant_domain.error()); }
            if (constant_domain.value() == Tribool::True)
                { return PeriodEvidence{Tribool::True, std::nullopt, true}; }
            return PeriodEvidence{};
        }
        if (nonzero != Tribool::True) { return PeriodEvidence{}; }
        auto positive = detail::query_positive_value(detail::node(*slope), facts_, context_);
        if (!positive) { return Result<PeriodEvidence>::failure(positive.error()); }
        std::shared_ptr<SymbolicExpr> magnitude;
        if (positive.value() == Tribool::True) {
            magnitude = slope;
        } else if (positive.value() == Tribool::False) {
            magnitude = SymbolicExpr::multiply(SymbolicExpr::number(-1), slope);
        } else {
            magnitude = detail::make_expression_ptr(detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs, std::vector<Node>{detail::node(*slope)}));
        }
        auto exact_pi = detail::constant_expression("pi");
        if (!exact_pi) { return Result<PeriodEvidence>::failure(exact_pi.error()); }
        auto numerator = type == FunctionNode::FuncType::Tan ? exact_pi.value() :
            SymbolicExpr::multiply(SymbolicExpr::number(2), exact_pi.value());
        auto period = detail::simplify_expression(SymbolicExpr::divide(numerator, magnitude), context_);
        if (!period) { return Result<PeriodEvidence>::failure(period.error()); }
        return PeriodEvidence{Tribool::True, *period.value(), true};
    }

    bool independent(const Node& node) const {
        detail::RewriteBudget budget(context_, context_.limits().max_recursion_depth,
            std::numeric_limits<std::size_t>::max(), periodic_operation);
        const auto names = free_variables(node, &budget);
        for (const auto& name : names) {
            auto step = context_.consume_steps(1, periodic_operation);
            if (!step) throw step.error();
            if (name == variable_ || assumptions_.has_period_declarations(name)) return false;
        }
        return true;
    }

    Result<Tribool> real_defined(const Node& node) const {
        auto defined = detail::query_definedness(node, facts_, Domain::Real, context_);
        if (!defined || defined.value() != Tribool::True) return defined;
        return detail::query_real_value(node, facts_, context_);
    }

    const AssumptionContext& assumptions_;
    detail::AssumptionFacts facts_;
    const std::string& variable_;
    ComputationContext& context_;
    const int max_depth_;
    int depth_ = 0;
};

Result<PeriodEvidence> periodicity_checked(
    const SymbolicExpr& expression, const std::string& variable,
    const AssumptionContext& assumptions, const InferenceEngine& engine,
    ComputationContext& context) {
    if (variable.empty()) {
        return Result<PeriodEvidence>::failure(CasErrc::InvalidArgument,
            "independent variable must not be empty", periodic_operation);
    }
    return checked_inference_result<PeriodEvidence>(expression, periodic_operation, context,
        [&]() -> Result<PeriodEvidence> {
            PeriodInference inference(assumptions, engine, variable, context);
            return inference.infer(detail::node(expression));
        });
}
}

InferenceTriboolResult InferenceEngine::query_periodic_checked(
    const SymbolicExpr& expression, const std::string& variable) const {
    ComputationContext context;
    return query_periodic_checked(expression, variable, context);
}

InferenceTriboolResult InferenceEngine::query_periodic_checked(
    const SymbolicExpr& expression, const std::string& variable, ComputationContext& context) const {
    auto proof = periodicity_checked(expression, variable, impl_->ctx, *this, context);
    if (!proof) return InferenceTriboolResult::failure(proof.error());
    return proof.value().periodic;
}

InferencePeriodResult InferenceEngine::infer_period_checked(
    const SymbolicExpr& expression, const std::string& variable) const {
    ComputationContext context;
    return infer_period_checked(expression, variable, context);
}

InferencePeriodResult InferenceEngine::infer_period_checked(
    const SymbolicExpr& expression, const std::string& variable, ComputationContext& context) const {
    auto proof = periodicity_checked(expression, variable, impl_->ctx, *this, context);
    if (!proof) return InferencePeriodResult::failure(proof.error());
    return std::move(proof.value().period);
}

}
