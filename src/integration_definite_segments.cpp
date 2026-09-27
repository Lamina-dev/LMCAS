#include "internal/integration_definite_support.hpp"
#include "limit_result.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "residual_verification.hpp"
#include "root_of_identity.hpp"

#include <unordered_map>

namespace LMCAS::detail::definite {
namespace {
class SegmentFacts final : public FactsQuery {
public:
    SegmentFacts(const AssumptionContext& assumptions, const std::string& variable,
        const Boundary& lower, const Boundary& upper)
        : base_(assumptions), var_(variable), lower_(lower), upper_(upper) {}

    mutable std::optional<CasError> error;
    Result<std::optional<int>> value_sign(const Node& node, ComputationContext& context) const {
        auto result = value_sign_impl(node, context);
        if (!result) error = result.error();
        return result;
    }

    Result<std::optional<int>> value_sign_impl(const Node& node, ComputationContext& context) const {
        auto access = context.consume_steps(0, "facts.segment");
        if (!access) { return Result<std::optional<int>>::failure(access.error()); }
        const auto found = signs_.find(node);
        if (found != signs_.end()) return found->second;
        auto value = rational_checked(node, var_, context);
        if (!value) return Result<std::optional<int>>::failure(value.error());
        if (!value.value()) return std::optional<int>{};
        auto denominator = polynomial_sign_checked(value.value()->denominator, lower_, upper_, context);
        if (!denominator) return denominator;
        if (!denominator.value() || *denominator.value() == 0) return std::optional<int>{};
        auto numerator = polynomial_sign_checked(value.value()->numerator, lower_, upper_, context);
        if (!numerator) return numerator;
        std::optional<int> result;
        if (numerator.value()) result = *numerator.value() * *denominator.value();
        signs_.emplace(node, result);
        return result;
    }

    Result<Tribool> is_positive(const Node& node, ComputationContext& context) const override { return property(node, 0, context); }
    Result<Tribool> is_negative(const Node& node, ComputationContext& context) const override { return property(node, 1, context); }
    Result<Tribool> is_nonnegative(const Node& node, ComputationContext& context) const override { return property(node, 2, context); }
    Result<Tribool> is_nonzero(const Node& node, ComputationContext& context) const override { return property(node, 3, context); }
    Result<Tribool> is_real(const Node& node, ComputationContext& context) const override {
        auto value = value_sign(node, context);
        if (!value) return Result<Tribool>::failure(value.error());
        if (value.value()) return Tribool::True;
        return remember(base_.is_real(node, context));
    }
private:
    Result<Tribool> remember(Result<Tribool> result) const {
        if (!result) error = result.error();
        return result;
    }
    Result<Tribool> property(const Node& node, int property, ComputationContext& context) const {
        auto value = value_sign(node, context);
        if (!value) { return Result<Tribool>::failure(value.error()); }
        if (value.value()) {
            int s = *value.value();
            bool truth = property == 0 ? s > 0 : property == 1 ? s < 0 : property == 2 ? s >= 0 : s != 0;
            return truth ? Tribool::True : Tribool::False;
        }
        switch (property) {
        case 0: return remember(base_.is_positive(node, context));
        case 1: return remember(base_.is_negative(node, context));
        case 2: return remember(base_.is_nonnegative(node, context));
        default: return remember(base_.is_nonzero(node, context));
        }
    }
    detail::AssumptionFacts base_;
    const std::string& var_;
    const Boundary& lower_;
    const Boundary& upper_;
    mutable std::unordered_map<Node, std::optional<int>> signs_;
};

static Result<Tribool> predicate_truth(const Node& node, const SegmentFacts& facts, ComputationContext& context);

static Tribool binary_logical_truth(LogicalNode::Op op, Tribool a, Tribool b) {
    if (op == LogicalNode::Op::And) {
        if (a == Tribool::False || b == Tribool::False) { return Tribool::False; }
        return a == Tribool::True && b == Tribool::True ? Tribool::True : Tribool::Unknown;
    }
    if (a == Tribool::True || b == Tribool::True) { return Tribool::True; }
    return a == Tribool::False && b == Tribool::False ? Tribool::False : Tribool::Unknown;
}

static Result<Tribool> logical_predicate_truth(const LogicalNode& logical,
    const SegmentFacts& facts, ComputationContext& context) {
    auto left = predicate_truth(logical.left(), facts, context);
    if (!left) { return left; }
    const auto negate = [](Tribool truth) {
        return truth == Tribool::True ? Tribool::False :
            truth == Tribool::False ? Tribool::True : Tribool::Unknown;
    };
    if (logical.op() == LogicalNode::Op::Not) { return negate(left.value()); }
    auto right = predicate_truth(logical.right(), facts, context);
    if (!right) { return right; }
    auto a = left.value(), b = right.value();
    if (logical.op() == LogicalNode::Op::Implies) { a = negate(a); }
    return binary_logical_truth(logical.op(), a, b);
}

static Result<Tribool> relational_predicate_truth(const RelationalNode& relation,
    const SegmentFacts& facts, ComputationContext& context) {
    auto difference = sym_sub(detail::expression_from_node(relation.left()),
        detail::expression_from_node(relation.right()));
    auto value = facts.value_sign(detail::node(difference), context);
    if (!value) { return Result<Tribool>::failure(value.error()); }
    if (!value.value()) { return Tribool::Unknown; }
    const int s = *value.value();
    bool truth;
    switch (relation.op()) {
    case RelationOp::EQ: { truth = s == 0; break; }
    case RelationOp::NEQ: { truth = s != 0; break; }
    case RelationOp::LT: { truth = s < 0; break; }
    case RelationOp::LEQ: { truth = s <= 0; break; }
    case RelationOp::GT: { truth = s > 0; break; }
    case RelationOp::GEQ: { truth = s >= 0; break; }
    default: { return Tribool::Unknown; }
    }
    return truth ? Tribool::True : Tribool::False;
}

static Result<Tribool> predicate_truth(const Node& node, const SegmentFacts& facts, ComputationContext& context) {
    if (const auto* logical = dynamic_cast<const LogicalNode*>(node.get())) {
        return logical_predicate_truth(*logical, facts, context);
    }
    if (const auto* relation = dynamic_cast<const RelationalNode*>(node.get())) {
        return relational_predicate_truth(*relation, facts, context);
    }
    return Tribool::Unknown;
}

/**
 * @brief 在已认证开区间上进行恒等改写，保留对数分支。
 * 对数分支变更仅由独立的原函数候选生成器处理。
 */
class SegmentExpression final : public detail::SymbolicRewriter {
public:
    SegmentExpression(const SegmentFacts& facts, ComputationContext& context) : facts_(facts), context_(context) {}
    std::optional<CasError> error;
    void visit(const PiecewiseNode& piecewise) override {
        for (const auto& branch : piecewise.branches()) {
            auto truth = predicate_truth(branch.condition, facts_, context_);
            if (!truth) { error = truth.error(); set_result(current()); return; }
            if (truth.value() == Tribool::Unknown) {
                error = CasError{CasErrc::Inconclusive, "piecewise branch is undecided", operation};
                set_result(current());
                return;
            }
            if (truth.value() == Tribool::True) {
                set_result(rewrite(branch.expression));
                return;
            }
        }
        if (piecewise.default_expr()) set_result(rewrite(piecewise.default_expr()));
        else {
            error = CasError{CasErrc::DomainError, "piecewise function has no value on the segment",
                "integrate.definite.convergence"};
            set_result(current());
        }
    }
    void visit(const PowerNode& node) override {
        bool changed = false;
        auto base = rewrite_child(node.base(), changed);
        auto exponent = rewrite_child(node.exponent(), changed);
        auto rewritten = changed ? detail::make_node<PowerNode>(base, exponent) : current();
        if (!specialize_absolute(absolute_argument(*rewritten))) set_result(rewritten);
    }
    void visit(const FunctionNode& node) override {
        bool changed = false;
        auto arguments = rewrite_children(node.arguments(), changed);
        auto rewritten = changed ? detail::make_node<FunctionNode>(node.type(), std::move(arguments)) : current();
        if (!specialize_absolute(absolute_argument(*rewritten))) set_result(rewritten);
    }
private:
    bool specialize_absolute(const Node& argument) {
        if (!argument) return false;
        auto real = detail::query_real_value(argument, facts_, context_);
        if (!real) { error = real.error(); set_result(current()); return true; }
        auto value = facts_.value_sign(argument, context_);
        if (!value) { error = value.error(); set_result(current()); return true; }
        if (!value.value()) {
            auto nonnegative = detail::query_nonnegative_value(argument, facts_, context_);
            if (!nonnegative) { error = nonnegative.error(); set_result(current()); return true; }
            if (nonnegative.value() == Tribool::True) value.value() = 1;
            else if (nonnegative.value() == Tribool::False) value.value() = -1;
        }
        if (real.value() != Tribool::True || !value.value()) return false;
        set_result(*value.value() < 0
            ? detail::node(SymbolicExpr::multiply(SymbolicExpr::number(-1),
                detail::make_expression_ptr(argument))) : argument);
        return true;
    }
    const SegmentFacts& facts_;
    ComputationContext& context_;
};

}

bool needs_segment_expression(const Node& node) {
    class Detector final : public detail::RecursiveSymbolicVisitor {
    public:
        bool found = false;
        void visit(const PiecewiseNode&) override { found = true; }
        void visit(const FunctionNode& function) override {
            if (absolute_argument(function)) found = true;
            else detail::RecursiveSymbolicVisitor::visit(function);
        }
        void visit(const PowerNode& power) override {
            if (absolute_argument(power)) found = true;
            else detail::RecursiveSymbolicVisitor::visit(power);
        }
    } detector;
    node->accept(detector);
    return detector.found;
}

static Result<void> collect_predicate_roots(const Node& condition, const std::string& var,
    std::vector<Algebraic>& roots, ComputationContext& context) {
    auto step = context.consume_steps(1, operation);
    if (!step) { return step; }
    if (const auto* logical = dynamic_cast<const LogicalNode*>(condition.get())) {
        auto left = collect_predicate_roots(logical->left(), var, roots, context);
        if (!left) { return left; }
        if (logical->right()) { return collect_predicate_roots(logical->right(), var, roots, context); }
        return Result<void>::success();
    }
    const auto* relation = dynamic_cast<const RelationalNode*>(condition.get());
    if (!relation) { return undecided<void>(); }
    auto difference = sym_sub(detail::expression_from_node(relation->left()),
        detail::expression_from_node(relation->right()));
    auto rational = rational_checked(detail::node(difference), var, context);
    if (!rational) { return Result<void>::failure(rational.error()); }
    if (!rational.value()) { return undecided<void>(); }
    for (const auto* polynomial : {&rational.value()->numerator, &rational.value()->denominator}) {
        if (polynomial->degree() <= 0) { continue; }
        auto isolated = roots_checked(*polynomial, context);
        if (!isolated) { return Result<void>::failure(isolated.error()); }
        roots.insert(roots.end(), isolated.value().begin(), isolated.value().end());
    }
    return Result<void>::success();
}

Result<AssumptionContext> real_assumptions(const std::string& var, ComputationContext& context) {
    AssumptionContext assumptions = context.assumptions() ? *context.assumptions() : AssumptionContext{};
    assumptions.push();
    auto real = assumptions.assume_domain_checked(var, Domain::Real);
    if (!real) return Result<AssumptionContext>::failure(real.error());
    return assumptions;
}

static Result<void> original_conditions_checked(const Node& node, const FactsQuery& facts,
    std::vector<std::shared_ptr<SymbolicExpr>>& conditions, ComputationContext& context) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (const auto* piecewise = dynamic_cast<const PiecewiseNode*>(node.get())) {
        for (const auto& branch : piecewise->branches()) {
            conditions.push_back(detail::make_expression_ptr(branch.condition));
            auto child = original_conditions_checked(branch.expression, facts, conditions, context);
            if (!child) return child;
        }
        if (piecewise->default_expr())
            return original_conditions_checked(piecewise->default_expr(), facts, conditions, context);
        return Result<void>::success();
    }
    auto domain = detail::domain_constraints(node, facts, Domain::Real, context);
    if (!domain) return Result<void>::failure(domain.error());
    if (!domain.value()) return undecided<void>();
    conditions.insert(conditions.end(), domain.value()->begin(), domain.value()->end());
    class AbsoluteBoundaries final : public detail::RecursiveSymbolicVisitor {
    public:
        AbsoluteBoundaries(std::vector<std::shared_ptr<SymbolicExpr>>& out, const FactsQuery& facts,
                           ComputationContext& context)
            : out_(out), facts_(facts), context_(context) {}
        std::optional<CasError> error;
        void visit(const FunctionNode& function) override {
            collect(function);
            detail::RecursiveSymbolicVisitor::visit(function);
        }
        void visit(const PowerNode& power) override {
            collect(power);
            detail::RecursiveSymbolicVisitor::visit(power);
        }
    private:
        void collect(const SymbolicNode& node) {
            if (auto argument = absolute_argument(node)) {
                auto nonnegative = detail::query_nonnegative_value(argument, facts_, context_);
                if (!nonnegative) { error = nonnegative.error(); return; }
                auto nonzero = detail::query_nonzero_value(argument, facts_, Domain::Real, context_);
                if (!nonzero) { error = nonzero.error(); return; }
                if (nonnegative.value() == Tribool::Unknown && nonzero.value() == Tribool::Unknown)
                    out_.push_back(detail::make_expression_ptr(detail::make_node<RelationalNode>(
                        argument, detail::node(SymbolicExpr::number(0)), RelationOp::EQ)));
            }
        }
        std::vector<std::shared_ptr<SymbolicExpr>>& out_;
        const FactsQuery& facts_;
        ComputationContext& context_;
    } boundaries(conditions, facts, context);
    node->accept(boundaries);
    if (boundaries.error) return Result<void>::failure(*boundaries.error);
    return Result<void>::success();
}

Result<SymbolicExpr> segment_expression_checked(const SymbolicExpr& expression,
    const std::string& var, const SymbolicExpr& lower, const SymbolicExpr& upper,
    ComputationContext& context) {
    if (!needs_segment_expression(detail::node(expression))) return expression;
    auto a = boundary_checked(lower, context);
    auto b = boundary_checked(upper, context);
    if (!a) return Result<SymbolicExpr>::failure(a.error());
    if (!b) return Result<SymbolicExpr>::failure(b.error());
    auto assumptions = real_assumptions(var, context);
    if (!assumptions) return Result<SymbolicExpr>::failure(assumptions.error());
    SegmentFacts facts(assumptions.value(), var, a.value(), b.value());
    SegmentExpression rewrite(facts, context);
    auto selected = rewrite.rewrite(detail::node(expression));
    if (rewrite.error) return Result<SymbolicExpr>::failure(*rewrite.error);
    return detail::expression_from_node(selected);
}

static Result<std::vector<Algebraic>> ordered_partition_roots(
    std::vector<Algebraic>& roots, const Boundary& a, const Boundary& b,
    ComputationContext& context) {
    std::vector<Algebraic> ordered;
    for (auto& root : roots) {
        auto left = compare_boundaries(a, Boundary{0, root}, context);
        auto right = compare_boundaries(Boundary{0, root}, b, context);
        if (!left) { return Result<std::vector<Algebraic>>::failure(left.error()); }
        if (!right) { return Result<std::vector<Algebraic>>::failure(right.error()); }
        if (left.value() >= 0 || right.value() >= 0) { continue; }
        auto position = ordered.begin();
        bool duplicate = false;
        for (; position != ordered.end(); ++position) {
            auto comparison = compare_boundaries(Boundary{0, root}, Boundary{0, *position}, context);
            if (!comparison) { return Result<std::vector<Algebraic>>::failure(comparison.error()); }
            if (comparison.value() == 0) { duplicate = true; break; }
            if (comparison.value() < 0) { break; }
        }
        if (!duplicate) { ordered.insert(position, std::move(root)); }
    }
    return ordered;
}

Result<std::vector<SymbolicExpr>> definite_partition_checked(const SymbolicExpr& original,
    const std::string& var, const SymbolicExpr& lower, const SymbolicExpr& upper,
    ComputationContext& context) {
    auto assumptions = real_assumptions(var, context);
    if (!assumptions) { return Result<std::vector<SymbolicExpr>>::failure(assumptions.error()); }
    detail::AssumptionFacts facts(assumptions.value());
    std::vector<std::shared_ptr<SymbolicExpr>> conditions;
    auto captured = original_conditions_checked(detail::node(original), facts, conditions, context);
    if (!captured) { return Result<std::vector<SymbolicExpr>>::failure(captured.error()); }
    if (conditions.empty()) { return std::vector<SymbolicExpr>{lower, upper}; }
    auto a = boundary_checked(lower, context);
    auto b = boundary_checked(upper, context);
    if (!a) { return Result<std::vector<SymbolicExpr>>::failure(a.error()); }
    if (!b) { return Result<std::vector<SymbolicExpr>>::failure(b.error()); }
    std::vector<Algebraic> roots;
    for (const auto& condition : conditions) {
        auto collected = collect_predicate_roots(detail::node(condition), var, roots, context);
        if (!collected) { return Result<std::vector<SymbolicExpr>>::failure(collected.error()); }
    }
    auto ordered = ordered_partition_roots(roots, a.value(), b.value(), context);
    if (!ordered) { return Result<std::vector<SymbolicExpr>>::failure(ordered.error()); }
    std::vector<SymbolicExpr> result{lower};
    for (const auto& root : ordered.value()) {
        if (root.is_rational()) { result.push_back(*SymbolicExpr::number(root.lower)); }
        else {
            auto expression = make_rootof_checked(poly_to_symbolic(root.polynomial), var, root.root_index, context);
            if (!expression) { return Result<std::vector<SymbolicExpr>>::failure(expression.error()); }
            result.push_back(*expression.value());
        }
    }
    result.push_back(upper);
    return result;
}

static Result<Tribool> defined_on_interval(const Node& node, const FactsQuery& facts,
    ComputationContext& context) {
    auto defined = detail::query_definedness(node, facts, Domain::Real, context);
    if (!defined || defined.value() != Tribool::Unknown) { return defined; }
    const auto* interval = dynamic_cast<const SegmentFacts*>(&facts);
    if (!interval) { return defined; }
    auto conditions = detail::domain_constraints(node, facts, Domain::Real, context);
    if (!conditions) { return Result<Tribool>::failure(conditions.error()); }
    if (!conditions.value()) { return Tribool::Unknown; }
    bool unknown = false;
    for (const auto& condition : *conditions.value()) {
        auto truth = predicate_truth(detail::node(condition), *interval, context);
        if (!truth) { return truth; }
        if (truth.value() == Tribool::False) { return Tribool::False; }
        unknown |= truth.value() == Tribool::Unknown;
    }
    return unknown ? Tribool::Unknown : Tribool::True;
}

static bool continuous_function_family(FunctionNode::FuncType type) {
    using F = FunctionNode::FuncType;
    switch (type) {
    case F::Sin: { return true; }
    case F::Cos: { return true; }
    case F::Exp: { return true; }
    case F::Abs: { return true; }
    case F::ArcTan: { return true; }
    case F::Ln: { return true; }
    case F::Log: { return true; }
    case F::Sqrt: { return true; }
    case F::Sinh: { return true; }
    case F::Cosh: { return true; }
    case F::Tan: { return true; }
    case F::ArcSin: { return true; }
    case F::ArcCos: { return true; }
    default: { return false; }
    }
}

static Result<Tribool> structurally_continuous(const Node& node, const FactsQuery& facts,
    ComputationContext& context);

static Result<Tribool> continuous_children(const std::vector<Node>& children,
    const FactsQuery& facts, ComputationContext& context) {
    for (const auto& child : children) {
        auto continuous = structurally_continuous(child, facts, context);
        if (!continuous || continuous.value() != Tribool::True) { return continuous; }
    }
    return Tribool::True;
}

static Result<Tribool> structurally_continuous(const Node& node, const FactsQuery& facts,
    ComputationContext& context) {
    auto step = context.consume_steps(1, operation);
    if (!step) { return Result<Tribool>::failure(step.error()); }
    auto defined = defined_on_interval(node, facts, context);
    if (!defined || defined.value() != Tribool::True) { return defined; }
    if (dynamic_cast<const NumberNode*>(node.get()) || dynamic_cast<const VariableNode*>(node.get()) ||
        dynamic_cast<const RootOfNode*>(node.get())) { return Tribool::True; }
    const auto* sum = dynamic_cast<const AddNode*>(node.get());
    const auto* product = dynamic_cast<const MultiplyNode*>(node.get());
    const auto* function = dynamic_cast<const FunctionNode*>(node.get());
    const auto* power = dynamic_cast<const PowerNode*>(node.get());
    if (power) {
        auto base = structurally_continuous(power->base(), facts, context);
        if (!base || base.value() != Tribool::True) { return base; }
        return structurally_continuous(power->exponent(), facts, context);
    }
    if (sum) { return continuous_children(sum->operands(), facts, context); }
    if (product) { return continuous_children(product->operands(), facts, context); }
    if (function) {
        if (!continuous_function_family(function->type())) { return Tribool::Unknown; }
        return continuous_children(function->arguments(), facts, context);
    }
    return Tribool::Unknown;
}

Result<Tribool> certify_continuity_checked(const SymbolicExpr& expression,
    const std::string& var, const SymbolicExpr& lower, const SymbolicExpr& upper,
    ComputationContext& context) {
    auto effective = segment_expression_checked(expression, var, lower, upper, context);
    if (!effective) { return is_undecided(effective.error()) ? Result<Tribool>(Tribool::Unknown) :
        Result<Tribool>::failure(effective.error()); }
    auto assumptions = real_assumptions(var, context);
    if (!assumptions) { return Result<Tribool>::failure(assumptions.error()); }
    detail::AssumptionFacts global(assumptions.value());
    auto everywhere = structurally_continuous(detail::node(effective.value()), global, context);
    if (!everywhere || everywhere.value() != Tribool::Unknown) { return everywhere; }
    auto a = boundary_checked(lower, context);
    auto b = boundary_checked(upper, context);
    if (!a) { return is_undecided(a.error()) ? Result<Tribool>(Tribool::Unknown) : Result<Tribool>::failure(a.error()); }
    if (!b) { return is_undecided(b.error()) ? Result<Tribool>(Tribool::Unknown) : Result<Tribool>::failure(b.error()); }
    SegmentFacts facts(assumptions.value(), var, a.value(), b.value());
    return structurally_continuous(detail::node(effective.value()), facts, context);
}

namespace {
/**
 * @brief 生成实原函数候选，与被积函数的恒等改写分离。
 * 证明 u 恒为正或恒为负后，才将 ln(abs(u)) 分别化为 ln(u) 或 ln(-u)。
 */
class RealPrimitiveCandidate final : public detail::SymbolicRewriter {
public:
    RealPrimitiveCandidate(const SegmentFacts& facts, ComputationContext& context) : facts_(facts), context_(context) {}
    std::optional<CasError> error;
    void visit(const FunctionNode& node) override {
        bool changed = false;
        auto arguments = rewrite_children(node.arguments(), changed);
        if (arguments.size() == 1 && (node.type() == FunctionNode::FuncType::Ln ||
            node.type() == FunctionNode::FuncType::Abs)) {
            auto nonzero = detail::query_nonzero_value(arguments[0], facts_, Domain::Real, context_);
            if (!nonzero) { error = nonzero.error(); set_result(current()); return; }
            auto real = detail::query_real_value(arguments[0], facts_, context_);
            if (!real) { error = real.error(); set_result(current()); return; }
            auto sign_value = facts_.value_sign(arguments[0], context_);
            if (!sign_value) { error = sign_value.error(); set_result(current()); return; }
            if (real.value() == Tribool::True && nonzero.value() == Tribool::True && sign_value.value()) {
                if (*sign_value.value() < 0) {
                    arguments[0] = detail::node(SymbolicExpr::multiply(SymbolicExpr::number(-1),
                        detail::make_expression_ptr(arguments[0])));
                    changed = true;
                }
                if (node.type() == FunctionNode::FuncType::Abs) {
                    set_result(arguments[0]);
                    return;
                }
            }
        }
        set_result(changed ? detail::make_node<FunctionNode>(node.type(), std::move(arguments)) : current());
    }
private:
    const SegmentFacts& facts_;
    ComputationContext& context_;
};

}

static Result<SymbolicExpr> boundary_value_checked(const SymbolicExpr& primitive,
    const std::string& var, const SymbolicExpr& point, LimitDirection direction,
    ComputationContext& context) {
    auto result = limit_checked(detail::make_expression_ptr(primitive), var, detail::make_expression_ptr(point), direction, context);
    if (!result) return Result<SymbolicExpr>::failure(result.error());
    const auto& outcome = result.value().value;
    if (const auto* finite = std::get_if<FiniteLimit>(&outcome)) return *finite->value;
    if (std::holds_alternative<PositiveInfinityLimit>(outcome)) return *SymbolicExpr::infinity(1);
    if (std::holds_alternative<NegativeInfinityLimit>(outcome)) return *SymbolicExpr::infinity(-1);
    return divergent<SymbolicExpr>();
}

static ExpressionResult certify_primitive_derivative(SymbolicExpr& candidate,
    const std::string& var, const SymbolicExpr& a, const SymbolicExpr& b,
    const Result<Boundary>& lower, const Result<Boundary>& upper,
    ComputationContext& context) {
    if (!lower && !is_undecided(lower.error())) {
        return ExpressionResult::failure(lower.error());
    }
    if (!upper && !is_undecided(upper.error())) {
        return ExpressionResult::failure(upper.error());
    }
    if (lower && upper) {
        auto assumptions = real_assumptions(var, context);
        if (!assumptions) { return ExpressionResult::failure(assumptions.error()); }
        SegmentFacts facts(assumptions.value(), var, lower.value(), upper.value());
        RealPrimitiveCandidate rewrite(facts, context);
        candidate = detail::expression_from_node(rewrite.rewrite(detail::node(candidate)));
        if (rewrite.error) { return ExpressionResult::failure(*rewrite.error); }
    }
    auto primitive_domain = certify_continuity_checked(candidate, var, a, b, context);
    if (!primitive_domain) { return ExpressionResult::failure(primitive_domain.error()); }
    if (primitive_domain.value() != Tribool::True) { return undecided<std::shared_ptr<SymbolicExpr>>(); }
    auto derivative = candidate.differentiate(var);
    if (!derivative) { return ExpressionResult::failure(CasErrc::InternalInvariant,
        "primitive differentiation returned no expression", operation); }
    auto derivative_domain = certify_continuity_checked(*derivative, var, a, b, context);
    if (!derivative_domain) { return ExpressionResult::failure(derivative_domain.error()); }
    if (derivative_domain.value() != Tribool::True) { return undecided<std::shared_ptr<SymbolicExpr>>(); }
    return derivative;
}

static Result<void> prove_rational_derivative(const SymbolicExpr& effective,
    const std::shared_ptr<SymbolicExpr>& derivative, const std::string& var,
    ComputationContext& context) {
    auto f = rational_checked(detail::node(effective), var, context);
    auto d = rational_checked(detail::node(derivative), var, context);
    if (!f) { return Result<void>::failure(f.error()); }
    if (!d) { return Result<void>::failure(d.error()); }
    if (!f.value() || !d.value()) { return undecided<void>(); }
    auto left = product_checked(f.value()->numerator, d.value()->denominator, context);
    auto right = product_checked(d.value()->numerator, f.value()->denominator, context);
    if (!left) { return Result<void>::failure(left.error()); }
    if (!right) { return Result<void>::failure(right.error()); }
    auto proof = check_zero_residual(poly_to_symbolic(left.value() - right.value()), context);
    if (!proof) { return Result<void>::failure(proof.error()); }
    if (!std::holds_alternative<ProvedZeroResidual>(proof.value())) { return undecided<void>(); }
    return Result<void>::success();
}

static Result<void> prove_segment_derivative(SymbolicExpr& effective,
    std::shared_ptr<SymbolicExpr>& derivative, const std::string& var,
    const Result<Boundary>& lower, const Result<Boundary>& upper,
    ComputationContext& context) {
    if (lower && upper) {
        auto assumptions = real_assumptions(var, context);
        if (!assumptions) { return Result<void>::failure(assumptions.error()); }
        SegmentFacts facts(assumptions.value(), var, lower.value(), upper.value());
        NormalizationVisitor normalize_derivative(context, facts);
        detail::node(derivative)->accept(normalize_derivative);
        if (facts.error) { return Result<void>::failure(*facts.error); }
        derivative = detail::make_expression_ptr(normalize_derivative.get_result());
        NormalizationVisitor normalize_original(context, facts);
        detail::node(effective)->accept(normalize_original);
        if (facts.error) { return Result<void>::failure(*facts.error); }
        effective = detail::expression_from_node(normalize_original.get_result());
    }
    auto residual = sym_sub(*derivative, effective);
    if (lower && upper) {
        auto assumptions = real_assumptions(var, context);
        if (!assumptions) { return Result<void>::failure(assumptions.error()); }
        SegmentFacts facts(assumptions.value(), var, lower.value(), upper.value());
        NormalizationVisitor normalize_residual(context, facts);
        detail::node(residual)->accept(normalize_residual);
        if (facts.error) { return Result<void>::failure(*facts.error); }
        residual = detail::make_expression_ptr(normalize_residual.get_result());
    }
    auto proof = check_zero_residual(residual, context);
    if (!proof) { return Result<void>::failure(proof.error()); }
    if (!std::holds_alternative<ProvedZeroResidual>(proof.value())) {
        return prove_rational_derivative(effective, derivative, var, context);
    }
    return Result<void>::success();
}

static Result<SymbolicExpr> evaluate_segment_endpoints(const SymbolicExpr& candidate,
    const std::string& var, const SymbolicExpr& a, const SymbolicExpr& b,
    bool concrete_endpoints, ComputationContext& context) {
    if (!concrete_endpoints) {
        auto high = candidate.substitute(var, detail::make_expression_ptr(b));
        auto low = candidate.substitute(var, detail::make_expression_ptr(a));
        if (!high || !low) { return Result<SymbolicExpr>::failure(CasErrc::InternalInvariant,
            "primitive substitution returned no expression", operation); }
        return *sym_sub(*high, *low)->simplify();
    }
    auto high = boundary_value_checked(candidate, var, b, LimitDirection::FromBelow, context);
    auto low = boundary_value_checked(candidate, var, a, LimitDirection::FromAbove, context);
    if (!high) { return high; }
    if (!low) { return low; }
    const int hs = infinity_sign(detail::node(high.value()));
    const int ls = -infinity_sign(detail::node(low.value()));
    if (hs && ls && hs != ls) { return divergent<SymbolicExpr>(); }
    if (hs || ls) { return *SymbolicExpr::infinity(hs ? hs : ls); }
    return *sym_sub(high.value(), low.value())->simplify();
}

Result<SymbolicExpr> definite_segment_checked(const SymbolicExpr& original,
    const SymbolicExpr& primitive, const std::string& var,
    const SymbolicExpr& a, const SymbolicExpr& b, ComputationContext& context) {
    auto continuous = certify_continuity_checked(original, var, a, b, context);
    if (!continuous) { return Result<SymbolicExpr>::failure(continuous.error()); }
    if (continuous.value() == Tribool::False) { return divergent<SymbolicExpr>(); }
    if (continuous.value() != Tribool::True) { return undecided<SymbolicExpr>(); }
    auto effective = segment_expression_checked(original, var, a, b, context);
    if (!effective) { return effective; }
    SymbolicExpr candidate = primitive;
    auto lower = boundary_checked(a, context);
    auto upper = boundary_checked(b, context);
    auto derivative = certify_primitive_derivative(candidate, var, a, b, lower, upper, context);
    if (!derivative) { return Result<SymbolicExpr>::failure(derivative.error()); }
    auto proof = prove_segment_derivative(effective.value(), derivative.value(), var, lower, upper, context);
    if (!proof) { return Result<SymbolicExpr>::failure(proof.error()); }
    return evaluate_segment_endpoints(candidate, var, a, b, lower && upper, context);
}
}
