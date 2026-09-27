#include "matcher.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS {

void RewriteEngine::add_rule(const Rule& rule) {
    rules.push_back(rule);
}

class RewriteVisitor : public LMCAS::detail::SymbolicVisitor {
public:
    const RewriteEngine& engine;
    ComputationContext& context;
    std::shared_ptr<const SymbolicNode> result;
    bool changed = false;

    RewriteVisitor(const RewriteEngine& e, ComputationContext& ctx)
        : engine(e), context(ctx) {}

    std::shared_ptr<const SymbolicNode> get_result() const { return result; }

    std::shared_ptr<const SymbolicNode> try_match(std::shared_ptr<const SymbolicNode> node) {
        auto current_expr = LMCAS::detail::expression_from_node(node);
        const auto& rules = engine.get_rules();
        const AssumptionContext* ctx = context.assumptions().get();
        for (const auto& rule : rules) {
            MatchMap bindings;

            if (Matcher::match(rule.pattern, current_expr, rule.wildcards, bindings)) {

                if (rule.condition && !rule.condition(bindings)) {
                    continue;
                }
                if (rule.assumption_condition &&
                    (!ctx || !rule.assumption_condition(bindings, ctx))) {
                    continue;
                }

                SymbolicExpr new_expr = Matcher::replace(rule.replacement, bindings, true);

                changed = true;
                return LMCAS::detail::node(new_expr);
            }
        }
        return node;
    }

    std::shared_ptr<const SymbolicNode> visit_child(const std::shared_ptr<const SymbolicNode>& child,
                                             bool& child_changed) {
        bool saved_changed = changed;
        changed = false;
        child->accept(*this);
        auto rewritten = result;
        if (changed || !rewritten->equals(*child)) {
            child_changed = true;
        }
        changed = saved_changed;
        return rewritten;
    }

    void finish_rewrite(std::shared_ptr<const SymbolicNode> node, bool child_changed,
                        bool original_changed) {
        changed = child_changed;
        auto matched = try_match(node);
        if (matched != node && !matched->equals(*node)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const NumberNode& node) override {

        result = try_match(node.clone());
    }

    void visit(const VariableNode& node) override {
        result = try_match(node.clone());
    }

    void visit(const AddNode& node) override {

        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool child_changed = false;

        bool original_changed = changed;
        changed = false;

        for (auto& op : node.operands()) {
            bool current_changed = changed;
            changed = false;
            op->accept(*this);
            new_ops.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        std::shared_ptr<const SymbolicNode> node_to_match;
        if (child_changed) {
            node_to_match = SymbolicFactory::create_add(new_ops);
        } else {
            node_to_match = LMCAS::detail::make_node<AddNode>(node.operands());
        }

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const MultiplyNode& node) override {
         std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool child_changed = false;
        bool original_changed = changed;
        changed = false;

        for (auto& op : node.operands()) {
            bool current_changed = changed;
            changed = false;
            op->accept(*this);
            new_ops.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        std::shared_ptr<const SymbolicNode> node_to_match;
        if (child_changed) {
            node_to_match = SymbolicFactory::create_multiply(new_ops);
        } else {
             node_to_match = LMCAS::detail::make_node<MultiplyNode>(node.operands());
        }

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const PowerNode& node) override {
        bool original_changed = changed;
        changed = false;

        node.base()->accept(*this);
        auto new_base = result;
        bool base_changed = changed;

        changed = false;
        node.exponent()->accept(*this);
        auto new_exp = result;
        bool exp_changed = changed;

        bool child_changed = base_changed || exp_changed;
        auto node_to_match = LMCAS::detail::make_node<PowerNode>(new_base, new_exp);

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const FunctionNode& node) override {
        bool original_changed = changed;
        changed = false;
        bool child_changed = false;

        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (auto& arg : node.arguments()) {
            bool current_changed = changed;
            changed = false;
            arg->accept(*this);
            new_args.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        auto node_to_match = LMCAS::detail::make_node<FunctionNode>(node.type(), new_args);

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const MatrixNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        if (const auto* dense =
                std::get_if<MatrixNode::DenseStorage>(&node.storage())) {
            MatrixNode::DenseStorage entries;
            entries.reserve(dense->size());
            for (const auto& entry : *dense) {
                entries.push_back(visit_child(entry, child_changed));
            }
            finish_rewrite(LMCAS::detail::make_node<MatrixNode>(
                               node.rows(), node.cols(), std::move(entries)),
                           child_changed, original_changed);
            return;
        }
        MatrixNode::SparseStorage entries;
        for (const auto& [index, entry] :
             std::get<MatrixNode::SparseStorage>(node.storage())) {
            entries.emplace(index, visit_child(entry, child_changed));
        }
        finish_rewrite(LMCAS::detail::make_node<MatrixNode>(
                           node.rows(), node.cols(), std::move(entries)),
                       child_changed, original_changed);
    }

    void visit(const RelationalNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto left = visit_child(node.left(), child_changed);
        auto right = visit_child(node.right(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<RelationalNode>(left, right, node.op()),
                       child_changed, original_changed);
    }

    void visit(const LogicalNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto left = visit_child(node.left(), child_changed);
        std::shared_ptr<const SymbolicNode> right = nullptr;
        if (node.right()) right = visit_child(node.right(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<LogicalNode>(left, right, node.op()),
                       child_changed, original_changed);
    }

    void visit(const PiecewiseNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(node.branches().size());
        for (const auto& branch : node.branches()) {
            auto expression = visit_child(branch.expression, child_changed);
            auto condition = visit_child(branch.condition, child_changed);
            branches.push_back({expression, condition});
        }
        std::shared_ptr<const SymbolicNode> default_expr = nullptr;
        if (node.default_expr()) default_expr = visit_child(node.default_expr(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<PiecewiseNode>(std::move(branches), default_expr),
                       child_changed, original_changed);
    }

    void visit(const SummationNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto body = visit_child(node.body(), child_changed);
        auto lower = visit_child(node.lower_bound(), child_changed);
        auto upper = visit_child(node.upper_bound(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<SummationNode>(body, node.index_var(), lower, upper),
                       child_changed, original_changed);
    }

    void visit(const ProductNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto body = visit_child(node.body(), child_changed);
        auto lower = visit_child(node.lower_bound(), child_changed);
        auto upper = visit_child(node.upper_bound(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<ProductNode>(body, node.index_var(), lower, upper),
                       child_changed, original_changed);
    }

    void visit(const TransformNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto body = visit_child(node.body(), child_changed);
        auto target = visit_child(node.target(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<TransformNode>(
                           node.transform_type(), body, node.source_var(), target),
                       child_changed, original_changed);
    }

    void visit(const QuantifierNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto domain = visit_child(node.domain(), child_changed);
        auto predicate = visit_child(node.predicate(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<QuantifierNode>(
                           node.quantifier_type(), node.bound_var(), domain, predicate),
                       child_changed, original_changed);
    }

    void visit(const SetBuilderNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto domain = visit_child(node.domain(), child_changed);
        auto predicate = visit_child(node.predicate(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<SetBuilderNode>(
                           node.element_var(), domain, predicate),
                       child_changed, original_changed);
    }

    void visit(const ComplexNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto real = visit_child(node.real(), child_changed);
        auto imag = visit_child(node.imag(), child_changed);
        finish_rewrite(SymbolicFactory::create_complex(real, imag),
                       child_changed, original_changed);
    }
    void visit(const FiniteSetNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        std::vector<std::shared_ptr<const SymbolicNode>> elements;
        for (const auto& element : node.elements()) {
            elements.push_back(visit_child(element, child_changed));
        }
        finish_rewrite(LMCAS::detail::make_node<FiniteSetNode>(std::move(elements)),
                       child_changed, original_changed);
    }
    void visit(const IntervalNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto lower = visit_child(node.lower(), child_changed);
        auto upper = visit_child(node.upper(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<IntervalNode>(
                           lower, upper, node.lower_closed(), node.upper_closed()),
                       child_changed, original_changed);
    }
    void visit(const MembershipNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto element = visit_child(node.element(), child_changed);
        auto set = visit_child(node.set(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<MembershipNode>(element, set),
                       child_changed, original_changed);
    }
    void visit(const UninterpretedFunctionNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        std::vector<std::shared_ptr<const SymbolicNode>> arguments;
        arguments.reserve(node.arguments().size());
        for (const auto& argument : node.arguments()) {
            arguments.push_back(visit_child(argument, child_changed));
        }
        finish_rewrite(LMCAS::detail::make_node<UninterpretedFunctionNode>(
                           node.name(), std::move(arguments)),
                       child_changed, original_changed);
    }
    void visit(const QuantityNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto value = visit_child(node.value(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<QuantityNode>(
                           value, node.dimension(), node.scale_to_base(), node.display_unit()),
                       child_changed, original_changed);
    }
    void visit(const IntegralNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto body = visit_child(node.body(), child_changed);
        auto lower = node.lower()
            ? visit_child(node.lower(), child_changed) : nullptr;
        auto upper = node.upper()
            ? visit_child(node.upper(), child_changed) : nullptr;
        finish_rewrite(LMCAS::detail::make_node<IntegralNode>(
                           body, node.variable(), lower, upper),
                       child_changed, original_changed);
    }
    void visit(const LimitNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto body = visit_child(node.body(), child_changed);
        auto point = visit_child(node.point(), child_changed);
        finish_rewrite(LMCAS::detail::make_node<LimitNode>(
                           body, node.variable(), point, node.direction()),
                       child_changed, original_changed);
    }
    void visit(const RootOfNode& node) override {
        const bool original_changed = changed;
        finish_rewrite(node.clone(), false, original_changed);
    }
};

Result<SymbolicExpr> RewriteEngine::apply_step_checked(
    const SymbolicExpr& expr,
    ComputationContext& context) const {
    auto budget = context.consume_steps(1, "rewrite.apply_step");
    if (!budget) return Result<SymbolicExpr>::failure(budget.error());
    if (!LMCAS::detail::node(expr)) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InvalidArgument, "rewrite expression must not be null",
            "rewrite.apply_step");
    }
    RewriteVisitor v(*this, context);
    LMCAS::detail::node(expr)->accept(v);
    return Result<SymbolicExpr>::success(
        LMCAS::detail::expression_from_node(v.get_result()));
}

Result<SymbolicExpr> RewriteEngine::apply_checked(
    const SymbolicExpr& expr,
    ComputationContext& context,
    int max_iterations) const {
    if (max_iterations < 0) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InvalidArgument, "rewrite iteration count must be non-negative",
            "rewrite.apply");
    }
    SymbolicExpr current = expr;
    for (int i = 0; i < max_iterations; ++i) {
        auto next_result = apply_step_checked(current, context);
        if (!next_result) return next_result;
        SymbolicExpr next = std::move(next_result.value());
        if (LMCAS::detail::node(next)->equals(*LMCAS::detail::node(current))) {
            return Result<SymbolicExpr>::success(std::move(current));
        }
        current = next;
    }
    return Result<SymbolicExpr>::success(std::move(current));
}

Result<SymbolicExpr> RewriteEngine::apply_step(const SymbolicExpr& expr) const {
    ComputationContext context;
    return apply_step_checked(expr, context);
}

Result<SymbolicExpr> RewriteEngine::apply(
    const SymbolicExpr& expr,
    int max_iterations) const {
    ComputationContext context;
    return apply_checked(expr, context, max_iterations);
}

}
