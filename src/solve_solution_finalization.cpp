#include "internal/transcendental_solver_support.hpp"
#include "assumption_context.hpp"
#include "root_of_utils.hpp"
#include "residual_verification.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/expression_analysis.hpp"

namespace LMCAS {

ExprPtr detail::equation_predicate(const ExprPtr& expression) {
    if (std::dynamic_pointer_cast<const RelationalNode>(detail::node(expression))) {
        return expression;
    }
    return detail::make_expression_ptr(detail::make_node<RelationalNode>(
        detail::node(expression), detail::node(SymbolicExpr::number(0)), RelationOp::EQ));
}

ExprPtr detail::substitute_raw(const ExprPtr& expression, const std::string& variable,
                               const ExprPtr& value) {
    return detail::make_expression_ptr(substitute_free(
        detail::node(expression), variable, detail::node(value)));
}

namespace {
using detail::equation_predicate;
using detail::substitute_raw;

Result<void> apply_root_representation_policy(ExprPtr& root, const SolveOptions& options) {
    if (!options.return_rootof) {
        if (std::dynamic_pointer_cast<const RootOfNode>(detail::node(root))) {
            root = rootof_simplify(root);
        }
        if (!root || detail::contains_node_type<RootOfNode>(detail::node(root))) {
            return Result<void>::failure(CasErrc::Inconclusive,
                "complete roots require RootOf representation", "solve.representation");
        }
    }
    return Result<void>::success();
}

ExprPtr equation_residual(const ExprPtr& expression) {
    if (auto relation = std::dynamic_pointer_cast<const RelationalNode>(detail::node(expression))) {
        return SymbolicExpr::add(detail::make_expression_ptr(relation->left()),
            SymbolicExpr::multiply(SymbolicExpr::number(-1),
                detail::make_expression_ptr(relation->right())));
    }
    return expression;
}

void append_condition(std::vector<ExprPtr>& conditions, const ExprPtr& condition) {
    for (const auto& existing : conditions) {
        if (detail::node(existing)->equals(*detail::node(condition))) { return; }
    }
    conditions.push_back(condition);
}

Tribool negate_truth(Tribool value) {
    return value == Tribool::True ? Tribool::False :
        value == Tribool::False ? Tribool::True : Tribool::Unknown;
}

Result<Tribool> solve_condition_truth(const ExprPtr& expression,
    const FactsQuery& facts, ComputationContext& context);

Result<Tribool> logical_condition_truth(const LogicalNode& logic,
    const FactsQuery& facts, ComputationContext& context) {
    auto left = solve_condition_truth(detail::make_expression_ptr(logic.left()), facts, context);
    if (!left) { return left; }
    if (logic.op() == LogicalNode::Op::Not) { return negate_truth(left.value()); }
    auto right = solve_condition_truth(detail::make_expression_ptr(logic.right()), facts, context);
    if (!right) { return right; }
    auto l = left.value();
    auto r = right.value();
    if (logic.op() == LogicalNode::Op::Implies) { l = negate_truth(l); }
    if (logic.op() == LogicalNode::Op::And) {
        if (l == Tribool::False || r == Tribool::False) { return Tribool::False; }
        return l == Tribool::True && r == Tribool::True ? Tribool::True : Tribool::Unknown;
    }
    if (l == Tribool::True || r == Tribool::True) { return Tribool::True; }
    return l == Tribool::False && r == Tribool::False ? Tribool::False : Tribool::Unknown;
}

Result<Tribool> equality_condition_truth(const ExprPtr& difference, Tribool nonzero,
    RelationOp op, ComputationContext& context) {
    if (nonzero == Tribool::Unknown) {
        auto residual = check_zero_residual(difference, context);
        if (!residual) { return Result<Tribool>::failure(residual.error()); }
        if (std::holds_alternative<ProvedZeroResidual>(residual.value())) { nonzero = Tribool::False; }
    }
    if (op == RelationOp::NEQ) { return nonzero; }
    return negate_truth(nonzero);
}

Result<Tribool> relational_condition_truth(const ExprPtr& expression,
    const RelationalNode& relation, const FactsQuery& facts, ComputationContext& context) {
    auto defined = detail::query_definedness(detail::node(expression), facts, Domain::Complex, context);
    if (!defined) { return Result<Tribool>::failure(defined.error()); }
    if (defined.value() != Tribool::True) { return Tribool::Unknown; }
    auto difference = equation_residual(expression)->simplify();
    auto nz = detail::query_nonzero_value(detail::node(difference), facts, Domain::Complex, context);
    if (!nz) { return nz; }
    if (relation.op() == RelationOp::EQ || relation.op() == RelationOp::NEQ) {
        return equality_condition_truth(difference, nz.value(), relation.op(), context);
    }
    if (relation.op() == RelationOp::LT || relation.op() == RelationOp::LEQ) {
        difference = SymbolicExpr::multiply(SymbolicExpr::number(-1), difference)->simplify();
    }
    if (relation.op() == RelationOp::GEQ || relation.op() == RelationOp::LEQ) {
        return detail::query_nonnegative_value(detail::node(difference), facts, context);
    }
    return detail::query_positive_value(detail::node(difference), facts, context);
}
Result<Tribool> solve_condition_truth(const ExprPtr& expression,
    const FactsQuery& facts, ComputationContext& context) {
    auto step = context.consume_steps(1, "solve.conditions");
    if (!step) { return Result<Tribool>::failure(step.error()); }
    const auto node = detail::node(expression);
    if (auto logic = std::dynamic_pointer_cast<const LogicalNode>(node)) {
        return logical_condition_truth(*logic, facts, context);
    }
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(node);
    if (!relation) { return Tribool::Unknown; }
    return relational_condition_truth(expression, *relation, facts, context);
}

class SolutionFinalizer {
    const ExprPtr& original;
    const std::string& variable;
    const FactsQuery& facts;
    ComputationContext& context;
    Domain domain;
    const SolveOptions& options;

    Result<bool> finish_conditions(std::vector<ExprPtr>& conditions) const {
        std::vector<ExprPtr> remaining;
        for (const auto& condition : conditions) {
            auto truth = solve_condition_truth(condition, facts, context);
            if (!truth) { return Result<bool>::failure(truth.error()); }
            if (truth.value() == Tribool::False) { return false; }
            if (truth.value() == Tribool::Unknown) { append_condition(remaining, condition); }
        }
        conditions = std::move(remaining);
        return true;
    }

    Result<bool> finalize_value(const ExprPtr& value, std::vector<ExprPtr>& conditions) const {
        for (auto& condition : conditions) {
            condition = substitute_raw(condition, variable, value);
        }
        auto bound = substitute_raw(original, variable, value);
        auto defined = detail::query_definedness(detail::node(bound), facts, domain, context);
        if (!defined) { return Result<bool>::failure(defined.error()); }
        if (defined.value() == Tribool::False) { return false; }
        auto domain_conditions = detail::domain_constraints(detail::node(bound), facts, domain, context);
        if (!domain_conditions) { return Result<bool>::failure(domain_conditions.error()); }
        if (!domain_conditions.value()) {
            return Result<bool>::failure(CasErrc::Inconclusive,
                "candidate domain cannot be completely expressed", "solve.finalize");
        }
        for (const auto& condition : *domain_conditions.value()) { append_condition(conditions, condition); }
        auto valid = finish_conditions(conditions);
        if (!valid || !valid.value()) { return valid; }
        auto residual = check_zero_residual(equation_residual(bound), context);
        if (!residual) { return Result<bool>::failure(residual.error()); }
        if (!std::holds_alternative<ProvedZeroResidual>(residual.value())) {
            auto pointwise_nonzero = detail::query_nonzero_value(detail::node(equation_residual(bound)->simplify()), facts, domain, context);
            if (!pointwise_nonzero) { return Result<bool>::failure(pointwise_nonzero.error()); }
            if (pointwise_nonzero.value() == Tribool::True) { return false; }
            if (pointwise_nonzero.value() == Tribool::Unknown) {
                append_condition(conditions, equation_predicate(bound));
            }
        }
        return true;
    }

    template<class Solutions>
    SolveResult filter_values(Solutions& solutions) const {
        Solutions kept;
        for (auto& solution : solutions.values) {
            auto valid = finalize_value(solution.value, solution.conditions);
            if (!valid) { return SolveResult::failure(valid.error()); }
            if (valid.value()) {
                auto represented = apply_root_representation_policy(solution.value, options);
                if (!represented) { return SolveResult::failure(represented.error()); }
                kept.values.push_back(std::move(solution));
            }
        }
        if (kept.values.empty()) { return SolutionSet{EmptySolutions{}}; }
        return SolutionSet{std::move(kept)};
    }

public:
    SolutionFinalizer(const ExprPtr& original, const std::string& variable,
        const FactsQuery& facts, ComputationContext& context, Domain domain,
        const SolveOptions& options)
        : original(original), variable(variable), facts(facts), context(context),
          domain(domain), options(options) {}

    SolveResult finish(SolutionSet solutions) const {
        if (auto finite = std::get_if<FiniteSolutions>(&solutions)) { return filter_values(*finite); }
        if (auto parametric = std::get_if<ParametricSolutions>(&solutions)) { return filter_values(*parametric); }
        auto obligations = detail::domain_constraints(detail::node(original), facts, domain, context);
        if (!obligations) { return SolveResult::failure(obligations.error()); }
        if (!obligations.value()) {
            return SolveResult::failure(CasErrc::Inconclusive,
                "original equation domain cannot be completely expressed", "solve.finalize");
        }
        if (auto conditional = std::get_if<ConditionalSolutions>(&solutions)) {
            for (const auto& condition : *obligations.value()) {
                append_condition(conditional->value.conditions, condition);
            }
            auto valid = finish_conditions(conditional->value.conditions);
            if (!valid) { return SolveResult::failure(valid.error()); }
            if (!valid.value()) { return SolutionSet{EmptySolutions{}}; }
            return solutions;
        }
        auto conditions = std::move(*obligations.value());
        auto valid = finish_conditions(conditions);
        if (!valid) { return SolveResult::failure(valid.error()); }
        if (!valid.value()) { return SolutionSet{EmptySolutions{}}; }
        if (conditions.empty()) { return solutions; }
        return SolutionSet{ConditionalSolutions{ConditionSet{variable,
            equation_predicate(original), std::move(conditions)}}};
    }
};
}

SolveResult detail::finalize_solution_set(const ExprPtr& original,
    const std::string& variable, SolutionSet solutions, ComputationContext& context,
    Domain domain, const SolveOptions& options) {
    if (std::holds_alternative<EmptySolutions>(solutions)) { return solutions; }
    AssumptionContext local_facts = context.assumptions() ?
        *context.assumptions() : AssumptionContext{};
    if (domain == Domain::Real) {
        auto declared = local_facts.assume_domain_checked(variable, Domain::Real);
        if (!declared) { return SolveResult::failure(declared.error()); }
    }
    if (const auto* families = std::get_if<ParametricSolutions>(&solutions)) {
        for (const auto& family : families->values) {
            for (const auto& parameter : family.integer_parameters) {
                auto declared = local_facts.assume_domain_checked(parameter, Domain::Integer);
                if (!declared) { return SolveResult::failure(declared.error()); }
            }
        }
    }
    detail::AssumptionFacts adapter(local_facts);
    return SolutionFinalizer(original, variable, adapter, context, domain, options).finish(std::move(solutions));
}

}
