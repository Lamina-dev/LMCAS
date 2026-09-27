#include "limit_result.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/facts_query.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/limit_visitor.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/pointwise_comparison.hpp"
#include "polynomial_conversion.hpp"
#include "symbolic.hpp"
#include <exception>

namespace LMCAS {
namespace {
using NodePtr = std::shared_ptr<const SymbolicNode>;
LimitResult outcome(LimitOutcome value) {
    return LimitResult::success(Verified<LimitOutcome>{std::move(value), ByConstructionProof{}});
}
LimitResult unknown() {
    return LimitResult::failure(CasErrc::Inconclusive,
        "limit could not be proved in the supported domain", "limit");
}
LimitResult assumption_limit_failure(const CasError& error) {
    return error.code == CasErrc::InvalidArgument ? unknown() : LimitResult::failure(error);
}
int infinity_sign(const LimitOutcome& value) {
    if (std::holds_alternative<PositiveInfinityLimit>(value)) return 1;
    if (std::holds_alternative<NegativeInfinityLimit>(value)) return -1;
    return 0;
}
bool has_infinity(const NodePtr& node) {
    class Detector final : public detail::RecursiveSymbolicVisitor {
    public:
        bool found = false;
        void visit(const FunctionNode& f) override {
            if (f.type() == FunctionNode::FuncType::Infinity) found = true;
            else detail::RecursiveSymbolicVisitor::visit(f);
        }
    } detector;
    if (node) node->accept(detector);
    return detector.found;
}
Result<int> finite_sign(const NodePtr& node, const FactsQuery& facts, ComputationContext& context) {
    auto nz = detail::query_nonzero_value(node, facts, Domain::Real, context);
    if (!nz) { return Result<int>::failure(nz.error()); }
    if (nz.value() == Tribool::False) { return 0; }
    auto positive = detail::query_positive_value(node, facts, context);
    if (!positive) { return Result<int>::failure(positive.error()); }
    if (positive.value() == Tribool::True) { return 1; }
    auto nonnegative = detail::query_nonnegative_value(node, facts, context);
    if (!nonnegative) { return Result<int>::failure(nonnegative.error()); }
    if (nonnegative.value() == Tribool::False) { return -1; }
    auto negative = facts.is_negative(node, context);
    if (!negative) { return Result<int>::failure(negative.error()); }
    if (negative.value() == Tribool::True) { return -1; }
    return Result<int>::failure(CasErrc::Inconclusive, "finite sign is unknown", "limit");
}
LimitResult classify_limit_result(const NodePtr& node, const std::string& variable,
    ComputationContext& context, const FactsQuery& facts, Domain domain = Domain::Real);

LimitResult classify_infinite_sum(const AddNode& sum, const std::string& variable,
    ComputationContext& context, const FactsQuery& facts, Domain domain) {
    int sign = 0;
    std::vector<NodePtr> finite_terms;
    for (const auto& term : sum.operands()) {
        auto value = classify_limit_result(term, variable, context, facts, domain);
        if (!value) { return value; }
        int next = infinity_sign(value.value().value);
        if (sign && next && sign != next) { return unknown(); }
        if (next) { sign = next; }
        else { finite_terms.push_back(detail::node(std::get<FiniteLimit>(value.value().value).value)); }
    }
    if (sign) { return sign > 0 ? outcome(PositiveInfinityLimit{}) : outcome(NegativeInfinityLimit{}); }
    auto finite = detail::make_expression_ptr(detail::make_node<AddNode>(std::move(finite_terms)))->simplify();
    return finite ? classify_limit_result(detail::node(finite), variable, context, facts, domain) : unknown();
}

LimitResult classify_infinite_product(const MultiplyNode& product, const std::string& variable,
    ComputationContext& context, const FactsQuery& facts, Domain domain) {
    int sign = 1;
    bool infinite = false;
    std::vector<NodePtr> finite_factors;
    for (const auto& factor : product.operands()) {
        auto value = classify_limit_result(factor, variable, context, facts, domain);
        if (!value) { return value; }
        int next = infinity_sign(value.value().value);
        if (next) { infinite = true; sign *= next; }
        else { finite_factors.push_back(detail::node(std::get<FiniteLimit>(value.value().value).value)); }
    }
    if (infinite) {
        for (const auto& factor : finite_factors) {
            auto s = finite_sign(factor, facts, context);
            if (!s) { return LimitResult::failure(s.error()); }
            if (!s.value()) { return unknown(); }
            sign *= s.value();
        }
        return sign > 0 ? outcome(PositiveInfinityLimit{}) : outcome(NegativeInfinityLimit{});
    }
    auto finite = detail::make_expression_ptr(detail::make_node<MultiplyNode>(std::move(finite_factors)))->simplify();
    return finite ? classify_limit_result(detail::node(finite), variable, context, facts, domain) : unknown();
}

LimitResult classify_infinite_power(const PowerNode& power, const std::string& variable,
    ComputationContext& context, const FactsQuery& facts, Domain domain) {
    auto base = classify_limit_result(power.base(), variable, context, facts, domain);
    if (!base) { return base; }
    const int sign = infinity_sign(base.value().value);
    BigInt integer;
    auto number = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
    auto exponent = extract_coeff_value<Rational>(detail::make_expression_ptr(power.exponent()));
    bool integral = number && try_get_integer_value(number, integer);
    if (!sign || !exponent || (sign <= 0 && !integral)) { return unknown(); }
    if (exponent.value() < Rational(0)) { return outcome(FiniteLimit{SymbolicExpr::number(0)}); }
    if (exponent.value() > Rational(0)) { return sign < 0 && integer.is_odd()
        ? outcome(NegativeInfinityLimit{}) : outcome(PositiveInfinityLimit{}); }
    return unknown();
}

LimitResult classify_infinite_function(const FunctionNode& function, const std::string& variable,
    ComputationContext& context, const FactsQuery& facts, Domain domain) {
    if (function.arguments().size() != 1) { return unknown(); }
    auto argument = classify_limit_result(function.arguments()[0], variable, context, facts, domain);
    if (!argument) { return argument; }
    int sign = infinity_sign(argument.value().value);
    if (!sign) {
        auto applied = detail::make_node<FunctionNode>(function.type(), std::vector<NodePtr>{
            detail::node(std::get<FiniteLimit>(argument.value().value).value)});
        return classify_limit_result(applied, variable, context, facts, domain);
    }
    switch (function.type()) {
        case FunctionNode::FuncType::Exp: {
            return sign > 0 ? outcome(PositiveInfinityLimit{}) : outcome(FiniteLimit{SymbolicExpr::number(0)});
        }
        case FunctionNode::FuncType::Abs: { return outcome(PositiveInfinityLimit{}); }
        case FunctionNode::FuncType::Ln: case FunctionNode::FuncType::Sqrt: {
            return sign > 0 ? outcome(PositiveInfinityLimit{}) : unknown();
        }
        case FunctionNode::FuncType::ArcTan: {
            return outcome(FiniteLimit{SymbolicExpr::multiply(SymbolicExpr::number(Rational(sign,2)),
                SymbolicExpr::variable("pi"))});
        }
        default: { return unknown(); }
    }
}

LimitResult classify_limit_result(const NodePtr& node, const std::string& variable,
    ComputationContext& context, const FactsQuery& facts, Domain domain) {
    auto step = context.consume_steps(1, "limit.classify");
    if (!step) { return LimitResult::failure(step.error()); }
    if (!node || detail::contains_node_type<LimitNode>(node) || expression_depends_on_variable(node, variable))
        { return unknown(); }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() == FunctionNode::FuncType::Infinity) { return outcome(PositiveInfinityLimit{}); }
    }
    if (has_infinity(node)) {
        if (auto sum = std::dynamic_pointer_cast<const AddNode>(node)) {
            return classify_infinite_sum(*sum, variable, context, facts, domain);
        }
        if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
            return classify_infinite_product(*product, variable, context, facts, domain);
        }
        if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
            return classify_infinite_power(*power, variable, context, facts, domain);
        }
        if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
            return classify_infinite_function(*function, variable, context, facts, domain);
        }
        return unknown();
    }
    auto defined = detail::query_definedness(node, facts, domain, context);
    if (!defined) { return LimitResult::failure(defined.error()); }
    if (defined.value() != Tribool::True) { return unknown(); }
    return outcome(FiniteLimit{detail::make_expression_ptr(node)});
}

struct ComplexContinuityPoint {
    const std::string& variable;
    const NodePtr& point;
    const FactsQuery& facts;
    ComputationContext& context;

    bool positive(const NodePtr& argument) const {
        auto substituted = substitute_free(argument, variable, point);
        auto value = detail::query_positive_value(substituted, facts, context);
        if (!value) throw value.error();
        return value.value() == Tribool::True;
    }
};

bool elementary_continuous_structure(
    const NodePtr& node, const ComplexContinuityPoint* complex_point = nullptr);

bool power_continuous_structure(const PowerNode& power, const ComplexContinuityPoint* point) {
    if (point) {
        BigInt integer;
        auto exponent = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
        if ((!exponent || !try_get_integer_value(exponent, integer)) &&
            !point->positive(power.base())) { return false; }
    }
    return elementary_continuous_structure(power.base(), point) &&
        elementary_continuous_structure(power.exponent(), point);
}

bool continuous_function_branch(FunctionNode::FuncType type, const NodePtr& argument,
    const ComplexContinuityPoint* point) {
    switch (type) {
        case FunctionNode::FuncType::Sin: { return true; }
        case FunctionNode::FuncType::Cos: { return true; }
        case FunctionNode::FuncType::Tan: { return true; }
        case FunctionNode::FuncType::Exp: { return true; }
        case FunctionNode::FuncType::Abs: { return true; }
        case FunctionNode::FuncType::Ln: { return !point || point->positive(argument); }
        case FunctionNode::FuncType::Sqrt: { return !point || point->positive(argument); }
        case FunctionNode::FuncType::ArcTan: { return point == nullptr; }
        case FunctionNode::FuncType::ArcSin: { return point == nullptr; }
        case FunctionNode::FuncType::ArcCos: { return point == nullptr; }
        default: { return false; }
    }
}
bool elementary_continuous_structure(const NodePtr& node, const ComplexContinuityPoint* point) {
    if (std::dynamic_pointer_cast<const NumberNode>(node) ||
        std::dynamic_pointer_cast<const VariableNode>(node)) { return true; }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return power_continuous_structure(*power, point);
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return point && elementary_continuous_structure(complex->real(), point) &&
            elementary_continuous_structure(complex->imag(), point);
    }
    const std::vector<NodePtr>* operands = nullptr;
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(node)) {
        operands = &sum->operands();
    } else if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        operands = &product->operands();
    }
    if (operands) {
        for (const auto& operand : *operands) {
            if (!elementary_continuous_structure(operand, point)) { return false; }
        }
        return true;
    }
    auto function = std::dynamic_pointer_cast<const FunctionNode>(node);
    return function && function->arguments().size() == 1 &&
        continuous_function_branch(function->type(), function->arguments()[0], point) &&
        elementary_continuous_structure(function->arguments()[0], point);
}

class SideBranchRewriter final : public detail::SymbolicRewriter {
public:
    SideBranchRewriter(LimitVisitor& visitor, const FactsQuery& facts, ComputationContext& context)
        : visitor_(visitor), facts_(facts), context_(context) {}
    void visit(const PiecewiseNode& node) override {
        auto selected = visitor_.select_piecewise_near_point(node);
        set_result(selected ? rewrite(selected) : current());
    }
    void visit(const PowerNode& node) override {
        bool changed = false;
        auto base = rewrite_child(node.base(), changed);
        auto exponent = rewrite_child(node.exponent(), changed);
        const auto square = std::dynamic_pointer_cast<const PowerNode>(base);
        if (square) {
            auto outer = extract_coeff_value<Rational>(detail::make_expression_ptr(exponent));
            auto inner = extract_coeff_value<Rational>(detail::make_expression_ptr(square->exponent()));
            if (outer && inner && outer.value() == Rational(1, 2) &&
                inner.value() == Rational(2) && rewrite_absolute(square->base())) return;
        }
        set_result(changed ? detail::make_node<PowerNode>(base, exponent) : current());
    }
    void visit(const FunctionNode& node) override {
        bool changed = false;
        auto args = rewrite_children(node.arguments(), changed);
        if (args.size() == 1) {
            if (node.type() == FunctionNode::FuncType::Abs && rewrite_absolute(args[0])) { return; }
            if (node.type() == FunctionNode::FuncType::Sqrt && rewrite_square_absolute(args[0])) { return; }
        }
        set_result(changed ? detail::make_node<FunctionNode>(node.type(), std::move(args)) : current());
    }
private:
    bool rewrite_square_absolute(const NodePtr& argument) {
        auto square = std::dynamic_pointer_cast<const PowerNode>(argument);
        if (!square) { return false; }
        auto exponent = extract_coeff_value<Rational>(detail::make_expression_ptr(square->exponent()));
        return exponent && exponent.value() == Rational(2) && rewrite_absolute(square->base());
    }
    bool rewrite_absolute(const NodePtr& argument) {
        auto real = detail::query_real_value(argument, facts_, context_);
        if (!real) throw real.error();
        auto sign = real.value() == Tribool::True ? visitor_.sign_near_point(argument) : std::nullopt;
        if (!sign) return false;
        set_result(*sign < 0 ? detail::node(SymbolicExpr::multiply(
            SymbolicExpr::number(-1), detail::make_expression_ptr(argument))) : argument);
        return true;
    }
    LimitVisitor& visitor_;
    const FactsQuery& facts_;
    ComputationContext& context_;
};

std::optional<LimitResult> scaled_one_sided_limit(const MultiplyNode& product,
    const std::string& variable, const LimitExprPtr& point, LimitDirection direction,
    ComputationContext& context, const FactsQuery& facts) {
    std::shared_ptr<const SymbolicNode> varying;
    bool invertible_scale = true;
    for (const auto& factor : product.operands()) {
        if (expression_depends_on_variable(factor, variable)) {
            if (varying) { invertible_scale = false; break; }
            varying = factor;
        } else {
            auto nonzero = detail::query_nonzero_value(factor, facts, Domain::Real, context);
            if (!nonzero) { return LimitResult::failure(nonzero.error()); }
            if (nonzero.value() != Tribool::True) { invertible_scale = false; break; }
        }
    }
    if (invertible_scale && varying) {
        auto inner = limit_checked(detail::make_expression_ptr(varying), variable, point, direction, context);
        if (!inner && inner.error().code != CasErrc::Inconclusive) { return inner; }
        if (inner && std::holds_alternative<LimitDoesNotExist>(inner.value().value)) { return inner; }
    }
    return std::nullopt;
}

std::optional<LimitResult> oscillatory_one_sided_limit(const NodePtr& node,
    const std::string& variable, const LimitExprPtr& point, LimitDirection direction,
    ComputationContext& context, const FactsQuery& facts, LimitVisitor& visitor) {
    auto f = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!f || f->arguments().size() != 1) { return std::nullopt; }
    if (f->type() != FunctionNode::FuncType::Sin && f->type() != FunctionNode::FuncType::Cos) {
        return std::nullopt;
    }
    auto argument = detail::make_expression_ptr(f->arguments()[0]);
    auto inner = limit_checked(argument, variable, point, direction, context);
    if (!inner && inner.error().code != CasErrc::Inconclusive) { return inner; }
    if (!inner || !infinity_sign(inner.value().value) || !elementary_continuous_structure(f->arguments()[0])) {
        return std::nullopt;
    }
    auto domain = detail::domain_constraints(f->arguments()[0], facts, Domain::Real, context);
    if (!domain) { return LimitResult::failure(domain.error()); }
    if (!domain.value()) { return std::nullopt; }
    for (const auto& condition : *domain.value()) {
        if (!visitor.condition_holds_near_point(detail::node(condition))) { return std::nullopt; }
    }
    return outcome(LimitDoesNotExist{});
}

LimitResult evaluate_one_sided_node(const NodePtr& node, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, ComputationContext& context,
    const AssumptionContext& local, const FactsQuery& facts, LimitVisitor& visitor) {
    if (!expression_depends_on_variable(node, variable))
        { return classify_limit_result(node, variable, context, facts); }
    auto original_domain = detail::domain_constraints(node, facts, Domain::Real, context);
    if (!original_domain) { return LimitResult::failure(original_domain.error()); }
    if (!original_domain.value()) { return unknown(); }
    for (const auto& condition : *original_domain.value())
        { if (!visitor.condition_holds_near_point(detail::node(condition))) { return unknown(); } }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        auto scaled = scaled_one_sided_limit(*product, variable, point, direction, context, facts);
        if (scaled) { return std::move(*scaled); }
    }
    auto oscillatory = oscillatory_one_sided_limit(node, variable, point, direction, context, facts, visitor);
    if (oscillatory) { return std::move(*oscillatory); }
    auto prepared = local.simplify(detail::expression_from_node(node));
    if (!prepared) { return unknown(); }
    detail::node(prepared)->accept(visitor);
    return classify_limit_result(visitor.get_result(), variable, context, facts);
}


LimitResult limit_one_sided_checked(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, ComputationContext& context) {
    AssumptionContext local = context.assumptions() ? *context.assumptions() : AssumptionContext{};
    auto real_variable = local.assume_domain_checked(variable, Domain::Real);
    if (!real_variable) { return assumption_limit_failure(real_variable.error()); }
    if (has_infinity(detail::node(point))) {
        const int sign = std::dynamic_pointer_cast<const FunctionNode>(detail::node(point)) ? 1 : -1;
        auto tail = local.assume_sign_checked(variable, sign < 0 ? Sign::Negative : Sign::Positive);
        if (!tail) { return assumption_limit_failure(tail.error()); }
    }
    detail::AssumptionFacts facts(local);
    const char* token = direction == LimitDirection::FromBelow ? "-" : "+";
    LimitVisitor visitor(variable, detail::node(point), token, &local, &context);
    if (!has_infinity(detail::node(point))) {
        auto sign = visitor.sign_near_point(detail::make_node<VariableNode>(variable));
        if (sign && *sign != 0) {
            auto side = local.assume_sign_checked(variable, *sign < 0 ? Sign::Negative : Sign::Positive);
            if (!side) { return assumption_limit_failure(side.error()); }
        }
    }
    SideBranchRewriter branches(visitor, facts, context);
    auto node = branches.rewrite(detail::node(expression));
    return evaluate_one_sided_node(node, variable, point, direction, context, local, facts, visitor);
}

LimitResult combine_one_sided_limits(const LimitOutcome& left, const LimitOutcome& right,
    ComputationContext& context) {
    if (std::holds_alternative<LimitDoesNotExist>(left) || std::holds_alternative<LimitDoesNotExist>(right))
        return outcome(LimitDoesNotExist{});
    int ls = infinity_sign(left), rs = infinity_sign(right);
    if (ls || rs) {
        if (ls == rs) return outcome(left);
        return outcome(LimitDoesNotExist{});
    }
    const auto& lv = std::get<FiniteLimit>(left).value;
    const auto& rv = std::get<FiniteLimit>(right).value;
    auto equal = detail::compare_pointwise_values(lv, rv, context, Domain::Real);
    if (!equal) return LimitResult::failure(equal.error());
    if (equal.value() == Tribool::True) return outcome(left);
    if (equal.value() == Tribool::False) return outcome(LimitDoesNotExist{});
    return unknown();
}

class ComplexPuncturedFacts final : public FactsQuery {
public:
    ComplexPuncturedFacts(const FactsQuery& base, const std::vector<NodePtr>& nonzero)
        : base_(base), nonzero_(nonzero) {}
    Result<Tribool> is_positive(const NodePtr& node, ComputationContext& context) const override { return base_.is_positive(node, context); }
    Result<Tribool> is_negative(const NodePtr& node, ComputationContext& context) const override { return base_.is_negative(node, context); }
    Result<Tribool> is_nonnegative(const NodePtr& node, ComputationContext& context) const override { return base_.is_nonnegative(node, context); }
    Result<Tribool> is_real(const NodePtr& node, ComputationContext& context) const override { return base_.is_real(node, context); }
    Result<Tribool> is_nonzero(const NodePtr& node, ComputationContext& context) const override {
        auto access = context.consume_steps(0, "facts.complex_punctured");
        if (!access) { return Result<Tribool>::failure(access.error()); }
        for (const auto& proven : nonzero_) {
            auto step = context.consume_steps(1, "facts.complex_punctured");
            if (!step) { return Result<Tribool>::failure(step.error()); }
            if (node->equals(*proven)) { return Tribool::True; }
        }
        return base_.is_nonzero(node, context);
    }
private:
    const FactsQuery& base_;
    const std::vector<NodePtr>& nonzero_;
};

Result<bool> complex_polynomial_nonzero(const RelationalNode& relation, const std::string& variable,
    const FactsQuery& facts, ComputationContext& context, std::vector<NodePtr>& nonzero_bases) {
    auto difference = SymbolicExpr::add(detail::make_expression_ptr(relation.left()),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), detail::make_expression_ptr(relation.right())));
    auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(difference, variable);
    if (!polynomial) {
        if (polynomial.error().code == CasErrc::UnsupportedExpression) { return false; }
        return Result<bool>::failure(polynomial.error());
    }
    bool nonzero = false;
    for (const auto& coefficient : polynomial.value().coeffs) {
        auto step = context.consume_steps(1, "limit.complex_domain");
        if (!step) { return Result<bool>::failure(step.error()); }
        auto defined = detail::query_definedness(detail::node(coefficient.val), facts, Domain::Complex, context);
        if (!defined) { return Result<bool>::failure(defined.error()); }
        if (defined.value() != Tribool::True) { return false; }
        auto value = detail::query_nonzero_value(detail::node(coefficient.val), facts, Domain::Complex, context);
        if (!value) { return Result<bool>::failure(value.error()); }
        nonzero |= value.value() == Tribool::True;
    }
    if (!nonzero) { return false; }
    auto normalized = difference->simplify();
    if (!normalized) { return false; }
    nonzero_bases.push_back(detail::node(normalized));
    return true;
}

Result<bool> complex_punctured_condition(const NodePtr& condition, const std::string& variable,
    const FactsQuery& facts, ComputationContext& context, std::vector<NodePtr>& nonzero_bases) {
    auto step = context.consume_steps(1, "limit.complex_condition");
    if (!step) { return Result<bool>::failure(step.error()); }
    auto logical = std::dynamic_pointer_cast<const LogicalNode>(condition);
    if (logical && logical->op() == LogicalNode::Op::Or) {
        auto left = complex_punctured_condition(logical->left(), variable, facts, context, nonzero_bases);
        if (!left || left.value()) { return left; }
        return complex_punctured_condition(logical->right(), variable, facts, context, nonzero_bases);
    }
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(condition);
    if (!relation || relation->op() != RelationalNode::Op::NEQ) { return false; }
    return complex_polynomial_nonzero(*relation, variable, facts, context, nonzero_bases);
}

Result<bool> complex_punctured_domain(const NodePtr& node, const std::string& variable,
    const FactsQuery& facts, ComputationContext& context, std::vector<NodePtr>& nonzero_bases) {
    auto domain = detail::domain_constraints(node, facts, Domain::Complex, context);
    if (!domain) { return Result<bool>::failure(domain.error()); }
    if (!domain.value()) { return false; }
    for (const auto& condition : *domain.value()) {
        auto proven = complex_punctured_condition(detail::node(condition), variable, facts, context, nonzero_bases);
        if (!proven || !proven.value()) { return proven; }
    }
    return true;
}

Result<void> validate_limit_input(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, Domain domain) {
    if (!expression || !point || !detail::node(expression) || !detail::node(point) || variable.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
            "limit requires an expression, named variable, and point", "limit");
    }
    if (expression_depends_on_variable(detail::node(point), variable)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
            "limit point depends on the approaching variable", "limit");
    }
    if (domain != Domain::Real && domain != Domain::Complex) {
        return Result<void>::failure(CasErrc::InvalidArgument, "limit domain must be Real or Complex", "limit");
    }
    if (domain == Domain::Complex && direction != LimitDirection::Both) {
        return Result<void>::failure(CasErrc::InvalidArgument, "complex limits have no ordered sides", "limit");
    }
    return Result<void>::success();
}

LimitResult limit_complex_at_point(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, ComputationContext& context, const FactsQuery& facts) {
    if (has_infinity(detail::node(expression))) { return unknown(); }
    if (!expression_depends_on_variable(detail::node(expression), variable)) {
        return classify_limit_result(detail::node(expression), variable, context, facts, Domain::Complex);
    }
    ComplexContinuityPoint continuity_point{variable, detail::node(point), facts, context};
    if (!elementary_continuous_structure(detail::node(expression), &continuity_point)) { return unknown(); }
    std::vector<NodePtr> nonzero_bases;
    auto domain_proof = complex_punctured_domain(detail::node(expression), variable, facts, context, nonzero_bases);
    if (!domain_proof) { return LimitResult::failure(domain_proof.error()); }
    if (!domain_proof.value()) { return unknown(); }
    ComplexPuncturedFacts punctured(facts, nonzero_bases);
    NormalizationVisitor normalizer(context, punctured, Domain::Complex);
    detail::node(expression)->accept(normalizer);
    auto prepared = normalizer.get_result();
    if (!prepared) { return unknown(); }
    auto substituted = substitute_free(prepared, variable, detail::node(point));
    if (!substituted) { return unknown(); }
    auto defined = detail::query_definedness(substituted, facts, Domain::Complex, context);
    if (!defined) { return LimitResult::failure(defined.error()); }
    if (defined.value() != Tribool::True) { return unknown(); }
    substituted->accept(normalizer);
    auto value = normalizer.get_result();
    return value ? classify_limit_result(value, variable, context, facts, Domain::Complex) : unknown();
}

LimitResult limit_real_at_point(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, int endpoint_sign, ComputationContext& context) {
    if (direction != LimitDirection::Both || endpoint_sign) {
        return limit_one_sided_checked(expression, variable, point, direction, context);
    }
    auto left = limit_one_sided_checked(expression, variable, point, LimitDirection::FromBelow, context);
    if (!left && left.error().code != CasErrc::Inconclusive) { return left; }
    auto right = limit_one_sided_checked(expression, variable, point, LimitDirection::FromAbove, context);
    if (!right && right.error().code != CasErrc::Inconclusive) { return right; }
    if ((left && std::holds_alternative<LimitDoesNotExist>(left.value().value)) ||
        (right && std::holds_alternative<LimitDoesNotExist>(right.value().value))) {
        return outcome(LimitDoesNotExist{});
    }
    if (!left || !right) { return unknown(); }
    return combine_one_sided_limits(left.value().value, right.value().value, context);
}

LimitResult limit_at_checked_point(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, ComputationContext& context,
    const FactsQuery& facts, Domain domain) {
    auto endpoint = classify_limit_result(detail::node(point), variable, context, facts, domain);
    if (!endpoint) { return endpoint; }
    const int endpoint_sign = infinity_sign(endpoint.value().value);
    if (domain == Domain::Complex && endpoint_sign) { return unknown(); }
    auto canonical_point = endpoint_sign ? SymbolicExpr::infinity(endpoint_sign)
        : std::get<FiniteLimit>(endpoint.value().value).value->simplify();
    if (!canonical_point) { return unknown(); }
    if (domain == Domain::Complex) {
        return limit_complex_at_point(expression, variable, canonical_point, context, facts);
    }
    return limit_real_at_point(expression, variable, canonical_point, direction, endpoint_sign, context);
}
}

LimitResult limit_checked(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, ComputationContext& context, Domain domain) {
    auto valid = validate_limit_input(expression, variable, point, direction, domain);
    if (!valid) { return LimitResult::failure(valid.error()); }
    auto step = context.consume_steps(1, "limit");
    if (!step) { return LimitResult::failure(step.error()); }
    try {
        std::optional<detail::AssumptionFacts> assumed;
        if (context.assumptions()) { assumed.emplace(*context.assumptions()); }
        const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
        return limit_at_checked_point(expression, variable, point, direction, context, facts, domain);
    } catch (const CasError& error) {
        return LimitResult::failure(error);
    } catch (const std::bad_alloc&) {
        return LimitResult::failure(CasErrc::ResourceLimit, "limit allocation failed", "limit");
    } catch (const std::exception& error) {
        return LimitResult::failure(CasErrc::InternalInvariant, error.what(), "limit");
    }
}
LimitResult limit_checked(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, Domain domain) {
    ComputationContext context;
    return limit_checked(expression, variable, point, direction, context, domain);
}
LimitExpressionResult limit_expression_checked(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, ComputationContext& context, Domain domain) {
    auto result = limit_checked(expression, variable, point, direction, context, domain);
    if (!result) return LimitExpressionResult::failure(result.error());
    const auto& value = result.value().value;
    if (auto finite = std::get_if<FiniteLimit>(&value)) return LimitExpressionResult::success(finite->value);
    if (int sign = infinity_sign(value)) return LimitExpressionResult::success(SymbolicExpr::infinity(sign));
    return LimitExpressionResult::failure(CasErrc::Inconclusive, "limit does not exist", "limit");
}
LimitExpressionResult limit_expression_checked(const LimitExprPtr& expression, const std::string& variable,
    const LimitExprPtr& point, LimitDirection direction, Domain domain) {
    ComputationContext context;
    return limit_expression_checked(expression, variable, point, direction, context, domain);
}
}
