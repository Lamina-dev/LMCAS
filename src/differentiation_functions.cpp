#include "internal/visitors/differentiation_visitor.hpp"

namespace LMCAS {

namespace {

std::shared_ptr<const SymbolicNode> trigonometric_outer(
    const FunctionNode& node, detail::RewriteBudget* budget) {
    std::shared_ptr<const SymbolicNode> d_outer;
    switch (node.type()) {
    case FunctionNode::FuncType::Sin:
        d_outer = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cos, node.arguments());
        break;
    case FunctionNode::FuncType::Cos:
        d_outer = SymbolicFactory::create_multiply({
            SymbolicFactory::create_number(BigInt(-1)),
            LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sin, node.arguments())
        }, budget);
        break;
    case FunctionNode::FuncType::Tan:
        {
            auto sec = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sec, node.arguments());
            d_outer = LMCAS::detail::make_node<PowerNode>(sec, SymbolicFactory::create_number(BigInt(2)));
        }
        break;
    case FunctionNode::FuncType::Cot:
        {
            auto csc = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Csc, node.arguments());
            auto csc_sq = LMCAS::detail::make_node<PowerNode>(csc, SymbolicFactory::create_number(BigInt(2)));
            d_outer = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)),
                csc_sq
            }, budget);
        }
        break;
    case FunctionNode::FuncType::Sec:
        {
            auto sec = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sec, node.arguments());
            auto tan = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Tan, node.arguments());
            d_outer = SymbolicFactory::create_multiply({sec, tan}, budget);
        }
        break;
    case FunctionNode::FuncType::Csc:
        {
            auto csc = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Csc, node.arguments());
            auto cot = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cot, node.arguments());
            d_outer = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)),
                csc, cot
            }, budget);
        }
        break;
    default:
        return nullptr;
    }
    budget->measure(d_outer);
    return d_outer;
}

std::shared_ptr<const SymbolicNode> inverse_trigonometric_outer(
    const FunctionNode& node, detail::RewriteBudget* budget) {
    const auto& arg = node.arguments()[0];
    std::shared_ptr<const SymbolicNode> d_outer;
    switch (node.type()) {
    case FunctionNode::FuncType::ArcSin:
        {
            auto arg_sq = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
            auto neg_arg_sq = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), arg_sq
            }, budget);
            auto one_minus_sq = SymbolicFactory::create_add({
                SymbolicFactory::create_number(BigInt(1)), neg_arg_sq
            }, budget);
            d_outer = LMCAS::detail::make_node<PowerNode>(
                one_minus_sq, SymbolicFactory::create_number(Rational(-1, 2)));
        }
        break;
    case FunctionNode::FuncType::ArcCos:
        {
            auto arg_sq = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
            auto neg_arg_sq = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), arg_sq
            }, budget);
            auto one_minus_sq = SymbolicFactory::create_add({
                SymbolicFactory::create_number(BigInt(1)), neg_arg_sq
            }, budget);
            auto inv_sqrt = LMCAS::detail::make_node<PowerNode>(
                one_minus_sq, SymbolicFactory::create_number(Rational(-1, 2)));
            d_outer = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), inv_sqrt
            }, budget);
        }
        break;
    case FunctionNode::FuncType::ArcTan:
        {
            auto arg_sq = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
            auto one_plus_sq = SymbolicFactory::create_add({
                SymbolicFactory::create_number(BigInt(1)), arg_sq
            }, budget);
            d_outer = LMCAS::detail::make_node<PowerNode>(
                one_plus_sq, SymbolicFactory::create_number(BigInt(-1)));
        }
        break;
    default:
        return nullptr;
    }
    budget->measure(d_outer);
    return d_outer;
}

std::shared_ptr<const SymbolicNode> hyperbolic_outer(
    const FunctionNode& node, detail::RewriteBudget* budget) {
    std::shared_ptr<const SymbolicNode> d_outer;
    switch (node.type()) {
    case FunctionNode::FuncType::Sinh:
        d_outer = LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Cosh, node.arguments());
        break;
    case FunctionNode::FuncType::Cosh:
        d_outer = LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sinh, node.arguments());
        break;
    case FunctionNode::FuncType::Tanh:
        {
            auto cosh = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Cosh, node.arguments());
            d_outer = LMCAS::detail::make_node<PowerNode>(
                cosh, SymbolicFactory::create_number(BigInt(-2)));
        }
        break;
    default:
        return nullptr;
    }
    budget->measure(d_outer);
    return d_outer;
}

std::shared_ptr<const SymbolicNode> elementary_outer(
    const FunctionNode& node, detail::RewriteBudget* budget) {
    const auto& arg = node.arguments()[0];
    std::shared_ptr<const SymbolicNode> d_outer;
    switch (node.type()) {
    case FunctionNode::FuncType::Exp:
         d_outer = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Exp, node.arguments());
         break;
    case FunctionNode::FuncType::Ln:
         d_outer = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
         break;
    case FunctionNode::FuncType::Sqrt:
         d_outer = SymbolicFactory::create_multiply({
            SymbolicFactory::create_number(Rational(1, 2)),
            LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(Rational(-1, 2)))
         }, budget);
         break;
    case FunctionNode::FuncType::Abs:
        {
            auto abs_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs, node.arguments());
            auto abs_inv = LMCAS::detail::make_node<PowerNode>(abs_arg, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({arg, abs_inv}, budget);
        }
        break;
    default:
        return nullptr;
    }
    budget->measure(d_outer);
    return d_outer;
}

std::shared_ptr<const SymbolicNode> special_outer(
    const FunctionNode& node, detail::RewriteBudget* budget) {
    const auto& arg = node.arguments()[0];
    std::shared_ptr<const SymbolicNode> d_outer;
    switch (node.type()) {
    case FunctionNode::FuncType::LambertW:
        {
            auto W = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::LambertW, node.arguments());
            auto one = SymbolicFactory::create_number(BigInt(1));
            auto one_plus_W = SymbolicFactory::create_add({one, W}, budget);
            auto minus_W = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), W}, budget);
            auto numerator = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Exp,
                std::vector<std::shared_ptr<const SymbolicNode>>{minus_W});
            auto denominator_inverse = LMCAS::detail::make_node<PowerNode>(
                one_plus_W, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({numerator, denominator_inverse}, budget);
        }
        break;
    case FunctionNode::FuncType::Erf:
        {
            auto pi = LMCAS::detail::make_node<VariableNode>("pi", true);
            auto sqrt_pi = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Sqrt,
                std::vector<std::shared_ptr<const SymbolicNode>>{pi});
            auto sqrt_pi_inv = LMCAS::detail::make_node<PowerNode>(sqrt_pi, SymbolicFactory::create_number(BigInt(-1)));
            auto two = SymbolicFactory::create_number(BigInt(2));
            auto neg_one = SymbolicFactory::create_number(BigInt(-1));
            auto arg_sq = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
            auto neg_arg_sq = SymbolicFactory::create_multiply({neg_one, arg_sq}, budget);
            auto exp_term = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Exp,
                std::vector<std::shared_ptr<const SymbolicNode>>{neg_arg_sq});
            d_outer = SymbolicFactory::create_multiply({two, sqrt_pi_inv, exp_term}, budget);
        }
        break;
    case FunctionNode::FuncType::Ei:
        {
            auto exp_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Exp, node.arguments());
            auto arg_inv = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({exp_arg, arg_inv}, budget);
        }
        break;
    case FunctionNode::FuncType::Si:
        {
            auto sin_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Sin, node.arguments());
            auto arg_inv = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({sin_arg, arg_inv}, budget);
        }
        break;
    case FunctionNode::FuncType::Ci:
        {
            auto cos_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Cos, node.arguments());
            auto arg_inv = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({cos_arg, arg_inv}, budget);
        }
        break;
    case FunctionNode::FuncType::Li:
        {
            auto ln_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Ln, node.arguments());
            d_outer = LMCAS::detail::make_node<PowerNode>(ln_arg, SymbolicFactory::create_number(BigInt(-1)));
        }
        break;
    default:
        return nullptr;
    }
    budget->measure(d_outer);
    return d_outer;
}

bool valid_log_base(const std::shared_ptr<const SymbolicNode>& base) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(base);
    return number && number->is_positive() && !number->is_one();
}

std::shared_ptr<const SymbolicNode> logarithmic_outer(
    const FunctionNode& node, detail::RewriteBudget* budget) {
    const auto& arg = node.arguments()[0];
    const auto& base = node.arguments()[1];
    auto ln_base = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Ln,
        std::vector<std::shared_ptr<const SymbolicNode>>{base});
    return LMCAS::detail::make_node<PowerNode>(
        SymbolicFactory::create_multiply({arg, ln_base}, budget),
        SymbolicFactory::create_number(BigInt(-1)));
}

}

void DifferentiationVisitor::differentiate_atan2(const FunctionNode& node) {
    auto* budget = rewrite_budget();
    const auto& y = node.arguments()[0];
    const auto& x = node.arguments()[1];
    y->accept(*this);
    auto dy = result;
    x->accept(*this);
    auto dx = result;
    auto numerator = SymbolicFactory::create_add({
        SymbolicFactory::create_multiply({x, dy}, budget),
        SymbolicFactory::create_multiply({
            SymbolicFactory::create_number(BigInt(-1)), y, dx}, budget)
    }, budget);
    auto two = SymbolicFactory::create_number(BigInt(2));
    auto denominator = SymbolicFactory::create_add({
        LMCAS::detail::make_node<PowerNode>(x, two),
        LMCAS::detail::make_node<PowerNode>(y, two)
    }, budget);
    auto inverse = LMCAS::detail::make_node<PowerNode>(
        denominator, SymbolicFactory::create_number(BigInt(-1)));
    set_result(SymbolicFactory::create_multiply({numerator, inverse}, budget));
}

void DifferentiationVisitor::visit(const FunctionNode& node) {
    auto* budget = rewrite_budget();
    if (node.type() == FunctionNode::FuncType::Atan2 && node.arguments().size() == 2) {
        differentiate_atan2(node);
        return;
    }
    if (node.type() == FunctionNode::FuncType::Log && node.arguments().size() == 2) {
        const auto& arg = node.arguments()[0];
        if (!valid_log_base(node.arguments()[1])) {
            unsupported("FunctionNode");
        }
        arg->accept(*this);
        auto d_arg = result;
        if (d_arg->is_zero()) return;
        auto inverse = logarithmic_outer(node, budget);
        set_result(SymbolicFactory::create_multiply({std::move(d_arg), inverse}, budget));
        return;
    }
    if (node.arguments().size() != 1) unsupported("FunctionNode");
    node.arguments()[0]->accept(*this);
    auto d_arg = result;
    if (d_arg->is_zero()) {
        return;
    }

    auto d_outer = trigonometric_outer(node, budget);
    if (!d_outer) d_outer = inverse_trigonometric_outer(node, budget);
    if (!d_outer) d_outer = hyperbolic_outer(node, budget);
    if (!d_outer) d_outer = elementary_outer(node, budget);
    if (!d_outer) d_outer = special_outer(node, budget);
    if (!d_outer) unsupported("FunctionNode");
    set_result(SymbolicFactory::create_multiply({std::move(d_outer), std::move(d_arg)}, budget));
}

}
