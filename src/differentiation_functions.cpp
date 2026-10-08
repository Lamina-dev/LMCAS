#include "internal/visitors/differentiation_visitor.hpp"

namespace LMCAS {

namespace {

std::shared_ptr<const SymbolicNode> trigonometric_outer(const FunctionNode& node) {
    std::shared_ptr<const SymbolicNode> d_outer;
    switch (node.type()) {
    case FunctionNode::FuncType::Sin:
        d_outer = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cos, node.arguments());
        break;
    case FunctionNode::FuncType::Cos:
        d_outer = SymbolicFactory::create_multiply({
            SymbolicFactory::create_number(BigInt(-1)),
            LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sin, node.arguments())
        });
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
            });
        }
        break;
    case FunctionNode::FuncType::Sec:
        {
            auto sec = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sec, node.arguments());
            auto tan = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Tan, node.arguments());
            d_outer = SymbolicFactory::create_multiply({sec, tan});
        }
        break;
    case FunctionNode::FuncType::Csc:
        {
            auto csc = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Csc, node.arguments());
            auto cot = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cot, node.arguments());
            d_outer = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)),
                csc, cot
            });
        }
        break;
    default:
        return nullptr;
    }
    return d_outer;
}

std::shared_ptr<const SymbolicNode> inverse_trigonometric_outer(const FunctionNode& node) {
    const auto& arg = node.arguments()[0];
    std::shared_ptr<const SymbolicNode> d_outer;
    switch (node.type()) {
    case FunctionNode::FuncType::ArcSin:
        {
            auto arg_sq = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
            auto neg_arg_sq = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), arg_sq
            });
            auto one_minus_sq = SymbolicFactory::create_add({
                SymbolicFactory::create_number(BigInt(1)), neg_arg_sq
            });
            d_outer = LMCAS::detail::make_node<PowerNode>(
                one_minus_sq, SymbolicFactory::create_number(Rational(-1, 2)));
        }
        break;
    case FunctionNode::FuncType::ArcCos:
        {
            auto arg_sq = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
            auto neg_arg_sq = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), arg_sq
            });
            auto one_minus_sq = SymbolicFactory::create_add({
                SymbolicFactory::create_number(BigInt(1)), neg_arg_sq
            });
            auto inv_sqrt = LMCAS::detail::make_node<PowerNode>(
                one_minus_sq, SymbolicFactory::create_number(Rational(-1, 2)));
            d_outer = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), inv_sqrt
            });
        }
        break;
    case FunctionNode::FuncType::ArcTan:
        {
            auto arg_sq = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
            auto one_plus_sq = SymbolicFactory::create_add({
                SymbolicFactory::create_number(BigInt(1)), arg_sq
            });
            d_outer = LMCAS::detail::make_node<PowerNode>(
                one_plus_sq, SymbolicFactory::create_number(BigInt(-1)));
        }
        break;
    default:
        return nullptr;
    }
    return d_outer;
}

std::shared_ptr<const SymbolicNode> hyperbolic_outer(const FunctionNode& node) {
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
    return d_outer;
}

std::shared_ptr<const SymbolicNode> elementary_outer(const FunctionNode& node) {
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
         });
         break;
    case FunctionNode::FuncType::Abs:
        {
            auto abs_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs, node.arguments());
            auto abs_inv = LMCAS::detail::make_node<PowerNode>(abs_arg, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({arg, abs_inv});
        }
        break;
    default:
        return nullptr;
    }
    return d_outer;
}

std::shared_ptr<const SymbolicNode> special_outer(const FunctionNode& node) {
    const auto& arg = node.arguments()[0];
    std::shared_ptr<const SymbolicNode> d_outer;
    switch (node.type()) {
    case FunctionNode::FuncType::LambertW:
        {
            auto W = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::LambertW, node.arguments());
            auto one = SymbolicFactory::create_number(BigInt(1));
            auto one_plus_W = SymbolicFactory::create_add({one, W});
            auto minus_W = SymbolicFactory::create_multiply({
                SymbolicFactory::create_number(BigInt(-1)), W});
            auto numerator = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Exp,
                std::vector<std::shared_ptr<const SymbolicNode>>{minus_W});
            auto denominator_inverse = LMCAS::detail::make_node<PowerNode>(
                one_plus_W, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({numerator, denominator_inverse});
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
            auto neg_arg_sq = SymbolicFactory::create_multiply({neg_one, arg_sq});
            auto exp_term = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Exp,
                std::vector<std::shared_ptr<const SymbolicNode>>{neg_arg_sq});
            d_outer = SymbolicFactory::create_multiply({two, sqrt_pi_inv, exp_term});
        }
        break;
    case FunctionNode::FuncType::Ei:
        {
            auto exp_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Exp, node.arguments());
            auto arg_inv = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({exp_arg, arg_inv});
        }
        break;
    case FunctionNode::FuncType::Si:
        {
            auto sin_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Sin, node.arguments());
            auto arg_inv = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({sin_arg, arg_inv});
        }
        break;
    case FunctionNode::FuncType::Ci:
        {
            auto cos_arg = LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Cos, node.arguments());
            auto arg_inv = LMCAS::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
            d_outer = SymbolicFactory::create_multiply({cos_arg, arg_inv});
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
    return d_outer;
}

bool valid_log_base(const std::shared_ptr<const SymbolicNode>& base) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(base);
    return number && number->is_positive() && !number->is_one();
}

std::shared_ptr<const SymbolicNode> logarithmic_outer(const FunctionNode& node) {
    const auto& arg = node.arguments()[0];
    const auto& base = node.arguments()[1];
    auto ln_base = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Ln,
        std::vector<std::shared_ptr<const SymbolicNode>>{base});
    return LMCAS::detail::make_node<PowerNode>(
        SymbolicFactory::create_multiply({arg, ln_base}),
        SymbolicFactory::create_number(BigInt(-1)));
}

}

void DifferentiationVisitor::differentiate_atan2(const FunctionNode& node) {
    const auto& y = node.arguments()[0];
    const auto& x = node.arguments()[1];
    y->accept(*this);
    auto dy = result;
    x->accept(*this);
    auto dx = result;
    auto numerator = SymbolicFactory::create_add({
        SymbolicFactory::create_multiply({x, dy}),
        SymbolicFactory::create_multiply({
            SymbolicFactory::create_number(BigInt(-1)), y, dx})
    });
    auto two = SymbolicFactory::create_number(BigInt(2));
    auto denominator = SymbolicFactory::create_add({
        LMCAS::detail::make_node<PowerNode>(x, two),
        LMCAS::detail::make_node<PowerNode>(y, two)
    });
    auto inverse = LMCAS::detail::make_node<PowerNode>(
        denominator, SymbolicFactory::create_number(BigInt(-1)));
    result = SymbolicFactory::create_multiply({numerator, inverse});
}

void DifferentiationVisitor::visit(const FunctionNode& node) {
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
        auto inverse = logarithmic_outer(node);
        result = SymbolicFactory::create_multiply({std::move(d_arg), inverse});
        return;
    }
    if (node.arguments().size() != 1) unsupported("FunctionNode");
    node.arguments()[0]->accept(*this);
    auto d_arg = result;
    if (d_arg->is_zero()) {
        return;
    }

    auto d_outer = trigonometric_outer(node);
    if (!d_outer) d_outer = inverse_trigonometric_outer(node);
    if (!d_outer) d_outer = hyperbolic_outer(node);
    if (!d_outer) d_outer = elementary_outer(node);
    if (!d_outer) d_outer = special_outer(node);
    if (!d_outer) unsupported("FunctionNode");
    result = SymbolicFactory::create_multiply({std::move(d_outer), std::move(d_arg)});
}

}
