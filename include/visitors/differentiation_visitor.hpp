/**
 * @file differentiation_visitor.hpp
 * @brief 微分访问器，对 AST 执行符号求导（含隐函数模式）。
 */
#pragma once
#include "../symbolic_ast.hpp"
#include <stdexcept>
#include <string>

/** @brief 微分访问器，对符号 AST 执行求导运算，支持普通求导和隐函数求导 */
class DifferentiationVisitor : public lamina::detail::SymbolicVisitor {
    std::string var;

    [[noreturn]] void unsupported(const char* node_type) const {
        throw std::runtime_error(std::string(node_type) +
                                 " differentiation is outside the current supported domain");
    }

public:
    std::shared_ptr<const SymbolicNode> result;  ///< 求导结果节点

    std::string implicit_var;    ///< 隐函数因变量名
    bool implicit_mode = false;  ///< 是否为隐函数求导模式

    /**
     * @brief 构造普通求导访问器
     * @param v 求导变量名
     */
    DifferentiationVisitor(const std::string& v) : var(v), result(nullptr) {}

    /**
     * @brief 构造隐函数求导访问器
     * @param v 自变量名
     * @param implicit_v 隐函数因变量名
     */
    DifferentiationVisitor(const std::string& v, const std::string& implicit_v)
        : var(v), result(nullptr), implicit_var(implicit_v), implicit_mode(true) {}

    /**
     * @brief 获取求导结果
     * @return 求导后的 AST 节点
     */
    std::shared_ptr<const SymbolicNode> get_result() const {
        return result;
    }
    void visit(const NumberNode&) override {
        result = SymbolicFactory::create_number(BigInt(0));
    }

    void visit(const VariableNode& node) override {
        if (node.name() == var) {
            result = SymbolicFactory::create_number(BigInt(1));
        } else {
            result = SymbolicFactory::create_number(BigInt(0));
        }
    }

    void visit(const AddNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> diff_ops;
        for (const auto& op : node.operands()) {
            op->accept(*this);
            diff_ops.push_back(result);
        }
        result = SymbolicFactory::create_add(diff_ops);
    }

    void visit(const MultiplyNode& node) override {
        if (node.operands().empty()) {
             result = SymbolicFactory::create_number(BigInt(0));
             return;
        }

        std::vector<std::shared_ptr<const SymbolicNode>> sum_terms;

        for (size_t i = 0; i < node.operands().size(); ++i) {

            node.operands()[i]->accept(*this);
            auto d_term = result;

            if (d_term->is_zero()) continue;

            std::vector<std::shared_ptr<const SymbolicNode>> prod_terms;

            for (size_t j = 0; j < node.operands().size(); ++j) {
                if (i == j) {
                    prod_terms.push_back(d_term);
                } else {
                    prod_terms.push_back(node.operands()[j]);
                }
            }
            sum_terms.push_back(SymbolicFactory::create_multiply(prod_terms));
        }

        if (sum_terms.empty()) {
             result = SymbolicFactory::create_number(BigInt(0));
        } else {
             result = SymbolicFactory::create_add(sum_terms);
        }
    }

    void visit(const PowerNode& node) override {

        node.base()->accept(*this);
        auto du = result;
        node.exponent()->accept(*this);
        auto dv = result;

        auto u = node.base();
        auto v = node.exponent();

        if (dv->is_zero()) {

            auto n = v;

            auto n_minus_1 = SymbolicFactory::create_add({n, SymbolicFactory::create_number(BigInt(-1))});

            auto u_pow = lamina::detail::make_node<PowerNode>(u, n_minus_1);

            result = SymbolicFactory::create_multiply({n, u_pow, du});
        } else {

            auto ln_u = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, std::vector<std::shared_ptr<const SymbolicNode>>{u});

            auto t1 = SymbolicFactory::create_multiply({dv, ln_u});

            auto u_inv = lamina::detail::make_node<PowerNode>(u, SymbolicFactory::create_number(BigInt(-1)));
            auto t2 = SymbolicFactory::create_multiply({v, du, u_inv});

            auto sum = SymbolicFactory::create_add({t1, t2});

            auto u_pow_v = lamina::detail::make_node<PowerNode>(u, v);
            result = SymbolicFactory::create_multiply({u_pow_v, sum});
        }
    }

    void visit(const FunctionNode& node) override {
        if (node.arguments().size() != 1) {
             result = SymbolicFactory::create_number(BigInt(0));
             return;
        }

        auto& arg = node.arguments()[0];
        arg->accept(*this);
        auto d_arg = result;

        if (d_arg->is_zero()) {
            result = SymbolicFactory::create_number(BigInt(0));
            return;
        }

        std::shared_ptr<const SymbolicNode> d_outer;

        switch (node.type()) {
            case FunctionNode::FuncType::Sin:
                d_outer = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cos, node.arguments());
                break;
            case FunctionNode::FuncType::Cos:
                d_outer = SymbolicFactory::create_multiply({
                    SymbolicFactory::create_number(BigInt(-1)),
                    lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sin, node.arguments())
                });
                break;
            case FunctionNode::FuncType::Tan:
                {
                    auto sec = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sec, node.arguments());
                    d_outer = lamina::detail::make_node<PowerNode>(sec, SymbolicFactory::create_number(BigInt(2)));
                }
                break;
            case FunctionNode::FuncType::Cot:
                {
                    // d/dx cot(x) = -csc(x)^2
                    auto csc = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Csc, node.arguments());
                    auto csc_sq = lamina::detail::make_node<PowerNode>(csc, SymbolicFactory::create_number(BigInt(2)));
                    d_outer = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(BigInt(-1)),
                        csc_sq
                    });
                }
                break;
            case FunctionNode::FuncType::Sec:
                {
                    // d/dx sec(x) = sec(x) * tan(x)
                    auto sec = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sec, node.arguments());
                    auto tan = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Tan, node.arguments());
                    d_outer = SymbolicFactory::create_multiply({sec, tan});
                }
                break;
            case FunctionNode::FuncType::Csc:
                {
                    // d/dx csc(x) = -csc(x) * cot(x)
                    auto csc = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Csc, node.arguments());
                    auto cot = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cot, node.arguments());
                    d_outer = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(BigInt(-1)),
                        csc, cot
                    });
                }
                break;
            case FunctionNode::FuncType::ArcSin:
                {
                    // d/dx arcsin(x) = 1/sqrt(1 - x^2) = (1 - x^2)^(-1/2)
                    auto arg_sq = lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
                    auto neg_arg_sq = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(BigInt(-1)), arg_sq
                    });
                    auto one_minus_sq = SymbolicFactory::create_add({
                        SymbolicFactory::create_number(BigInt(1)), neg_arg_sq
                    });
                    d_outer = lamina::detail::make_node<PowerNode>(
                        one_minus_sq, SymbolicFactory::create_number(Rational(-1, 2)));
                }
                break;
            case FunctionNode::FuncType::ArcCos:
                {
                    // d/dx arccos(x) = -1/sqrt(1 - x^2) = -(1 - x^2)^(-1/2)
                    auto arg_sq = lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
                    auto neg_arg_sq = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(BigInt(-1)), arg_sq
                    });
                    auto one_minus_sq = SymbolicFactory::create_add({
                        SymbolicFactory::create_number(BigInt(1)), neg_arg_sq
                    });
                    auto inv_sqrt = lamina::detail::make_node<PowerNode>(
                        one_minus_sq, SymbolicFactory::create_number(Rational(-1, 2)));
                    d_outer = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(BigInt(-1)), inv_sqrt
                    });
                }
                break;
            case FunctionNode::FuncType::ArcTan:
                {
                    // d/dx arctan(x) = 1/(1 + x^2) = (1 + x^2)^(-1)
                    auto arg_sq = lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
                    auto one_plus_sq = SymbolicFactory::create_add({
                        SymbolicFactory::create_number(BigInt(1)), arg_sq
                    });
                    d_outer = lamina::detail::make_node<PowerNode>(
                        one_plus_sq, SymbolicFactory::create_number(BigInt(-1)));
                }
                break;
            case FunctionNode::FuncType::Sinh:
                d_outer = lamina::detail::make_node<FunctionNode>(
                    FunctionNode::FuncType::Cosh, node.arguments());
                break;
            case FunctionNode::FuncType::Cosh:
                d_outer = lamina::detail::make_node<FunctionNode>(
                    FunctionNode::FuncType::Sinh, node.arguments());
                break;
            case FunctionNode::FuncType::Tanh:
                {
                    // d/dx tanh(x) = sech(x)^2 = 1/cosh(x)^2 = cosh(x)^(-2)
                    auto cosh = lamina::detail::make_node<FunctionNode>(
                        FunctionNode::FuncType::Cosh, node.arguments());
                    d_outer = lamina::detail::make_node<PowerNode>(
                        cosh, SymbolicFactory::create_number(BigInt(-2)));
                }
                break;
            case FunctionNode::FuncType::Exp:
                 d_outer = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Exp, node.arguments());
                 break;
            case FunctionNode::FuncType::Ln:
                 d_outer = lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
                 break;
            case FunctionNode::FuncType::Sqrt:
                 d_outer = SymbolicFactory::create_multiply({
                    SymbolicFactory::create_number(Rational(1, 2)),
                    lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(Rational(-1, 2)))
                 });
                 break;
            case FunctionNode::FuncType::Abs:
                {
                    /// d/dx |u| = u / |u| = sgn(u)  (undefined at u=0)
                    auto abs_arg = lamina::detail::make_node<FunctionNode>(
                        FunctionNode::FuncType::Abs, node.arguments());
                    auto abs_inv = lamina::detail::make_node<PowerNode>(abs_arg, SymbolicFactory::create_number(BigInt(-1)));
                    d_outer = SymbolicFactory::create_multiply({arg, abs_inv});
                }
                break;
            case FunctionNode::FuncType::LambertW:
                {
                    auto W = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::LambertW, node.arguments());
                    auto one = SymbolicFactory::create_number(BigInt(1));
                    auto one_plus_W = SymbolicFactory::create_add({one, W});
                    auto denom = SymbolicFactory::create_multiply({arg, one_plus_W});
                    auto denom_inv = lamina::detail::make_node<PowerNode>(denom, SymbolicFactory::create_number(BigInt(-1)));
                    d_outer = SymbolicFactory::create_multiply({W, denom_inv});
                }
                break;
            case FunctionNode::FuncType::RootOf:
                throw std::runtime_error("RootOf differentiation is not mathematically supported yet.");
            case FunctionNode::FuncType::Erf:
                {
                    // d/dx erf(x) = (2/sqrt(pi)) * exp(-x^2)
                    auto pi = lamina::detail::make_node<VariableNode>("pi");
                    auto sqrt_pi = lamina::detail::make_node<FunctionNode>(
                        FunctionNode::FuncType::Sqrt,
                        std::vector<std::shared_ptr<const SymbolicNode>>{pi});
                    auto sqrt_pi_inv = lamina::detail::make_node<PowerNode>(sqrt_pi, SymbolicFactory::create_number(BigInt(-1)));
                    auto two = SymbolicFactory::create_number(BigInt(2));
                    auto neg_one = SymbolicFactory::create_number(BigInt(-1));
                    auto arg_sq = lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(2)));
                    auto neg_arg_sq = SymbolicFactory::create_multiply({neg_one, arg_sq});
                    auto exp_term = lamina::detail::make_node<FunctionNode>(
                        FunctionNode::FuncType::Exp,
                        std::vector<std::shared_ptr<const SymbolicNode>>{neg_arg_sq});
                    d_outer = SymbolicFactory::create_multiply({two, sqrt_pi_inv, exp_term});
                }
                break;
            case FunctionNode::FuncType::Ei:
                {
                    // d/dx Ei(x) = exp(x) / x
                    auto exp_arg = lamina::detail::make_node<FunctionNode>(
                        FunctionNode::FuncType::Exp, node.arguments());
                    auto arg_inv = lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
                    d_outer = SymbolicFactory::create_multiply({exp_arg, arg_inv});
                }
                break;
            case FunctionNode::FuncType::Si:
                {
                    // d/dx Si(x) = sin(x) / x
                    auto sin_arg = lamina::detail::make_node<FunctionNode>(
                        FunctionNode::FuncType::Sin, node.arguments());
                    auto arg_inv = lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
                    d_outer = SymbolicFactory::create_multiply({sin_arg, arg_inv});
                }
                break;
            case FunctionNode::FuncType::Ci:
                {
                    // d/dx Ci(x) = cos(x) / x
                    auto cos_arg = lamina::detail::make_node<FunctionNode>(
                        FunctionNode::FuncType::Cos, node.arguments());
                    auto arg_inv = lamina::detail::make_node<PowerNode>(arg, SymbolicFactory::create_number(BigInt(-1)));
                    d_outer = SymbolicFactory::create_multiply({cos_arg, arg_inv});
                }
                break;
            case FunctionNode::FuncType::Li:
                {
                    // d/dx Li(x) = 1 / ln(x)
                    auto ln_arg = lamina::detail::make_node<FunctionNode>(
                        FunctionNode::FuncType::Ln, node.arguments());
                    d_outer = lamina::detail::make_node<PowerNode>(ln_arg, SymbolicFactory::create_number(BigInt(-1)));
                }
                break;
            default:
                d_outer = SymbolicFactory::create_number(BigInt(0));
        }

        result = SymbolicFactory::create_multiply(std::vector<std::shared_ptr<const SymbolicNode>>{d_outer, d_arg});
    }

    void visit(const MatrixNode& node) override {

        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
             const auto& dense = std::get<MatrixNode::DenseStorage>(node.storage());
             MatrixNode::DenseStorage new_dense;
             for(const auto& item : dense) {
                 if (item) {
                    item->accept(*this);
                    new_dense.push_back(result);
                 } else {
                    new_dense.push_back(nullptr);
                 }
             }
             result = lamina::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_dense);
        } else {
             const auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage());
             MatrixNode::SparseStorage new_sparse;
             for(const auto& [idx, item] : sparse) {
                 item->accept(*this);
                 if (!result->is_zero()) {
                     new_sparse[idx] = result;
                 }
             }
             result = lamina::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_sparse);
        }
    }

    void visit(const RelationalNode&) override {
        unsupported("RelationalNode");
    }

    void visit(const LogicalNode&) override {
        unsupported("LogicalNode");
    }

    void visit(const ComplexNode& node) override {
        node.real()->accept(*this);
        auto dr = result;
        node.imag()->accept(*this);
        auto di = result;
        result = SymbolicFactory::create_complex(dr, di);
    }

    void visit(const PiecewiseNode& node) override {
        std::vector<PiecewiseNode::Branch> new_branches;
        for (const auto& br : node.branches()) {
            PiecewiseNode::Branch new_br;
            br.expression->accept(*this);
            new_br.expression = result;
            new_br.condition = br.condition;
            new_branches.push_back(new_br);
        }
        std::shared_ptr<const SymbolicNode> new_def = nullptr;
        if (node.default_expr()) {
            node.default_expr()->accept(*this);
            new_def = result;
        }
        result = lamina::detail::make_node<PiecewiseNode>(std::move(new_branches), new_def);
    }

    void visit(const SummationNode&) override {
        unsupported("SummationNode");
    }

    void visit(const ProductNode&) override {
        unsupported("ProductNode");
    }

    void visit(const TransformNode&) override {
        unsupported("TransformNode");
    }

    void visit(const QuantifierNode&) override {
        unsupported("QuantifierNode");
    }

    void visit(const SetBuilderNode&) override {
        unsupported("SetBuilderNode");
    }
    void visit(const FiniteSetNode&) override { unsupported("FiniteSetNode"); }
    void visit(const IntervalNode&) override { unsupported("IntervalNode"); }
    void visit(const MembershipNode&) override { unsupported("MembershipNode"); }
    void visit(const QuantityNode& node) override {
        node.value()->accept(*this);
        result = lamina::detail::make_node<QuantityNode>(
            result, node.dimension(), node.scale_to_base(), node.display_unit());
    }
};
