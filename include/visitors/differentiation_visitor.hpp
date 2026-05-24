/**
 * @file differentiation_visitor.hpp
 * @brief 微分访问器，对 AST 执行符号求导（含隐函数模式）。
 */
#pragma once
#include "../symbolic_ast.hpp"

/** @brief 微分访问器，对符号 AST 执行求导运算，支持普通求导和隐函数求导 */
class DifferentiationVisitor : public SymbolicVisitor {
    std::string var;
public:
    std::shared_ptr<SymbolicNode> result;  ///< 求导结果节点

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
        : var(v), implicit_var(implicit_v), result(nullptr), implicit_mode(true) {}

    /**
     * @brief 获取求导结果
     * @return 求导后的 AST 节点
     */
    std::shared_ptr<SymbolicNode> get_result() const {
        return result;
    }
    void visit(NumberNode& node) override {
        result = SymbolicFactory::create_number(0.0);
    }

    void visit(VariableNode& node) override {
        if (node.name == var) {
            result = SymbolicFactory::create_number(1.0);
        } else {
            result = SymbolicFactory::create_number(0.0);
        }
    }

    void visit(AddNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> diff_ops;
        for (const auto& op : node.operands) {
            op->accept(*this);
            diff_ops.push_back(result);
        }
        result = SymbolicFactory::create_add(diff_ops);
    }

    void visit(MultiplyNode& node) override {
        if (node.operands.empty()) {
             result = SymbolicFactory::create_number(0.0);
             return;
        }

        std::vector<std::shared_ptr<SymbolicNode>> sum_terms;

        for (size_t i = 0; i < node.operands.size(); ++i) {

            node.operands[i]->accept(*this);
            auto d_term = result;

            if (d_term->is_zero()) continue;

            std::vector<std::shared_ptr<SymbolicNode>> prod_terms;

            for (size_t j = 0; j < node.operands.size(); ++j) {
                if (i == j) {
                    prod_terms.push_back(d_term);
                } else {
                    prod_terms.push_back(node.operands[j]);
                }
            }
            sum_terms.push_back(SymbolicFactory::create_multiply(prod_terms));
        }

        if (sum_terms.empty()) {
             result = SymbolicFactory::create_number(0.0);
        } else {
             result = SymbolicFactory::create_add(sum_terms);
        }
    }

    void visit(PowerNode& node) override {

        node.base->accept(*this);
        auto du = result;
        node.exponent->accept(*this);
        auto dv = result;

        auto u = node.base;
        auto v = node.exponent;

        if (dv->is_zero()) {

            auto n = v;

            auto n_minus_1 = SymbolicFactory::create_add({n, SymbolicFactory::create_number(-1.0)});

            auto u_pow = std::make_shared<PowerNode>(u, n_minus_1);

            result = SymbolicFactory::create_multiply({n, u_pow, du});
        } else {

            auto ln_u = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, std::vector<std::shared_ptr<SymbolicNode>>{u});

            auto t1 = SymbolicFactory::create_multiply({dv, ln_u});

            auto u_inv = std::make_shared<PowerNode>(u, SymbolicFactory::create_number(-1.0));
            auto t2 = SymbolicFactory::create_multiply({v, du, u_inv});

            auto sum = SymbolicFactory::create_add({t1, t2});

            auto u_pow_v = std::make_shared<PowerNode>(u, v);
            result = SymbolicFactory::create_multiply({u_pow_v, sum});
        }
    }

    void visit(FunctionNode& node) override {
        if (node.arguments.size() != 1) {
             result = SymbolicFactory::create_number(0.0);
             return;
        }

        auto& arg = node.arguments[0];
        arg->accept(*this);
        auto d_arg = result;

        if (d_arg->is_zero()) {
            result = SymbolicFactory::create_number(0.0);
            return;
        }

        std::shared_ptr<SymbolicNode> d_outer;

        switch (node.type) {
            case FunctionNode::FuncType::Sin:
                d_outer = std::make_shared<FunctionNode>(FunctionNode::FuncType::Cos, node.arguments);
                break;
            case FunctionNode::FuncType::Cos:
                d_outer = SymbolicFactory::create_multiply({
                    SymbolicFactory::create_number(-1.0),
                    std::make_shared<FunctionNode>(FunctionNode::FuncType::Sin, node.arguments)
                });
                break;
            case FunctionNode::FuncType::Tan:
                {
                    auto sec = std::make_shared<FunctionNode>(FunctionNode::FuncType::Sec, node.arguments);
                    d_outer = std::make_shared<PowerNode>(sec, SymbolicFactory::create_number(2.0));
                }
                break;
            case FunctionNode::FuncType::Cot:
                {
                    // d/dx cot(x) = -csc(x)^2
                    auto csc = std::make_shared<FunctionNode>(FunctionNode::FuncType::Csc, node.arguments);
                    auto csc_sq = std::make_shared<PowerNode>(csc, SymbolicFactory::create_number(2.0));
                    d_outer = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(-1.0),
                        csc_sq
                    });
                }
                break;
            case FunctionNode::FuncType::Sec:
                {
                    // d/dx sec(x) = sec(x) * tan(x)
                    auto sec = std::make_shared<FunctionNode>(FunctionNode::FuncType::Sec, node.arguments);
                    auto tan = std::make_shared<FunctionNode>(FunctionNode::FuncType::Tan, node.arguments);
                    d_outer = SymbolicFactory::create_multiply({sec, tan});
                }
                break;
            case FunctionNode::FuncType::Csc:
                {
                    // d/dx csc(x) = -csc(x) * cot(x)
                    auto csc = std::make_shared<FunctionNode>(FunctionNode::FuncType::Csc, node.arguments);
                    auto cot = std::make_shared<FunctionNode>(FunctionNode::FuncType::Cot, node.arguments);
                    d_outer = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(-1.0),
                        csc, cot
                    });
                }
                break;
            case FunctionNode::FuncType::ArcSin:
                {
                    // d/dx arcsin(x) = 1/sqrt(1 - x^2) = (1 - x^2)^(-1/2)
                    auto arg_sq = std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(2.0));
                    auto neg_arg_sq = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(-1.0), arg_sq
                    });
                    auto one_minus_sq = SymbolicFactory::create_add({
                        SymbolicFactory::create_number(1.0), neg_arg_sq
                    });
                    d_outer = std::make_shared<PowerNode>(
                        one_minus_sq, SymbolicFactory::create_number(-0.5));
                }
                break;
            case FunctionNode::FuncType::ArcCos:
                {
                    // d/dx arccos(x) = -1/sqrt(1 - x^2) = -(1 - x^2)^(-1/2)
                    auto arg_sq = std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(2.0));
                    auto neg_arg_sq = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(-1.0), arg_sq
                    });
                    auto one_minus_sq = SymbolicFactory::create_add({
                        SymbolicFactory::create_number(1.0), neg_arg_sq
                    });
                    auto inv_sqrt = std::make_shared<PowerNode>(
                        one_minus_sq, SymbolicFactory::create_number(-0.5));
                    d_outer = SymbolicFactory::create_multiply({
                        SymbolicFactory::create_number(-1.0), inv_sqrt
                    });
                }
                break;
            case FunctionNode::FuncType::ArcTan:
                {
                    // d/dx arctan(x) = 1/(1 + x^2) = (1 + x^2)^(-1)
                    auto arg_sq = std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(2.0));
                    auto one_plus_sq = SymbolicFactory::create_add({
                        SymbolicFactory::create_number(1.0), arg_sq
                    });
                    d_outer = std::make_shared<PowerNode>(
                        one_plus_sq, SymbolicFactory::create_number(-1.0));
                }
                break;
            case FunctionNode::FuncType::Sinh:
                d_outer = std::make_shared<FunctionNode>(
                    FunctionNode::FuncType::Cosh, node.arguments);
                break;
            case FunctionNode::FuncType::Cosh:
                d_outer = std::make_shared<FunctionNode>(
                    FunctionNode::FuncType::Sinh, node.arguments);
                break;
            case FunctionNode::FuncType::Tanh:
                {
                    // d/dx tanh(x) = sech(x)^2 = 1/cosh(x)^2 = cosh(x)^(-2)
                    auto cosh = std::make_shared<FunctionNode>(
                        FunctionNode::FuncType::Cosh, node.arguments);
                    d_outer = std::make_shared<PowerNode>(
                        cosh, SymbolicFactory::create_number(-2.0));
                }
                break;
            case FunctionNode::FuncType::Exp:
                 d_outer = std::make_shared<FunctionNode>(FunctionNode::FuncType::Exp, node.arguments);
                 break;
            case FunctionNode::FuncType::Ln:
                 d_outer = std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(-1.0));
                 break;
            case FunctionNode::FuncType::Sqrt:
                 d_outer = SymbolicFactory::create_multiply({
                    SymbolicFactory::create_number(0.5),
                    std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(-0.5))
                 });
                 break;
            case FunctionNode::FuncType::LambertW:
                {
                    auto W = std::make_shared<FunctionNode>(FunctionNode::FuncType::LambertW, node.arguments);
                    auto one = SymbolicFactory::create_number(1.0);
                    auto one_plus_W = SymbolicFactory::create_add({one, W});
                    auto denom = SymbolicFactory::create_multiply({arg, one_plus_W});
                    auto denom_inv = std::make_shared<PowerNode>(denom, SymbolicFactory::create_number(-1.0));
                    d_outer = SymbolicFactory::create_multiply({W, denom_inv});
                }
                break;
            case FunctionNode::FuncType::RootOf:
                throw std::runtime_error("RootOf differentiation is not mathematically supported yet.");
            case FunctionNode::FuncType::Erf:
                {
                    // d/dx erf(x) = (2/sqrt(pi)) * exp(-x^2)
                    auto pi = std::make_shared<VariableNode>("pi");
                    auto sqrt_pi = std::make_shared<FunctionNode>(
                        FunctionNode::FuncType::Sqrt,
                        std::vector<std::shared_ptr<SymbolicNode>>{pi});
                    auto sqrt_pi_inv = std::make_shared<PowerNode>(sqrt_pi, SymbolicFactory::create_number(-1.0));
                    auto two = SymbolicFactory::create_number(2.0);
                    auto neg_one = SymbolicFactory::create_number(-1.0);
                    auto arg_sq = std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(2.0));
                    auto neg_arg_sq = SymbolicFactory::create_multiply({neg_one, arg_sq});
                    auto exp_term = std::make_shared<FunctionNode>(
                        FunctionNode::FuncType::Exp,
                        std::vector<std::shared_ptr<SymbolicNode>>{neg_arg_sq});
                    d_outer = SymbolicFactory::create_multiply({two, sqrt_pi_inv, exp_term});
                }
                break;
            case FunctionNode::FuncType::Ei:
                {
                    // d/dx Ei(x) = exp(x) / x
                    auto exp_arg = std::make_shared<FunctionNode>(
                        FunctionNode::FuncType::Exp, node.arguments);
                    auto arg_inv = std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(-1.0));
                    d_outer = SymbolicFactory::create_multiply({exp_arg, arg_inv});
                }
                break;
            case FunctionNode::FuncType::Si:
                {
                    // d/dx Si(x) = sin(x) / x
                    auto sin_arg = std::make_shared<FunctionNode>(
                        FunctionNode::FuncType::Sin, node.arguments);
                    auto arg_inv = std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(-1.0));
                    d_outer = SymbolicFactory::create_multiply({sin_arg, arg_inv});
                }
                break;
            case FunctionNode::FuncType::Ci:
                {
                    // d/dx Ci(x) = cos(x) / x
                    auto cos_arg = std::make_shared<FunctionNode>(
                        FunctionNode::FuncType::Cos, node.arguments);
                    auto arg_inv = std::make_shared<PowerNode>(arg, SymbolicFactory::create_number(-1.0));
                    d_outer = SymbolicFactory::create_multiply({cos_arg, arg_inv});
                }
                break;
            case FunctionNode::FuncType::Li:
                {
                    // d/dx Li(x) = 1 / ln(x)
                    auto ln_arg = std::make_shared<FunctionNode>(
                        FunctionNode::FuncType::Ln, node.arguments);
                    d_outer = std::make_shared<PowerNode>(ln_arg, SymbolicFactory::create_number(-1.0));
                }
                break;
            default:
                d_outer = SymbolicFactory::create_number(0.0);
        }

        result = SymbolicFactory::create_multiply({d_outer, d_arg});
    }

    void visit(MatrixNode& node) override {

        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage)) {
             const auto& dense = std::get<MatrixNode::DenseStorage>(node.storage);
             MatrixNode::DenseStorage new_dense;
             for(const auto& item : dense) {
                 if (item) {
                    item->accept(*this);
                    new_dense.push_back(result);
                 } else {
                    new_dense.push_back(nullptr);
                 }
             }
             result = std::make_shared<MatrixNode>(node.rows, node.cols, new_dense);
        } else {
             const auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage);
             MatrixNode::SparseStorage new_sparse;
             for(const auto& [idx, item] : sparse) {
                 item->accept(*this);
                 if (!result->is_zero()) {
                     new_sparse[idx] = result;
                 }
             }
             result = std::make_shared<MatrixNode>(node.rows, node.cols, new_sparse);
        }
    }
};
