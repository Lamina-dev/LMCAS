#pragma once
#include "../symbolic_ast.hpp"


class DifferentiationVisitor : public SymbolicVisitor {
    std::string var;
public:
    std::shared_ptr<SymbolicNode> result;
    // 隐函数微分相关
    std::string implicit_var; // 例如 y
    bool implicit_mode = false;

    DifferentiationVisitor(const std::string& v) : var(v), result(nullptr) {}
    DifferentiationVisitor(const std::string& v, const std::string& implicit_v)
        : var(v), implicit_var(implicit_v), result(nullptr), implicit_mode(true) {}

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
