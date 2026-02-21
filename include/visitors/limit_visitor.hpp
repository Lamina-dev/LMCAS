#pragma once

#include "../symbolic_ast.hpp"
#include "normalization_visitor.hpp"
#include <iostream>
#include <cmath>

class LimitVisitor : public SymbolicVisitor {
    std::string var;
    std::shared_ptr<SymbolicNode> point;
    std::string direction; 
public:
    std::shared_ptr<SymbolicNode> result;

    LimitVisitor(std::string v, std::shared_ptr<SymbolicNode> p, std::string dir = "") 
        : var(std::move(v)), point(std::move(p)), direction(std::move(dir)) {}

    std::shared_ptr<SymbolicNode> get_result() const { return result; }

    void visit(NumberNode& node) override { result = node.clone(); }
    void visit(VariableNode& node) override {
        if (node.name == var) result = point->clone();
        else result = node.clone();
    }

    void visit(AddNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for(auto& op : node.operands) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        NormalizationVisitor norm;
        std::make_shared<AddNode>(new_ops)->accept(norm);
        result = norm.get_result();
    }

    void visit(MultiplyNode& node) override {
        // Try simple substitution first
        std::vector<std::shared_ptr<SymbolicNode>> subs_ops;
        for(auto& op : node.operands) {
            op->accept(*this);
            subs_ops.push_back(result);
        }
        
        NormalizationVisitor norm;
        auto subs_res = std::make_shared<MultiplyNode>(subs_ops);
        subs_res->accept(norm);
        auto final_subs = norm.get_result();
        
        // If it's a valid number/expression (not NaN or Inf), return it?
        // But 0 * Inf is Indeterminate.
        
        // Let's identify Indeterminate forms 0*Inf
        bool has_zero = false;
        bool has_inf = false;
        
        // Check for 0/0 or Inf/Inf forms by splitting num/den
        std::vector<std::shared_ptr<SymbolicNode>> num_nodes, den_nodes;
        for(auto& op : node.operands) {
            if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                if (auto num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                    double e = 0;
                    if (std::holds_alternative<double>(num->value)) e = std::get<double>(num->value);
                    else if (std::holds_alternative<BigInt>(num->value)) e = std::get<BigInt>(num->value).to_double();
                    else if (std::holds_alternative<Rational>(num->value)) e = std::get<Rational>(num->value).to_double();
                    
                    if (e < 0) {
                         // op is 1/base^{-e}
                         auto new_exp = std::make_shared<NumberNode>(-e);
                         den_nodes.push_back(std::make_shared<PowerNode>(pow->base, new_exp));
                         continue;
                    }
                }
            }
            num_nodes.push_back(op);
        }

        if (den_nodes.empty()) {
            result = final_subs;
            return;
        }

        auto N = num_nodes.size() == 1 ? num_nodes[0] : std::make_shared<MultiplyNode>(num_nodes);
        auto D = den_nodes.size() == 1 ? den_nodes[0] : std::make_shared<MultiplyNode>(den_nodes);
        
        // Evaluate N(point) and D(point)
        LimitVisitor sub_vis(var, point, direction);
        N->accept(sub_vis);
        auto val_n = sub_vis.get_result();
        D->accept(sub_vis);
        auto val_d = sub_vis.get_result();
        
        NormalizationVisitor norm2;
        val_n->accept(norm2); val_n = norm2.get_result();
        val_d->accept(norm2); val_d = norm2.get_result();
        
        bool n_zero = val_n->is_zero();
        bool d_zero = val_d->is_zero();
        bool n_inf = is_inf(val_n);
        bool d_inf = is_inf(val_d);

        if ((n_zero && d_zero) || (n_inf && d_inf)) {
            // L'Hopital
            DifferentiationVisitor diff_vis(var);
            N->accept(diff_vis); auto dN = diff_vis.get_result();
            D->accept(diff_vis); auto dD = diff_vis.get_result();
            
            // Limit(dN/dD, var, point)
            std::vector<std::shared_ptr<SymbolicNode>> m_ops = { dN, std::make_shared<PowerNode>(dD, std::make_shared<NumberNode>(BigInt(-1))) };
            auto ratio = std::make_shared<MultiplyNode>(m_ops);
            ratio->accept(*this);
            return;
        }
        
        // k / 0 -> Inf
        if (!n_zero && d_zero) {
            // Sign analysis...
            result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, std::vector<std::shared_ptr<SymbolicNode>>{});
            return;
        }

        result = final_subs;
    }

    void visit(PowerNode& node) override {
        node.base->accept(*this);
        auto b = result;
        node.exponent->accept(*this);
        auto e = result;
        
        // 1^Inf is Indeterminate
        // Handling p^q as exp(q * ln(p))
        if (is_inf(e) && b->is_one()) {
             // ...
        }
        
        NormalizationVisitor norm;
        std::make_shared<PowerNode>(b, e)->accept(norm);
        result = norm.get_result();
    }

    void visit(FunctionNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_args;
        for(auto& a : node.arguments) {
            a->accept(*this);
            new_args.push_back(result);
        }
        NormalizationVisitor norm;
        std::make_shared<FunctionNode>(node.type, new_args)->accept(norm);
        result = norm.get_result();
    }

    void visit(MatrixNode& node) override { result = node.clone(); }

private:
    bool is_inf(const std::shared_ptr<SymbolicNode>& node) {
        if (!node) return false;
        if (auto f = std::dynamic_pointer_cast<FunctionNode>(node)) {
            return f->type == FunctionNode::FuncType::Infinity;
        }
        if (auto m = std::dynamic_pointer_cast<MultiplyNode>(node)) {
            for(auto& op : m->operands) if (is_inf(op)) return true;
        }
        return false;
    }
};
