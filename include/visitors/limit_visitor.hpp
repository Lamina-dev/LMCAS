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

        std::vector<std::shared_ptr<SymbolicNode>> subs_ops;
        for(auto& op : node.operands) {
            op->accept(*this);
            subs_ops.push_back(result);
        }

        NormalizationVisitor norm;
        auto subs_res = std::make_shared<MultiplyNode>(subs_ops);
        subs_res->accept(norm);
        auto final_subs = norm.get_result();

        bool has_zero = false;
        bool has_inf = false;

        std::vector<std::shared_ptr<SymbolicNode>> num_nodes, den_nodes;
        for(auto& op : node.operands) {
            if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                if (auto num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                    double e = 0;
                    if (std::holds_alternative<double>(num->value)) e = std::get<double>(num->value);
                    else if (std::holds_alternative<BigInt>(num->value)) e = std::get<BigInt>(num->value).to_double();
                    else if (std::holds_alternative<Rational>(num->value)) e = std::get<Rational>(num->value).to_double();

                    if (e < 0) {

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

            DifferentiationVisitor diff_vis(var);
            N->accept(diff_vis); auto dN = diff_vis.get_result();
            D->accept(diff_vis); auto dD = diff_vis.get_result();

            std::vector<std::shared_ptr<SymbolicNode>> m_ops = { dN, std::make_shared<PowerNode>(dD, std::make_shared<NumberNode>(BigInt(-1))) };
            auto ratio = std::make_shared<MultiplyNode>(m_ops);
            ratio->accept(*this);
            return;
        }

        if (!n_zero && d_zero) {

            int sign_n = 1;
            if (auto num = std::dynamic_pointer_cast<NumberNode>(val_n)) {
                 if (std::holds_alternative<double>(num->value)) sign_n = std::get<double>(num->value) > 0 ? 1 : -1;
                 else if (std::holds_alternative<BigInt>(num->value)) sign_n = std::get<BigInt>(num->value) > BigInt(0) ? 1 : -1;
                 else if (std::holds_alternative<Rational>(num->value)) sign_n = std::get<Rational>(num->value) > Rational(0) ? 1 : -1;
            }

            DifferentiationVisitor diff_vis(var);
            std::shared_ptr<SymbolicNode> curr_d = D;
            int sign_d = 0;

            for(int i=1; i<=3; ++i) {
                 curr_d->accept(diff_vis);
                 auto deriv = diff_vis.get_result();
                 if(!deriv) break;

                 LimitVisitor sub_vis(var, point, direction);
                 deriv->accept(sub_vis);
                 auto val_deriv = sub_vis.get_result();

                 NormalizationVisitor norm;
                 val_deriv->accept(norm);
                 val_deriv = norm.get_result();

                 if (!val_deriv->is_zero()) {
                      int s_deriv = 1;
                      if (auto num = std::dynamic_pointer_cast<NumberNode>(val_deriv)) {
                           if (std::holds_alternative<double>(num->value)) s_deriv = std::get<double>(num->value) > 0 ? 1 : -1;
                           else if (std::holds_alternative<BigInt>(num->value)) s_deriv = std::get<BigInt>(num->value) > BigInt(0) ? 1 : -1;
                           else if (std::holds_alternative<Rational>(num->value)) s_deriv = std::get<Rational>(num->value) > Rational(0) ? 1 : -1;
                      }

                      if (i % 2 != 0) {
                          if (direction == "-") s_deriv *= -1;
                      }

                      sign_d = s_deriv;
                      break;
                 }
                 curr_d = deriv;
            }
            if (sign_d == 0) sign_d = (direction == "-" ? -1 : 1);

            int final_sign = sign_n * sign_d;
            std::vector<std::shared_ptr<SymbolicNode>> inf_args;
            auto inf_node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);

            if (final_sign < 0) {
                std::vector<std::shared_ptr<SymbolicNode>> m_args = { std::make_shared<NumberNode>(BigInt(-1)), inf_node };
                result = std::make_shared<MultiplyNode>(m_args);
            } else {
                result = inf_node;
            }
            return;
        }

        result = final_subs;
    }

    void visit(PowerNode& node) override {
        node.base->accept(*this);
        auto b = result;
        node.exponent->accept(*this);
        auto e = result;

        if (b->is_zero()) {
             bool neg_exp = false;
             bool odd_exp = false;
             if (auto num = std::dynamic_pointer_cast<NumberNode>(e)) {
                 if (std::holds_alternative<BigInt>(num->value)) {
                      BigInt v = std::get<BigInt>(num->value);
                      if (v < BigInt(0)) {
                           neg_exp = true;
                           if (v.is_odd()) odd_exp = true;
                      }
                 } else if (std::holds_alternative<double>(num->value)) {
                      if (std::get<double>(num->value) < 0) neg_exp = true;
                 }
             }

             if (neg_exp) {
                  DifferentiationVisitor diff_vis(var);
                  std::shared_ptr<SymbolicNode> curr_base = node.base;
                  int sign_base = 0;

                  for(int i=1; i<=3; ++i) {
                       curr_base->accept(diff_vis);
                       auto deriv = diff_vis.get_result();
                       if (!deriv) break;

                       LimitVisitor sub(var, point, direction);
                       deriv->accept(sub);
                       auto val = sub.get_result();

                       NormalizationVisitor n; val->accept(n); val = n.get_result();

                       if (!val->is_zero()) {
                            int s = 1;
                            if (auto nval = std::dynamic_pointer_cast<NumberNode>(val)) {
                                 if (std::holds_alternative<double>(nval->value)) s = std::get<double>(nval->value) > 0 ? 1 : -1;
                                 else if (std::holds_alternative<BigInt>(nval->value)) s = std::get<BigInt>(nval->value) > BigInt(0) ? 1 : -1;
                                 else if (std::holds_alternative<Rational>(nval->value)) s = std::get<Rational>(nval->value) > Rational(0) ? 1 : -1;
                            }

                            if (i % 2 != 0 && direction == "-") s *= -1;
                            sign_base = s;
                            break;
                       }
                       curr_base = deriv;
                  }
                  if (sign_base == 0) sign_base = (direction == "-" ? -1 : 1);

                  int final_sign = 1;
                  if (odd_exp) final_sign = sign_base;

                  std::vector<std::shared_ptr<SymbolicNode>> inf_args;
                  auto inf_node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);

                  if (final_sign < 0) {
                        std::vector<std::shared_ptr<SymbolicNode>> m = {std::make_shared<NumberNode>(BigInt(-1)), inf_node};
                        result = std::make_shared<MultiplyNode>(m);
                  } else {
                        result = inf_node;
                  }
                  return;
             }
        }

        if (is_inf(e) && b->is_one()) {

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
