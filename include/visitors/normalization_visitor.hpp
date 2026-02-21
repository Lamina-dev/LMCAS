#pragma once

#include "../symbolic_ast.hpp"
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>


inline int get_node_degree_helper(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return 0;
    if (std::dynamic_pointer_cast<VariableNode>(node)) return 1;
    if (auto p = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (auto e = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
             if (std::holds_alternative<BigInt>(e->value)) return (int)std::get<BigInt>(e->value).to_int();
             if (std::holds_alternative<double>(e->value)) return (int)std::get<double>(e->value);
        }
        return 1;
    }
    if (auto m = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        int d = 0;
        for (auto& op : m->operands) d += get_node_degree_helper(op);
        return d;
    }
    return 0;
}

struct NodeCompare {
    bool operator()(const std::shared_ptr<SymbolicNode>& lhs, const std::shared_ptr<SymbolicNode>& rhs) const {
        if (!lhs && !rhs) return false;
        if (!lhs) return true;
        if (!rhs) return false;

        // Numbers always come last in additions (degree 0)
        int d1 = get_node_degree_helper(lhs);
        int d2 = get_node_degree_helper(rhs);
        if (d1 != d2) return d1 > d2; 

        // If same degree, numbers first or something
        bool isNum1 = std::dynamic_pointer_cast<NumberNode>(lhs) != nullptr;
        bool isNum2 = std::dynamic_pointer_cast<NumberNode>(rhs) != nullptr;
        if (isNum1 != isNum2) return isNum1; 

        if (lhs->type_priority() != rhs->type_priority()) {
            return lhs->type_priority() < rhs->type_priority();
        }

        return lhs->compare(*rhs) < 0;
    }
};


inline std::shared_ptr<NumberNode> add_numbers(const std::shared_ptr<NumberNode>& a, const std::shared_ptr<NumberNode>& b) {
     if (std::holds_alternative<double>(a->value) || std::holds_alternative<double>(b->value)) {
         double v1 = std::holds_alternative<double>(a->value) ? std::get<double>(a->value) : 
                     (std::holds_alternative<Rational>(a->value) ? std::get<Rational>(a->value).to_double() : std::get<BigInt>(a->value).to_double());
         double v2 = std::holds_alternative<double>(b->value) ? std::get<double>(b->value) : 
                     (std::holds_alternative<Rational>(b->value) ? std::get<Rational>(b->value).to_double() : std::get<BigInt>(b->value).to_double());
         return std::make_shared<NumberNode>(v1 + v2);
     }
     
     if (std::holds_alternative<Rational>(a->value) || std::holds_alternative<Rational>(b->value)) {
         Rational r1 = std::holds_alternative<Rational>(a->value) ? std::get<Rational>(a->value) : 
                       (std::holds_alternative<BigInt>(a->value) ? Rational(std::get<BigInt>(a->value)) : Rational(0));
         Rational r2 = std::holds_alternative<Rational>(b->value) ? std::get<Rational>(b->value) : 
                       (std::holds_alternative<BigInt>(b->value) ? Rational(std::get<BigInt>(b->value)) : Rational(0));
         return std::make_shared<NumberNode>(r1 + r2);
     }
     
     BigInt i1 = std::get<BigInt>(a->value);
     BigInt i2 = std::get<BigInt>(b->value);
     return std::make_shared<NumberNode>(i1 + i2);
}

inline std::shared_ptr<NumberNode> multiply_numbers(const std::shared_ptr<NumberNode>& a, const std::shared_ptr<NumberNode>& b) {
     if (std::holds_alternative<double>(a->value) || std::holds_alternative<double>(b->value)) {
         double v1 = std::holds_alternative<double>(a->value) ? std::get<double>(a->value) : 
                     (std::holds_alternative<Rational>(a->value) ? std::get<Rational>(a->value).to_double() : std::get<BigInt>(a->value).to_double());
         double v2 = std::holds_alternative<double>(b->value) ? std::get<double>(b->value) : 
                     (std::holds_alternative<Rational>(b->value) ? std::get<Rational>(b->value).to_double() : std::get<BigInt>(b->value).to_double());
         return std::make_shared<NumberNode>(v1 * v2);
     }
     
     if (std::holds_alternative<Rational>(a->value) || std::holds_alternative<Rational>(b->value)) {
         Rational r1 = std::holds_alternative<Rational>(a->value) ? std::get<Rational>(a->value) : 
                       (std::holds_alternative<BigInt>(a->value) ? Rational(std::get<BigInt>(a->value)) : Rational(1));
         Rational r2 = std::holds_alternative<Rational>(b->value) ? std::get<Rational>(b->value) : 
                       (std::holds_alternative<BigInt>(b->value) ? Rational(std::get<BigInt>(b->value)) : Rational(1));
         return std::make_shared<NumberNode>(r1 * r2);
     }
     
     BigInt i1 = std::get<BigInt>(a->value);
     BigInt i2 = std::get<BigInt>(b->value);
     return std::make_shared<NumberNode>(i1 * i2);
}

class NormalizationVisitor : public SymbolicVisitor {
public:
    std::shared_ptr<SymbolicNode> result;

    std::shared_ptr<SymbolicNode> get_result() const {
        return result;
    }

    std::shared_ptr<SymbolicNode> expand_product(const std::shared_ptr<SymbolicNode>& lhs, const std::shared_ptr<SymbolicNode>& rhs) {
        
        auto add_lhs = std::dynamic_pointer_cast<AddNode>(lhs);
        auto add_rhs = std::dynamic_pointer_cast<AddNode>(rhs);

        if (add_lhs && add_rhs) {
            
            std::vector<std::shared_ptr<SymbolicNode>> new_terms;
            for (const auto& op1 : add_lhs->operands) {
                for (const auto& op2 : add_rhs->operands) {
                    auto prod = expand_product(op1, op2);
                    new_terms.push_back(prod);
                }
            }
            
            return std::make_shared<AddNode>(new_terms);
        } else if (add_lhs) {
            
            std::vector<std::shared_ptr<SymbolicNode>> new_terms;
            for (const auto& op : add_lhs->operands) {
                 auto prod = expand_product(op, rhs);
                 new_terms.push_back(prod);
            }
            return std::make_shared<AddNode>(new_terms);
        } else if (add_rhs) {
             
            std::vector<std::shared_ptr<SymbolicNode>> new_terms;
            for (const auto& op : add_rhs->operands) {
                 auto prod = expand_product(lhs, op);
                 new_terms.push_back(prod);
            }
            return std::make_shared<AddNode>(new_terms);
        }

        
        
        
        
        std::shared_ptr<NumberNode> const_acc = std::make_shared<NumberNode>(BigInt(1));
        std::map<std::shared_ptr<SymbolicNode>, std::shared_ptr<NumberNode>, NodeCompare> bases;
        
        
        auto process_factor = [&](const std::shared_ptr<SymbolicNode>& factor) {
             if (auto num = std::dynamic_pointer_cast<NumberNode>(factor)) {
                 if (num->is_zero()) return false; 
                 const_acc = multiply_numbers(const_acc, num);
             } else {
                 std::shared_ptr<SymbolicNode> base = factor;
                 std::shared_ptr<NumberNode> exp = std::make_shared<NumberNode>(BigInt(1));

                 if (auto pow = std::dynamic_pointer_cast<PowerNode>(factor)) {
                     base = pow->base;
                     if (auto e_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                         exp = e_num;
                     } 
                 }
                 
                 
                 
                 
                 
                 
                 
                 
                 auto it = bases.find(base);
                 if (it == bases.end()) {
                     bases[base] = exp;
                 } else {
                     bases[base] = add_numbers(it->second, exp);
                 }
             }
             return true;
        };

        
        
        
        auto flatten_and_process = [&](const std::shared_ptr<SymbolicNode>& node) {
            if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
                for(const auto& op : mul->operands) {
                    if (!process_factor(op)) return false;
                }
            } else {
                if (!process_factor(node)) return false;
            }
            return true;
        };

        if (!flatten_and_process(lhs)) return std::make_shared<NumberNode>(BigInt(0));
        if (!flatten_and_process(rhs)) return std::make_shared<NumberNode>(BigInt(0));

        std::vector<std::shared_ptr<SymbolicNode>> final_ops;
        if (!const_acc->is_one()) {
             final_ops.push_back(const_acc);
        }
        
        std::vector<std::shared_ptr<SymbolicNode>> var_ops;
        for (auto const& [base, exp] : bases) {
            if (exp->is_zero()) {
            } else if (exp->is_one()) {
                var_ops.push_back(base);
            } else {
                var_ops.push_back(std::make_shared<PowerNode>(base, exp));
            }
        }
        
        // Sort variable ops
        std::sort(var_ops.begin(), var_ops.end(), [](const std::shared_ptr<SymbolicNode>& l, const std::shared_ptr<SymbolicNode>& r) {
            int d1 = get_node_degree_helper(l);
            int d2 = get_node_degree_helper(r);
            if (d1 != d2) return d1 > d2;
            return l->compare(*r) < 0;
        });

        final_ops.insert(final_ops.end(), var_ops.begin(), var_ops.end());
        
        if (final_ops.empty()) return std::make_shared<NumberNode>(BigInt(1));
        if (final_ops.size() == 1) return final_ops[0];
        return std::make_shared<MultiplyNode>(final_ops);
    }

    void visit(NumberNode& node) override {
        result = node.clone();
    }

    void visit(VariableNode& node) override {
        result = node.clone();
    }

    void visit(AddNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> simplified_ops;
        
        
        for (const auto& op : node.operands) {
            op->accept(*this);
            if (auto add = std::dynamic_pointer_cast<AddNode>(result)) {
                simplified_ops.insert(simplified_ops.end(), add->operands.begin(), add->operands.end());
            } else {
                simplified_ops.push_back(result);
            }
        }

        std::shared_ptr<NumberNode> constant_acc = std::make_shared<NumberNode>(BigInt(0));
        std::map<std::shared_ptr<SymbolicNode>, std::shared_ptr<NumberNode>, NodeCompare> terms;

        for (const auto& op : simplified_ops) {
            if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
                constant_acc = add_numbers(constant_acc, num);
            } else {
                std::shared_ptr<SymbolicNode> term_part = op;
                std::shared_ptr<NumberNode> coeff_part = std::make_shared<NumberNode>(BigInt(1));

                if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {
                    
                    std::shared_ptr<NumberNode> coeff = nullptr;
                    int coeff_idx = -1;
                    
                    if (!mul->operands.empty()) {
                        if (auto n_first = std::dynamic_pointer_cast<NumberNode>(mul->operands.front())) {
                            coeff = n_first;
                            coeff_idx = 0;
                        } else if (auto n_last = std::dynamic_pointer_cast<NumberNode>(mul->operands.back())) {
                            coeff = n_last;
                            coeff_idx = (int)mul->operands.size() - 1;
                        }
                    }

                    if (coeff) {
                         coeff_part = coeff;
                         if (mul->operands.size() == 2) {
                             term_part = mul->operands[coeff_idx == 0 ? 1 : 0];
                         } else {
                             std::vector<std::shared_ptr<SymbolicNode>> rest;
                             rest.reserve(mul->operands.size() - 1);
                             for(int k=0; k<(int)mul->operands.size(); ++k) {
                                 if (k != coeff_idx) rest.push_back(mul->operands[k]);
                             }
                             
                             term_part = std::make_shared<MultiplyNode>(rest);
                         }
                    }
                }

                auto it = terms.find(term_part);
                if (it == terms.end()) {
                    terms[term_part] = coeff_part;
                } else {
                    terms[term_part] = add_numbers(it->second, coeff_part);
                }
            }
        }

        std::vector<std::shared_ptr<SymbolicNode>> final_ops;
        
        if (!constant_acc->is_zero()) {
            final_ops.push_back(constant_acc);
        }

        for (auto const& [term, coeff] : terms) {
            if (coeff->is_zero()) continue;
            
            if (coeff->is_one()) {
                final_ops.push_back(term);
            } else {
                std::vector<std::shared_ptr<SymbolicNode>> m_ops;
                m_ops.push_back(coeff);
                if (auto m = std::dynamic_pointer_cast<MultiplyNode>(term)) {
                     m_ops.insert(m_ops.end(), m->operands.begin(), m->operands.end());
                } else {
                     m_ops.push_back(term);
                }
                final_ops.push_back(std::make_shared<MultiplyNode>(m_ops));
            }
        }

        if (final_ops.empty()) {
            result = std::make_shared<NumberNode>(BigInt(0));
        } else if (final_ops.size() == 1) {
            result = final_ops[0];
        } else {
            // Sort final_ops: Higher degree first, constants LAST
            std::sort(final_ops.begin(), final_ops.end(), [](const std::shared_ptr<SymbolicNode>& l, const std::shared_ptr<SymbolicNode>& r) {
                int d1 = get_node_degree_helper(l);
                int d2 = get_node_degree_helper(r);
                if (d1 != d2) return d1 > d2;
                
                bool isNum1 = std::dynamic_pointer_cast<NumberNode>(l) != nullptr;
                bool isNum2 = std::dynamic_pointer_cast<NumberNode>(r) != nullptr;
                if (isNum1 != isNum2) return isNum2; // Numbers last
                
                return l->compare(*r) < 0;
            });
            result = std::make_shared<AddNode>(final_ops);
        }
    }

    void visit(MultiplyNode& node) override {
        
        std::vector<std::shared_ptr<SymbolicNode>> sc;
        
        for (const auto& op : node.operands) {
            op->accept(*this);
            auto res = result;
            if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(res)) {
                
                sc.insert(sc.end(), mul->operands.begin(), mul->operands.end());
            } else {
                sc.push_back(res);
            }
        }

        
        bool has_add = false;
        
        
        
        
        
        
        
        
        
        
        for(const auto& op : sc) {
            if(std::dynamic_pointer_cast<AddNode>(op)) {
                has_add = true;
                break;
            }
        }

        if (has_add) {
            if (sc.empty()) {
                result = std::make_shared<NumberNode>(BigInt(1));
                return;
            }
            std::shared_ptr<SymbolicNode> current = sc[0];
            for(size_t i=1; i<sc.size(); ++i) {
                
                
                current = expand_product(current, sc[i]);
            }
            
            result = current;
            return;
        }

        
        bool has_matrix = false;
        for(const auto& op : sc) {
            if (std::dynamic_pointer_cast<MatrixNode>(op)) { has_matrix = true; break; }
            if (auto p = std::dynamic_pointer_cast<PowerNode>(op)) {
                if (std::dynamic_pointer_cast<MatrixNode>(p->base)) { has_matrix = true; break; }
            }
        }

        if (has_matrix) {
            std::vector<std::shared_ptr<SymbolicNode>> new_ops;
            std::shared_ptr<NumberNode> scalar_part = std::make_shared<NumberNode>(BigInt(1));
            
            for(const auto& op : sc) {
                if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
                    scalar_part = multiply_numbers(scalar_part, num);
                } else {
                    new_ops.push_back(op);
                }
            }
            
            
            
            
            std::vector<std::shared_ptr<SymbolicNode>> fused_ops;
            if (!new_ops.empty()) fused_ops.push_back(new_ops[0]);
            
            for(size_t i=1; i<new_ops.size(); ++i) {
                auto left = fused_ops.back();
                auto right = new_ops[i];
                
                auto m_left = std::dynamic_pointer_cast<MatrixNode>(left);
                auto m_right = std::dynamic_pointer_cast<MatrixNode>(right);
                
                if (m_left && m_right) {
                    
                    if (m_left->cols == m_right->rows) {
                         
                         if (std::holds_alternative<MatrixNode::DenseStorage>(m_left->storage) && 
                             std::holds_alternative<MatrixNode::DenseStorage>(m_right->storage)) {
                             
                             const auto& d_l = std::get<MatrixNode::DenseStorage>(m_left->storage);
                             const auto& d_r = std::get<MatrixNode::DenseStorage>(m_right->storage);
                             
                             size_t R = m_left->rows;
                             size_t C = m_right->cols;
                             size_t K = m_left->cols;

                             
                             MatrixNode::DenseStorage res_data;
                             res_data.reserve(R*C);
                             
                             for(size_t r=0; r<R; ++r) {
                                 for(size_t c=0; c<C; ++c) {
                                     
                                     std::vector<std::shared_ptr<SymbolicNode>> sum_ops;
                                     for(size_t k=0; k<K; ++k) {
                                         std::vector<std::shared_ptr<SymbolicNode>> prod_ops = {
                                             d_l[r*K + k], d_r[k*C + c]
                                         };
                                         sum_ops.push_back(std::make_shared<MultiplyNode>(prod_ops));
                                     }
                                     
                                     NormalizationVisitor elem_vis;
                                     auto elem_node = std::make_shared<AddNode>(sum_ops);
                                     elem_node->accept(elem_vis);
                                     
                                     auto res_val = elem_vis.get_result();
                                     if (auto num = std::dynamic_pointer_cast<NumberNode>(res_val)) {
                                         if (std::holds_alternative<double>(num->value)) {
                                             double v = std::get<double>(num->value);
                                             if (std::abs(v) < 1e-10) { 
                                                 res_val = std::make_shared<NumberNode>(BigInt(0));
                                             } else if (std::abs(v - 1.0) < 1e-10) {
                                                 res_val = std::make_shared<NumberNode>(BigInt(1));
                                             }
                                         }
                                     }
                                     res_data.push_back(res_val);
                                 }
                             }
                             
                             fused_ops.pop_back();
                             fused_ops.push_back(std::make_shared<MatrixNode>(R, C, res_data));
                             continue;
                         }
                    }
                }
                fused_ops.push_back(right);
            }
            
            if (!scalar_part->is_one()) {
                fused_ops.insert(fused_ops.begin(), scalar_part);
            }
            
            if (fused_ops.empty()) result = std::make_shared<NumberNode>(BigInt(1));
            else if (fused_ops.size() == 1) result = fused_ops[0];
            else result = std::make_shared<MultiplyNode>(fused_ops);
            
            return;
        }

        std::shared_ptr<NumberNode> const_acc = std::make_shared<NumberNode>(BigInt(1));
        std::map<std::shared_ptr<SymbolicNode>, std::shared_ptr<NumberNode>, NodeCompare> bases;

        for (const auto& op : sc) {
             if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
                 if (num->is_zero()) {
                     result = std::make_shared<NumberNode>(BigInt(0));
                     return;
                 }
                 const_acc = multiply_numbers(const_acc, num);
             } else {
                 std::shared_ptr<SymbolicNode> base = op;
                 std::shared_ptr<NumberNode> exp = std::make_shared<NumberNode>(BigInt(1));
                 bool is_number_power = false;

                 if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                     base = pow->base;
                     if (auto e_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                         exp = e_num;
                         
                         if (auto b_num = std::dynamic_pointer_cast<NumberNode>(base)) {
                             long long exp_val = 0;
                             bool exp_ok = false;
                             bool exp_is_half = false;
                             
                             if (std::holds_alternative<BigInt>(e_num->value)) {
                                 exp_val = (long long)std::get<BigInt>(e_num->value).to_double();
                                 exp_ok = true;
                             } else if (std::holds_alternative<double>(e_num->value)) {
                                 double d = std::get<double>(e_num->value);
                                 if (d == std::floor(d)) {
                                     exp_val = (long long)d;
                                     exp_ok = true;
                                 } else if (std::abs(d - 0.5) < 1e-9) {
                                     exp_is_half = true;
                                 }
                             } else if (std::holds_alternative<Rational>(e_num->value)) {
                                 Rational r = std::get<Rational>(e_num->value);
                                 if (r.get_denominator() == BigInt(1)) {
                                     exp_val = (long long)r.get_numerator().to_double();
                                     exp_ok = true;
                                 } else if (r.get_numerator() == BigInt(1) && r.get_denominator() == BigInt(2)) {
                                     exp_is_half = true;
                                 }
                             }

                             if (exp_ok) {
                                 std::shared_ptr<NumberNode> pow_val = nullptr;
                                 
                                 if (exp_val == -1) {
                                      if (std::holds_alternative<BigInt>(b_num->value)) {
                                          pow_val = std::make_shared<NumberNode>(Rational(BigInt(1), std::get<BigInt>(b_num->value)));
                                      } else if (std::holds_alternative<Rational>(b_num->value)) {
                                          Rational r = std::get<Rational>(b_num->value);
                                          if (!r.get_numerator().is_zero()) pow_val = std::make_shared<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                      } else if (std::holds_alternative<double>(b_num->value)) {
                                          pow_val = std::make_shared<NumberNode>(1.0 / std::get<double>(b_num->value));
                                      }
                                 } else if (exp_val == 0) {
                                      pow_val = std::make_shared<NumberNode>(BigInt(1));
                                 } else if (std::abs(exp_val) > 0 && std::abs(exp_val) < 64) { 
                                      
                                      long long abs_exp = std::abs(exp_val);
                                      std::shared_ptr<NumberNode> base_pow_val = nullptr;

                                      if (std::holds_alternative<BigInt>(b_num->value)) {
                                          BigInt b = std::get<BigInt>(b_num->value);
                                          BigInt res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      } else if (std::holds_alternative<Rational>(b_num->value)) {
                                          Rational b = std::get<Rational>(b_num->value);
                                          Rational res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      } else if (std::holds_alternative<double>(b_num->value)) {
                                          double b = std::get<double>(b_num->value);
                                          double res = 1.0;
                                          for(int k=0;k<abs_exp;++k) res *= b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      }

                                      if (base_pow_val) {
                                          if (exp_val > 0) {
                                              pow_val = base_pow_val;
                                          } else {
                                              
                                              if (std::holds_alternative<BigInt>(base_pow_val->value)) {
                                                  pow_val = std::make_shared<NumberNode>(Rational(BigInt(1), std::get<BigInt>(base_pow_val->value)));
                                              } else if (std::holds_alternative<Rational>(base_pow_val->value)) {
                                                  Rational r = std::get<Rational>(base_pow_val->value);
                                                  pow_val = std::make_shared<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                              } else if (std::holds_alternative<double>(base_pow_val->value)) {
                                                  pow_val = std::make_shared<NumberNode>(1.0 / std::get<double>(base_pow_val->value));
                                              }
                                          }
                                      }
                                 }
                                 
                                 if (pow_val) {
                                     const_acc = multiply_numbers(const_acc, pow_val);
                                     is_number_power = true;
                                 }
                             } else if (exp_is_half) {
                                 
                                 std::shared_ptr<NumberNode> root_val = nullptr;
                                 if (std::holds_alternative<BigInt>(b_num->value)) {
                                     
                                     double d = std::get<BigInt>(b_num->value).to_double();
                                     if (d >= 0) {
                                         double r = std::sqrt(d);
                                         if (std::abs(r - std::round(r)) < 1e-9) {
                                             BigInt bi((long long)std::round(r));
                                             if (bi * bi == std::get<BigInt>(b_num->value)) {
                                                 root_val = std::make_shared<NumberNode>(bi);
                                             }
                                         }
                                     }
                                 } else if (std::holds_alternative<double>(b_num->value)) {
                                     double d = std::get<double>(b_num->value);
                                     if (d >= 0) root_val = std::make_shared<NumberNode>(std::sqrt(d));
                                 }
                                 
                                 if (root_val) {
                                     const_acc = multiply_numbers(const_acc, root_val);
                                     is_number_power = true;
                                 }
                             }
                         }
                     } 
                 }

                 if (!is_number_power) {
                     auto it = bases.find(base);
                     if (it == bases.end()) {
                         bases[base] = exp;
                     } else {
                         bases[base] = add_numbers(it->second, exp);
                     }
                 }
             }
        }

        std::vector<std::shared_ptr<SymbolicNode>> final_ops;
        if (!const_acc->is_one()) {
             final_ops.push_back(const_acc);
        }
        
        std::vector<std::shared_ptr<SymbolicNode>> var_ops;
        for (auto const& [base, exp] : bases) {
            if (exp->is_zero()) {
            } else if (exp->is_one()) {
                var_ops.push_back(base);
            } else {
                var_ops.push_back(std::make_shared<PowerNode>(base, exp));
            }
        }
        
        // Sort variable ops to be consistent
        std::sort(var_ops.begin(), var_ops.end(), [](const std::shared_ptr<SymbolicNode>& l, const std::shared_ptr<SymbolicNode>& r) {
            int d1 = get_node_degree_helper(l);
            int d2 = get_node_degree_helper(r);
            if (d1 != d2) return d1 > d2;
            return l->compare(*r) < 0;
        });

        final_ops.insert(final_ops.end(), var_ops.begin(), var_ops.end());
        
        if (final_ops.empty()) result = std::make_shared<NumberNode>(BigInt(1));
        else if (final_ops.size() == 1) result = final_ops[0];
        else result = std::make_shared<MultiplyNode>(final_ops);
    }

    void visit(PowerNode& node) override {
        node.base->accept(*this);
        auto s_base = result;
        node.exponent->accept(*this);
        auto s_exp = result;
         
        if (s_exp->is_zero()) {
            result = std::make_shared<NumberNode>(BigInt(1));
            return;
        }
        if (s_exp->is_one()) {
            result = s_base;
            return;
        }

        if (s_base->is_zero()) {
             result = std::make_shared<NumberNode>(BigInt(0));
             return;
        }
        if (s_base->is_one()) {
            result = std::make_shared<NumberNode>(BigInt(1));
            return;
        }

        
        if (auto b_num = std::dynamic_pointer_cast<NumberNode>(s_base)) {
            if (auto e_num = std::dynamic_pointer_cast<NumberNode>(s_exp)) {
                 long long exp_val = 0;
                 bool exp_ok = false;
                 bool exp_is_half = false;
                 
                 if (std::holds_alternative<BigInt>(e_num->value)) {
                     exp_val = (long long)std::get<BigInt>(e_num->value).to_double();
                     exp_ok = true;
                 } else if (std::holds_alternative<double>(e_num->value)) {
                     double d = std::get<double>(e_num->value);
                     if (d == std::floor(d)) {
                         exp_val = (long long)d;
                         exp_ok = true;
                     } else if (std::abs(d - 0.5) < 1e-9) {
                         exp_is_half = true;
                     }
                 } else if (std::holds_alternative<Rational>(e_num->value)) {
                     Rational r = std::get<Rational>(e_num->value);
                     if (r.get_denominator() == BigInt(1)) {
                         exp_val = (long long)r.get_numerator().to_double();
                         exp_ok = true;
                     } else if (r.get_numerator() == BigInt(1) && r.get_denominator() == BigInt(2)) {
                         exp_is_half = true;
                     }
                 }

                 if (exp_ok) {
                     std::shared_ptr<NumberNode> pow_val = nullptr;
                     if (exp_val == -1) {
                          if (std::holds_alternative<BigInt>(b_num->value)) {
                              pow_val = std::make_shared<NumberNode>(Rational(BigInt(1), std::get<BigInt>(b_num->value)));
                          } else if (std::holds_alternative<Rational>(b_num->value)) {
                              Rational r = std::get<Rational>(b_num->value);
                              if (!r.get_numerator().is_zero()) pow_val = std::make_shared<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                          } else if (std::holds_alternative<double>(b_num->value)) {
                              pow_val = std::make_shared<NumberNode>(1.0 / std::get<double>(b_num->value));
                          }
                     } else if (exp_val == 0) {
                          pow_val = std::make_shared<NumberNode>(BigInt(1));
                     } else if (std::abs(exp_val) > 0 && std::abs(exp_val) < 64) {
                          
                                      long long abs_exp = std::abs(exp_val);
                                      std::shared_ptr<NumberNode> base_pow_val = nullptr;

                                      if (std::holds_alternative<BigInt>(b_num->value)) {
                                          BigInt b = std::get<BigInt>(b_num->value);
                                          BigInt res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      } else if (std::holds_alternative<Rational>(b_num->value)) {
                                          Rational b = std::get<Rational>(b_num->value);
                                          Rational res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      } else if (std::holds_alternative<double>(b_num->value)) {
                                          double b = std::get<double>(b_num->value);
                                          double res = 1.0;
                                          for(int k=0;k<abs_exp;++k) res *= b;
                                          base_pow_val = std::make_shared<NumberNode>(res);
                                      }

                                      if (base_pow_val) {
                                          if (exp_val > 0) {
                                              pow_val = base_pow_val;
                                          } else {
                                              
                                              if (std::holds_alternative<BigInt>(base_pow_val->value)) {
                                                  pow_val = std::make_shared<NumberNode>(Rational(BigInt(1), std::get<BigInt>(base_pow_val->value)));
                                              } else if (std::holds_alternative<Rational>(base_pow_val->value)) {
                                                  Rational r = std::get<Rational>(base_pow_val->value);
                                                  pow_val = std::make_shared<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                              } else if (std::holds_alternative<double>(base_pow_val->value)) {
                                                  pow_val = std::make_shared<NumberNode>(1.0 / std::get<double>(base_pow_val->value));
                                              }
                                          }
                                      }
                     }
                     
                     if (pow_val) {
                         result = pow_val;
                         return;
                     }
                 } else if (exp_is_half) {
                     if (std::holds_alternative<BigInt>(b_num->value)) {
                         double d = std::get<BigInt>(b_num->value).to_double();
                         if (d >= 0) {
                             double r = std::sqrt(d);
                             if (std::abs(r - std::round(r)) < 1e-9) {
                                 BigInt bi((long long)std::round(r));
                                 if (bi * bi == std::get<BigInt>(b_num->value)) {
                                     result = std::make_shared<NumberNode>(bi);
                                     return;
                                 }
                             }
                         }
                     } else if (std::holds_alternative<double>(b_num->value)) {
                         double d = std::get<double>(b_num->value);
                         if (d >= 0) {
                             result = std::make_shared<NumberNode>(std::sqrt(d));
                             return;
                         }
                     }
                 }
            }
        }
        
        
        if (auto m_base = std::dynamic_pointer_cast<MultiplyNode>(s_base)) {
            std::vector<std::shared_ptr<SymbolicNode>> new_ops;
            for(auto& op : m_base->operands) {
                
                auto term_pow = std::make_shared<PowerNode>(op, s_exp);
                term_pow->accept(*this); 
                
                
                if (auto mul_res = std::dynamic_pointer_cast<MultiplyNode>(result)) {
                     new_ops.insert(new_ops.end(), mul_res->operands.begin(), mul_res->operands.end());
                } else {
                     new_ops.push_back(result);
                }
            }
            
            auto final_mul = std::make_shared<MultiplyNode>(new_ops);
            final_mul->accept(*this);
            return;
        }

        
        if (auto p_base = std::dynamic_pointer_cast<PowerNode>(s_base)) {
             
             std::vector<std::shared_ptr<SymbolicNode>> exp_ops;
             exp_ops.push_back(p_base->exponent);
             exp_ops.push_back(s_exp);
             
             auto mul_exp = std::make_shared<MultiplyNode>(exp_ops);
             mul_exp->accept(*this); 
             
             auto new_pow = std::make_shared<PowerNode>(p_base->base, result);
             new_pow->accept(*this); 
             return;
        }

        result = std::make_shared<PowerNode>(s_base, s_exp);
    }

    void visit(FunctionNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> s_args;
        for(auto& a : node.arguments) {
            a->accept(*this);
            s_args.push_back(result);
        }
        
        // Basic constant evaluation
        if (s_args.size() == 1) {
            if (auto num = std::dynamic_pointer_cast<NumberNode>(s_args[0])) {
                double val = 0.0;
                if (std::holds_alternative<double>(num->value)) val = std::get<double>(num->value);
                else if (std::holds_alternative<BigInt>(num->value)) val = std::get<BigInt>(num->value).to_double();
                else if (std::holds_alternative<Rational>(num->value)) val = std::get<Rational>(num->value).to_double();

                switch (node.type) {
                    case FunctionNode::FuncType::Sin:
                        if (std::abs(val) < 1e-12) { result = std::make_shared<NumberNode>(BigInt(0)); return; }
                        break;
                    case FunctionNode::FuncType::Cos:
                        if (std::abs(val) < 1e-12) { result = std::make_shared<NumberNode>(BigInt(1)); return; }
                        break;
                    case FunctionNode::FuncType::Tan:
                        if (std::abs(val) < 1e-12) { result = std::make_shared<NumberNode>(BigInt(0)); return; }
                        break;
                    case FunctionNode::FuncType::Exp:
                        if (std::abs(val) < 1e-12) { result = std::make_shared<NumberNode>(BigInt(1)); return; }
                        break;
                    case FunctionNode::FuncType::Ln:
                        if (std::abs(val - 1.0) < 1e-12) { result = std::make_shared<NumberNode>(BigInt(0)); return; }
                        break;
                    case FunctionNode::FuncType::Abs:
                        result = std::make_shared<NumberNode>(std::abs(val));
                        return;
                    case FunctionNode::FuncType::Sqrt:
                        if (val >= 0) {
                             double sq = std::sqrt(val);
                             if (std::abs(sq - std::round(sq)) < 1e-12) {
                                 result = std::make_shared<NumberNode>(BigInt((long long)std::round(sq)));
                             } else {
                                 result = std::make_shared<NumberNode>(sq);
                             }
                             return;
                        }
                        break;
                    default: break;
                }
            }
        }
        
        result = std::make_shared<FunctionNode>(node.type, s_args);
    }

    void visit(MatrixNode& node) override {
        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage)) {
             auto& dense = std::get<MatrixNode::DenseStorage>(node.storage);
             MatrixNode::DenseStorage new_dense;
             for(auto& item : dense) {
                 if(item) {
                    item->accept(*this);
                    new_dense.push_back(result);
                 } else {
                    new_dense.push_back(nullptr);
                 }
             }
             result = std::make_shared<MatrixNode>(node.rows, node.cols, new_dense);
        } else {
             auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage);
             MatrixNode::SparseStorage new_sparse;
             for(auto& [idx, item] : sparse) {
                 item->accept(*this);
                 if (!result->is_zero()) {
                     new_sparse[idx] = result;
                 }
             }
             result = std::make_shared<MatrixNode>(node.rows, node.cols, new_sparse);
        }
    }

    void visit(RelationalNode& node) override {
        
        std::shared_ptr<SymbolicNode> new_left = nullptr;
        std::shared_ptr<SymbolicNode> new_right = nullptr;

        if (node.left) {
            node.left->accept(*this);
            new_left = result;
        }
        if (node.right) {
            
            
            
            
            
            node.right->accept(*this);
            new_right = result;
        }
        
        if (!new_left) new_left = node.left;
        if (!new_right) new_right = node.right;

        result = std::make_shared<RelationalNode>(new_left, new_right, node.op);
    }
};
