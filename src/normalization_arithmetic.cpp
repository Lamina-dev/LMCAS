#include "visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"

std::shared_ptr<const SymbolicNode> NormalizationVisitor::get_result() const {
        return result;
    }

std::shared_ptr<const SymbolicNode> NormalizationVisitor::expand_product(const std::shared_ptr<const SymbolicNode>& lhs, const std::shared_ptr<const SymbolicNode>& rhs) {

        auto add_lhs = std::dynamic_pointer_cast<const AddNode>(lhs);
        auto add_rhs = std::dynamic_pointer_cast<const AddNode>(rhs);

        auto c_lhs = std::dynamic_pointer_cast<const ComplexNode>(lhs);
        auto c_rhs = std::dynamic_pointer_cast<const ComplexNode>(rhs);

        if (c_lhs && c_rhs) {
            auto ac = expand_product(c_lhs->real(), c_rhs->real());
            auto bd = expand_product(c_lhs->imag(), c_rhs->imag());
            auto ad = expand_product(c_lhs->real(), c_rhs->imag());
            auto bc = expand_product(c_lhs->imag(), c_rhs->real());
            auto neg_one = lamina::detail::make_node<NumberNode>(BigInt(-1));
            auto neg_bd = expand_product(neg_one, bd);
            auto real_part = SymbolicFactory::create_add({ac, neg_bd});
            auto imag_part = SymbolicFactory::create_add({ad, bc});
            NormalizationVisitor norm(assumptions_);
            real_part->accept(norm); auto norm_r = norm.get_result();
            imag_part->accept(norm); auto norm_i = norm.get_result();
            return SymbolicFactory::create_complex(norm_r, norm_i);
        } else if (c_lhs) {
            auto nr = expand_product(c_lhs->real(), rhs);
            auto ni = expand_product(c_lhs->imag(), rhs);
            NormalizationVisitor norm(assumptions_);
            nr->accept(norm); auto norm_r = norm.get_result();
            ni->accept(norm); auto norm_i = norm.get_result();
            return SymbolicFactory::create_complex(norm_r, norm_i);
        } else if (c_rhs) {
            auto nr = expand_product(lhs, c_rhs->real());
            auto ni = expand_product(lhs, c_rhs->imag());
            NormalizationVisitor norm(assumptions_);
            nr->accept(norm); auto norm_r = norm.get_result();
            ni->accept(norm); auto norm_i = norm.get_result();
            return SymbolicFactory::create_complex(norm_r, norm_i);
        }

        if (add_lhs && add_rhs) {

            std::vector<std::shared_ptr<const SymbolicNode>> new_terms;
            for (const auto& op1 : add_lhs->operands()) {
                for (const auto& op2 : add_rhs->operands()) {
                    auto prod = expand_product(op1, op2);
                    new_terms.push_back(prod);
                }
            }

            return lamina::detail::make_node<AddNode>(new_terms);
        } else if (add_lhs) {

            std::vector<std::shared_ptr<const SymbolicNode>> new_terms;
            for (const auto& op : add_lhs->operands()) {
                 auto prod = expand_product(op, rhs);
                 new_terms.push_back(prod);
            }
            return lamina::detail::make_node<AddNode>(new_terms);
        } else if (add_rhs) {

            std::vector<std::shared_ptr<const SymbolicNode>> new_terms;
            for (const auto& op : add_rhs->operands()) {
                 auto prod = expand_product(lhs, op);
                 new_terms.push_back(prod);
            }
            return lamina::detail::make_node<AddNode>(new_terms);
        }

        std::shared_ptr<const NumberNode> const_acc = lamina::detail::make_node<NumberNode>(BigInt(1));
        struct FactorAccum {
            std::shared_ptr<const NumberNode> exponent;
            std::vector<std::shared_ptr<const SymbolicNode>> originals;
        };
        std::map<std::shared_ptr<const SymbolicNode>, FactorAccum, NodeCompare> bases;

        auto process_factor = [&](const std::shared_ptr<const SymbolicNode>& factor) {
             if (auto num = std::dynamic_pointer_cast<const NumberNode>(factor)) {
                 if (num->is_zero()) return false;
                 const_acc = multiply_numbers(const_acc, num);
             } else {
                 std::shared_ptr<const SymbolicNode> base = factor;
                 std::shared_ptr<const NumberNode> exp = lamina::detail::make_node<NumberNode>(BigInt(1));

                 if (auto pow = std::dynamic_pointer_cast<const PowerNode>(factor)) {
                     if (auto e_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                         if (is_positive_integer_number(e_num)) {
                             base = pow->base();
                             exp = e_num;
                         }
                     }
                 }

                 auto it = bases.find(base);
                 if (it == bases.end()) {
                     bases.emplace(base, FactorAccum{exp, {factor}});
                 } else {
                     it->second.exponent = add_numbers(it->second.exponent, exp);
                     it->second.originals.push_back(factor);
                 }
             }
             return true;
        };

        auto flatten_and_process = [&](const std::shared_ptr<const SymbolicNode>& node) {
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
                for(const auto& op : mul->operands()) {
                    if (!process_factor(op)) return false;
                }
            } else {
                if (!process_factor(node)) return false;
            }
            return true;
        };

        if (!flatten_and_process(lhs)) return lamina::detail::make_node<NumberNode>(BigInt(0));
        if (!flatten_and_process(rhs)) return lamina::detail::make_node<NumberNode>(BigInt(0));

        std::vector<std::shared_ptr<const SymbolicNode>> final_ops;
        if (!const_acc->is_one()) {
             final_ops.push_back(const_acc);
        }

        std::vector<std::shared_ptr<const SymbolicNode>> var_ops;
        for (auto const& [base, acc] : bases) {
            const auto& exp = acc.exponent;
            if (exp->is_zero()) {
                var_ops.insert(var_ops.end(), acc.originals.begin(), acc.originals.end());
            } else if (exp->is_one()) {
                var_ops.push_back(base);
            } else {
                var_ops.push_back(lamina::detail::make_node<PowerNode>(base, exp));
            }
        }

        std::sort(var_ops.begin(), var_ops.end(), [](const std::shared_ptr<const SymbolicNode>& l, const std::shared_ptr<const SymbolicNode>& r) {
            int d1 = get_node_degree_helper(l);
            int d2 = get_node_degree_helper(r);
            if (d1 != d2) return d1 > d2;
            return l->compare(*r) < 0;
        });

        final_ops.insert(final_ops.end(), var_ops.begin(), var_ops.end());

        if (final_ops.empty()) return lamina::detail::make_node<NumberNode>(BigInt(1));
        if (final_ops.size() == 1) return final_ops[0];

        if (final_ops.size() > 1 && std::dynamic_pointer_cast<const NumberNode>(final_ops.back())) {
             std::rotate(final_ops.begin(), final_ops.end() - 1, final_ops.end());
        }
        return make_normalized_multiply_node(final_ops);
    }
void NormalizationVisitor::visit(const NumberNode& node) {
        result = node.clone();
    }
void NormalizationVisitor::visit(const VariableNode& node) {
        result = node.clone();
    }
void NormalizationVisitor::visit(const AddNode& node) {
        std::vector<std::shared_ptr<const SymbolicNode>> simplified_ops;

        for (const auto& op : node.operands()) {
            op->accept(*this);
            if (auto add = std::dynamic_pointer_cast<const AddNode>(result)) {
                simplified_ops.insert(simplified_ops.end(), add->operands().begin(), add->operands().end());
            } else {
                simplified_ops.push_back(result);
            }
        }

        /// Merge ComplexNodes and NumberNodes
        std::vector<std::shared_ptr<const SymbolicNode>> real_parts;
        std::vector<std::shared_ptr<const SymbolicNode>> imag_parts;
        std::vector<std::shared_ptr<const SymbolicNode>> non_complex_ops;
        bool has_complex = false;

        for (const auto& op : simplified_ops) {
            if (auto c = std::dynamic_pointer_cast<const ComplexNode>(op)) {
                has_complex = true;
                real_parts.push_back(c->real());
                imag_parts.push_back(c->imag());
            } else if (std::dynamic_pointer_cast<const NumberNode>(op)) {
                real_parts.push_back(op);
                imag_parts.push_back(SymbolicFactory::create_number(BigInt(0)));
            } else {
                non_complex_ops.push_back(op);
            }
        }

        if (has_complex) {
            auto real_sum = SymbolicFactory::create_add(real_parts);
            auto imag_sum = SymbolicFactory::create_add(imag_parts);
            NormalizationVisitor sub_norm(assumptions_);
            real_sum->accept(sub_norm); auto norm_real = sub_norm.get_result();
            imag_sum->accept(sub_norm); auto norm_imag = sub_norm.get_result();
            auto merged_complex = SymbolicFactory::create_complex(norm_real, norm_imag);
            if (!merged_complex->is_zero()) {
                non_complex_ops.push_back(merged_complex);
            }
            simplified_ops = non_complex_ops;
        }

        if (!simplified_ops.empty() && std::dynamic_pointer_cast<const MatrixNode>(simplified_ops[0])) {
            auto first_mat = std::dynamic_pointer_cast<const MatrixNode>(simplified_ops[0]);
            size_t rows = first_mat->rows();
            size_t cols = first_mat->cols();
            bool all_matrices = true;
            for (const auto& op : simplified_ops) {
                 auto m = std::dynamic_pointer_cast<const MatrixNode>(op);
                 if (!m || m->rows() != rows || m->cols() != cols) {
                     all_matrices = false; break;
                 }
            }

            if (all_matrices) {
                 std::vector<std::shared_ptr<const SymbolicNode>> new_elements;
                 new_elements.reserve(rows * cols);

                 for (size_t i = 0; i < rows * cols; ++i) {
                     std::vector<std::shared_ptr<const SymbolicNode>> elem_ops;
                     for (const auto& op : simplified_ops) {
                         auto m = std::dynamic_pointer_cast<const MatrixNode>(op);
                         std::shared_ptr<const SymbolicNode> val;
                         if (std::holds_alternative<MatrixNode::DenseStorage>(m->storage())) {
                             const auto& dense = std::get<MatrixNode::DenseStorage>(m->storage());
                             if (i < dense.size()) val = dense[i];
                             else val = lamina::detail::make_node<NumberNode>(BigInt(0));
                         } else {
                             const auto& sparse = std::get<MatrixNode::SparseStorage>(m->storage());
                             auto it = sparse.find(i);
                             if (it != sparse.end()) val = it->second;
                             else val = lamina::detail::make_node<NumberNode>(BigInt(0));
                         }
                         if (!val) val = lamina::detail::make_node<NumberNode>(BigInt(0));
                         elem_ops.push_back(val);
                     }

                     auto elem_add = lamina::detail::make_node<AddNode>(elem_ops);
                     elem_add->accept(*this);
                     new_elements.push_back(result);
                 }

                 result = lamina::detail::make_node<MatrixNode>(rows, cols, new_elements);
                 return;
            }
        }

        std::shared_ptr<const NumberNode> constant_acc = lamina::detail::make_node<NumberNode>(BigInt(0));
        std::map<std::shared_ptr<const SymbolicNode>, std::shared_ptr<const NumberNode>, NodeCompare> terms;

        for (const auto& op : simplified_ops) {
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(op)) {
                constant_acc = add_numbers(constant_acc, num);
            } else {
                std::shared_ptr<const SymbolicNode> term_part = op;
                std::shared_ptr<const NumberNode> coeff_part = lamina::detail::make_node<NumberNode>(BigInt(1));

                if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(op)) {

                    std::shared_ptr<const NumberNode> coeff = nullptr;
                    int coeff_idx = -1;

                    if (!mul->operands().empty()) {
                        if (auto n_first = std::dynamic_pointer_cast<const NumberNode>(mul->operands().front())) {
                            coeff = n_first;
                            coeff_idx = 0;
                        } else if (auto n_last = std::dynamic_pointer_cast<const NumberNode>(mul->operands().back())) {
                            coeff = n_last;
                            coeff_idx = (int)mul->operands().size() - 1;
                        }
                    }

                    if (coeff) {
                         coeff_part = coeff;
                         if (mul->operands().size() == 2) {
                             term_part = mul->operands()[coeff_idx == 0 ? 1 : 0];
                         } else {
                             std::vector<std::shared_ptr<const SymbolicNode>> rest;
                             rest.reserve(mul->operands().size() - 1);
                             for(int k=0; k<(int)mul->operands().size(); ++k) {
                                 if (k != coeff_idx) rest.push_back(mul->operands()[k]);
                             }

                             auto rest_mul = make_normalized_multiply_node(rest);
                             NormalizationVisitor rest_norm(assumptions_);
                             rest_mul->accept(rest_norm);
                             term_part = rest_norm.get_result();
                         }
                    }
                }

                auto it = terms.end();
                for (auto scan = terms.begin(); scan != terms.end(); ++scan) {
                    if (term_part->compare(*scan->first) == 0) {
                        it = scan;
                        break;
                    }
                }
                if (it == terms.end()) {
                    terms[term_part] = coeff_part;
                } else {
                    it->second = add_numbers(it->second, coeff_part);
                }
            }
        }

        std::vector<std::shared_ptr<const SymbolicNode>> final_ops;

        if (!constant_acc->is_zero()) {
            final_ops.push_back(constant_acc);
        }

        for (auto const& [term, coeff] : terms) {
            if (coeff->is_zero()) continue;

            if (coeff->is_one()) {
                final_ops.push_back(term);
            } else {
                std::vector<std::shared_ptr<const SymbolicNode>> m_ops;
                m_ops.push_back(coeff);
                if (auto m = std::dynamic_pointer_cast<const MultiplyNode>(term)) {
                     m_ops.insert(m_ops.end(), m->operands().begin(), m->operands().end());
                } else {
                     m_ops.push_back(term);
                }
                final_ops.push_back(make_normalized_multiply_node(m_ops));
            }
        }

        if (final_ops.empty()) {
            result = lamina::detail::make_node<NumberNode>(BigInt(0));
        } else if (final_ops.size() == 1) {
            result = final_ops[0];
        } else {

            std::sort(final_ops.begin(), final_ops.end(), [](const std::shared_ptr<const SymbolicNode>& l, const std::shared_ptr<const SymbolicNode>& r) {
                int d1 = get_node_degree_helper(l);
                int d2 = get_node_degree_helper(r);
                if (d1 != d2) return d1 > d2;

                bool isNum1 = std::dynamic_pointer_cast<const NumberNode>(l) != nullptr;
                bool isNum2 = std::dynamic_pointer_cast<const NumberNode>(r) != nullptr;
                if (isNum1 != isNum2) return isNum2;

                return l->compare(*r) < 0;
            });
            result = lamina::detail::make_node<AddNode>(final_ops);
        }
    }
void NormalizationVisitor::visit(const ComplexNode& node) {
        node.real()->accept(*this);
        auto norm_r = result;
        node.imag()->accept(*this);
        auto norm_i = result;
        result = SymbolicFactory::create_complex(norm_r, norm_i);
    }
void NormalizationVisitor::visit(const MultiplyNode& node) {

        std::vector<std::shared_ptr<const SymbolicNode>> sc;

        for (const auto& op : node.operands()) {
            op->accept(*this);
            auto res = result;
            if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(res)) {

                sc.insert(sc.end(), mul->operands().begin(), mul->operands().end());
            } else {
                sc.push_back(res);
            }
        }

        // `simplify()` must not apply distributivity. Expansion changes the
        // expression shape, can trigger term explosion, and may hide domain
        // conditions introduced by later cancellation. Explicit `expand()` is
        // the only path that may call expand_product().

        bool has_matrix = false;
        for(const auto& op : sc) {
            if (std::dynamic_pointer_cast<const MatrixNode>(op)) { has_matrix = true; break; }
            if (auto p = std::dynamic_pointer_cast<const PowerNode>(op)) {
                if (std::dynamic_pointer_cast<const MatrixNode>(p->base())) { has_matrix = true; break; }
            }
        }

        if (has_matrix) {
            std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
            std::shared_ptr<const NumberNode> scalar_part = lamina::detail::make_node<NumberNode>(BigInt(1));

            for(const auto& op : sc) {
                if (auto num = std::dynamic_pointer_cast<const NumberNode>(op)) {
                    scalar_part = multiply_numbers(scalar_part, num);
                } else {
                    new_ops.push_back(op);
                }
            }

            std::vector<std::shared_ptr<const SymbolicNode>> fused_ops;
            if (!new_ops.empty()) fused_ops.push_back(new_ops[0]);

            for(size_t i=1; i<new_ops.size(); ++i) {
                auto left = fused_ops.back();
                auto right = new_ops[i];

                auto m_left = std::dynamic_pointer_cast<const MatrixNode>(left);
                auto m_right = std::dynamic_pointer_cast<const MatrixNode>(right);

                if (m_left && m_right) {

                    if (m_left->cols() == m_right->rows()) {

                         if (std::holds_alternative<MatrixNode::DenseStorage>(m_left->storage()) &&
                             std::holds_alternative<MatrixNode::DenseStorage>(m_right->storage())) {

                             const auto& d_l = std::get<MatrixNode::DenseStorage>(m_left->storage());
                             const auto& d_r = std::get<MatrixNode::DenseStorage>(m_right->storage());

                             size_t R = m_left->rows();
                             size_t C = m_right->cols();
                             size_t K = m_left->cols();

                             MatrixNode::DenseStorage res_data;
                             res_data.reserve(R*C);

                             for(size_t r=0; r<R; ++r) {
                                 for(size_t c=0; c<C; ++c) {

                                     std::vector<std::shared_ptr<const SymbolicNode>> sum_ops;
                                     for(size_t k=0; k<K; ++k) {
                                         std::vector<std::shared_ptr<const SymbolicNode>> prod_ops = {
                                             d_l[r*K + k], d_r[k*C + c]
                                         };
                                         sum_ops.push_back(make_normalized_multiply_node(prod_ops));
                                     }

                                     NormalizationVisitor elem_vis(assumptions_);
                                     auto elem_node = lamina::detail::make_node<AddNode>(sum_ops);
                                     elem_node->accept(elem_vis);

                                     auto res_val = elem_vis.get_result();
                                     if (auto num = std::dynamic_pointer_cast<const NumberNode>(res_val)) {
                                         if (std::holds_alternative<lmmc_real_t>(num->value())) {
                                              lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                                             int eq_0, eq_1;
                                             lmmc_double_nearly_equal_tol(v, 0.0, 1e-10, 1e-10, &eq_0);
                                             lmmc_double_nearly_equal_tol(v, 1.0, 1e-10, 1e-10, &eq_1);
                                             if (eq_0) {
                                                 res_val = lamina::detail::make_node<NumberNode>(BigInt(0));
                                             } else if (eq_1) {
                                                 res_val = lamina::detail::make_node<NumberNode>(BigInt(1));
                                             }
                                         }
                                     }
                                     res_data.push_back(res_val);
                                 }
                             }

                             fused_ops.pop_back();
                             fused_ops.push_back(lamina::detail::make_node<MatrixNode>(R, C, res_data));
                             continue;
                         }
                    }
                }
                fused_ops.push_back(right);
            }

            if (!scalar_part->is_one()) {
                fused_ops.insert(fused_ops.begin(), scalar_part);
            }

            if (fused_ops.empty()) result = lamina::detail::make_node<NumberNode>(BigInt(1));
            else if (fused_ops.size() == 1) result = fused_ops[0];
            else result = make_normalized_multiply_node(fused_ops);

            return;
        }

        std::shared_ptr<const NumberNode> const_acc = lamina::detail::make_node<NumberNode>(BigInt(1));
        struct FactorAccum {
            std::shared_ptr<const NumberNode> exponent;
            std::vector<std::shared_ptr<const SymbolicNode>> originals;
        };
        std::map<std::shared_ptr<const SymbolicNode>, FactorAccum, NodeCompare> bases;

        for (const auto& op : sc) {
             if (auto num = std::dynamic_pointer_cast<const NumberNode>(op)) {
                 if (num->is_zero()) {
                     result = lamina::detail::make_node<NumberNode>(BigInt(0));
                     return;
                 }
                 const_acc = multiply_numbers(const_acc, num);
             } else {
                 std::shared_ptr<const SymbolicNode> base = op;
                 std::shared_ptr<const NumberNode> exp = lamina::detail::make_node<NumberNode>(BigInt(1));
                 bool is_number_power = false;

                 if (auto pow = std::dynamic_pointer_cast<const PowerNode>(op)) {
                     // Only split base/exponent of a PowerNode if the exponent is a NumberNode.
                     // For symbolic exponents (e.g. 2^x), keep the PowerNode atomic so we
                     // don't silently drop the exponent when accumulating into `bases`.
                     auto e_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
                     if (e_num) {
                         auto pow_base = pow->base();

                         if (auto b_num = std::dynamic_pointer_cast<const NumberNode>(pow_base)) {
                             long long exp_val = 0;
                             bool exp_ok = false;
                             bool exp_is_half = false;

                             if (std::holds_alternative<BigInt>(e_num->value())) {
                                 auto converted =
                                     std::get<BigInt>(e_num->value()).try_to_int64();
                                 if (converted) {
                                     exp_val = static_cast<long long>(*converted);
                                     exp_ok = true;
                                 }
                                 } else if (std::holds_alternative<lmmc_real_t>(e_num->value())) {
                                      lmmc_real_t d = std::get<lmmc_real_t>(e_num->value());
                                     int eq_half;
                                     lmmc_double_nearly_equal_tol(d, 0.5, 1e-9, 1e-9, &eq_half);
                                     if (std::isfinite(d) && d == std::floor(d) &&
                                         d >= static_cast<lmmc_real_t>(
                                             std::numeric_limits<long long>::min()) &&
                                         d <= static_cast<lmmc_real_t>(
                                             std::numeric_limits<long long>::max())) {
                                         exp_val = (long long)d;
                                         exp_ok = true;
                                     } else if (eq_half) {
                                         exp_is_half = true;
                                     }
                                 } else if (std::holds_alternative<Rational>(e_num->value())) {
                                 Rational r = std::get<Rational>(e_num->value());
                                 if (r.get_denominator() == BigInt(1)) {
                                     auto converted = r.get_numerator().try_to_int64();
                                     if (converted) {
                                         exp_val = static_cast<long long>(*converted);
                                         exp_ok = true;
                                     }
                                 } else if (r.get_numerator() == BigInt(1) && r.get_denominator() == BigInt(2)) {
                                     exp_is_half = true;
                                 }
                             }

                             if (exp_ok) {
                                 std::shared_ptr<const NumberNode> pow_val = nullptr;

                                 if (exp_val == -1) {
                                      if (std::holds_alternative<BigInt>(b_num->value())) {
                                          const auto& bi = std::get<BigInt>(b_num->value());
                                          if (!bi.is_zero()) pow_val = lamina::detail::make_node<NumberNode>(Rational(BigInt(1), bi));
                                      } else if (std::holds_alternative<Rational>(b_num->value())) {
                                          Rational r = std::get<Rational>(b_num->value());
                                          if (!r.get_numerator().is_zero()) pow_val = lamina::detail::make_node<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                      } else if (std::holds_alternative<lmmc_real_t>(b_num->value())) {
                                          lmmc_real_t bv = std::get<lmmc_real_t>(b_num->value());
                                          if (bv != 0.0) pow_val = lamina::detail::make_node<NumberNode>(1.0 / bv);
                                      }
                                 } else if (exp_val == 0) {
                                      pow_val = lamina::detail::make_node<NumberNode>(BigInt(1));
                                 } else if (std::abs(exp_val) > 0 && std::abs(exp_val) < 64) {

                                      long long abs_exp = std::abs(exp_val);
                                      std::shared_ptr<const NumberNode> base_pow_val = nullptr;

                                      if (std::holds_alternative<BigInt>(b_num->value())) {
                                          BigInt b = std::get<BigInt>(b_num->value());
                                          BigInt res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = lamina::detail::make_node<NumberNode>(res);
                                      } else if (std::holds_alternative<Rational>(b_num->value())) {
                                          Rational b = std::get<Rational>(b_num->value());
                                          Rational res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = lamina::detail::make_node<NumberNode>(res);
                                      } else if (std::holds_alternative<lmmc_real_t>(b_num->value())) {
                                          lmmc_real_t b = std::get<lmmc_real_t>(b_num->value());
                                          lmmc_real_t res = 1.0;
                                          for(int k=0;k<abs_exp;++k) res *= b;
                                           base_pow_val = lamina::detail::make_node<NumberNode>(res);
                                      }

                                      if (base_pow_val) {
                                          if (exp_val > 0) {
                                              pow_val = base_pow_val;
                                          } else {

                                              if (std::holds_alternative<BigInt>(base_pow_val->value())) {
                                                  const auto& bi = std::get<BigInt>(base_pow_val->value());
                                                  if (!bi.is_zero()) pow_val = lamina::detail::make_node<NumberNode>(Rational(BigInt(1), bi));
                                              } else if (std::holds_alternative<Rational>(base_pow_val->value())) {
                                                  Rational r = std::get<Rational>(base_pow_val->value());
                                                  if (!r.get_numerator().is_zero()) pow_val = lamina::detail::make_node<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                              } else if (std::holds_alternative<lmmc_real_t>(base_pow_val->value())) {
                                                  lmmc_real_t bv = std::get<lmmc_real_t>(base_pow_val->value());
                                                  if (bv != 0.0) pow_val = lamina::detail::make_node<NumberNode>(1.0 / bv);
                                              }
                                          }
                                      }
                                 }

                                 if (pow_val) {
                                     const_acc = multiply_numbers(const_acc, pow_val);
                                     is_number_power = true;
                                 }
                             } else if (exp_is_half) {

                                 std::shared_ptr<const NumberNode> root_val = nullptr;
                                 if (std::holds_alternative<BigInt>(b_num->value())) {
                                     const BigInt value = std::get<BigInt>(b_num->value());
                                     if (value >= BigInt(0)) {
                                         const BigInt root = value.sqrt();
                                         if (root * root == value) {
                                             root_val = lamina::detail::make_node<NumberNode>(root);
                                         }
                                     }
                                 } else if (std::holds_alternative<Rational>(b_num->value())) {
                                     const Rational value = std::get<Rational>(b_num->value());
                                     if (value >= Rational(0)) {
                                         const BigInt numerator_root =
                                             value.get_numerator().sqrt();
                                         const BigInt denominator_root =
                                             value.get_denominator().sqrt();
                                         if (numerator_root * numerator_root ==
                                                 value.get_numerator() &&
                                             denominator_root * denominator_root ==
                                                 value.get_denominator()) {
                                             root_val = lamina::detail::make_node<NumberNode>(Rational(
                                                 numerator_root, denominator_root));
                                         }
                                     }
                                 } else if (std::holds_alternative<lmmc_real_t>(b_num->value())) {
                                     double d = std::get<lmmc_real_t>(b_num->value());
                                     if (d >= 0) root_val = lamina::detail::make_node<NumberNode>(std::sqrt(d));
                                 }

                                 if (root_val) {
                                     const_acc = multiply_numbers(const_acc, root_val);
                                     is_number_power = true;
                                 }
                             }
                         }

                         if (!is_number_power && is_positive_integer_number(e_num)) {
                             base = pow_base;
                             exp = e_num;
                         }
                     }
                 }

                 if (!is_number_power) {
                     auto it = bases.find(base);
                     if (it == bases.end()) {
                         bases.emplace(base, FactorAccum{exp, {op}});
                     } else {
                         it->second.exponent = add_numbers(it->second.exponent, exp);
                         it->second.originals.push_back(op);
                     }
                 }
             }
        }

        std::vector<std::shared_ptr<const SymbolicNode>> final_ops;
        if (!const_acc->is_one()) {
             final_ops.push_back(const_acc);
        }

        std::vector<std::shared_ptr<const SymbolicNode>> var_ops;
        for (auto const& [base, acc] : bases) {
            const auto& exp = acc.exponent;
            if (exp->is_zero()) {
                var_ops.insert(var_ops.end(), acc.originals.begin(), acc.originals.end());
            } else if (exp->is_one()) {
                var_ops.push_back(base);
            } else {
                var_ops.push_back(lamina::detail::make_node<PowerNode>(base, exp));
            }
        }

        std::sort(var_ops.begin(), var_ops.end(), [](const std::shared_ptr<const SymbolicNode>& l, const std::shared_ptr<const SymbolicNode>& r) {
            int d1 = get_node_degree_helper(l);
            int d2 = get_node_degree_helper(r);
            if (d1 != d2) return d1 > d2;
            return l->compare(*r) < 0;
        });

        final_ops.insert(final_ops.end(), var_ops.begin(), var_ops.end());

        if (final_ops.empty()) result = lamina::detail::make_node<NumberNode>(BigInt(1));
        else if (final_ops.size() == 1) result = final_ops[0];
        else result = make_normalized_multiply_node(final_ops);
    }
