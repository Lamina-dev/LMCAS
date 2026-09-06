#include "visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/squared_norm.hpp"

namespace LMCAS {

bool NormalizationVisitor::try_normalize_squared_norm(
    const SymbolicNode& node, std::shared_ptr<const SymbolicNode>& argument) {
    const auto squares = LMCAS::detail::squared_norm_terms(node);
    if (!squares[0]) return false;
    squares[0]->base()->accept(*this);
    auto first = result;
    std::shared_ptr<const SymbolicNode> second;
    if (squares[1]) {
        squares[1]->base()->accept(*this);
        second = result;
    }
    const auto* a = dynamic_cast<const NumberNode*>(first.get());
    const auto* b = dynamic_cast<const NumberNode*>(second.get());
    const bool numeric = a && (!second || b);
    const bool approximate = numeric &&
        (std::holds_alternative<lmmc_real_t>(a->value()) ||
         (b && std::holds_alternative<lmmc_real_t>(b->value())));
    if (approximate) {
        const auto value = [](const NumberNode& number) {
            const auto& v = number.value();
            if (const auto* real = std::get_if<lmmc_real_t>(&v)) return *real;
            if (const auto* rational = std::get_if<Rational>(&v)) return rational->to_double();
            return std::get<BigInt>(v).to_double();
        };
        const double magnitude = b ? std::hypot(value(*a), value(*b)) : std::abs(value(*a));
        if (std::isfinite(magnitude)) {
            result = LMCAS::detail::make_node<NumberNode>(magnitude);
            return true;
        }
    }
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    terms.reserve(second ? 2 : 1);
    if (!first->is_zero() || !second || second->is_zero()) {
        terms.push_back(LMCAS::detail::make_node<PowerNode>(first, squares[0]->exponent()));
    }
    if (second && !second->is_zero()) {
        terms.push_back(LMCAS::detail::make_node<PowerNode>(second, squares[1]->exponent()));
    }
    argument = terms.size() == 1 ? terms.front()
                                : LMCAS::detail::make_node<AddNode>(std::move(terms));
    if ((numeric && !approximate) ||
        dynamic_cast<const ComplexNode*>(first.get()) ||
        dynamic_cast<const ComplexNode*>(second.get())) {
        argument->accept(*this);
        argument = result;
        return false;
    }
    if (const auto* power = dynamic_cast<const PowerNode*>(&node)) {
        result = LMCAS::detail::make_node<PowerNode>(argument, power->exponent());
    } else {
        auto norm = LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sqrt,
            std::vector<std::shared_ptr<const SymbolicNode>>{argument});
        auto simplified = try_assumption_simplify(norm);
        result = simplified ? simplified : norm;
    }
    return true;
}

void NormalizationVisitor::visit(const PowerNode& node) {
        std::shared_ptr<const SymbolicNode> s_base;
        if (try_normalize_squared_norm(node, s_base)) return;
        if (!s_base) {
            node.base()->accept(*this);
            s_base = result;
        }
        node.exponent()->accept(*this);
        auto s_exp = result;

        if (s_exp->is_zero()) {
            if (s_base->is_zero()) {
                result = LMCAS::detail::make_node<PowerNode>(s_base, s_exp);
                return;
            }
            if (is_provably_nonzero(s_base)) {
                result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                return;
            }
            result = LMCAS::detail::make_node<PowerNode>(s_base, s_exp);
            return;
        }
        if (s_exp->is_one()) {
            result = s_base;
            return;
        }

        if (s_base->is_zero()) {
             // 0^x = 0 only when x is provably positive; otherwise keep the
             // PowerNode so domain issues (e.g. 0^(-1/2)) are not silently
             // turned into 0. This mirrors SymbolicFactory::create_power.
             if (s_exp->is_positive()) {
                 result = LMCAS::detail::make_node<NumberNode>(BigInt(0));
                 return;
             }
             result = LMCAS::detail::make_node<PowerNode>(s_base, s_exp);
             return;
        }
        if (s_base->is_one()) {
            result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
            return;
        }

        if (auto b_num = std::dynamic_pointer_cast<const NumberNode>(s_base)) {
            if (auto e_num = std::dynamic_pointer_cast<const NumberNode>(s_exp)) {
                 long long exp_val = 0;
                 const bool exp_ok = try_get_integer_value(e_num, exp_val);
                 bool exp_is_half = false;
                 if (!exp_ok) {
                     if (const auto* real = std::get_if<lmmc_real_t>(&e_num->value())) {
                         exp_is_half = *real == 0.5;
                     } else if (const auto* rational = std::get_if<Rational>(&e_num->value())) {
                         exp_is_half = *rational == Rational(1, 2);
                     }
                 }

                 if (exp_ok) {
                     if (const auto* real = std::get_if<lmmc_real_t>(&b_num->value());
                         real && exp_val > -64 && exp_val < 64) {
                         // Floating powers retain their final range; exact powers
                         // continue through integer and rational arithmetic below.
                         const double value = std::pow(*real, static_cast<double>(exp_val));
                         if (std::isfinite(value)) {
                             result = LMCAS::detail::make_node<NumberNode>(value);
                         } else {
                             result = LMCAS::detail::make_node<PowerNode>(s_base, s_exp);
                         }
                         return;
                     }
                     std::shared_ptr<const NumberNode> pow_val = nullptr;
                     if (exp_val == -1) {
                          if (std::holds_alternative<BigInt>(b_num->value())) {
                              pow_val = LMCAS::detail::make_node<NumberNode>(Rational(BigInt(1), std::get<BigInt>(b_num->value())));
                          } else if (std::holds_alternative<Rational>(b_num->value())) {
                              Rational r = std::get<Rational>(b_num->value());
                              if (!r.get_numerator().is_zero()) pow_val = LMCAS::detail::make_node<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                          }
                     } else if (exp_val == 0) {
                          pow_val = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                     } else if (exp_val > -64 && exp_val < 64) {

                                      long long abs_exp = std::abs(exp_val);
                                      std::shared_ptr<const NumberNode> base_pow_val = nullptr;

                                      if (std::holds_alternative<BigInt>(b_num->value())) {
                                          BigInt b = std::get<BigInt>(b_num->value());
                                          BigInt res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = LMCAS::detail::make_node<NumberNode>(res);
                                      } else if (std::holds_alternative<Rational>(b_num->value())) {
                                          Rational b = std::get<Rational>(b_num->value());
                                          Rational res(1);
                                          for(int k=0;k<abs_exp;++k) res = res * b;
                                          base_pow_val = LMCAS::detail::make_node<NumberNode>(res);
                                      }

                                      if (base_pow_val) {
                                          if (exp_val > 0) {
                                              pow_val = base_pow_val;
                                          } else {

                                              if (std::holds_alternative<BigInt>(base_pow_val->value())) {
                                                  pow_val = LMCAS::detail::make_node<NumberNode>(Rational(BigInt(1), std::get<BigInt>(base_pow_val->value())));
                                              } else if (std::holds_alternative<Rational>(base_pow_val->value())) {
                                                  Rational r = std::get<Rational>(base_pow_val->value());
                                                  pow_val = LMCAS::detail::make_node<NumberNode>(Rational(r.get_denominator(), r.get_numerator()));
                                              }
                                          }
                                      }
                     }

                     if (pow_val) {
                         result = pow_val;
                         return;
                     }
                 } else if (exp_is_half) {
                     if (std::holds_alternative<BigInt>(b_num->value())) {
                         const BigInt val = std::get<BigInt>(b_num->value());
                         if (val >= BigInt(0)) {
                             const BigInt root = val.sqrt();
                             if (root * root == val) {
                                 result = LMCAS::detail::make_node<NumberNode>(root);
                                 return;
                             }
                             if (val > BigInt(0) && val < BigInt(1000000)) {
                                 const long long v_ll = val.to_int();
                                 const long long approximate_root =
                                     static_cast<long long>(std::sqrt(
                                         static_cast<double>(v_ll)));
                                 for (long long i = approximate_root; i >= 2; --i) {
                                     if (v_ll % (i*i) == 0) {
                                         long long s_ll = i;
                                         long long k_ll = v_ll / (i*i);

                                         auto s_node = LMCAS::detail::make_node<NumberNode>(BigInt(s_ll));
                                         auto k_node = LMCAS::detail::make_node<NumberNode>(BigInt(k_ll));
                                         auto half_node = LMCAS::detail::make_node<NumberNode>(Rational(1, 2));
                                         auto pow_node = LMCAS::detail::make_node<PowerNode>(k_node, half_node);

                                         result = SymbolicFactory::create_multiply({s_node, pow_node});
                                         return;
                                     }
                                 }
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
                                 result = LMCAS::detail::make_node<NumberNode>(Rational(
                                     numerator_root, denominator_root));
                                 return;
                             }
                         }
                     } else if (std::holds_alternative<lmmc_real_t>(b_num->value())) {
                         double d = std::get<lmmc_real_t>(b_num->value());
                         if (d >= 0) {
                             result = LMCAS::detail::make_node<NumberNode>(std::sqrt(d));
                             return;
                         }
                     }
                 }

                 /// Simplify powers of imaginary unit: (-1)^(n/2) for odd n
                 /// i² = -1, i³ = -i, i⁴ = 1 (handled via (-1)^(n/2) reduction)
                 if (std::holds_alternative<Rational>(e_num->value())) {
                     Rational exp_r = std::get<Rational>(e_num->value());
                     bool base_is_neg_one = false;
                     if (std::holds_alternative<BigInt>(b_num->value())) {
                         base_is_neg_one = (std::get<BigInt>(b_num->value()) == BigInt(-1));
                     } else if (std::holds_alternative<Rational>(b_num->value())) {
                         base_is_neg_one = (std::get<Rational>(b_num->value()) == Rational(-1));
                     } else if (std::holds_alternative<lmmc_real_t>(b_num->value())) {
                         base_is_neg_one = std::get<lmmc_real_t>(b_num->value()) == -1.0;
                     }

                     if (base_is_neg_one && exp_r.get_denominator() == BigInt(2)) {
                         /// (-1)^(n/2) where n is the numerator
                         BigInt remainder = exp_r.get_numerator() % BigInt(4);
                         if (remainder.IsNegative()) remainder = remainder + BigInt(4);
                         const long long r_mod = remainder.to_int();
                         /// i = (-1)^(1/2)
                         auto i_node = LMCAS::detail::make_node<PowerNode>(
                             LMCAS::detail::make_node<NumberNode>(BigInt(-1)),
                             LMCAS::detail::make_node<NumberNode>(Rational(1, 2)));
                         if (r_mod == 0) {
                             /// i⁴ = 1
                             result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                             return;
                         } else if (r_mod == 1) {
                             /// i¹ = i (keep as (-1)^(1/2))
                             result = i_node;
                             return;
                         } else if (r_mod == 2) {
                             /// i² = -1
                             result = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
                             return;
                         } else { // r_mod == 3
                             /// i³ = -i
                             std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {
                                 LMCAS::detail::make_node<NumberNode>(BigInt(-1)), i_node
                             };
                             result = make_normalized_multiply_node(mul_ops);
                             return;
                         }
                     }
                 }
            }
        }

        if (auto m_base = std::dynamic_pointer_cast<const MultiplyNode>(s_base)) {
            std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
            for(auto& op : m_base->operands()) {

                auto term_pow = LMCAS::detail::make_node<PowerNode>(op, s_exp);
                term_pow->accept(*this);

                if (auto mul_res = std::dynamic_pointer_cast<const MultiplyNode>(result)) {
                     new_ops.insert(new_ops.end(), mul_res->operands().begin(), mul_res->operands().end());
                } else {
                     new_ops.push_back(result);
                }
            }

            auto final_mul = make_normalized_multiply_node(new_ops);
            final_mul->accept(*this);
            return;
        }

        if (auto p_base = std::dynamic_pointer_cast<const PowerNode>(s_base)) {
             auto inner_exp = std::dynamic_pointer_cast<const NumberNode>(p_base->exponent());
             auto outer_exp = std::dynamic_pointer_cast<const NumberNode>(s_exp);
             long long outer_integer = 0;
             bool outer_is_integer = try_get_integer_value(outer_exp, outer_integer);
             if (outer_is_integer) {
                 auto inner_base_num = std::dynamic_pointer_cast<const NumberNode>(p_base->base());
                 bool inner_base_is_neg_one = false;
                 if (inner_base_num) {
                     if (std::holds_alternative<BigInt>(inner_base_num->value())) {
                         inner_base_is_neg_one = (std::get<BigInt>(inner_base_num->value()) == BigInt(-1));
                     } else if (std::holds_alternative<Rational>(inner_base_num->value())) {
                         inner_base_is_neg_one = (std::get<Rational>(inner_base_num->value()) == Rational(-1));
                     } else if (std::holds_alternative<lmmc_real_t>(inner_base_num->value())) {
                         inner_base_is_neg_one =
                             std::get<lmmc_real_t>(inner_base_num->value()) == -1.0;
                     }
                 }

                 bool inner_exp_is_half = false;
                 if (inner_exp && std::holds_alternative<Rational>(inner_exp->value())) {
                     const auto& inner_r = std::get<Rational>(inner_exp->value());
                     inner_exp_is_half = (inner_r.get_numerator() == BigInt(1) &&
                                          inner_r.get_denominator() == BigInt(2));
                 }

                 if (inner_base_is_neg_one && inner_exp_is_half) {
                     long long r_mod = ((outer_integer % 4) + 4) % 4;
                     auto i_node = LMCAS::detail::make_node<PowerNode>(
                         LMCAS::detail::make_node<NumberNode>(BigInt(-1)),
                         LMCAS::detail::make_node<NumberNode>(Rational(1, 2)));
                     if (r_mod == 0) {
                         result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                         return;
                     }
                     if (r_mod == 1) {
                         result = i_node;
                         return;
                     }
                     if (r_mod == 2) {
                         result = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
                         return;
                     }
                     result = make_normalized_multiply_node(
                         std::vector<std::shared_ptr<const SymbolicNode>>{
                             LMCAS::detail::make_node<NumberNode>(BigInt(-1)), i_node});
                     return;
                 }
             }

             if (is_positive_integer_number(inner_exp) && is_positive_integer_number(outer_exp)) {
                 std::vector<std::shared_ptr<const SymbolicNode>> exp_ops;
                 exp_ops.push_back(p_base->exponent());
                 exp_ops.push_back(s_exp);

                 auto mul_exp = make_normalized_multiply_node(exp_ops);
                 mul_exp->accept(*this);

                 auto new_pow = LMCAS::detail::make_node<PowerNode>(p_base->base(), result);
                 new_pow->accept(*this);
                 return;
             }

             result = LMCAS::detail::make_node<PowerNode>(s_base, s_exp);
             return;
        }

        result = LMCAS::detail::make_node<PowerNode>(s_base, s_exp);
    }
void NormalizationVisitor::visit(const FunctionNode& node) {
        std::shared_ptr<const SymbolicNode> argument;
        if (try_normalize_squared_norm(node, argument)) return;
        std::vector<std::shared_ptr<const SymbolicNode>> s_args;
        if (argument) {
            s_args.push_back(std::move(argument));
        } else {
            for (const auto& a : node.arguments()) {
                a->accept(*this);
                s_args.push_back(result);
            }
        }

        if (s_args.size() == 1) {
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(s_args[0])) {
                const bool approximate =
                    std::holds_alternative<lmmc_real_t>(num->value());
                const lmmc_real_t val = approximate
                    ? std::get<lmmc_real_t>(num->value())
                    : 0.0;
                const bool exact_zero =
                    (std::holds_alternative<BigInt>(num->value()) &&
                     std::get<BigInt>(num->value()).is_zero()) ||
                    (std::holds_alternative<Rational>(num->value()) &&
                     std::get<Rational>(num->value()) == Rational(0));
                const bool exact_one =
                    (std::holds_alternative<BigInt>(num->value()) &&
                     std::get<BigInt>(num->value()) == BigInt(1)) ||
                    (std::holds_alternative<Rational>(num->value()) &&
                     std::get<Rational>(num->value()) == Rational(1));

                switch (node.type()) {
                    case FunctionNode::FuncType::Sin:
                    {
                        if (exact_zero || (approximate && val == 0.0)) {
                            result = LMCAS::detail::make_node<NumberNode>(BigInt(0));
                            return;
                        }
                        break;
                    }
                    case FunctionNode::FuncType::Cos:
                    {
                        if (exact_zero || (approximate && val == 0.0)) {
                            result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                            return;
                        }
                        break;
                    }
                    case FunctionNode::FuncType::Tan:
                    {
                        if (exact_zero || (approximate && val == 0.0)) {
                            result = LMCAS::detail::make_node<NumberNode>(BigInt(0));
                            return;
                        }
                        break;
                    }
                    case FunctionNode::FuncType::Exp:
                    {
                        if (exact_zero || (approximate && val == 0.0)) {
                            result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                            return;
                        }
                        break;
                    }
                    case FunctionNode::FuncType::Ln:
                    {
                        if (exact_one || (approximate && val == 1.0)) {
                            result = LMCAS::detail::make_node<NumberNode>(BigInt(0));
                            return;
                        }
                        if (exact_zero || (approximate && val == 0.0)) {
                             std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
                             auto inf = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
                             std::vector<std::shared_ptr<const SymbolicNode>> m_args = {LMCAS::detail::make_node<NumberNode>(BigInt(-1)), inf};
                             result = make_normalized_multiply_node(m_args);
                             return;
                        }
                        break;
                    }
                    case FunctionNode::FuncType::Log:
                        if (s_args.size() == 2) {

                             std::vector<std::shared_ptr<const SymbolicNode>> args_x = { s_args[0] };
                             std::vector<std::shared_ptr<const SymbolicNode>> args_b = { s_args[1] };
                             auto ln_x = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, args_x);
                             auto ln_b = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, args_b);
                             auto ln_b_inv = LMCAS::detail::make_node<PowerNode>(ln_b, LMCAS::detail::make_node<NumberNode>(BigInt(-1)));
                             std::vector<std::shared_ptr<const SymbolicNode>> m_args = { ln_x, ln_b_inv };
                             auto prod = make_normalized_multiply_node(m_args);

                             NormalizationVisitor v(assumptions_);
                             prod->accept(v);
                             result = v.get_result();
                             return;
                        }
                        break;
                    case FunctionNode::FuncType::Abs:
                        if (std::holds_alternative<BigInt>(num->value())) {
                            result = LMCAS::detail::make_node<NumberNode>(
                                std::get<BigInt>(num->value()).Abs());
                        } else if (std::holds_alternative<Rational>(num->value())) {
                            const Rational value = std::get<Rational>(num->value());
                            result = LMCAS::detail::make_node<NumberNode>(
                                value.get_numerator().IsNegative()
                                    ? Rational(value.get_numerator().Abs(),
                                               value.get_denominator())
                                    : value);
                        } else {
                            result = LMCAS::detail::make_node<NumberNode>(std::abs(val));
                        }
                        return;
                    case FunctionNode::FuncType::Sqrt:
                    {
                        if (std::holds_alternative<BigInt>(num->value())) {
                            const BigInt n = std::get<BigInt>(num->value());
                            if (n >= BigInt(0)) {
                                const BigInt root = n.sqrt();
                                if (root * root == n) {
                                    result = LMCAS::detail::make_node<NumberNode>(root);
                                    return;
                                }
                            }
                            result = LMCAS::detail::make_node<FunctionNode>(node.type(), s_args);
                            return;
                        } else if (std::holds_alternative<Rational>(num->value())) {
                            const Rational value = std::get<Rational>(num->value());
                            if (value >= Rational(0)) {
                                const BigInt numerator_root =
                                    value.get_numerator().sqrt();
                                const BigInt denominator_root =
                                    value.get_denominator().sqrt();
                                if (numerator_root * numerator_root ==
                                        value.get_numerator() &&
                                    denominator_root * denominator_root ==
                                        value.get_denominator()) {
                                    result = LMCAS::detail::make_node<NumberNode>(Rational(
                                        numerator_root, denominator_root));
                                    return;
                                }
                            }
                            result = LMCAS::detail::make_node<FunctionNode>(node.type(), s_args);
                            return;
                        }

                        if (std::holds_alternative<lmmc_real_t>(num->value()) && val >= 0) {
                             result = LMCAS::detail::make_node<NumberNode>(std::sqrt(val));
                             return;
                        }
                    }
                        break;
                     case FunctionNode::FuncType::LambertW:
                     {
                         if (exact_zero) {
                             result = LMCAS::detail::make_node<NumberNode>(BigInt(0));
                             return;
                         }
                         if (approximate) {
                             lmmc_real_t w_res;
                             if (lmmc_lambertw(val, &w_res) == LMMC_STATUS_OK) {
                                 result = LMCAS::detail::make_node<NumberNode>(w_res);
                                 return;
                             }
                         }

                         break;
                     }
                     default: break;
                }
            }
        }

        if (node.type() == FunctionNode::FuncType::Ln && s_args.size() == 1) {
             if (auto pow = std::dynamic_pointer_cast<const PowerNode>(s_args[0])) {
                  auto y = pow->exponent();
                  auto x = pow->base();

                  std::vector<std::shared_ptr<const SymbolicNode>> ln_args = { x };
                  auto ln_x = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, ln_args);

                  std::vector<std::shared_ptr<const SymbolicNode>> m_args = { y, ln_x };
                  auto prod = make_normalized_multiply_node(m_args);

                  NormalizationVisitor v(assumptions_);
                  prod->accept(v);
                  result = v.get_result();
                  return;
             }
             if (auto func = std::dynamic_pointer_cast<const FunctionNode>(s_args[0])) {
                 if (func->type() == FunctionNode::FuncType::Exp && func->arguments().size() == 1) {
                     result = func->arguments()[0];
                     return;
                 }
             }
        }

        if (node.type() == FunctionNode::FuncType::Log && s_args.size() == 2) {
             if (auto pow = std::dynamic_pointer_cast<const PowerNode>(s_args[0])) {
                 if (pow->base()->compare(*s_args[1]) == 0) {
                     result = pow->exponent();
                     return;
                 }
             }
             std::vector<std::shared_ptr<const SymbolicNode>> args_x = { s_args[0] };
             std::vector<std::shared_ptr<const SymbolicNode>> args_b = { s_args[1] };
             auto ln_x = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, args_x);
             auto ln_b = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, args_b);
             auto ln_b_inv = LMCAS::detail::make_node<PowerNode>(ln_b, LMCAS::detail::make_node<NumberNode>(BigInt(-1)));
             std::vector<std::shared_ptr<const SymbolicNode>> m_args = { ln_x, ln_b_inv };
             auto prod = make_normalized_multiply_node(m_args);
             NormalizationVisitor v(assumptions_);
             prod->accept(v);
             result = v.get_result();
             return;
        }

        if (s_args.size() == 1) {
            Rational k_pi_val;
            if (get_pi_coeff(s_args[0], k_pi_val)) {
                BigInt n = k_pi_val.get_numerator();
                BigInt d = k_pi_val.get_denominator();

                BigInt two_d = d * BigInt(2);
                BigInt reduced_n = n % two_d;
                if (reduced_n.IsNegative()) reduced_n = reduced_n + two_d;

                auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
                auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                auto minus_one = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
                auto half = LMCAS::detail::make_node<NumberNode>(Rational(1, 2));
                auto minus_half = LMCAS::detail::make_node<NumberNode>(Rational(-1, 2));

                auto root2 = LMCAS::detail::make_node<PowerNode>(LMCAS::detail::make_node<NumberNode>(BigInt(2)), LMCAS::detail::make_node<NumberNode>(Rational(1, 2)));
                std::vector<std::shared_ptr<const SymbolicNode>> half_root2_args = { half, root2 };
                auto half_root2 = make_normalized_multiply_node(half_root2_args);
                std::vector<std::shared_ptr<const SymbolicNode>> minus_half_root2_args = { minus_half, root2 };
                auto minus_half_root2 = make_normalized_multiply_node(minus_half_root2_args);

                auto root3 = LMCAS::detail::make_node<PowerNode>(LMCAS::detail::make_node<NumberNode>(BigInt(3)), LMCAS::detail::make_node<NumberNode>(Rational(1, 2)));
                std::vector<std::shared_ptr<const SymbolicNode>> half_root3_args = { half, root3 };
                auto half_root3 = make_normalized_multiply_node(half_root3_args);
                std::vector<std::shared_ptr<const SymbolicNode>> minus_half_root3_args = { minus_half, root3 };
                auto minus_half_root3 = make_normalized_multiply_node(minus_half_root3_args);

                std::vector<std::shared_ptr<const SymbolicNode>> minus_root3_args = { minus_one, root3 };
                auto minus_root3 = make_normalized_multiply_node(minus_root3_args);

                auto third = LMCAS::detail::make_node<NumberNode>(Rational(1, 3));
                std::vector<std::shared_ptr<const SymbolicNode>> third_root3_args = { third, root3 };
                auto third_root3 = make_normalized_multiply_node(third_root3_args);
                auto minus_third = LMCAS::detail::make_node<NumberNode>(Rational(-1, 3));
                std::vector<std::shared_ptr<const SymbolicNode>> minus_third_root3_args = { minus_third, root3 };
                auto minus_third_root3 = make_normalized_multiply_node(minus_third_root3_args);

                if (d == BigInt(1)) {

                    if (node.type() == FunctionNode::FuncType::Sin || node.type() == FunctionNode::FuncType::Tan) {
                        result = zero; return;
                    } else if (node.type() == FunctionNode::FuncType::Cos) {
                        if (reduced_n == BigInt(0)) result = one;
                        else result = minus_one;
                        return;
                    }
                } else if (d == BigInt(2)) {

                    if (node.type() == FunctionNode::FuncType::Sin) {
                        if (reduced_n == BigInt(1)) result = one;
                        else if (reduced_n == BigInt(3)) result = minus_one;
                        return;
                    } else if (node.type() == FunctionNode::FuncType::Cos) {
                        if (reduced_n == BigInt(1) || reduced_n == BigInt(3)) {
                            result = zero; return;
                        }
                    }

                } else if (d == BigInt(3)) {

                    if (node.type() == FunctionNode::FuncType::Sin) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(2)) result = half_root3;
                         else if (reduced_n == BigInt(4) || reduced_n == BigInt(5)) result = minus_half_root3;
                         return;
                    } else if (node.type() == FunctionNode::FuncType::Cos) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(5)) result = half;
                         else if (reduced_n == BigInt(2) || reduced_n == BigInt(4)) result = minus_half;
                         return;
                    } else if (node.type() == FunctionNode::FuncType::Tan) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(4)) result = root3;
                         else if (reduced_n == BigInt(2) || reduced_n == BigInt(5)) result = minus_root3;
                         return;
                    }
                } else if (d == BigInt(4)) {

                     if (node.type() == FunctionNode::FuncType::Sin) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(3)) result = half_root2;
                         else if (reduced_n == BigInt(5) || reduced_n == BigInt(7)) result = minus_half_root2;
                         return;
                    } else if (node.type() == FunctionNode::FuncType::Cos) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(7)) result = half_root2;
                         else if (reduced_n == BigInt(3) || reduced_n == BigInt(5)) result = minus_half_root2;
                         return;
                    } else if (node.type() == FunctionNode::FuncType::Tan) {
                         if (reduced_n == BigInt(1) || reduced_n == BigInt(5)) result = one;
                         else if (reduced_n == BigInt(3) || reduced_n == BigInt(7)) result = minus_one;
                         return;
                    }
                } else if (d == BigInt(6)) {

                    if (node.type() == FunctionNode::FuncType::Sin) {
                        if (reduced_n == BigInt(1) || reduced_n == BigInt(5)) result = half;
                        else if (reduced_n == BigInt(7) || reduced_n == BigInt(11)) result = minus_half;
                        return;
                    } else if (node.type() == FunctionNode::FuncType::Cos) {
                        if (reduced_n == BigInt(1) || reduced_n == BigInt(11)) result = half_root3;
                        else if (reduced_n == BigInt(5) || reduced_n == BigInt(7)) result = minus_half_root3;
                        return;
                    } else if (node.type() == FunctionNode::FuncType::Tan) {
                        if (reduced_n == BigInt(1) || reduced_n == BigInt(7)) result = third_root3;
                        else if (reduced_n == BigInt(5) || reduced_n == BigInt(11)) result = minus_third_root3;
                        return;
                    }
                }
            }
        }

        if (s_args.size() == 1) {
             std::shared_ptr<const SymbolicNode> pos_arg = nullptr;

             if (check_negative_arg(s_args[0], pos_arg)) {
                 if (node.type() == FunctionNode::FuncType::Sin ||
                     node.type() == FunctionNode::FuncType::Tan ||
                     node.type() == FunctionNode::FuncType::ArcTan ||
                     node.type() == FunctionNode::FuncType::ArcSin) {

                     std::vector<std::shared_ptr<const SymbolicNode>> new_args = { pos_arg };
                     auto new_func = LMCAS::detail::make_node<FunctionNode>(node.type(), new_args);
                     std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = { LMCAS::detail::make_node<NumberNode>(BigInt(-1)), new_func };
                     result = make_normalized_multiply_node(mul_ops);
                     return;
                 } else if (node.type() == FunctionNode::FuncType::Cos) {

                     std::vector<std::shared_ptr<const SymbolicNode>> new_args = { pos_arg };
                     result = LMCAS::detail::make_node<FunctionNode>(node.type(), new_args);
                     return;
                 }
             }
        }

        if (node.type() == FunctionNode::FuncType::Log && s_args.size() == 2) {
             auto val = s_args[0];
             auto base = s_args[1];

             if (val->equals(*base)) {
                 result = LMCAS::detail::make_node<NumberNode>(BigInt(1));
                 return;
             }

             if (val->is_one()) {
                 result = LMCAS::detail::make_node<NumberNode>(BigInt(0));
                 return;
             }

             if (auto pow = std::dynamic_pointer_cast<const PowerNode>(val)) {
                 if (pow->base()->equals(*base)) {
                     result = pow->exponent();
                     return;
                 }
             }
        }

        // Attempt assumption-based simplification before falling through
        auto func_node = LMCAS::detail::make_node<FunctionNode>(node.type(), s_args);
        if (auto simplified = try_assumption_simplify(func_node)) {
            // Recursively normalize the simplified result
            simplified->accept(*this);
            return;
        }
        result = func_node;
    }
void NormalizationVisitor::visit(const UninterpretedFunctionNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> arguments;
    arguments.reserve(node.arguments().size());
    for (const auto& argument : node.arguments()) {
        argument->accept(*this);
        arguments.push_back(result);
    }
    result = LMCAS::detail::make_node<UninterpretedFunctionNode>(
        node.name(), std::move(arguments));
}
bool NormalizationVisitor::try_get_integer_value(const std::shared_ptr<const NumberNode>& node, long long& value) {
        if (!node) return false;
        if (std::holds_alternative<BigInt>(node->value())) {
            auto converted = std::get<BigInt>(node->value()).try_to_int64();
            if (!converted) return false;
            value = static_cast<long long>(*converted);
            return true;
        }
        if (std::holds_alternative<Rational>(node->value())) {
            const auto& rational = std::get<Rational>(node->value());
            if (!(rational.get_denominator() == BigInt(1))) return false;
            auto converted = rational.get_numerator().try_to_int64();
            if (!converted) return false;
            value = static_cast<long long>(*converted);
            return true;
        }
        const auto real = std::get<lmmc_real_t>(node->value());
        if (!std::isfinite(real) || std::floor(real) != real ||
            real < static_cast<lmmc_real_t>(
                std::numeric_limits<long long>::min()) ||
            real >= -static_cast<lmmc_real_t>(
                std::numeric_limits<long long>::min())) return false;
        value = static_cast<long long>(real);
        return true;
    }
bool NormalizationVisitor::is_positive_integer_number(const std::shared_ptr<const NumberNode>& node) {
        if (!node) return false;
        if (std::holds_alternative<BigInt>(node->value())) {
            const auto& value = std::get<BigInt>(node->value());
            return !value.IsNegative() && !(value == BigInt(0));
        }
        if (std::holds_alternative<Rational>(node->value())) {
            const auto& value = std::get<Rational>(node->value());
            return value.get_denominator() == BigInt(1) && value > Rational(0);
        }
        const auto value = std::get<lmmc_real_t>(node->value());
        return std::isfinite(value) && value > 0.0 && std::floor(value) == value;
    }
bool NormalizationVisitor::is_provably_nonzero(const std::shared_ptr<const SymbolicNode>& node) const {
        if (!node) return false;
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
            return !num->is_zero();
        }
        if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
            return is_provably_nonzero(complex->real()) || is_provably_nonzero(complex->imag());
        }
        if (!assumptions_) return false;
        auto expr = LMCAS::detail::expression_from_node(node);
        auto nonzero = assumptions_->is_nonzero(expr);
        return nonzero &&
               nonzero.value() == LMCAS::Tribool::True;
    }
std::shared_ptr<const SymbolicNode> NormalizationVisitor::try_assumption_simplify(const std::shared_ptr<const SymbolicNode>& node) {
        if (!assumptions_) return nullptr;

        auto func = std::dynamic_pointer_cast<const FunctionNode>(node);
        if (!func || func->arguments().size() != 1) return nullptr;

        const auto& arg = func->arguments()[0];

        // Rule: sqrt(x²) → x when x is NonNegative
        // Rule: sqrt(x²) → abs(x) when x is Real (but not NonNegative)
        if (func->type() == FunctionNode::FuncType::Sqrt) {
            // Check if argument is a PowerNode with exponent 2
            auto pow = std::dynamic_pointer_cast<const PowerNode>(arg);
            if (pow) {
                auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
                if (exp_num) {
                    bool is_exp_two = false;
                    if (std::holds_alternative<BigInt>(exp_num->value())) {
                        is_exp_two = (std::get<BigInt>(exp_num->value()) == BigInt(2));
                    } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
                        is_exp_two = std::get<lmmc_real_t>(exp_num->value()) == 2.0;
                    } else if (std::holds_alternative<Rational>(exp_num->value())) {
                        is_exp_two = (std::get<Rational>(exp_num->value()) == Rational(2));
                    }

                    if (is_exp_two) {
                        // We have sqrt(base²) — query the base's properties
                        auto base_expr = LMCAS::detail::expression_from_node(pow->base());
                        auto nonnegative =
                            assumptions_->is_nonnegative(base_expr);
                        if (nonnegative &&
                            nonnegative.value() == LMCAS::Tribool::True) {
                            // sqrt(x²) → x when x is NonNegative
                            return pow->base();
                        }

                        auto real = assumptions_->is_real(base_expr);
                        if (real &&
                            real.value() == LMCAS::Tribool::True) {
                            // sqrt(x²) → abs(x) when x is Real (but not NonNegative)
                            std::vector<std::shared_ptr<const SymbolicNode>> abs_args = { pow->base() };
                            return LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Abs, abs_args);
                        }
                    }
                }
            }
        }

        // Rule: abs(x) → x when x is Positive
        // Rule: abs(x) → -x when x is Negative
        if (func->type() == FunctionNode::FuncType::Abs) {
            auto arg_expr = LMCAS::detail::expression_from_node(arg);
            auto positive = assumptions_->is_positive(arg_expr);
            if (positive &&
                positive.value() == LMCAS::Tribool::True) {
                // abs(x) → x when x is Positive
                return arg;
            }

            auto negative = assumptions_->is_negative(arg_expr);
            if (negative &&
                negative.value() == LMCAS::Tribool::True) {
                // abs(x) → -x when x is Negative
                std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {
                    LMCAS::detail::make_node<NumberNode>(BigInt(-1)), arg
                };
                return make_normalized_multiply_node(mul_ops);
            }
        }

        return nullptr;
    }

} // namespace LMCAS
