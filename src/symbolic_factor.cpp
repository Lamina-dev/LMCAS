#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include "../include/symbolic.hpp"
#include "../include/symbolic_ast.hpp"
#include "../include/poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "../include/transcendental_factor.hpp"
#include "../include/solve_polynomial.hpp"
#include "../include/solve_strategies.hpp"
#include "../include/multivariate_factor.hpp"

namespace {


/// 多元因式分解辅助函数

/**
 * @internal
 * @brief 判断符号节点是否为多项式表达式(不含超越函数,负指数等)
 */
static bool is_poly_expr_node(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return false;
    if (std::dynamic_pointer_cast<const NumberNode>(node)) return true;
    if (std::dynamic_pointer_cast<const VariableNode>(node)) return true;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (!is_poly_expr_node(op)) return false;
        }
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (!is_poly_expr_node(op)) return false;
        }
        return true;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
            long long exp_val = 0;
            if (std::holds_alternative<BigInt>(exp_num->value())) {
                auto bi = std::get<BigInt>(exp_num->value());
                if (bi.IsNegative()) return false;
                exp_val = bi.to_int();
            } else if (std::holds_alternative<Rational>(exp_num->value())) {
                auto r = std::get<Rational>(exp_num->value());
                if (!r.is_integer()) return false;
                auto bi = r.to_BigInt();
                if (bi.IsNegative()) return false;
                exp_val = bi.to_int();
            } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
                lmmc_real_t d = std::get<lmmc_real_t>(exp_num->value());
                if (!std::isfinite(d) || d < 0 || d != std::floor(d)) return false;
                exp_val = static_cast<long long>(d);
            } else {
                return false;
            }
            if (exp_val < 0 || exp_val > 100) return false;
            return is_poly_expr_node(pow->base());
        }
        return false;
    }

    return false;
}

} // namespace

/**
 * @internal
 * @brief 递归将符号节点转换为 MultiPoly
 *
 * 假设节点已通过 is_poly_expr_node 验证为多项式.
 * @param[in] node 符号节点
 * @param[in] vars 变量名列表(确定单项式各分量的含义)
 * @return 对应的 MultiPoly
 */
static lamina::MultiPoly symbolic_node_to_multipoly(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::vector<std::string>& vars)
{
    if (!node) return lamina::MultiPoly(Rational(0), vars);

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        Rational coeff(0);
        if (std::holds_alternative<BigInt>(num->value())) {
            coeff = Rational(std::get<BigInt>(num->value()));
        } else if (std::holds_alternative<Rational>(num->value())) {
            coeff = std::get<Rational>(num->value());
        } else if (std::holds_alternative<lmmc_real_t>(num->value())) {
            coeff = Rational::from_double(std::get<lmmc_real_t>(num->value()));
        }
        return lamina::MultiPoly(coeff, vars);
    }

    if (auto var_node = std::dynamic_pointer_cast<const VariableNode>(node)) {
        lamina::Monomial mono(vars.size(), 0);
        for (size_t i = 0; i < vars.size(); ++i) {
            if (vars[i] == var_node->name()) {
                mono[i] = 1;
                break;
            }
        }
        std::vector<lamina::MultiPoly::Term> terms;
        terms.push_back({mono, Rational(1)});
        return lamina::MultiPoly(std::move(terms), vars);
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        lamina::MultiPoly result(Rational(0), vars);
        for (const auto& op : add->operands()) {
            result = result + symbolic_node_to_multipoly(op, vars);
        }
        return result;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        lamina::MultiPoly result(Rational(1), vars);
        for (const auto& op : mul->operands()) {
            result = result * symbolic_node_to_multipoly(op, vars);
        }
        return result;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base_poly = symbolic_node_to_multipoly(pow->base(), vars);
        int exp_val = 0;
        if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
            if (std::holds_alternative<BigInt>(exp_num->value())) {
                exp_val = std::get<BigInt>(exp_num->value()).to_int();
            } else if (std::holds_alternative<Rational>(exp_num->value())) {
                exp_val = std::get<Rational>(exp_num->value()).to_BigInt().to_int();
            } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
                exp_val = static_cast<int>(std::get<lmmc_real_t>(exp_num->value()));
            }
        }
        if (exp_val == 0) return lamina::MultiPoly(Rational(1), vars);
        lamina::MultiPoly result(Rational(1), vars);
        for (int i = 0; i < exp_val; ++i) {
            result = result * base_poly;
        }
        return result;
    }

    return lamina::MultiPoly(Rational(0), vars);
}

/**
 * @internal
 * @brief 将 MultiPoly 转换为符号表达式
 * @param[in] poly 多元多项式
 * @return 对应的符号表达式
 */
static std::shared_ptr<SymbolicExpr> multipoly_to_symbolic(const lamina::MultiPoly& poly) {
    if (poly.is_zero()) return SymbolicExpr::number(0);

    const auto& vars = poly.variables();
    const auto& terms = poly.terms();

    std::vector<std::shared_ptr<SymbolicExpr>> term_exprs;
    term_exprs.reserve(terms.size());

    for (const auto& [mono, coeff] : terms) {
        std::vector<std::shared_ptr<SymbolicExpr>> factors;

        /// 系数部分
        if (!(coeff == Rational(1)) || lamina::total_degree(mono) == 0) {
            if (coeff == Rational(-1) && lamina::total_degree(mono) > 0) {
                factors.push_back(SymbolicExpr::number(-1));
            } else {
                factors.push_back(SymbolicExpr::number(coeff));
            }
        }

        /// 变量部分
        for (size_t i = 0; i < vars.size() && i < mono.size(); ++i) {
            if (mono[i] == 0) continue;
            auto var_expr = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<VariableNode>(vars[i]));
            if (mono[i] == 1) {
                factors.push_back(var_expr);
            } else {
                factors.push_back(SymbolicExpr::power(var_expr,
                    SymbolicExpr::number(mono[i])));
            }
        }

        std::shared_ptr<SymbolicExpr> term_expr;
        if (factors.empty()) {
            term_expr = SymbolicExpr::number(1);
        } else if (factors.size() == 1) {
            term_expr = factors[0];
        } else {
            auto result = factors[0];
            for (size_t i = 1; i < factors.size(); ++i) {
                result = SymbolicExpr::multiply(result, factors[i]);
            }
            term_expr = result;
        }
        term_exprs.push_back(term_expr);
    }

    if (term_exprs.empty()) return SymbolicExpr::number(0);
    if (term_exprs.size() == 1) return term_exprs[0]->simplify();

    auto result = term_exprs[0];
    for (size_t i = 1; i < term_exprs.size(); ++i) {
        result = SymbolicExpr::add(result, term_exprs[i]);
    }
    return result->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::factor_impl(
    lamina::ComputationContext& context) const {
    auto simp = simplify();
    if (!simp || !lamina::detail::node(simp)) return simp;

    if (auto add_node = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(simp))) {

        const Rational content_value = lamina::detail::propagate_result(
            lamina::symbolic_polynomial_content(*simp, context));
        if (content_value != Rational(0) && content_value != Rational(1)) {
            auto content = content_value.is_integer()
                ? number(content_value.to_BigInt())
                : number(content_value);
            const Rational inverse_value = Rational(1) / content_value;
            auto inverse_content = inverse_value.is_integer()
                ? number(inverse_value.to_BigInt())
                : number(inverse_value);
            std::vector<std::shared_ptr<const SymbolicNode>> primitive_terms;
            primitive_terms.reserve(add_node->operands().size());
            for (const auto& operand : add_node->operands()) {
                auto term = lamina::detail::make_expression_ptr(operand);
                auto primitive_term = multiply(term, inverse_content)->simplify();
                primitive_terms.push_back(lamina::detail::node(primitive_term));
            }
            auto primitive_sum = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<AddNode>(primitive_terms));
            auto factored_primitive = lamina::detail::propagate_result(
                primitive_sum->factor_checked(context));
            return multiply(content, factored_primitive)->simplify();
        }

        /// 步骤 1:提取各项的公因式(GCD)
        std::shared_ptr<SymbolicExpr> common = nullptr;
        for (const auto& op : add_node->operands()) {
             auto expr_op = lamina::detail::make_expression_ptr(op);
             if (!common) common = expr_op;
             else {
                 common = lamina::detail::propagate_result(
                     lamina::symbolic_polynomial_gcd(
                         *common, *expr_op, context));
             }
        }

        if (common && !common->is_one() && !common->is_zero()) {
             std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
             for (const auto& op : add_node->operands()) {
                  auto term = lamina::detail::make_expression_ptr(op);

                  auto inv_common = power(common, number(-1));
                  auto quot = multiply(term, inv_common);
                  quot = quot->simplify();
                  new_ops.push_back(lamina::detail::node(quot));
             }
             auto new_sum = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(new_ops));

             /// 递归分解余下的和式
             auto factored_sum = lamina::detail::propagate_result(
                 new_sum->factor_checked(context));

             std::vector<std::shared_ptr<const SymbolicNode>> final_ops = {lamina::detail::node(common), lamina::detail::node(factored_sum)};
             return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(final_ops));
        }

        /// 步骤 2:一元多项式分解(支持任意次数)
        auto factor_variables = lamina::free_variables(lamina::detail::node(simp));
        if (factor_variables.size() == 1) {
             std::string var = *factor_variables.begin();
             try {
                 auto poly = lamina::symbolic_to_poly<Rational>(simp, var);
                 int deg = poly.degree();

                 if (deg >= 2) {
                      /// 使用有理根定理逐步分解
                      auto roots = lamina::find_rational_roots(poly);

                      if (!roots.empty()) {
                           auto x_node = lamina::detail::make_expression_ptr(lamina::detail::make_node<VariableNode>(var));
                           auto leading = poly.lead_coeff();

                           std::vector<std::shared_ptr<const SymbolicNode>> factors;

                           /// 首项系数
                           if (!(leading == Rational(1))) {
                               factors.push_back(lamina::detail::node(number(leading)));
                           }

                           /// 从根构建线性因子 (x - r)
                           for (const auto& r : roots) {
                               std::shared_ptr<SymbolicExpr> linear_factor;
                               if (r == Rational(0)) {
                                   linear_factor = x_node;
                               } else {
                                   auto neg_r = number(Rational(0) - r);
                                   linear_factor = SymbolicExpr::add(x_node, neg_r)->simplify();
                               }
                               factors.push_back(lamina::detail::node(linear_factor));
                           }

                           /// 计算余下的不可约因子:原多项式 / 已分解因子的乘积
                           lamina::Polynomial<Rational> factored_product({leading}, var);
                           for (const auto& r : roots) {
                               lamina::Polynomial<Rational> lin({Rational(0) - r, Rational(1)}, var);
                               factored_product = factored_product * lin;
                           }

                           auto [quotient, remainder] = poly.div_mod(factored_product);

                           if (remainder.is_zero() && quotient.degree() >= 1) {
                               /// 余下的不可约因子
                               if (quotient.degree() == 1 && quotient.lead_coeff() == Rational(1)) {
                                   /// 一次因子直接加入
                                   auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                   factors.push_back(lamina::detail::node(q_expr));
                               } else if (quotient.degree() >= 2) {
                                   /// 尝试递归分解余下部分
                                   auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                   auto q_factored = lamina::detail::propagate_result(
                                       q_expr->factor_checked(context));
                                   factors.push_back(lamina::detail::node(q_factored));
                               } else {
                                   auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                   factors.push_back(lamina::detail::node(q_expr));
                               }
                           } else if (!remainder.is_zero()) {
                               /// 除法有余数 -> 回退到原始方法
                               goto try_solve_quadratic;
                           }

                           if (factors.size() > 1) {
                               return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(factors));
                           } else if (factors.size() == 1) {
                               return lamina::detail::make_expression_ptr(factors[0]);
                           }
                      }
                 }
             } catch (const std::invalid_argument&) {
             } catch (const std::out_of_range&) {
             }

             try_solve_quadratic:
             /// 后备:二次多项式通过求解方程分解
             try {
                 auto poly = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(
                     simp, *factor_variables.begin());
                 if (poly.degree() == 2) {
                      std::string var = *factor_variables.begin();
                      auto solutions =
                          lamina::detail::propagate_result(
                              lamina::solve_finite_checked(simp, var, context));
                      if (solutions.size() == 2) {
                           auto leading = poly.coeffs[2].val;
                           if (!leading) leading = number(1);
                           leading = leading->simplify();

                           auto x_node = lamina::detail::make_expression_ptr(lamina::detail::make_node<VariableNode>(var));

                           auto t1 = SymbolicExpr::add(x_node, multiply(solutions[0], number(-1)))->simplify();
                           auto t2 = SymbolicExpr::add(x_node, multiply(solutions[1], number(-1)))->simplify();

                           std::vector<std::shared_ptr<const SymbolicNode>> factors;
                           if (!leading->is_one()) factors.push_back(lamina::detail::node(leading));
                           factors.push_back(lamina::detail::node(t1));
                           factors.push_back(lamina::detail::node(t2));

                           return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(factors));
                      }
                 }
             } catch (const std::invalid_argument&) {
             } catch (const std::out_of_range&) {
             }
        }

        /// 步骤 3:多元多项式 - 先尝试 factor_multivariate,再回退逐变量分解
        if (factor_variables.size() > 1) {
            /// 3a: 检测是否为多项式表达式,若是则使用 MultiPoly 路径
            if (is_poly_expr_node(lamina::detail::node(simp))) {
                try {
                    std::vector<std::string> var_list(
                        factor_variables.begin(), factor_variables.end());
                    auto mpoly = symbolic_node_to_multipoly(lamina::detail::node(simp), var_list);

                    if (!mpoly.is_zero() && !mpoly.is_constant()) {
                        const auto& checked =
                            lamina::detail::propagate_result(
                                lamina::factor_multivariate_checked(
                                    mpoly, context));
                        if (checked.completeness !=
                            lamina::Completeness::Complete) {
                            return simp;
                        }
                        const auto& result = checked.value;

                        /// 检查分解是否产生了多个因子
                        if (!result.factors.empty()) {
                            std::vector<std::shared_ptr<const SymbolicNode>> factor_nodes;

                            /// 常数因子
                            if (!(result.constant == Rational(1))) {
                                factor_nodes.push_back(
                                    lamina::detail::node(SymbolicExpr::number(result.constant)));
                            }

                            /// 各不可约因子
                            for (size_t i = 0; i < result.factors.size(); ++i) {
                                auto factor_expr = multipoly_to_symbolic(result.factors[i]);
                                if (!factor_expr || factor_expr->is_zero()) continue;

                                int mult = (i < result.multiplicities.size())
                                    ? result.multiplicities[i] : 1;

                                if (mult == 1) {
                                    factor_nodes.push_back(lamina::detail::node(factor_expr));
                                } else {
                                    auto pow_expr = SymbolicExpr::power(
                                        factor_expr, SymbolicExpr::number(mult));
                                    factor_nodes.push_back(lamina::detail::node(pow_expr));
                                }
                            }

                            if (factor_nodes.size() > 1) {
                                return lamina::detail::make_expression_ptr(
                                    lamina::detail::make_node<MultiplyNode>(std::move(factor_nodes)));
                            } else if (factor_nodes.size() == 1) {
                                return lamina::detail::make_expression_ptr(factor_nodes[0]);
                            }
                        }
                    }
                } catch (const std::invalid_argument&) {
                    /// MultiPoly conversion failed; use the variable fallback.
                } catch (const std::out_of_range&) {
                    /// MultiPoly conversion failed; use the variable fallback.
                }
            }

            /// 3b: 回退 - 逐变量尝试有理根分解
            for (const auto& var : factor_variables) {
                try {
                    auto poly = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(simp, var);
                    if (poly.degree() >= 2) {
                        /// 尝试用 Rational 系数做分解
                        auto poly_r = lamina::symbolic_to_poly<Rational>(simp, var);
                        auto roots = lamina::find_rational_roots(poly_r);
                        if (!roots.empty()) {
                            auto x_node = lamina::detail::make_expression_ptr(lamina::detail::make_node<VariableNode>(var));
                            std::vector<std::shared_ptr<const SymbolicNode>> factors;

                            auto leading = poly_r.lead_coeff();
                            if (!(leading == Rational(1))) {
                                factors.push_back(lamina::detail::node(number(leading)));
                            }

                            for (const auto& r : roots) {
                                std::shared_ptr<SymbolicExpr> linear_factor;
                                if (r == Rational(0)) {
                                    linear_factor = x_node;
                                } else {
                                    auto neg_r = number(Rational(0) - r);
                                    linear_factor = SymbolicExpr::add(x_node, neg_r)->simplify();
                                }
                                factors.push_back(lamina::detail::node(linear_factor));
                            }

                            /// 计算余下因子
                            lamina::Polynomial<Rational> factored_product({leading}, var);
                            for (const auto& r : roots) {
                                lamina::Polynomial<Rational> lin({Rational(0) - r, Rational(1)}, var);
                                factored_product = factored_product * lin;
                            }
                            auto [quotient, remainder] = poly_r.div_mod(factored_product);
                            if (remainder.is_zero() && quotient.degree() >= 1) {
                                auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                factors.push_back(lamina::detail::node(q_expr));
                            } else if (remainder.is_zero() && !quotient.is_zero() && quotient.degree() == 0) {
                                /// 常数商 -> 已完全分解
                                if (!(quotient.coeffs[0] == Rational(1))) {
                                    factors.push_back(
                                        lamina::detail::node(number(quotient.coeffs[0])));
                                }
                            }

                            if (factors.size() > 1) {
                                return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(factors));
                            }
                        }
                    }
                } catch (const std::invalid_argument&) {
                    continue;
                } catch (const std::out_of_range&) {
                    continue;
                }
            }
        }
    }

    /// 标准多项式分解未匹配时,后备路径检测超越函数并尝试超越因式分解.
    {
        const auto transform_variables =
            lamina::free_variables(lamina::detail::node(simp));
        std::string target_var;
        if (!transform_variables.empty()) {
            target_var = *transform_variables.begin();
        }

        if (!target_var.empty()) {
            try {
                auto trans_factors = lamina::factor_transcendental(simp, target_var);
                if (trans_factors.size() > 1) {
                    /// 超越分解成功:组装乘积表达式
                    std::vector<std::shared_ptr<const SymbolicNode>> factor_nodes;
                    factor_nodes.reserve(trans_factors.size());
                    for (const auto& f : trans_factors) {
                        if (f && lamina::detail::node(f)) {
                            factor_nodes.push_back(lamina::detail::node(f));
                        }
                    }
                    if (factor_nodes.size() > 1) {
                        return lamina::detail::make_expression_ptr(
                            lamina::detail::make_node<MultiplyNode>(std::move(factor_nodes)));
                    }
                }
            } catch (const std::invalid_argument&) {
                /// Transcendental conversion did not apply.
            } catch (const std::out_of_range&) {
                /// Transcendental conversion did not apply.
            }
        }
    }

    return simp;
}

lamina::ExpressionResult SymbolicExpr::factor_checked(
    lamina::ComputationContext& context) const
{
    constexpr const char* operation = "factor";
    auto budget = context.consume_steps(1, operation);
    if (!budget) return lamina::ExpressionResult::failure(budget.error());
    try {
        return lamina::ExpressionResult::success(factor_impl(context));
    } catch (const lamina::detail::ResultPropagation& propagation) {
        return lamina::ExpressionResult::failure(propagation.error());
    } catch (const std::bad_alloc&) {
        return lamina::ExpressionResult::failure(
            lamina::CasErrc::ResourceLimit,
            "allocation failed while factoring expression", operation);
    } catch (const std::exception& ex) {
        return lamina::ExpressionResult::failure(
            lamina::CasErrc::InternalInvariant, ex.what(), operation);
    }
}

lamina::ExpressionResult SymbolicExpr::factor_checked() const
{
    lamina::ComputationContext context;
    return factor_checked(context);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::cancel() const {
    auto simp = simplify();
    if (!simp || !lamina::detail::node(simp)) return simp;

    /// 辅助 lambda:从乘积节点中分离分子因子和分母因子.
    /// 分母因子 = 含负指数的 PowerNode.
    auto separate_num_den = [](const std::shared_ptr<const SymbolicNode>& node,
                               std::vector<std::shared_ptr<const SymbolicNode>>& num_out,
                               std::vector<std::shared_ptr<const SymbolicNode>>& den_out) {
        auto classify = [&](const std::shared_ptr<const SymbolicNode>& factor) {
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(factor)) {
                if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                    bool is_negative = false;
                    if (std::holds_alternative<BigInt>(exp_num->value())) {
                        is_negative = std::get<BigInt>(exp_num->value()).IsNegative();
                    } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
                        is_negative = std::get<lmmc_real_t>(exp_num->value()) < 0;
                    } else if (std::holds_alternative<Rational>(exp_num->value())) {
                        is_negative = std::get<Rational>(exp_num->value()).get_numerator().IsNegative();
                    }
                    if (is_negative) {
                        /// 取绝对值指数
                        std::shared_ptr<const SymbolicNode> pos_exp;
                        if (std::holds_alternative<BigInt>(exp_num->value())) {
                            pos_exp = lamina::detail::make_node<NumberNode>(BigInt(0) - std::get<BigInt>(exp_num->value()));
                        } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
                            pos_exp = lamina::detail::make_node<NumberNode>(-std::get<lmmc_real_t>(exp_num->value()));
                        } else {
                            auto r = std::get<Rational>(exp_num->value());
                            pos_exp = lamina::detail::make_node<NumberNode>(Rational(BigInt(0) - r.get_numerator(), r.get_denominator()));
                        }
                        auto pos_exp_num = std::dynamic_pointer_cast<const NumberNode>(pos_exp);
                        if (pos_exp_num && pos_exp_num->is_one()) {
                            den_out.push_back(pow->base());
                        } else {
                            den_out.push_back(lamina::detail::make_node<PowerNode>(pow->base(), pos_exp));
                        }
                        return;
                    }
                }
            }
            num_out.push_back(factor);
        };

        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
            for (const auto& op : mul->operands()) {
                classify(op);
            }
        } else {
            classify(node);
        }
    };

    /// 辅助 lambda:从因子列表构建乘积表达式
    auto build_product = [](const std::vector<std::shared_ptr<const SymbolicNode>>& factors) -> std::shared_ptr<SymbolicExpr> {
        if (factors.empty()) return SymbolicExpr::number(1);
        if (factors.size() == 1) return lamina::detail::make_expression_ptr(factors[0]);
        return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(factors));
    };

    /// 策略 1:表达式是 AddNode,各项可能含公共分母因子.
    /// 例如 simplify 后 (x^2-1)/(x-1) 变为 -1*(x-1)^-1 + x^2*(x-1)^-1
    /// 需要提取公共分母,重组为 (分子之和)/分母 再做 GCD 约分.
    if (auto add_node = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(simp))) {
        /// 对每个加法项分离分子/分母
        struct TermInfo {
            std::shared_ptr<SymbolicExpr> numerator;
            std::shared_ptr<SymbolicExpr> denominator;
        };
        std::vector<TermInfo> terms;
        bool has_denominator = false;

        for (const auto& term : add_node->operands()) {
            std::vector<std::shared_ptr<const SymbolicNode>> t_num, t_den;
            separate_num_den(term, t_num, t_den);
            auto num_expr = build_product(t_num)->simplify();
            auto den_expr = build_product(t_den)->simplify();
            if (!den_expr->is_one()) has_denominator = true;
            terms.push_back({num_expr, den_expr});
        }

        if (has_denominator) {
            /// 计算公共分母(所有项分母的 LCM,简化处理:乘积)
            /// 对于常见情况(所有项分母相同),直接提取.
            /// 检查是否所有分母相同
            bool all_same_den = true;
            auto first_den_str = terms[0].denominator->to_string();
            for (size_t i = 1; i < terms.size(); ++i) {
                if (terms[i].denominator->to_string() != first_den_str) {
                    all_same_den = false;
                    break;
                }
            }

            std::shared_ptr<SymbolicExpr> combined_num;
            std::shared_ptr<SymbolicExpr> combined_den;

            if (all_same_den) {
                /// 所有项分母相同:分子直接相加
                combined_den = terms[0].denominator;
                std::vector<std::shared_ptr<const SymbolicNode>> num_ops;
                for (const auto& t : terms) {
                    if (t.numerator && lamina::detail::node(t.numerator)) {
                        num_ops.push_back(lamina::detail::node(t.numerator));
                    }
                }
                if (num_ops.empty()) return SymbolicExpr::number(0);
                if (num_ops.size() == 1) combined_num = lamina::detail::make_expression_ptr(num_ops[0]);
                else combined_num = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(num_ops));
                combined_num = combined_num->simplify();
            } else {
                /// 分母不同:通分(乘以其他项的分母)
                /// 计算总分母 = 所有不同分母的乘积
                std::shared_ptr<SymbolicExpr> total_den = SymbolicExpr::number(1);
                /// 收集不同的分母
                std::vector<std::shared_ptr<SymbolicExpr>> unique_dens;
                for (const auto& t : terms) {
                    if (!t.denominator->is_one()) {
                        bool found = false;
                        for (const auto& ud : unique_dens) {
                            if (ud->to_string() == t.denominator->to_string()) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            unique_dens.push_back(t.denominator);
                        }
                    }
                }
                for (const auto& ud : unique_dens) {
                    total_den = SymbolicExpr::multiply(total_den, ud)->simplify();
                }
                combined_den = total_den;

                /// 每项分子乘以 (总分母 / 该项分母)
                std::vector<std::shared_ptr<const SymbolicNode>> num_ops;
                for (const auto& t : terms) {
                    auto factor = divide(total_den, t.denominator)->simplify();
                    auto adjusted_num = SymbolicExpr::multiply(t.numerator, factor)->simplify();
                    if (adjusted_num && lamina::detail::node(adjusted_num)) {
                        num_ops.push_back(lamina::detail::node(adjusted_num));
                    }
                }
                if (num_ops.empty()) return SymbolicExpr::number(0);
                if (num_ops.size() == 1) combined_num = lamina::detail::make_expression_ptr(num_ops[0]);
                else combined_num = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(num_ops));
                combined_num = combined_num->expand()->simplify();
            }

            /// 对 combined_num / combined_den 做多项式 GCD 约分
            combined_den = combined_den->expand()->simplify();

            auto rational_variables =
                lamina::free_variables(lamina::detail::node(combined_num));
            const auto denominator_variables =
                lamina::free_variables(lamina::detail::node(combined_den));
            rational_variables.insert(denominator_variables.begin(),
                                      denominator_variables.end());

            auto cur_num = combined_num;
            auto cur_den = combined_den;

            for (const auto& var : rational_variables) {
                try {
                    auto num_expanded = cur_num->expand()->simplify();
                    auto den_expanded = cur_den->expand()->simplify();

                    /// 使用 SymbolicPolyCoeff 保留其他变量作为符号系数
                    auto poly_num = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(num_expanded, var);
                    auto poly_den = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(den_expanded, var);

                    if (poly_num.is_zero()) return SymbolicExpr::number(0);
                    if (poly_den.is_zero()) return simp;

                    /// 尝试用 Rational 系数做 GCD(纯数值系数时更可靠)
                    auto poly_num_r = lamina::symbolic_to_poly<Rational>(num_expanded, var);
                    auto poly_den_r = lamina::symbolic_to_poly<Rational>(den_expanded, var);

                    /// 检查 Rational 转换是否丢失了信息(多元情况)
                    bool rational_ok = true;
                    auto reconstructed_num = lamina::poly_to_symbolic(poly_num_r)->expand()->simplify();
                    auto reconstructed_den = lamina::poly_to_symbolic(poly_den_r)->expand()->simplify();
                    auto diff_num = SymbolicExpr::add(num_expanded, SymbolicExpr::multiply(SymbolicExpr::number(-1), reconstructed_num))->simplify();
                    auto diff_den = SymbolicExpr::add(den_expanded, SymbolicExpr::multiply(SymbolicExpr::number(-1), reconstructed_den))->simplify();
                    if (!diff_num->is_zero() || !diff_den->is_zero()) {
                        rational_ok = false;
                    }

                    if (rational_ok) {
                        auto g = lamina::Polynomial<Rational>::gcd(poly_num_r, poly_den_r);
                        if (g.degree() >= 1) {
                            auto [q_num, r_num] = poly_num_r.div_mod(g);
                            auto [q_den, r_den] = poly_den_r.div_mod(g);
                            cur_num = lamina::poly_to_symbolic(q_num)->simplify();
                            cur_den = lamina::poly_to_symbolic(q_den)->simplify();
                        }
                    } else {
                        /// 多元情况:使用 SymbolicPolyCoeff 做 GCD
                        auto g = lamina::Polynomial<lamina::SymbolicPolyCoeff>::gcd(poly_num, poly_den);
                        if (g.degree() >= 1) {
                            auto [q_num, r_num] = poly_num.div_mod(g);
                            auto [q_den, r_den] = poly_den.div_mod(g);
                            cur_num = lamina::poly_to_symbolic(q_num)->simplify();
                            cur_den = lamina::poly_to_symbolic(q_den)->simplify();
                        }
                    }
                } catch (const std::invalid_argument&) {
                    continue;
                } catch (const std::out_of_range&) {
                    continue;
                } catch (const std::runtime_error&) {
                    continue;
                }
            }

            if (cur_den->is_one()) return cur_num;

            /// 分母为 -1 时取负
            if (lamina::detail::node(cur_den)) {
                auto diff = SymbolicExpr::add(lamina::detail::make_expression_ptr(lamina::detail::node(cur_den)), SymbolicExpr::number(-1))->simplify();
                if (diff->is_zero()) {
                    return SymbolicExpr::multiply(SymbolicExpr::number(-1), cur_num)->simplify();
                }
            }

            return divide(cur_num, cur_den)->simplify();
        }
    }

    /// 策略 2:表达式是 MultiplyNode,直接分离分子/分母.
    std::vector<std::shared_ptr<const SymbolicNode>> num_factors;
    std::vector<std::shared_ptr<const SymbolicNode>> den_factors;
    separate_num_den(lamina::detail::node(simp), num_factors, den_factors);

    if (den_factors.empty()) return simp;

    auto numerator = build_product(num_factors)->simplify();
    auto denominator = build_product(den_factors)->simplify();

    if (denominator->is_one()) return numerator;

    auto rational_variables = lamina::free_variables(lamina::detail::node(numerator));
    const auto denominator_variables =
        lamina::free_variables(lamina::detail::node(denominator));
    rational_variables.insert(denominator_variables.begin(),
                              denominator_variables.end());

    if (rational_variables.empty()) {
        return divide(numerator, denominator)->simplify();
    }

    auto cur_num = numerator;
    auto cur_den = denominator;

    for (const auto& var : rational_variables) {
        try {
            auto num_expanded = cur_num->expand()->simplify();
            auto den_expanded = cur_den->expand()->simplify();

            auto poly_num_r = lamina::symbolic_to_poly<Rational>(num_expanded, var);
            auto poly_den_r = lamina::symbolic_to_poly<Rational>(den_expanded, var);

            /// 检查 Rational 转换是否丢失了信息
            bool rational_ok = true;
            auto reconstructed_num = lamina::poly_to_symbolic(poly_num_r)->expand()->simplify();
            auto reconstructed_den = lamina::poly_to_symbolic(poly_den_r)->expand()->simplify();
            auto diff_num = SymbolicExpr::add(num_expanded, SymbolicExpr::multiply(SymbolicExpr::number(-1), reconstructed_num))->simplify();
            auto diff_den = SymbolicExpr::add(den_expanded, SymbolicExpr::multiply(SymbolicExpr::number(-1), reconstructed_den))->simplify();
            if (!diff_num->is_zero() || !diff_den->is_zero()) {
                rational_ok = false;
            }

            if (rational_ok) {
                if (poly_num_r.is_zero()) return SymbolicExpr::number(0);
                if (poly_den_r.is_zero()) return simp;

                auto g = lamina::Polynomial<Rational>::gcd(poly_num_r, poly_den_r);
                if (g.degree() >= 1) {
                    auto [q_num, r_num] = poly_num_r.div_mod(g);
                    auto [q_den, r_den] = poly_den_r.div_mod(g);
                    cur_num = lamina::poly_to_symbolic(q_num)->simplify();
                    cur_den = lamina::poly_to_symbolic(q_den)->simplify();
                }
            } else {
                auto poly_num = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(num_expanded, var);
                auto poly_den = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(den_expanded, var);

                if (poly_num.is_zero()) return SymbolicExpr::number(0);
                if (poly_den.is_zero()) return simp;

                auto g = lamina::Polynomial<lamina::SymbolicPolyCoeff>::gcd(poly_num, poly_den);
                if (g.degree() >= 1) {
                    auto [q_num, r_num] = poly_num.div_mod(g);
                    auto [q_den, r_den] = poly_den.div_mod(g);
                    cur_num = lamina::poly_to_symbolic(q_num)->simplify();
                    cur_den = lamina::poly_to_symbolic(q_den)->simplify();
                }
            }
        } catch (const std::invalid_argument&) {
            continue;
        } catch (const std::out_of_range&) {
            continue;
        } catch (const std::runtime_error&) {
            continue;
        }
    }

    if (cur_den->is_one()) return cur_num;

    if (lamina::detail::node(cur_den)) {
        auto diff = SymbolicExpr::add(lamina::detail::make_expression_ptr(lamina::detail::node(cur_den)), SymbolicExpr::number(-1))->simplify();
        if (diff->is_zero()) {
            return SymbolicExpr::multiply(SymbolicExpr::number(-1), cur_num)->simplify();
        }
    }

    return divide(cur_num, cur_den)->simplify();
}
