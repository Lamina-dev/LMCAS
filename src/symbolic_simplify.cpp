#include "symbolic.hpp"
#include "symbolic_internal.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

// Forward decl for internal use
static std::shared_ptr<SymbolicExpr> expand_power_integer(const std::shared_ptr<SymbolicExpr>& base, int n);

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_sqrt() const {
    if (operands.empty()) return std::make_shared<SymbolicExpr>(*this);
    
    auto simplified_operand = operands[0]->simplify();
	if (simplified_operand->type == SymbolicExpr::Type::Infinity) return simplified_operand;
    
    if (simplified_operand->is_number()) {
		
		auto scvrs = simplified_operand->convert_rational();
		
		if (simplified_operand->is_rational() && scvrs.get_denominator() == ::BigInt(1)) {
			::BigInt actual = scvrs.get_numerator();
			simplified_operand->number_value = actual;
		}
		
		// pair 格式：first 为系数，second 为根式下的值
		auto num_process = [](int n) -> std::pair<int, int> {
			if (n < 0) throw std::runtime_error("Square root of negative number");
			if (n == 0 || n == 1) return std::make_pair(n, 1);
			int sqrt_n = static_cast<int>(std::sqrt(n));
            if (sqrt_n * sqrt_n == n) return std::make_pair(sqrt_n, 1);
            int factor = 1, remaining = n;
            for (int i = 2; i * i <= remaining; ++i) {
                while (remaining % (i * i) == 0) {
                    factor *= i;
                    remaining /= (i * i);
                }
            }
            return std::make_pair(factor, remaining);
		};
		
		auto in_simplify_range = [](const ::BigInt& bi) -> bool {
			return bi <= BigInt(INT_MAX) && bi >= BigInt(INT_MIN);
		};
		
		auto generate_component = [](const ::Rational& rat) -> std::shared_ptr<SymbolicExpr> {
			if (rat.get_denominator() == ::BigInt(1)) return SymbolicExpr::number(rat.get_numerator());
			else return SymbolicExpr::number(rat);
		};
		
        auto num_val = simplified_operand->get_number();
        if (std::holds_alternative<int>(num_val)) {
            int n = std::get<int>(num_val);
			auto res = num_process(n);
			if (res.second == 1) return SymbolicExpr::number(res.first);
			else if (res.first == 1) return SymbolicExpr::sqrt(SymbolicExpr::number(res.second));
			else return SymbolicExpr::multiply(SymbolicExpr::number(res.first), SymbolicExpr::sqrt(SymbolicExpr::number(res.second)));
        }
        if (std::holds_alternative<::BigInt>(num_val)) {
            const auto& bi = std::get<::BigInt>(num_val);
            if (bi.negative) throw std::runtime_error("Square root of negative number");
            if (bi.is_zero() || bi == BigInt(1)) return SymbolicExpr::number(bi);
            if (bi.is_perfect_square()) return SymbolicExpr::number(bi.sqrt());

			// 尝试对大整数进行轻量级的平方因子提取：仅用一组小素数的平方因子来加速常见情形
			{
				BigInt factor(1), remaining(bi);
				const int small_primes[] = {2,3,5,7,11,13,17,19,23,29};
				for (int p : small_primes) {
					BigInt bp(p);
					BigInt psq = bp * bp;
					while ((remaining % psq).is_zero()) {
						factor = factor * bp;
						remaining = remaining / psq;
					}
				}
				if (factor > BigInt(1)) {
					if (remaining == BigInt(1)) return SymbolicExpr::number(factor);
					else return SymbolicExpr::multiply(SymbolicExpr::number(factor), SymbolicExpr::sqrt(SymbolicExpr::number(remaining)));
				}
			}
			// 暂时只判断可以转成int的
            if (in_simplify_range(bi)) {
                return SymbolicExpr::sqrt(SymbolicExpr::number(bi.to_int()))->simplify();
            }
        }
		// 分数化简中底数和指数分别判断
		if (std::holds_alternative<::Rational>(num_val)) {
			const auto &nobj = std::get<::Rational>(num_val);
			const auto &nume = nobj.get_numerator();
			const auto &deme = nobj.get_denominator();
			// 暂时只判断可以转成int的
			if (in_simplify_range(nume) && in_simplify_range(deme)) {
				auto numsimp = num_process(nume.to_int());
				auto demsimp = num_process(deme.to_int());
				::Rational numarea = ::Rational(numsimp.first, demsimp.first);
				::Rational sqarea = ::Rational(numsimp.second, demsimp.second);
				
				if (sqarea == ::Rational(1)) return SymbolicExpr::number(numarea);
				else if (numarea == ::Rational(1)) return SymbolicExpr::sqrt(SymbolicExpr::number(sqarea));
				return SymbolicExpr::multiply(generate_component(numarea), SymbolicExpr::sqrt(generate_component(sqarea)));
			}
		}
    }
    // sqrt(x*x) 或 sqrt(π*π) 直接返回 x 或 π
    if (simplified_operand->type == SymbolicExpr::Type::Multiply && simplified_operand->operands.size() == 2) {
        const auto& a = simplified_operand->operands[0];
        const auto& b = simplified_operand->operands[1];
        // sqrt(x*x) = x
        if (a->type == SymbolicExpr::Type::Variable && b->type == SymbolicExpr::Type::Variable && a->identifier == b->identifier) {
            return a;
        }
        // sqrt(π*π) = π
        if (a->type == SymbolicExpr::Type::Variable && b->type == SymbolicExpr::Type::Variable &&
            ((a->identifier == "π" && b->identifier == "π") || (a->identifier == "pi" && b->identifier == "pi"))) {
            return a;
        }
        if (a->to_string() == b->to_string()) {
            return a;
        }
        auto get_var = [](const std::shared_ptr<SymbolicExpr>& expr) -> std::string {
            if (expr->type == SymbolicExpr::Type::Variable) return expr->identifier;
            if (expr->type == SymbolicExpr::Type::Multiply && expr->operands.size() == 2) {
                if (expr->operands[1]->type == SymbolicExpr::Type::Variable) return expr->operands[1]->identifier;
            }
            return "";
        };
        std::string var_a = get_var(a);
        std::string var_b = get_var(b);
        if (!var_a.empty() && var_a == var_b) {
            // sqrt((c*π)*(d*π)) = sqrt((c*d)*π^2) = π*sqrt(c*d)
            auto pow2 = SymbolicExpr::power(SymbolicExpr::variable(var_a), SymbolicExpr::number(2));
            auto left_coeff = (a->type == SymbolicExpr::Type::Multiply && a->operands[0]->is_number()) ? a->operands[0] : SymbolicExpr::number(1);
            auto right_coeff = (b->type == SymbolicExpr::Type::Multiply && b->operands[0]->is_number()) ? b->operands[0] : SymbolicExpr::number(1);
            auto coeff_mul = SymbolicExpr::multiply(left_coeff, right_coeff)->simplify();
            // sqrt((c*π)*(d*π)) = π*sqrt(c*d)
            if (coeff_mul->is_number()) {
                auto sqrt_coeff = SymbolicExpr::sqrt(coeff_mul)->simplify();
                if (sqrt_coeff->is_number() && sqrt_coeff->to_string() == "1") {
                    return SymbolicExpr::variable(var_a);
                } else {
                    return SymbolicExpr::multiply(sqrt_coeff, SymbolicExpr::variable(var_a));
                }
            } else {
                // fallback: sqrt(π^2)
                return SymbolicExpr::sqrt(pow2)->simplify();
            }
        }
    }
    // sqrt(x^2) 或 sqrt(π^2) 直接返回 x 或 π
    if (simplified_operand->type == SymbolicExpr::Type::Power && simplified_operand->operands.size() == 2) {
        const auto& base = simplified_operand->operands[0];
        const auto& exp = simplified_operand->operands[1];
        if (exp->is_number()) {
            auto exp_val = exp->get_number();
            int n = 0;
            BigInt big_n;// 0 for default
            if (std::holds_alternative<int>(exp_val)) n = std::get<int>(exp_val);
            else if (std::holds_alternative<::Rational>(exp_val)) {
                ::Rational r = std::get<::Rational>(exp_val);
                if (r.is_integer()) big_n = r.get_numerator();
            }
            if (n == 2 || big_n.to_string() == "2") {
                // sqrt(x^2) = x
                if (base->type == SymbolicExpr::Type::Variable && (base->identifier == "π" || base->identifier == "pi")) {
                    return base;
                }
                return base;
            }
        }
    }
    return SymbolicExpr::sqrt(simplified_operand);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_multiply() const {
    if (operands.size() != 2) return std::make_shared<SymbolicExpr>(*this);

    auto left = operands[0]->simplify();
    auto right = operands[1]->simplify();

	// 0 and 1 checks
	if (left->is_number() && left->convert_rational() == ::Rational(0)) return left;
	if (right->is_number() && right->convert_rational() == ::Rational(0)) return right;
	if (left->is_number() && left->convert_rational() == ::Rational(1)) return right;
	if (right->is_number() && right->convert_rational() == ::Rational(1)) return left;

    // Check for Matrix multiplication (Matrix*Matrix or Scalar*Matrix)
    if (left->type == SymbolicExpr::Type::Matrix || right->type == SymbolicExpr::Type::Matrix) {
        auto res = multiply_matrices(left, right);
        if (res) return res;
    }

	/*
	if (left->type == SymbolicExpr::Type::Add || right->type == SymbolicExpr::Type::Add) {
		return SymbolicExpr::multiply(left, right)->expand();
	}
	*/

	std::vector<std::shared_ptr<SymbolicExpr>> terms;
	::Rational coeff(1);

	std::function<void(const std::shared_ptr<SymbolicExpr>&)> collect = [&](const std::shared_ptr<SymbolicExpr>& e) {
		if (e->type == SymbolicExpr::Type::Multiply) {
			collect(e->operands[0]);
			collect(e->operands[1]);
		} else if (e->is_number()) {
			coeff = coeff * e->convert_rational();
		} else {
			terms.push_back(e);
		}
	};
	collect(left);
	collect(right);

	if (coeff == ::Rational(0)) return SymbolicExpr::number(0);

	std::sort(terms.begin(), terms.end(), [](const auto& a, const auto& b) {
		return a->compare(b) < 0; 
	});

	std::vector<std::shared_ptr<SymbolicExpr>> merged;
	if (!terms.empty()) {
		auto get_base_exp = [](const std::shared_ptr<SymbolicExpr>& t) -> std::pair<std::shared_ptr<SymbolicExpr>, ::Rational> {
			if (t->type == SymbolicExpr::Type::Power && t->operands[1]->is_number()) {
				return {t->operands[0], t->operands[1]->convert_rational()};
			}
			return {t, ::Rational(1)};
		};

		auto make_exp_node = [](const ::Rational& r) -> std::shared_ptr<SymbolicExpr> {
			if (r.is_integer()) {
				::BigInt num = r.get_numerator();
				int val = num.to_int();
				if (::BigInt(val) == num) return SymbolicExpr::number(val);
			}
			return SymbolicExpr::number(r);
		};

		auto current = get_base_exp(terms[0]);
		auto curr_base = current.first;
		auto curr_exp = current.second;

		for (size_t i = 1; i < terms.size(); ++i) {
			auto next = get_base_exp(terms[i]);
			auto next_base = next.first;
			auto next_exp = next.second;
			
			if (curr_base->compare(next_base) == 0) {
				curr_exp = curr_exp + next_exp;
			} else {
				if (curr_exp != ::Rational(0)) {
					if (curr_exp == ::Rational(1)) merged.push_back(curr_base);
					else merged.push_back(SymbolicExpr::power(curr_base, make_exp_node(curr_exp)));
				}
				curr_base = next_base;
				curr_exp = next_exp;
			}
		}
		if (curr_exp != ::Rational(0)) {
			if (curr_exp == ::Rational(1)) merged.push_back(curr_base);
			else merged.push_back(SymbolicExpr::power(curr_base, make_exp_node(curr_exp)));
		}
	}
	
	std::shared_ptr<SymbolicExpr> result;
    
    if (merged.empty()) {
        result = SymbolicExpr::number(coeff);
    } else {
        // Build right-associative chain: t1 * (t2 * (...))
        result = merged.back();
        for (int i = (int)merged.size() - 2; i >= 0; --i) {
            result = SymbolicExpr::multiply(merged[i], result);
        }
        
        // Multiply by coeff at the top
        if (coeff != ::Rational(1)) {
            result = SymbolicExpr::multiply(SymbolicExpr::number(coeff), result);
        }
    }
	
    auto return_res = result ? result : SymbolicExpr::number(1);
    return return_res;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_add() const {
    if (operands.size() != 2) return std::make_shared<SymbolicExpr>(*this);

    auto left_s = operands[0]->simplify();
    auto right_s = operands[1]->simplify();

    if (left_s->type == SymbolicExpr::Type::Matrix && right_s->type == SymbolicExpr::Type::Matrix) {
        auto res = add_matrices(left_s, right_s);
        if (res) return res;
    }

	if (left_s->type == SymbolicExpr::Type::Infinity) return left_s;
	if (right_s->type == SymbolicExpr::Type::Infinity) return right_s;

    // Helper to get coeff and term
    auto get_coeff_term = [](const std::shared_ptr<SymbolicExpr>& e) -> std::pair<::Rational, std::shared_ptr<SymbolicExpr>> {
        if (e->is_number()) {
            return {e->convert_rational(), SymbolicExpr::number(1)};
        }
        if (e->type == SymbolicExpr::Type::Multiply && e->operands.size() == 2) {
            if (e->operands[0]->is_number())
                return {e->operands[0]->convert_rational(), e->operands[1]};
            if (e->operands[1]->is_number())
                return {e->operands[1]->convert_rational(), e->operands[0]};
        }
        return {::Rational(1), e};
    };

    // Flatten
    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    std::function<void(const std::shared_ptr<SymbolicExpr>&)> flatten = [&](const std::shared_ptr<SymbolicExpr>& node) {
        if (node->type == SymbolicExpr::Type::Add) {
            flatten(node->operands[0]);
            flatten(node->operands[1]);
        } else {
            terms.push_back(node);
        }
    };
    flatten(left_s);
    flatten(right_s);

    // Sort
    std::sort(terms.begin(), terms.end(), [&](const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b){
        auto [ca, ta] = get_coeff_term(a);
        auto [cb, tb] = get_coeff_term(b);
        int cmp = ta->compare(tb);
        if (cmp != 0) return cmp < 0; 
        return ca < cb;
    });

    // Merge adjacent
    std::vector<std::shared_ptr<SymbolicExpr>> merged;
    if (!terms.empty()) {
		auto make_num_node = [](const ::Rational& r) -> std::shared_ptr<SymbolicExpr> {
			if (r.is_integer()) {
				::BigInt num = r.get_numerator();
				int val = num.to_int();
				if (::BigInt(val) == num) return SymbolicExpr::number(val);
			}
			return SymbolicExpr::number(r);
		};

        auto [curr_c, curr_t] = get_coeff_term(terms[0]);
        
        for (size_t i = 1; i < terms.size(); ++i) {
             auto [next_c, next_t] = get_coeff_term(terms[i]);
             if (symbolic_equal(curr_t, next_t)) {
                 curr_c = curr_c + next_c;
             } else {
                 // Push current
                 if (curr_t->is_number() && curr_t->convert_rational() == ::Rational(1)) {
                     if (curr_c != ::Rational(0)) {
                         merged.push_back(make_num_node(curr_c));
                     }
                 } else {
                     if (curr_c == ::Rational(0)) {
                         // Drop
                     } else if (curr_c == ::Rational(1)) {
                         merged.push_back(curr_t);
                     } else {
                         merged.push_back(SymbolicExpr::multiply(make_num_node(curr_c), curr_t));
                     }
                 }
                 // Reset
                 curr_c = next_c;
                 curr_t = next_t;
             }
        }
        // Push last
        if (curr_c != ::Rational(0)) {
            if (curr_t->is_number() && curr_t->convert_rational() == ::Rational(1)) {
                merged.push_back(make_num_node(curr_c));
            } else if (curr_c == ::Rational(1)) {
                merged.push_back(curr_t);
            } else {
                merged.push_back(SymbolicExpr::multiply(make_num_node(curr_c), curr_t));
            }
        }
    }
    
    if (merged.empty()) return SymbolicExpr::number(0);
    if (merged.size() == 1) return merged[0];
    
    // Rebuild Right-Associative
    // [a, b, c] -> a + (b + c)
    auto res = merged.back();
    for (int i = (int)merged.size() - 2; i >= 0; --i) {
        res = SymbolicExpr::add(merged[i], res); 
    }
    return res;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_power() const {
    // 如果底数是数字，且指数是整数，用数字储存
    auto base = operands[0]->simplify();
    auto exponent = operands[1]->simplify();
	
	// 处理指数为 0 的情况，但避免将 0^0 或 inf^0 错误化简
	if (exponent->type == SymbolicExpr::Type::Number && exponent->convert_rational() == ::Rational(0)) {
		if (base->is_number() && base->convert_rational() == ::Rational(0)) {
			// 0^0 未定，保留原表达式
			return SymbolicExpr::power(base, exponent);
		}
		if (base->type == SymbolicExpr::Type::Infinity) {
			// inf^0 未定，保留原表达式
			return SymbolicExpr::power(base, exponent);
		}
		return SymbolicExpr::number(1);
	}

	// 处理底数为 0 的情况，避免 0^负数 被错误化简
	if (base->type == SymbolicExpr::Type::Number && base->convert_rational() == ::Rational(0)) {
		if (exponent->is_number() && exponent->convert_rational() < ::Rational(0)) {
			// 0^(-n) 为未定/无穷，保留原式以避免错误
			return SymbolicExpr::power(base, exponent);
		}
		// 非负指数则为 0
		return SymbolicExpr::number(0);
	}
	if (exponent->is_number() && exponent->convert_rational() == ::Rational(1)) return base;
	if (base->is_number() && base->convert_rational() == ::Rational(1)) return base;
	
	if (base->type == SymbolicExpr::Type::Infinity) return base;
	if (exponent->type == SymbolicExpr::Type::Infinity) return exponent;
	
    if (base->is_number() && (exponent->is_int() || exponent->is_big_int())) {
		auto banum = base->convert_rational();
		auto exnum = exponent->convert_rational();
		
		if (exnum < ::Rational(0)) {
			return SymbolicExpr::power(SymbolicExpr::number(banum.reciprocal()), SymbolicExpr::number(::Rational(0) - exnum))->simplify();
		}
		
        auto expr = std::make_shared<SymbolicExpr>(Type::Number);
        // 底数是分数，结果为分数
        if (base->is_rational() || exponent->is_rational()) {
            if (exponent->is_int()) {
				expr->number_value = (base->get_rational()).power(BigInt(exponent->get_int()));
            } else if (exponent->is_big_int()) {
				expr->number_value = (base->get_rational()).power(exponent->get_big_int());
            }
        } else {// 否则结果为大整数
            BigInt b;
            if (base->is_int()) {
                b = BigInt(std::get<int>(base->get_number()));
            } else {
                b = std::get<BigInt>(base->get_number());
            }
            BigInt e;
            if (exponent->is_int()) {
                e = BigInt(std::get<int>(exponent->get_number()));
            } else {
                e = std::get<BigInt>(exponent->get_number());
            }
            if (e.to_int() >= 0) {
                expr->number_value = b.power(e);
            } else {
                
                expr->number_value = Rational(BigInt(1), b.power(e.negate()));
            }
        }

        return expr;
    } else if (base->is_number() && exponent->is_rational()) {
        // 底数是数字，指数是分数
		auto bsr = base->convert_rational();
		auto expr = exponent->convert_rational();
		
		auto in_range = [](const ::Rational& val) -> bool {
			auto vn = val.get_numerator(), vd = val.get_denominator();
			const ::BigInt lower = ::BigInt(INT_MIN), upper = ::BigInt(INT_MAX);
			return (vn >= lower && vn <= upper) && (vd >= lower && vd <= upper);
		};
		
		if (in_range(bsr) && in_range(expr)) {
			if (expr == ::Rational(1)) return SymbolicExpr::number(bsr);
			
			int bs_n = bsr.get_numerator().to_int(), bs_d = bsr.get_denominator().to_int();
			int es_n = expr.get_numerator().to_int(), es_d = expr.get_denominator().to_int();
			
			std::function<int(int,int)> __int_gcd;
			__int_gcd = [&__int_gcd](int a, int b) -> int {
				if (b == 0) return a;
				else return __int_gcd(b, a%b);
			};
			
			auto simplify_inner = [&__int_gcd](int& origin, const int& denom) -> int {
				if (denom == 1) return 1;
				int ediv = denom, target = origin;
				for (int i = 2; 1ll * i * i <= target; i++) {
					int exphere = 0;
					while (target % i == 0) {
						exphere++;
						target /= i;
					}
					if (exphere) {
						ediv = __int_gcd(ediv, exphere);
					}
				}
				if (ediv <= 1) return 0;
				int answer = 1;
				target = origin;
				for (int i = 2; 1ll * i * i <= target; i++) {
					int exphere = 0;
					while (target % i == 0) {
						exphere++;
						target /= i;
					}
					if (exphere && (exphere % ediv == 0)) {
						int contb = exphere / ediv;
						for (int j = 0; j < contb; j++) answer *= i;
					} else return false;
				}
				if (target != 1) {
					if (ediv != 1) return 0;
					answer *= target;
				}
				origin = answer;
				return ediv;
			};
			
			int simp1 = 1, simp2 = 1;
			if ((simp1 = simplify_inner(bs_n, es_d)) >= 1 && (simp2 = simplify_inner(bs_d, es_d)) >= 1) {
				int simps = __int_gcd(simp1, simp2);
				if (simps >= 1) {
					es_d /= simps;
					auto current_new_base = SymbolicExpr::number((::Rational(bs_n, bs_d)).power(::BigInt(es_n)));
					if (es_d == 1) return current_new_base;
					return SymbolicExpr::power(current_new_base, SymbolicExpr::number(::Rational(::BigInt(1), ::BigInt(es_d))));
				}
			}
		}
		
		// 避免修改，重新获取
		auto rconv = exponent->convert_rational();
		if (rconv.get_denominator() == ::BigInt(2) && rconv.get_numerator() >= ::BigInt(-3) 
			&& rconv.get_numerator() <= ::BigInt(3)) {
			return SymbolicExpr::sqrt(SymbolicExpr::power(base, SymbolicExpr::number(rconv.get_numerator())))->simplify();
		}
    } else if (base->type == SymbolicExpr::Type::Power || base->type == SymbolicExpr::Type::Sqrt) {
		if (base->type == SymbolicExpr::Type::Sqrt) {
			base = SymbolicExpr::power(base->operands[0], SymbolicExpr::number(::Rational(1, 2)));
		}
		auto pwr = SymbolicExpr::multiply(base->operands[1], exponent)->simplify();
		if (pwr->type == SymbolicExpr::Type::Number && pwr->convert_rational() == ::Rational(1))
			return base->operands[0]->simplify();
		return SymbolicExpr::power(base->operands[0]->simplify(), pwr);
	}
	
	if (exponent->is_int() || exponent->is_big_int()) {
		// (a*b)^n = a^n * b^n
		if (base->type == SymbolicExpr::Type::Multiply) {
			auto left = SymbolicExpr::power(base->operands[0], exponent);
			auto right = SymbolicExpr::power(base->operands[1], exponent);
			return SymbolicExpr::multiply(left, right)->simplify();
		}

		auto rconv = exponent->convert_rational();
		if (rconv == ::Rational(0)) return SymbolicExpr::number(1);
		if (rconv == ::Rational(1)) return std::make_shared<SymbolicExpr>(*base);
		if (rconv == ::Rational(-1)) {
			// 这里一定不是整数，尝试分母有理化
			std::function<bool(const std::shared_ptr<SymbolicExpr> &)> processable;
			processable = [&processable](const std::shared_ptr<SymbolicExpr> &obj) -> bool {
				return obj->type == SymbolicExpr::Type::Number || obj->type == SymbolicExpr::Type::Sqrt
					|| (obj->type == SymbolicExpr::Type::Multiply && obj->operands.size() == 2
						&& processable(obj->operands[0]) && processable(obj->operands[1]));
			};
			if (base->type == SymbolicExpr::Type::Add && base->operands.size() == 2 &&
				processable(base->operands[0]) && processable(base->operands[1])) {
				auto new_term = SymbolicExpr::multiply(SymbolicExpr::number(-1), base->operands[1])->simplify();
				auto new_nume = SymbolicExpr::add(base->operands[0], new_term);
				auto new_denom = SymbolicExpr::multiply(base, new_nume)->simplify();
				if (new_denom->type == SymbolicExpr::Type::Number) {
					return SymbolicExpr::multiply(SymbolicExpr::number(new_denom->convert_rational().reciprocal()), 
							new_nume)->simplify();
				} 
			}
		}
		// 防止死循环
		if (rconv.get_denominator() == ::BigInt(1) && rconv.get_numerator() >= ::BigInt(-3) && rconv.get_numerator() < ::BigInt(-1)) {
			// 转为倒数的情况
			return SymbolicExpr::power(SymbolicExpr::power(base, SymbolicExpr::number(rconv.get_numerator().Abs())), SymbolicExpr::number(-1))->simplify();
		}
		if (rconv.get_denominator() == ::BigInt(1) && rconv.get_numerator() > ::BigInt(1) && rconv.get_numerator() <= ::BigInt(4)) {
			int exps = rconv.get_numerator().to_int();
			std::shared_ptr<SymbolicExpr> result = std::make_shared<SymbolicExpr>(*base);
			for (int i = 2; i <= exps; i++)
				result = SymbolicExpr::multiply(result, base)->simplify();
			return result;
		}
	}

    return SymbolicExpr::power(base, exponent);
}

// Helpers for expand
std::shared_ptr<SymbolicExpr> distribute_multiply(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    if (a->type == SymbolicExpr::Type::Add) {
        // (x + y) * b = xb + yb
        auto terms = a->operands; 
        return SymbolicExpr::add(
            distribute_multiply(terms[0], b),
            distribute_multiply(terms[1], b)
        )->expand(); 
    }
    if (b->type == SymbolicExpr::Type::Add) {
        // a * (x + y) = ax + ay
        auto terms = b->operands;
        return SymbolicExpr::add(
            distribute_multiply(a, terms[0]),
            distribute_multiply(a, terms[1])
        )->expand();
    }
    
    return SymbolicExpr::multiply(a, b);
}

static std::shared_ptr<SymbolicExpr> expand_power_integer(const std::shared_ptr<SymbolicExpr>& base, int n) {
    if (n == 0) return SymbolicExpr::number(1);
    if (n == 1) return base;
    if (n < 0) return SymbolicExpr::power(base, SymbolicExpr::number(n)); 
    
    // n > 1
    auto rest = expand_power_integer(base, n - 1);
    // Be careful with large powers; recursion depth.
    return distribute_multiply(base, rest);
}

// --- Trigonometric Simplifications ---

// Helper: matches c * pi or pi * c or pi
// Returns pair {bool matched, Rational c}
static std::pair<bool, Rational> match_pi_multiple(const std::shared_ptr<SymbolicExpr>& expr) {
    if (expr->type == SymbolicExpr::Type::Variable && (expr->identifier == "π" || expr->identifier == "pi")) {
        return {true, Rational(1)};
    }
    
    if (expr->type == SymbolicExpr::Type::Multiply && expr->operands.size() == 2) {
        auto left = expr->operands[0];
        auto right = expr->operands[1];
        
        bool left_is_pi = (left->type == SymbolicExpr::Type::Variable && (left->identifier == "π" || left->identifier == "pi"));
        bool right_is_pi = (right->type == SymbolicExpr::Type::Variable && (right->identifier == "π" || right->identifier == "pi"));
        
        if (left_is_pi && right->is_number()) return {true, right->convert_rational()};
        if (right_is_pi && left->is_number()) return {true, left->convert_rational()};
    }
    
    return {false, Rational(0)};
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_sin() const {
    auto arg = operands[0]->simplify();
    
    // sin(0) = 0
    if (arg->is_number() && arg->convert_rational() == Rational(0)) return SymbolicExpr::number(0);
    
    // sin(-x) = -sin(x)
    if (arg->type == SymbolicExpr::Type::Multiply && arg->operands.size() == 2) {
        if (arg->operands[0]->is_number() && arg->operands[0]->convert_rational() == Rational(-1)) {
            // sin(-u) -> -sin(u)
            return SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::sin(arg->operands[1]))->simplify();
        }
    }
    
    // Check for k * pi
    auto pi_match = match_pi_multiple(arg);
    if (pi_match.first) {
        Rational k = pi_match.second;
        // sin(n * pi) = 0 for integer n
        if (k.is_integer()) return SymbolicExpr::number(0);
        
        // sin( (2n+1)/2 * pi ) = (-1)^n
        // k = n + 1/2 -> 2k = 2n + 1 an odd integer
        Rational two_k = k * Rational(2);
        if (two_k.is_integer()) {
            BigInt bk = two_k.get_numerator(); // This is 2k
            // Check if 2k is odd (not even)
            bool is_even = (bk._size == 0) || !(bk._data[0] & 1);
            if (!is_even) {
                // n = (2k - 1) / 4 ?? No.
                // Let 2k = m (odd). k = m/2.
                // sin(m/2 * pi).
                // m = 1 (pi/2) -> 1
                // m = 3 (3pi/2) -> -1
                // m = 5 -> 1
                // Pattern: (m-1)/2 is even -> 1, odd -> -1 ?
                // m=1 -> 0 -> 1. m=3 -> 1 -> -1.
                // index = (m-1)/2. if index even -> 1.
                
                // We need m % 4.
                // 1 % 4 = 1 -> 1
                // 3 % 4 = 3 -> -1
                // 5 % 4 = 1 -> 1
                // -1 % 4 = -1 -> 3 -> -1
                
                BigInt m = bk;
                long long m_ll = 0;
                try { m_ll = std::stoll(m.to_string()); } catch(...) { return SymbolicExpr::sin(arg); } // Too big
                
                long long rem = m_ll % 4;
                if (rem < 0) rem += 4;
                
                if (rem == 1) return SymbolicExpr::number(1);
                if (rem == 3) return SymbolicExpr::number(-1);
            }
        }
        
        // Standard values: pi/6, pi/4, pi/3
        // sin(pi/6) = 1/2
        if (k == Rational(1, 6)) return SymbolicExpr::number(Rational(1, 2));
        // sin(pi/4) = sqrt(2)/2
        if (k == Rational(1, 4)) {
            auto sq2 = SymbolicExpr::sqrt(SymbolicExpr::number(2));
            return SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), sq2)->simplify();
        }
        // sin(pi/3) = sqrt(3)/2
        if (k == Rational(1, 3)) {
            auto sq3 = SymbolicExpr::sqrt(SymbolicExpr::number(3));
            return SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), sq3)->simplify();
        }
    }

    return SymbolicExpr::sin(arg);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_cos() const {
    auto arg = operands[0]->simplify();
    
    // cos(0) = 1
    if (arg->is_number() && arg->convert_rational() == Rational(0)) return SymbolicExpr::number(1);
    
    // cos(-x) = cos(x)
    if (arg->type == SymbolicExpr::Type::Multiply && arg->operands.size() == 2) {
        if (arg->operands[0]->is_number() && arg->operands[0]->convert_rational() == Rational(-1)) {
            return SymbolicExpr::cos(arg->operands[1])->simplify();
        }
    }
    
    auto pi_match = match_pi_multiple(arg);
    if (pi_match.first) {
        Rational k = pi_match.second;
        // cos(pi) = -1, cos(2pi) = 1
        if (k.is_integer()) {
            BigInt num = k.get_numerator();
            bool num_even = (num._size == 0 || !(num._data[0] & 1));
            if (num_even) return SymbolicExpr::number(1);
            else return SymbolicExpr::number(-1);
        }
        
        // cos(pi/2 + n*pi) = 0
        Rational two_k = k * Rational(2);
        if (two_k.is_integer()) {
            BigInt bk = two_k.get_numerator();
            bool bk_even = (bk._size == 0 || !(bk._data[0] & 1));
            if (!bk_even) return SymbolicExpr::number(0);
        }
        
        // cos(pi/6) = sqrt(3)/2
        if (k == Rational(1, 6)) {
             auto sq3 = SymbolicExpr::sqrt(SymbolicExpr::number(3));
             return SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), sq3)->simplify();
        }
        // cos(pi/4) = sqrt(2)/2
        if (k == Rational(1, 4)) {
             auto sq2 = SymbolicExpr::sqrt(SymbolicExpr::number(2));
             return SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), sq2)->simplify();
        }
        // cos(pi/3) = 1/2
        if (k == Rational(1, 3)) {
             return SymbolicExpr::number(Rational(1, 2));
        }
    }

    return SymbolicExpr::cos(arg);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_tan() const {
    auto arg = operands[0]->simplify();
    
    if (arg->is_number() && arg->convert_rational() == Rational(0)) return SymbolicExpr::number(0);
    
    if (arg->type == SymbolicExpr::Type::Multiply && arg->operands.size() == 2) {
        if (arg->operands[0]->is_number() && arg->operands[0]->convert_rational() == Rational(-1)) {
             return SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::tan(arg->operands[1]))->simplify();
        }
    }
    
    // tan(pi/4) = 1
    auto pi_match = match_pi_multiple(arg);
    if (pi_match.first) {
        Rational k = pi_match.second;
        if (k == Rational(1, 4)) return SymbolicExpr::number(1);
        if (k.is_integer()) return SymbolicExpr::number(0);
    }
    
    return SymbolicExpr::tan(arg);
}

/*
std::shared_ptr<SymbolicExpr> SymbolicExpr::expand() const {
    // Moved to symbolic_poly.cpp
    return std::const_pointer_cast<SymbolicExpr>(shared_from_this());
}
*/

