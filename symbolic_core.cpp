#include "symbolic.hpp"
#include "symbolic_internal.hpp"
#include <iostream>
#include <cmath>
#include <climits>
#include <algorithm>
#include <map>

// 局部调试开关（仅影响本文件的调试输出）
namespace {
	// 当需要大量调试信息时，可在编译时或者运行时修改此值
	static bool SYMBOLIC_DEBUG = false;
}

// 表达式结构等价比较（尽量保守）：
// - 数字使用有理数比较
// - 变量按名字比较
// - 幂、根、加、乘等以结构递归比较（对乘法/加法采用字符串回退以避免复杂置换）
bool symbolic_equal(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
	if (a == b) return true;
	if (!a || !b) return false;
	if (a->type != b->type) return false;
	switch (a->type) {
		case SymbolicExpr::Type::Number:
			return a->convert_rational() == b->convert_rational();
		case SymbolicExpr::Type::Variable:
			return a->identifier == b->identifier;
		case SymbolicExpr::Type::Sqrt:
			if (a->operands.empty() || b->operands.empty()) return false;
			return symbolic_equal(a->operands[0], b->operands[0]);
		case SymbolicExpr::Type::Power:
			if (a->operands.size() < 2 || b->operands.size() < 2) return false;
			return symbolic_equal(a->operands[0], b->operands[0]) && symbolic_equal(a->operands[1], b->operands[1]);
		case SymbolicExpr::Type::Multiply:
		case SymbolicExpr::Type::Add: {
			// 支持交换律：把同类项扁平化为项列表，比较多重集合（顺序无关）
			std::function<void(const std::shared_ptr<SymbolicExpr>&, std::vector<std::shared_ptr<SymbolicExpr>> &)> collect;
			collect = [&](const std::shared_ptr<SymbolicExpr>& obj, std::vector<std::shared_ptr<SymbolicExpr>> &out) {
				if (!obj) return;
				if (obj->type == a->type && obj->operands.size() == 2) {
					collect(obj->operands[0], out);
					collect(obj->operands[1], out);
				} else {
					out.push_back(obj->simplify());
				}
			};
			std::vector<std::shared_ptr<SymbolicExpr>> la, lb;
			collect(a, la);
			collect(b, lb);
			if (la.size() != lb.size()) return false;
			std::vector<char> used(lb.size(), 0);
			for (size_t i = 0; i < la.size(); ++i) {
				bool matched = false;
				for (size_t j = 0; j < lb.size(); ++j) {
					if (used[j]) continue;
					if (symbolic_equal(la[i], lb[j])) { used[j] = 1; matched = true; break; }
				}
				if (!matched) return false;
			}
			return true;
		}
		case SymbolicExpr::Type::Infinity:
			return a->number_value == b->number_value;
		default:
			return a->to_string() == b->to_string();
	}
}

// 符号表达式的化简实现（分发器）
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify() const {
	// 添加“化简”标记，避免 simplify 重复调用导致效率降低（似乎暂时不可用？）
	//if (already_simplified) return std::make_shared<SymbolicExpr>(*this);
	
	static int current_simplify_level = 0;
	const int max_simplify_level = 200; // Increased limit for complex calculus
	if (current_simplify_level > max_simplify_level) {
		// Stop recursion to prevent stack overflow and endless warnings
		// err_stream << "[Warning] SymbolicExpr: Reaching maximum simplifying depth (" << max_simplify_level << "). Return unsimplified.\n";
		return std::make_shared<SymbolicExpr>(*this);
	}
	current_simplify_level++;
	
    auto intcall = [&]() {
		switch (type) {
			case Type::Number:
			case Type::Variable:
				return std::make_shared<SymbolicExpr>(*this);
				
			case Type::Sqrt:
				return simplify_sqrt();
				
			case Type::Multiply:
				return simplify_multiply();
				
			case Type::Add:
				return simplify_add();

			case Type::Power:
				return simplify_power();
			
                        case Type::Log: {
                                auto val = operands[0]->simplify();
                                auto base = operands[1]->simplify();
                                // log_b(x) -> ln(x)/ln(b)
                                return SymbolicExpr::multiply(
                                    SymbolicExpr::ln(val),
                                    SymbolicExpr::power(SymbolicExpr::ln(base), SymbolicExpr::number(-1))
                                )->simplify();
                        }
                        
                        case Type::Ln:
                                return simplify_ln();

                        case Type::Sin: return simplify_sin();
                        case Type::Cos: return simplify_cos();
                        case Type::Tan: return simplify_tan();

                        default:
                                // TODO: Recursive simplify for other functions
                                return std::make_shared<SymbolicExpr>(*this);
                }
        };
	
	auto res = intcall();
	current_simplify_level--;
	res->already_simplified = true;
	return res;
}

std::string SymbolicExpr::to_string() const {
	
	auto get_output = [](std::shared_ptr<const SymbolicExpr> expr) -> std::string {
		const std::string lbrace = std::string("("), rbrace = std::string(")");
		if ((expr->type == SymbolicExpr::Type::Number && (expr->convert_rational().get_denominator() == ::BigInt(1))) || expr->type == SymbolicExpr::Type::Variable
			|| expr->type == SymbolicExpr::Type::Sqrt) return expr->to_string();
		else return lbrace + expr->to_string() + rbrace;
	};
	
    switch (type) {
        case Type::Number:
            if (std::holds_alternative<int>(number_value)) {
                return std::to_string(std::get<int>(number_value));
            } else if (std::holds_alternative<::BigInt>(number_value)) {
                return std::get<::BigInt>(number_value).to_string();
            } else if (std::holds_alternative<::Rational>(number_value)) {
                return std::get<::Rational>(number_value).to_string();
            }
            return "0";
            
        case Type::Variable:
            return identifier;
			
		case Type::Infinity:
			if (std::get<int>(number_value) > 0) return "inf";
			else return "-inf";
            
        case Type::Sqrt:
            if (operands.empty()) return "√()";
            return "√" + get_output(operands[0]);
            
        case Type::Multiply:
            if (operands.size() < 2) return "*(?)";

            if (operands[0]->is_number() && operands[1]->type == Type::Sqrt) {
                return operands[0]->to_string() + operands[1]->to_string();
            }
            return get_output(operands[0]) + "*" + get_output(operands[1]);
            
        case Type::Add: {
            if (operands.size() < 2) return "+(?)";
			
            std::vector<std::shared_ptr<SymbolicExpr>> terms;
            std::function<void(const std::shared_ptr<SymbolicExpr>&)> flatten_add;
            flatten_add = [&](const std::shared_ptr<SymbolicExpr>& expr) {
                if (expr->type == Type::Add && expr->operands.size() == 2) {
                    flatten_add(expr->operands[0]);
                    flatten_add(expr->operands[1]);
                } else {
                    terms.push_back(expr);
                }
            };
            flatten_add(std::make_shared<SymbolicExpr>(*this));

            std::vector<std::string> result_terms;
			result_terms.reserve(terms.size());
            for (auto &term : terms) {
                result_terms.push_back(get_output(term));
            }

            if (result_terms.empty()) return "0";
            std::string res = result_terms[0];
            for (size_t i = 1; i < result_terms.size(); ++i) {
                res += "+" + result_terms[i];
            }
            return res;
        }
            
        case Type::Power:
            if (operands.size() < 2) return "^(?)";
            return get_output(operands[0]) + "^" + get_output(operands[1]);
            
        case Type::Sin: return "sin(" + operands[0]->to_string() + ")";
        case Type::Cos: return "cos(" + operands[0]->to_string() + ")";
        case Type::Tan: return "tan(" + operands[0]->to_string() + ")";
        case Type::Ln:  return "ln(" + operands[0]->to_string() + ")";
        case Type::Log: return "log_" + operands[1]->to_string() + "(" + operands[0]->to_string() + ")";
        
        case Type::Diff: return "diff(" + operands[0]->to_string() + ", " + identifier + ")";
        case Type::Integral: return "int(" + operands[0]->to_string() + ", " + identifier + ")";
        case Type::Limit: return "lim(" + operands[0]->to_string() + ", " + identifier + "->" + operands[1]->to_string() + ")";
        
        case Type::Matrix: {
            std::string res = "[";
            for (size_t i = 0; i < operands.size(); ++i) {
                if (i > 0) res += ", ";
                res += operands[i]->to_string();
            }
            res += "]";
            return res;
        }
        case Type::Vector: {
            std::string res = "[";
            for (size_t i = 0; i < operands.size(); ++i) {
                if (i > 0) res += ", ";
                auto& op = operands[i];
                const std::string lbrace = std::string("("), rbrace = std::string(")");
                if ((op->type == SymbolicExpr::Type::Number && (op->convert_rational().get_denominator() == ::BigInt(1))) || op->type == SymbolicExpr::Type::Variable || op->type == SymbolicExpr::Type::Sqrt)
                     res += op->to_string();
                else res += op->to_string(); // Vector elements usually don't need outer parens unless complex
            }
            res += "]";
            return res;
        }

        default:
            return "Unknown";
    }
}

double SymbolicExpr::to_double() const {
    switch (type) {
        case Type::Number:
            if (std::holds_alternative<int>(number_value)) {
                return static_cast<double>(std::get<int>(number_value));
            } else if (std::holds_alternative<::BigInt>(number_value)) {
                return std::get<::BigInt>(number_value).to_double();
            } else if (std::holds_alternative<::Rational>(number_value)) {
                return std::get<::Rational>(number_value).to_double();
            }
            return 0.0;

        case Type::Variable:
            if (identifier == "π" || identifier == "pi") {
                #ifdef M_PI
                return M_PI;
                #else
                return 3.14159265358979323846;
                #endif
            }
			if (identifier == "e") {
				return 2.718281828459045;
			}
            // 其他变量仍抛异常
            throw std::runtime_error("Symbolic variable cannot be converted to double");

		case Type::Infinity:
			if (std::get<int>(number_value) > 0) return std::numeric_limits<double>::infinity();
			else return -std::numeric_limits<double>::infinity();

        case Type::Sqrt:
            if (!operands.empty()) {
                return std::sqrt(operands[0]->to_double());
            }
            return 0.0;

        case Type::Multiply:
            if (operands.size() >= 2) {
                return operands[0]->to_double() * operands[1]->to_double();
            }
            return 0.0;
			
		case Type::Add:
			if (operands.size() >= 2) {
                return operands[0]->to_double() + operands[1]->to_double();
            }
            return 0.0;
			
		case Type::Power:
			return std::pow(operands[0]->to_double(), operands[1]->to_double());

        default:
            return 0.0;
    }
}

int SymbolicExpr::compare(const std::shared_ptr<SymbolicExpr>& other) const {
    if (this == other.get()) return 0;
    if (!other) return 1; 

    // 1. Compare Type Hierarchy
    if (this->type != other->type) {
        return (static_cast<int>(this->type) < static_cast<int>(other->type)) ? -1 : 1;
    }

    // 2. Compare Content based on Type
    switch (this->type) {
        case Type::Number: {
             auto r1 = this->convert_rational();
             auto r2 = other->convert_rational();
             if (r1 < r2) return -1;
             if (r1 > r2) return 1;
             return 0;
        }
        case Type::Variable: {
            if (this->identifier < other->identifier) return -1;
            if (this->identifier > other->identifier) return 1;
            return 0;
        }
        case Type::Infinity: {
             if (this->number_value < other->number_value) return -1;
             if (this->number_value > other->number_value) return 1;
             return 0;
        }
        default: break; 
    }

    // 3. Compare Operands Lexicographically
    if (this->operands.size() != other->operands.size()) {
        return (this->operands.size() < other->operands.size()) ? -1 : 1;
    }

    for (size_t i = 0; i < this->operands.size(); ++i) {
        int cmp = this->operands[i]->compare(other->operands[i]);
        if (cmp != 0) return cmp;
    }

    return 0;
}
