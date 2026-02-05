#pragma once
#include "bigint.hpp"
#include "rational.hpp"
#include <memory>
#include <string>
#include <variant>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>
#include <functional>
#include <iostream>
#include <cstdlib>

#ifndef _SYMBOLIC_DEBUG
// Default to silent unless explicitly enabled
#define _SYMBOLIC_DEBUG 0
#endif
#ifdef _WIN32
#ifdef LAMINA_CORE_EXPORTS
#define LAMINA_API __declspec(dllexport)
#else
#define LAMINA_API __declspec(dllimport)
#endif
#else
#define LAMINA_API
#endif

// Unified debug stream: controlled at compile-time by _SYMBOLIC_DEBUG and at
// runtime by the environment variable `LAMINA_SYMBOLIC_DEBUG` (0 or 1).
// Usage in source remains the same: `err_stream << "msg" << std::endl;`

class _NullBuffer : public std::streambuf {
public:
	int overflow(int c) override { return c; }
};

static _NullBuffer _null_buffer;
static std::ostream _null_stream(&_null_buffer);

inline bool symbolic_debug_default_enabled() {
#if _SYMBOLIC_DEBUG
	return true;
#else
	return false;
#endif
}

inline bool symbolic_debug_runtime_enabled() {
	static int state = -1;
	if (state != -1) return state == 1;
	// If environment variable LAMINA_SYMBOLIC_DEBUG is set, it overrides default.
	const char *ev = std::getenv("LAMINA_SYMBOLIC_DEBUG");
	if (ev) {
		if (ev[0] == '1') state = 1;
		else state = 0;
	} else {
		state = symbolic_debug_default_enabled() ? 1 : 0;
	}
	return state == 1;
}

inline std::ostream &debug_stream() {
	return symbolic_debug_runtime_enabled() ? std::cerr : _null_stream;
}

#define err_stream debug_stream()

// 符号表达式系统
// 支持精确的数学表达式，不进行数值近似

class LAMINA_API SymbolicExpr {
public:
    enum class Type {
        Number,      // 数字 (BigInt, Rational, int)
        Sqrt,        // 平方根 √
        Root,        // n次方根 √[n]，未使用
        Power,       // 幂次 ^
        Multiply,    // 乘法 *
        Add,         // 加法 +
        Subtract,    // 减法 -，未使用
        Infinity,      // 无限大（number_value 中的 int 决定是正或负）
        Variable,     // 变量 (如 π, e)
        
        // Trigonometric
        Sin, Cos, Tan, Cot, Sec, Csc,
        // Inverse Trigonometric
        ArcSin, ArcCos, ArcTan,
        // Hyperbolic
        Sinh, Cosh, Tanh,
        // Logarithmic
        Ln, Log,
        // Other
        Abs, Fac,
        // Calculus operators (symbolic representation)
        Diff, Integral, Limit,

        // Linear Algebra
        Matrix, // Represents a matrix (list of rows)
        Vector,  // Represents a row (list of elements)
        
    };

    // Compare precedence for canonical ordering
    // Returns -1 if this < other, 1 if this > other, 0 if equal
    int compare(const std::shared_ptr<SymbolicExpr>& other) const;

    // Generic Substitution
    // Replace variable with an expression
    std::shared_ptr<SymbolicExpr> substitute(const std::string& var_name, const std::shared_ptr<SymbolicExpr>& value) const;

    // Expansion (distributive law, etc.)
    std::shared_ptr<SymbolicExpr> expand() const;

    // Polynomial GCD
    static std::shared_ptr<SymbolicExpr> poly_gcd(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

    // Polynomial Resultant
    static std::shared_ptr<SymbolicExpr> poly_resultant(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b, const std::string& var);
    
    // Matrix Determinant
    static std::shared_ptr<SymbolicExpr> determinant(const std::shared_ptr<SymbolicExpr>& mat);

    // Matrix Transpose
    static std::shared_ptr<SymbolicExpr> transpose(const std::shared_ptr<SymbolicExpr>& mat);

    // Matrix Inverse
    static std::shared_ptr<SymbolicExpr> inverse(const std::shared_ptr<SymbolicExpr>& mat);

    // Matrix RREF (Gaussian Elimination)
    static std::shared_ptr<SymbolicExpr> rref(const std::shared_ptr<SymbolicExpr>& mat);
    
    // Characteristic Polynomial
    static std::shared_ptr<SymbolicExpr> charpoly(const std::shared_ptr<SymbolicExpr>& mat, const std::string& lambda);

    // Eigenvalues (returns a list of values)
    static std::shared_ptr<SymbolicExpr> eigenvalues(const std::shared_ptr<SymbolicExpr>& mat);

    // Eigenvectors
    // Returns a list of pairs: {eigenvalue, {eigenvector1, eigenvector2, ...}}
    static std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> eigenvectors(const std::shared_ptr<SymbolicExpr>& mat);

    // System solver (Equations list, Variables list) -> List of solutions (Map: var -> val)
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solve_system(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations, 
        const std::vector<std::string>& vars);

    // Calculus and Limits
    // Calculate limit of this expression as var -> point
    std::shared_ptr<SymbolicExpr> limit(const std::string& var, const std::shared_ptr<SymbolicExpr>& point) const;
    
    // Calculate indefinite integral with respect to var
    std::shared_ptr<SymbolicExpr> integrate(const std::string& var) const;

    // Calculate Taylor/Maclaurin Series expansion
    // var: variable to expand around
    // point: point a in (x-a)^n (expression, e.g., number(0))
    // order: maximum power n
    std::shared_ptr<SymbolicExpr> series(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, int order) const;

    // Helper to create divide (since it is usually multiply by power -1, but having a static helper is nice)
    static std::shared_ptr<SymbolicExpr> divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den);

    // Helper to create integral node (symbolic representation)
    static std::shared_ptr<SymbolicExpr> make_integral(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);
    
    // Helper to create limit node (symbolic representation)
    static std::shared_ptr<SymbolicExpr> make_limit(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var, const std::shared_ptr<SymbolicExpr>& point);

    // Factorization (common factors, basic identities)
    std::shared_ptr<SymbolicExpr> factor() const;

	/*
	哈希要保证的一些因素：
	- 乘法交换律（乘法：奇数二进制位参与运算）
	- 加法交换律（加法：偶数二进制位参与运算）
	- 加法和乘法区分
	这个函数在加法化简时使用。
	*/
	struct HashData {
#define _HASH_PARAMS		ODDBIT, EVENBIT, SQRBIT, HALFBIT
		using HashType = unsigned long long;
		// TODO: 允许多个 HashData 对象之间的不同以减少哈希冲突概率
		// 下方变量的名字不重要。
#define ODDBIT_D 0x555555555555555ull
#define EVENBIT_D 0xAAAAAAAAAAAAAAAull
#define SQRBIT_D 0xBDEEBD77BDEEBD7ull
#define HALFBIT_D 0x969969669699696ull
#define EMPTY 0ull
#define INFINITY_D 0xFFF7FFFFDEADBEEFull
#define PI_H 0x1451419810C0000ull
#define E_H 0x9198101145C0000ull
#define UNKNOWN_H 0xAD0AA0BEEFC0000ull
		
		HashType ODDBIT;
		HashType EVENBIT;
		HashType SQRBIT;
		HashType HALFBIT;
		
		::Rational k = ::Rational(1), ksqrt = ::Rational(1);
		HashType hash = EMPTY;
		std::shared_ptr<SymbolicExpr> hash_obj = SymbolicExpr::number(1);	// 因为是乘法，1为默认状态。此处存储被 hash 的项目对应的值
		
		static HashType bigint_hash(const BigInt& rt) {
			// 直接哈希所有 digits
			HashType weight = 1ull, ans = 0ull;
			auto digits = rt.get_digits();
			for (auto &i : digits) {
				ans = ans * weight + (i + 3ull);
				weight *= 17ull;			// 不用 10 减少哈希冲突
			}
			if (rt.negative) return ~ans;
			return ans;
		}
		
		static HashType rational_hash(const Rational& rt) {
			return bigint_hash(rt.get_numerator()) ^ bigint_hash(rt.get_denominator());
		}
		
		HashType to_single_hash() {
			return (rational_hash(k) & HALFBIT) ^ (rational_hash(ksqrt) & SQRBIT) ^ hash;
		}
		
		// TODO: 考虑优化
		std::shared_ptr<SymbolicExpr> get_combined_k() {
			return SymbolicExpr::multiply(SymbolicExpr::number(k), SymbolicExpr::sqrt(SymbolicExpr::number(ksqrt)))->simplify();
		}
		
		HashData() {
			
		}
		
		HashData(std::shared_ptr<SymbolicExpr> obj, 
			HashType ODDBIT = ODDBIT_D, HashType EVENBIT = EVENBIT_D, HashType SQRBIT = SQRBIT_D, HashType HALFBIT = HALFBIT_D)
			: ODDBIT(ODDBIT), EVENBIT(EVENBIT), SQRBIT(SQRBIT), HALFBIT(HALFBIT) {
			// Evaluate hash
			HashData ld, rd;
			HashType prehash = 0, rterm = 0;
			switch (obj->type) {
				case Type::Number:
					this->k = obj->convert_rational();
					err_stream << "[HPP Debug] Return as value " << k.to_string() << "\n";
					break;
				case Type::Infinity:
					this->hash = INFINITY_D;
					break;
				
				case Type::Sqrt:
					ld = HashData(obj->operands[0], _HASH_PARAMS);
					// sqrt 里面还有 sqrt，取值异或哈希
					this->ksqrt = ld.k;
					ld.k = ::Rational(0);
					this->hash = ld.to_single_hash() * SQRBIT;	// 表明这是个 sqrt，里面没东西则恰好为 0
					this->hash_obj = SymbolicExpr::sqrt(ld.hash_obj);
					break;
				case Type::Multiply:
					ld = HashData(obj->operands[0], _HASH_PARAMS);
					rd = HashData(obj->operands[1], _HASH_PARAMS);
					this->k = ld.k * rd.k;
					this->ksqrt = ld.ksqrt * rd.ksqrt;
					this->hash = (obj->operands[0]->is_number() ? 1 : ld.hash) * (obj->operands[1]->is_number() ? 1 : rd.hash);
					err_stream << "[HPP Debug] LDHash: " << ld.hash_obj->to_string() << ", RDHash: " << rd.hash_obj->to_string() << std::endl;
					err_stream << "[HPP Debug] My hash value is " << this->hash << std::endl;
					err_stream << "[HPP Debug] L applied: " << (obj->operands[0]->is_number() ? 1 : ld.hash) <<
						", R applied: " << (obj->operands[1]->is_number() ? 1 : rd.hash) << std::endl;
					if (!(ld.hash | rd.hash)) this->hash = 0;	// 里面没有东西
					this->hash_obj = SymbolicExpr::multiply(ld.hash_obj, rd.hash_obj)->simplify();
					break;
				case Type::Add:
					ld = HashData(obj->operands[0], _HASH_PARAMS);
					rd = HashData(obj->operands[1], _HASH_PARAMS);
					this->hash = ld.to_single_hash() + rd.to_single_hash();
					this->hash_obj = obj;	// 没有做任何处理
					break;
				case Type::Power:
					// TODO: 此处引入类似根式化简的机制，暂时直接 hash（可能有问题）
					ld = HashData(obj->operands[0], _HASH_PARAMS);
					rd = HashData(obj->operands[1], _HASH_PARAMS);
					// 不是特别恰当，但可以先这样
					// 保证 1，2，-1 等常见数值
					rterm = rd.to_single_hash() - 1;
					this->hash = ld.to_single_hash() ^ rterm ^ (rterm << 8) ^ (rterm << 16) ^ (rterm << 32);
					this->hash_obj = obj;	// 没有做任何处理
					break;
				case Type::Variable:
					if (obj->identifier == "π" || obj->identifier == "pi") this->hash = PI_H;
					else if (obj->identifier == "e") this->hash = E_H;
					else this->hash = UNKNOWN_H;
					this->hash_obj = obj;	// 没有做任何处理
					break;
				default:
					// 如果某个 hash 不能用就重新计算
					this->hash = EMPTY;
					for (auto &i : obj->operands) {
						if (obj->type == Type::Add) {
							this->hash += HashData(i).to_single_hash();	// 令其自然溢出，同时避免异或消除
						} else {
							this->hash *= HashData(i).to_single_hash() + 1;	// 令其自然溢出，同时避免异或消除
						}
					}
					this->hash ^= prehash;
					this->hash_obj = obj;
			}
			// TODO: 可能考虑在这里做根式化简
		}
#undef ODDBIT_D
#undef EVENBIT_D
#undef SQRBIT_D
#undef HALFBIT_D
#undef EMPTY
#undef INFINITY_D
#undef PI_H
#undef E_H
#undef UNKNOWN_H		
	};

    Type type;

    // 数值存储
    std::variant<int, ::BigInt, ::Rational> number_value;

    // 表达式参数
    std::vector<std::shared_ptr<SymbolicExpr>> operands;

    // 字符串标识（用于变量名或操作符）
    std::string identifier;

	// 是否已经化简完成
	bool already_simplified = false;

    // 构造函数
    SymbolicExpr(Type t) : type(t) {}

    // 数字构造函数
    static std::shared_ptr<SymbolicExpr> number(int n) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Number);
        expr->number_value = n;
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> number(const ::BigInt& bi) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Number);
        expr->number_value = bi;
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> number(const ::Rational& r) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Number);
        expr->number_value = r;
        return expr;
    }

	static std::shared_ptr<SymbolicExpr> infinity(int k = 1) {
		auto expr = std::make_shared<SymbolicExpr>(Type::Infinity);
		expr->number_value = k;
		return expr;
	}

    // 平方根构造函数
    static std::shared_ptr<SymbolicExpr> sqrt(std::shared_ptr<SymbolicExpr> operands) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Sqrt);
        expr->operands.push_back(operands);
        return expr;
    }

    // 乘法构造函数
    static std::shared_ptr<SymbolicExpr> multiply(std::shared_ptr<SymbolicExpr> left, std::shared_ptr<SymbolicExpr> right) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Multiply);
        if (right->is_number()) {
            // 将数字移至前端
            expr->operands.push_back(right);
            expr->operands.push_back(left);
            return expr;
        }
        expr->operands.push_back(left);
        expr->operands.push_back(right);
        return expr;
    }

    // 加法构造函数
    static std::shared_ptr<SymbolicExpr> add(std::shared_ptr<SymbolicExpr> left, std::shared_ptr<SymbolicExpr> right) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Add);
        expr->operands.push_back(left);
        expr->operands.push_back(right);
        return expr;
    }

    // 幂次构造函数
    static std::shared_ptr<SymbolicExpr> power(std::shared_ptr<SymbolicExpr> base, std::shared_ptr<SymbolicExpr> exponent) {

        // 直接符号储存
        auto expr = std::make_shared<SymbolicExpr>(Type::Power);
        expr->operands.push_back(base);
        expr->operands.push_back(exponent);
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> sin(std::shared_ptr<SymbolicExpr> op) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Sin);
        expr->operands.push_back(op);
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> cos(std::shared_ptr<SymbolicExpr> op) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Cos);
        expr->operands.push_back(op);
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> tan(std::shared_ptr<SymbolicExpr> op) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Tan);
        expr->operands.push_back(op);
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> ln(std::shared_ptr<SymbolicExpr> op) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Ln);
        expr->operands.push_back(op);
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> log(std::shared_ptr<SymbolicExpr> val, std::shared_ptr<SymbolicExpr> base) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Log);
        expr->operands.push_back(val);
        expr->operands.push_back(base);
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> integral(std::shared_ptr<SymbolicExpr> op, const std::string& var) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Integral);
        expr->operands.push_back(op);
        expr->identifier = var;
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> limit_func(std::shared_ptr<SymbolicExpr> op, const std::string& var, std::shared_ptr<SymbolicExpr> target) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Limit);
        expr->operands.push_back(op);
        expr->operands.push_back(target);
        expr->identifier = var;
        return expr;
    }

    static std::shared_ptr<SymbolicExpr> matrix(const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& elements) {
        auto mat = std::make_shared<SymbolicExpr>(Type::Matrix);
        for (const auto& row : elements) {
            auto vec = std::make_shared<SymbolicExpr>(Type::Vector);
            for (const auto& elem : row) {
                vec->operands.push_back(elem);
            }
            mat->operands.push_back(vec);
        }
        return mat;
    }

    // 变量构造函数
    static std::shared_ptr<SymbolicExpr> variable(const std::string& name) {
        auto expr = std::make_shared<SymbolicExpr>(Type::Variable);
        expr->identifier = name;
        return expr;
    }

    // 化简表达式
    std::shared_ptr<SymbolicExpr> simplify() const;

    // 对变量 var_name 求导
    std::shared_ptr<SymbolicExpr> differentiate(const std::string& var_name) const;

    // 求解方程 eq == 0，返回解集
    static std::vector<std::shared_ptr<SymbolicExpr>> solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name);

    // 转换为字符串表示
    std::string to_string() const;

    // 检查是否为数字
    bool is_number() const { return type == Type::Number; }

    bool get_number_value_is_zero() const {
        if (!is_number()) return false;
        if (std::holds_alternative<int>(number_value)) return std::get<int>(number_value) == 0;
        if (std::holds_alternative<::BigInt>(number_value)) return std::get<::BigInt>(number_value).is_zero();
        if (std::holds_alternative<::Rational>(number_value)) return std::get<::Rational>(number_value).is_zero();
        return false;
    }

    // 检查是否为大整数
    bool is_big_int() const { return is_number() && std::holds_alternative<::BigInt>(number_value); }

    // 检查是否为分数（有理数）
    bool is_rational() const { return is_number() && std::holds_alternative<::Rational>(number_value); }

    // 检查是否为整数
    bool is_int() const { return is_number() && std::holds_alternative<int>(number_value); }

    // 获取数字值（如果是数字的话）
    std::variant<int, ::BigInt, ::Rational> get_number() const {
        if (is_number()) {
            return number_value;
        }
        throw std::runtime_error("Expression is not a number");
    }

    int get_int() const {
        if (is_int()) {
            return std::get<int>(get_number());
        }
        throw std::runtime_error("Expression is not a int");
    }
    ::BigInt get_big_int() const {
        if (is_big_int()) {
            return std::get<BigInt>(get_number());
        }
        throw std::runtime_error("Expression is not a BigInt");
    }
    ::Rational get_rational() const {
        if (is_rational()) {
            return std::get<Rational>(get_number());
        }
        throw std::runtime_error("Expression is not a Rational");
    }
    ::Rational convert_rational() const {
		if (!is_number()) {
			throw std::runtime_error("Expression cannot be converted into Rational");
		}
		if (is_rational()) return get_rational();
		else if (is_big_int()) return ::Rational(get_big_int());
		else return ::Rational(get_int());
	}
	

    // 尝试计算数值（如果可能的话）
    double to_double() const;

private:
    // 内部化简函数
    std::shared_ptr<SymbolicExpr> simplify_sqrt() const;
    std::shared_ptr<SymbolicExpr> simplify_multiply() const;
    std::shared_ptr<SymbolicExpr> simplify_add() const;
    std::shared_ptr<SymbolicExpr> simplify_power() const;
    std::shared_ptr<SymbolicExpr> simplify_sin() const;
    std::shared_ptr<SymbolicExpr> simplify_cos() const;
    std::shared_ptr<SymbolicExpr> simplify_tan() const;
};
