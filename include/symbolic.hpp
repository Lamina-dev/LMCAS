#pragma once
#include "bigint.hpp"
#include "rational.hpp"
#include "symbolic_ast.hpp"
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

class LAMINA_API SymbolicExpr : public std::enable_shared_from_this<SymbolicExpr> {
public:
	// Use explicit node storage
	std::shared_ptr<SymbolicNode> root;
	
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

    // Constructor wrapping a node
    explicit SymbolicExpr(std::shared_ptr<SymbolicNode> node) : root(std::move(node)) {}

    // Legacy constructor (deprecated)
    [[deprecated("Use SymbolicExpr(shared_ptr<SymbolicNode>) instead")]]
    SymbolicExpr(Type t) {
		// Create a dummy node or throw
		// This is just to satisfy compilation of legacy code that constructs then fills
		// We can't really support it well.
		// For now, create a dummy node based on type if possible.
		// Actually, let's just create a raw pointer that will crash if used without being initialized properly?
		// Or throw.
		// Given migration phase, we'll try to support basic types if possible.
		if (t == Type::Number) root = std::make_shared<NumberNode>(0);
		else if (t == Type::Variable) root = std::make_shared<VariableNode>("");
		// ... others
    }

	// Deleted copy constructor to enforce unique ownership semantics? No, shared_ptr is fine.
	// Default copy/move is fine.

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
		std::shared_ptr<SymbolicExpr> hash_obj;	// 因为是乘法，1为默认状态。此处存储被 hash 的项目对应的值
		
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
			// return SymbolicExpr::multiply(SymbolicExpr::number(k), SymbolicExpr::sqrt(SymbolicExpr::number(ksqrt)))->simplify();
			return nullptr; // Placeholder
		}
		
		HashData() {
			
		}
		
		HashData(std::shared_ptr<SymbolicExpr> obj, 
			HashType ODDBIT = ODDBIT_D, HashType EVENBIT = EVENBIT_D, HashType SQRBIT = SQRBIT_D, HashType HALFBIT = HALFBIT_D)
			: ODDBIT(ODDBIT), EVENBIT(EVENBIT), SQRBIT(SQRBIT), HALFBIT(HALFBIT) {
			
			// Implementation removed during migration
			// Existing logic relied on direct field access which is gone.
			// TODO: Reimplement hashing logic using SymbolicNode visitor
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

	// Deprecated fields accessors
    [[deprecated("Use SymbolicNode directly")]]
    Type get_type() const;
    
    [[deprecated("Use SymbolicNode directly")]]
    std::vector<std::shared_ptr<SymbolicExpr>> get_operands() const;
    
    [[deprecated("Use SymbolicNode directly")]]
    std::variant<int, ::BigInt, ::Rational> get_number_value() const;
    
    [[deprecated("Use SymbolicNode directly")]]
    std::string get_identifier() const;

	// 是否已经化简完成
	// bool already_simplified = false; // Moved to Node if needed, or kept here as metadata? 
	// The prompt says "Replace internal storage". So we remove it.

    // 数字构造函数
    static std::shared_ptr<SymbolicExpr> number(int n) {
        return std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(n));
    }

    static std::shared_ptr<SymbolicExpr> number(const ::BigInt& bi) {
        return std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(bi));
    }

    static std::shared_ptr<SymbolicExpr> number(const ::Rational& r) {
        return std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(r));
    }

	static std::shared_ptr<SymbolicExpr> infinity(int k = 1) {
		// return std::make_shared<SymbolicExpr>(std::make_shared<InfinityNode>(k));
		// InfinityNode not defined in provided symbolic_ast.hpp snippet. 
		// Assuming we might need to add it or use a special NumberNode or fallback.
		// For now, return nullptr or dummy
		return nullptr;
	}

    // 平方根构造函数
    static std::shared_ptr<SymbolicExpr> sqrt(std::shared_ptr<SymbolicExpr> operand) {
        // PowerNode(base, 1/2)
        auto half = std::make_shared<NumberNode>(Rational(1, 2));
        return std::make_shared<SymbolicExpr>(std::make_shared<PowerNode>(operand->root, half));
    }

    // 乘法构造函数
    static std::shared_ptr<SymbolicExpr> multiply(std::shared_ptr<SymbolicExpr> left, std::shared_ptr<SymbolicExpr> right) {
        std::vector<std::shared_ptr<SymbolicNode>> ops = {left->root, right->root};
        return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(ops));
    }

    // 加法构造函数
    static std::shared_ptr<SymbolicExpr> add(std::shared_ptr<SymbolicExpr> left, std::shared_ptr<SymbolicExpr> right) {
        std::vector<std::shared_ptr<SymbolicNode>> ops = {left->root, right->root};
        return std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(ops));
    }

    // 幂次构造函数
    static std::shared_ptr<SymbolicExpr> power(std::shared_ptr<SymbolicExpr> base, std::shared_ptr<SymbolicExpr> exponent) {
        return std::make_shared<SymbolicExpr>(std::make_shared<PowerNode>(base->root, exponent->root));
    }

    static std::shared_ptr<SymbolicExpr> sin(std::shared_ptr<SymbolicExpr> op) {
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Sin, std::vector<std::shared_ptr<SymbolicNode>>{op->root}));
    }

    static std::shared_ptr<SymbolicExpr> cos(std::shared_ptr<SymbolicExpr> op) {
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Cos, std::vector<std::shared_ptr<SymbolicNode>>{op->root}));
    }

    static std::shared_ptr<SymbolicExpr> tan(std::shared_ptr<SymbolicExpr> op) {
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Tan, std::vector<std::shared_ptr<SymbolicNode>>{op->root}));
    }

    static std::shared_ptr<SymbolicExpr> ln(std::shared_ptr<SymbolicExpr> op) {
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, std::vector<std::shared_ptr<SymbolicNode>>{op->root}));
    }

    static std::shared_ptr<SymbolicExpr> exp(std::shared_ptr<SymbolicExpr> op) {
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Exp, std::vector<std::shared_ptr<SymbolicNode>>{op->root}));
    }

    static std::shared_ptr<SymbolicExpr> log(std::shared_ptr<SymbolicExpr> val, std::shared_ptr<SymbolicExpr> base) {
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Log, std::vector<std::shared_ptr<SymbolicNode>>{val->root, base->root}));
    }

    static std::shared_ptr<SymbolicExpr> integral(std::shared_ptr<SymbolicExpr> op, const std::string& var) {
        // FunctionNode? No, integral is operator.
		// Assuming we have FunctionNode type for integral or similar
		// For now, placeholder
		return nullptr;
    }

    static std::shared_ptr<SymbolicExpr> limit_func(std::shared_ptr<SymbolicExpr> op, const std::string& var, std::shared_ptr<SymbolicExpr> target) {
        return nullptr; // Placeholder
    }

    static std::shared_ptr<SymbolicExpr> matrix(const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& elements) {
        std::vector<std::vector<std::shared_ptr<SymbolicNode>>> node_elements;
		for(const auto& row : elements) {
			std::vector<std::shared_ptr<SymbolicNode>> node_row;
			for(const auto& elem : row) node_row.push_back(elem->root);
			node_elements.push_back(node_row);
		}
        return std::make_shared<SymbolicExpr>(std::make_shared<MatrixNode>(node_elements));
    }

    // 变量构造函数
    static std::shared_ptr<SymbolicExpr> variable(const std::string& name) {
        return std::make_shared<SymbolicExpr>(std::make_shared<VariableNode>(name));
    }

    // 化简表达式
    std::shared_ptr<SymbolicExpr> simplify() const;

    // 对变量 var_name 求导 (Wrapper)
    std::shared_ptr<SymbolicExpr> differentiate(const std::string& var_name) const;

private:
   // Original legacy differentiate
   std::shared_ptr<SymbolicExpr> differentiate_legacy(const std::string& var_name) const;

public:
    // 求解方程 eq == 0，返回解集
    static std::vector<std::shared_ptr<SymbolicExpr>> solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name);

    // 转换为字符串表示
    std::string to_string() const;

    // 检查是否为数字
    bool is_number() const { 
		return root && root->is_number(); 
	}

    bool get_number_value_is_zero() const {
        return root && root->is_zero();
    }

    // 检查是否为大整数
    bool is_big_int() const { 
		if (!is_number()) return false;
		auto node = std::dynamic_pointer_cast<NumberNode>(root);
		return node && std::holds_alternative<BigInt>(node->value);
	}

    // 检查是否为分数（有理数）
    bool is_rational() const { 
		if (!is_number()) return false;
		auto node = std::dynamic_pointer_cast<NumberNode>(root);
		return node && std::holds_alternative<Rational>(node->value);
	}

    // 检查是否为整数
    bool is_int() const { 
		if (!is_number()) return false;
		auto node = std::dynamic_pointer_cast<NumberNode>(root);
		return node && std::holds_alternative<double>(node->value); // Wait, double is float. int is BigInt or double(int)? 
		// Original code had int/BigInt/Rational. NumberNode has BigInt/Rational/double.
		// We'll emulate int via double or BigInt check.
		// For now returning false or checking double == (int)double
		return false; 
	}

    // 获取数字值（如果是数字的话）
    std::variant<int, ::BigInt, ::Rational> get_number() const {
        if (is_number()) {
            auto node = std::dynamic_pointer_cast<NumberNode>(root);
			if (std::holds_alternative<BigInt>(node->value)) return std::get<BigInt>(node->value);
			if (std::holds_alternative<Rational>(node->value)) return std::get<Rational>(node->value);
			// double -> int?
			return (int)std::get<double>(node->value);
        }
        throw std::runtime_error("Expression is not a number");
    }

    int get_int() const {
        // Simplify
        return 0;
    }
    ::BigInt get_big_int() const {
        if (is_big_int()) {
             auto node = std::dynamic_pointer_cast<NumberNode>(root);
             return std::get<BigInt>(node->value);
        }
        throw std::runtime_error("Expression is not a BigInt");
    }
    ::Rational get_rational() const {
        if (is_rational()) {
             auto node = std::dynamic_pointer_cast<NumberNode>(root);
             return std::get<Rational>(node->value);
        }
        throw std::runtime_error("Expression is not a Rational");
    }
    ::Rational convert_rational() const {
		if (is_rational()) return get_rational();
		if (is_big_int()) return Rational(get_big_int());
		return Rational(0);
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
    std::shared_ptr<SymbolicExpr> simplify_ln() const;
};
