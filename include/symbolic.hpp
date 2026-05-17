#pragma once
#define _USE_MATH_DEFINES
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




class LAMINA_API SymbolicExpr : public std::enable_shared_from_this<SymbolicExpr> {
public:
	
	std::shared_ptr<SymbolicNode> root;
	
    enum class Type {
        Number,      
        Sqrt,        
        Root,        
        Power,       
        Multiply,    
        Add,         
        Subtract,    
        Infinity,      
        Variable,     
        
        
        Sin, Cos, Tan, Cot, Sec, Csc,
        
        ArcSin, ArcCos, ArcTan, Atan2,
        
        Sinh, Cosh, Tanh,
        
        Ln, Log,
        
        Abs, Fac,
        
        Diff, Integral, Limit,

        
        Matrix, 
        Vector,  
        
    };

    
    explicit SymbolicExpr(std::shared_ptr<SymbolicNode> node) : root(std::move(node)) {}

    SymbolicExpr() : root(nullptr) {}

    
    [[deprecated("Use SymbolicExpr(shared_ptr<SymbolicNode>) instead")]]
    SymbolicExpr(Type t) {
		
		
		
		
		
		
		
		if (t == Type::Number) root = std::make_shared<NumberNode>(0);
		else if (t == Type::Variable) root = std::make_shared<VariableNode>("");
		
    }

	
	

    
    
    int compare(const std::shared_ptr<SymbolicExpr>& other) const;

    
    
    std::shared_ptr<SymbolicExpr> substitute(const std::string& var_name, const std::shared_ptr<SymbolicExpr>& value) const;

    
    std::shared_ptr<SymbolicExpr> expand() const;

    bool is_zero() const;
    bool is_one() const;

    
    static std::shared_ptr<SymbolicExpr> poly_gcd(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b);

    
    static std::shared_ptr<SymbolicExpr> poly_resultant(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b, const std::string& var);
    
    
    static std::shared_ptr<SymbolicExpr> determinant(const std::shared_ptr<SymbolicExpr>& mat);

    
    static std::shared_ptr<SymbolicExpr> transpose(const std::shared_ptr<SymbolicExpr>& mat);

    
    static std::shared_ptr<SymbolicExpr> inverse(const std::shared_ptr<SymbolicExpr>& mat);

    
    static std::shared_ptr<SymbolicExpr> rref(const std::shared_ptr<SymbolicExpr>& mat);
    
    
    static std::shared_ptr<SymbolicExpr> charpoly(const std::shared_ptr<SymbolicExpr>& mat, const std::string& lambda);

    
    static std::shared_ptr<SymbolicExpr> eigenvalues(const std::shared_ptr<SymbolicExpr>& mat);

    
    
    static std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> eigenvectors(const std::shared_ptr<SymbolicExpr>& mat);

    
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solve_system(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations, 
        const std::vector<std::string>& vars);

    
    
    std::shared_ptr<SymbolicExpr> limit(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, const std::string& direction = "") const;
    
    
    std::shared_ptr<SymbolicExpr> integrate(const std::string& var) const;

    
    
    
    
    std::shared_ptr<SymbolicExpr> series(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, int order) const;

    
    static std::shared_ptr<SymbolicExpr> divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den);

    
    static std::shared_ptr<SymbolicExpr> make_integral(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);
    
    
    static std::shared_ptr<SymbolicExpr> make_limit(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var, const std::shared_ptr<SymbolicExpr>& point);

    
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
		std::shared_ptr<SymbolicExpr> hash_obj;	
		
		static HashType bigint_hash(const BigInt& rt) {
			
			HashType weight = 1ull, ans = 0ull;
			auto digits = rt.get_digits();
			for (auto &i : digits) {
				ans = ans * weight + (i + 3ull);
				weight *= 17ull;			
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
		
		
		std::shared_ptr<SymbolicExpr> get_combined_k() {
			
			return nullptr; 
		}
		
		HashData() {
			
		}
		
		HashData(std::shared_ptr<SymbolicExpr> obj, 
			HashType ODDBIT = ODDBIT_D, HashType EVENBIT = EVENBIT_D, HashType SQRBIT = SQRBIT_D, HashType HALFBIT = HALFBIT_D)
			: ODDBIT(ODDBIT), EVENBIT(EVENBIT), SQRBIT(SQRBIT), HALFBIT(HALFBIT) {
			
			
			
			
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

	
    [[deprecated("Use SymbolicNode directly")]]
    Type get_type() const;
    
    [[deprecated("Use SymbolicNode directly")]]
    std::vector<std::shared_ptr<SymbolicExpr>> get_operands() const;
    
    [[deprecated("Use SymbolicNode directly")]]
    std::variant<int, ::BigInt, ::Rational> get_number_value() const;
    
    [[deprecated("Use SymbolicNode directly")]]
    std::string get_identifier() const;

	
	
	

    
    static std::shared_ptr<SymbolicExpr> number(int n) {
        return std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(BigInt(n)));
    }

    static std::shared_ptr<SymbolicExpr> number(long long n) {
        return std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(BigInt(n)));
    }
    
    static std::shared_ptr<SymbolicExpr> number(double n) {
        return std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(static_cast<lmmc_real_t>(n)));
    }

    static std::shared_ptr<SymbolicExpr> number(const ::BigInt& bi) {
        return std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(bi));
    }

    static std::shared_ptr<SymbolicExpr> number(const ::Rational& r) {
        return std::make_shared<SymbolicExpr>(std::make_shared<NumberNode>(r));
    }

	static std::shared_ptr<SymbolicExpr> infinity(int k = 1) {
		auto inf_node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, std::vector<std::shared_ptr<SymbolicNode>>{});
		auto inf_expr = std::make_shared<SymbolicExpr>(inf_node);
		if (k < 0) {
			return SymbolicExpr::multiply(SymbolicExpr::number(-1), inf_expr);
		}
		return inf_expr;
	}

    
    static std::shared_ptr<SymbolicExpr> sqrt(std::shared_ptr<SymbolicExpr> operand) {
        
        auto half = std::make_shared<NumberNode>(Rational(1, 2));
        return std::make_shared<SymbolicExpr>(std::make_shared<PowerNode>(operand->root, half));
    }

    
    static std::shared_ptr<SymbolicExpr> multiply(std::shared_ptr<SymbolicExpr> left, std::shared_ptr<SymbolicExpr> right) {
        std::vector<std::shared_ptr<SymbolicNode>> ops = {left->root, right->root};
        return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(ops));
    }

    
    static std::shared_ptr<SymbolicExpr> add(std::shared_ptr<SymbolicExpr> left, std::shared_ptr<SymbolicExpr> right) {
        std::vector<std::shared_ptr<SymbolicNode>> ops = {left->root, right->root};
        return std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(ops));
    }

    
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

    static std::shared_ptr<SymbolicExpr> lambertw(std::shared_ptr<SymbolicExpr> op) {
        if (!op) return nullptr;
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::LambertW, std::vector<std::shared_ptr<SymbolicNode>>{op->root}));
    }

    static std::shared_ptr<SymbolicExpr> log(std::shared_ptr<SymbolicExpr> val, std::shared_ptr<SymbolicExpr> base) {
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Log, std::vector<std::shared_ptr<SymbolicNode>>{val->root, base->root}));
    }

    static std::shared_ptr<SymbolicExpr> atan2(std::shared_ptr<SymbolicExpr> y, std::shared_ptr<SymbolicExpr> x) {
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::Atan2, std::vector<std::shared_ptr<SymbolicNode>>{y->root, x->root}));
    }

    static std::shared_ptr<SymbolicExpr> root_of(std::shared_ptr<SymbolicExpr> poly, const std::string& var, int index) {
        if (!poly) return nullptr;
        auto v = SymbolicExpr::variable(var);
        auto k = SymbolicExpr::number(index);
        std::vector<std::shared_ptr<SymbolicNode>> args = {poly->root, v->root, k->root};
        return std::make_shared<SymbolicExpr>(std::make_shared<FunctionNode>(FunctionNode::FuncType::RootOf, args));
    }

	static std::shared_ptr<SymbolicExpr> eq(std::shared_ptr<SymbolicExpr> lhs, std::shared_ptr<SymbolicExpr> rhs) {
        return std::make_shared<SymbolicExpr>(std::make_shared<RelationalNode>(lhs->root, rhs->root, RelationalNode::Op::EQ));
    }

    static std::shared_ptr<SymbolicExpr> integral(std::shared_ptr<SymbolicExpr> op, const std::string& var) {
        if (!op) return nullptr;
        return op->integrate(var);
    }

    static std::shared_ptr<SymbolicExpr> limit_func(std::shared_ptr<SymbolicExpr> op, const std::string& var, std::shared_ptr<SymbolicExpr> target) {
        if (!op) return nullptr;
        return op->limit(var, target);
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

    
    static std::shared_ptr<SymbolicExpr> variable(const std::string& name) {
        return std::make_shared<SymbolicExpr>(std::make_shared<VariableNode>(name));
    }

    
    std::shared_ptr<SymbolicExpr> simplify() const;

    std::shared_ptr<SymbolicExpr> simplify_trig() const;

    std::shared_ptr<SymbolicExpr> differentiate(const std::string& var_name) const;

private:
   
   std::shared_ptr<SymbolicExpr> differentiate_legacy(const std::string& var_name) const;

public:
    
    static std::vector<std::shared_ptr<SymbolicExpr>> solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name);

    
    std::string to_string() const;

    
    bool is_number() const { 
		return root && root->is_number(); 
	}

    bool get_number_value_is_zero() const {
        return root && root->is_zero();
    }

    
    bool is_big_int() const { 
		if (!is_number()) return false;
		auto node = std::dynamic_pointer_cast<NumberNode>(root);
		return node && std::holds_alternative<BigInt>(node->value);
	}

    
    bool is_rational() const { 
		if (!is_number()) return false;
		auto node = std::dynamic_pointer_cast<NumberNode>(root);
		return node && std::holds_alternative<Rational>(node->value);
	}

    
    bool is_int() const { 
		if (!is_number()) return false;
		auto node = std::dynamic_pointer_cast<NumberNode>(root);
		if (!node) return false;
		
		// BigInt is always an integer
		if (std::holds_alternative<BigInt>(node->value)) return true;
		
		// Rational is integer when denominator is 1
		if (std::holds_alternative<Rational>(node->value)) {
		    return std::get<Rational>(node->value).is_integer();
		}
		
		// lmmc_real_t is integer when it equals its rounded value
		if (std::holds_alternative<lmmc_real_t>(node->value)) {
		    lmmc_real_t v = std::get<lmmc_real_t>(node->value);
		    lmmc_real_t rounded = std::round(v);
		    int eq;
		    lmmc_double_nearly_equal_tol(v, rounded, 1e-12, 1e-12, &eq);
		    return eq != 0;
		}
		
		return false;
	}

    
    std::variant<int, ::BigInt, ::Rational> get_number() const {
        if (is_number()) {
            auto node = std::dynamic_pointer_cast<NumberNode>(root);
			if (std::holds_alternative<BigInt>(node->value)) return std::get<BigInt>(node->value);
			if (std::holds_alternative<Rational>(node->value)) return std::get<Rational>(node->value);
			
			return (int)std::get<lmmc_real_t>(node->value);
        }
        throw std::runtime_error("Expression is not a number");
    }

    int get_int() const {
        
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
	

    
    lmmc_real_t to_numeric() const;

private:
    
    std::shared_ptr<SymbolicExpr> simplify_sqrt() const;
    std::shared_ptr<SymbolicExpr> simplify_multiply() const;
    std::shared_ptr<SymbolicExpr> simplify_add() const;
    std::shared_ptr<SymbolicExpr> simplify_power() const;
    std::shared_ptr<SymbolicExpr> simplify_sin() const;
    std::shared_ptr<SymbolicExpr> simplify_cos() const;
    std::shared_ptr<SymbolicExpr> simplify_tan() const;
    std::shared_ptr<SymbolicExpr> simplify_ln() const;
};
