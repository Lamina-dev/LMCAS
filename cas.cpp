/*
 * Lamina CAS (计算机代数系统) 扩展实现
 * 提供符号计算、微积分、方程求解等功能
 */
#include "cas.hpp"
#include "value.hpp"
#include "symbolic.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <map>

// 兼容旧接口的简单定义，如果需要导出给其他模块
extern "C" {
    // 这里可以放导出声明，如果通过 dlsym 调用
}

Value cas_parse(const std::vector<Value>& args) {
    if (args.size() != 1 || !args[0].is_string()) {
        std::cerr << "Error: cas_parse() requires one string argument" << std::endl;
        return Value();
    }
    std::string s = std::get<std::string>(args[0].data);
    return Value(SymbolicExpr::variable(s)); 
}

Value cas_simplify(const std::vector<Value>& args) {
    if (args.size() != 1) {
        std::cerr << "Error: cas_simplify() requires one argument" << std::endl;
        return Value();
    }

    try {
        if (!args[0].as_symbolic_compatible() && !args[0].is_symbolic()) {
             return args[0];
        }
        auto expr = args[0].as_symbolic();
        auto simplified = expr->simplify();
        return Value(simplified);
    } catch (const std::exception& e) {
        std::cerr << "CAS Simplify Error: " << e.what() << std::endl;
        return Value();
    }
}

Value cas_differentiate(const std::vector<Value>& args) {
    if (args.size() < 2) {
         std::cerr << "Usage: differentiate(expr, var)" << std::endl;
         return Value();
    }

    auto expr = args[0].as_symbolic();
    std::string var;
    
    if (args[1].is_string()) {
        var = std::get<std::string>(args[1].data);
    } else if (args[1].is_symbolic()) {
        auto s = args[1].as_symbolic();
        if (s->type == SymbolicExpr::Type::Variable) {
            var = s->identifier;
        } else {
             std::cerr << "Error: differentiation variable must be a symbol." << std::endl;
             return Value();
        }
    } else {
        std::cerr << "Error: differentiation variable invalid." << std::endl;
        return Value();
    }
    
    if (!expr) return Value();
    return Value(expr->differentiate(var));
}

Value cas_integrate(const std::vector<Value>& args) {
     if (args.size() < 2) return Value();
     auto expr = args[0].as_symbolic();
     std::string var;
     if (args[1].is_string()) var = std::get<std::string>(args[1].data);
     else if (args[1].is_symbolic() && args[1].as_symbolic()->type == SymbolicExpr::Type::Variable) var = args[1].as_symbolic()->identifier;
     else return Value();
     
     if (!expr) return Value();
     return Value(SymbolicExpr::integral(expr, var));
}

Value cas_limit(const std::vector<Value>& args) {
     if (args.size() < 3) return Value();
     auto expr = args[0].as_symbolic();
     std::string var;
     if (args[1].is_string()) var = std::get<std::string>(args[1].data);
     else if (args[1].is_symbolic() && args[1].as_symbolic()->type == SymbolicExpr::Type::Variable) var = args[1].as_symbolic()->identifier;
     else return Value();
     auto target = args[2].as_symbolic();
     
     if (!expr) return Value();
     return Value(SymbolicExpr::limit_func(expr, var, target));
}

Value cas_solve(const std::vector<Value>& args) {
     if (args.size() < 2) return Value();
     auto expr = args[0].as_symbolic();
     std::string var;
     if (args[1].is_string()) var = std::get<std::string>(args[1].data);
     else if (args[1].is_symbolic() && args[1].as_symbolic()->type == SymbolicExpr::Type::Variable) var = args[1].as_symbolic()->identifier;
     else return Value();
     
     if (!expr) return Value();
     auto solutions = SymbolicExpr::solve(expr, var);
     if (solutions.empty()) return Value("No solution");
     // Return first solution
     return Value(solutions[0]);
}

Value cas_evaluate(const std::vector<Value>& args) {
    if (args.size() < 1) {
        std::cerr << "Error: cas_evaluate() requires at least one argument" << std::endl;
        return Value();
    }

    try {
        auto expr = args[0].as_symbolic();
        // evaluate 目前 SymbolicExpr 只有 to_double()
        // 或者保留符号求值？
        // 如果要做变量代换，需要 traverse SymbolicExpr Tree.
        // 暂时只支持 to_double
        return Value(expr->to_double());
    } catch (const std::exception& e) {
        std::cerr << "CAS Evaluate Error: " << e.what() << std::endl;
        return Value();
    }
}

// 简单的存储映射
static std::map<std::string, std::shared_ptr<SymbolicExpr>> stored_expressions;

Value cas_store(const std::vector<Value>& args) {
    if (args.size() != 2 || !args[0].is_string()) {
        std::cerr << "Error: cas_store() requires name and expression" << std::endl;
        return Value();
    }
    std::string name = std::get<std::string>(args[0].data);
    auto expr = args[1].as_symbolic();
    stored_expressions[name] = expr;
    return Value("Stored " + name);
}

Value cas_load(const std::vector<Value>& args) {
    if (args.size() != 1 || !args[0].is_string()) {
        std::cerr << "Error: cas_load() requires name" << std::endl;
        return Value();
    }
    std::string name = std::get<std::string>(args[0].data);
    if (stored_expressions.count(name)) {
        return Value(stored_expressions[name]);
    }
    return Value();
}
