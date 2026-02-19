#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include "../include/symbolic.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include "../include/visitors/differentiation_visitor.hpp"

// =========================================================
// Helper Visitors
// =========================================================

class SubstituteVisitor : public SymbolicVisitor {
    std::string var_name;
    std::shared_ptr<SymbolicNode> new_val;
public:
    std::shared_ptr<SymbolicNode> result;

    SubstituteVisitor(std::string v, std::shared_ptr<SymbolicNode> val) 
        : var_name(std::move(v)), new_val(std::move(val)) {}

    std::shared_ptr<SymbolicNode> get_result() const { return result; }

    void visit(NumberNode& node) override {
        result = node.clone();
    }

    void visit(VariableNode& node) override {
        if (node.name == var_name) {
            result = new_val->clone();
        } else {
            result = node.clone();
        }
    }

    void visit(AddNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for (const auto& op : node.operands) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = std::make_shared<AddNode>(new_ops);
    }

    void visit(MultiplyNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for (const auto& op : node.operands) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = std::make_shared<MultiplyNode>(new_ops);
    }

    void visit(PowerNode& node) override {
        node.base->accept(*this);
        auto new_base = result;
        node.exponent->accept(*this);
        auto new_exp = result;
        result = std::make_shared<PowerNode>(new_base, new_exp);
    }

    void visit(FunctionNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_args;
        for (const auto& arg : node.arguments) {
            arg->accept(*this);
            new_args.push_back(result);
        }
        result = std::make_shared<FunctionNode>(node.type, new_args);
    }

    void visit(MatrixNode& node) override {
        // Simple matrix substitution
        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage)) {
            const auto& dense = std::get<MatrixNode::DenseStorage>(node.storage);
            MatrixNode::DenseStorage new_dense;
            for(const auto& e : dense) {
                e->accept(*this);
                new_dense.push_back(result);
            }
            result = std::make_shared<MatrixNode>(node.rows, node.cols, new_dense);
        } else {
            const auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage);
            MatrixNode::SparseStorage new_sparse;
            for(const auto& [idx, val] : sparse) {
                val->accept(*this);
                new_sparse[idx] = result;
            }
            result = std::make_shared<MatrixNode>(node.rows, node.cols, new_sparse);
        }
    }
};

class IntegrationVisitor : public SymbolicVisitor {
    std::string var;
public:
    std::shared_ptr<SymbolicNode> result;

    IntegrationVisitor(std::string v) : var(std::move(v)) {}

    std::shared_ptr<SymbolicNode> get_result() const { return result; }

    void visit(NumberNode& node) override {
        // int(c) dx = c*x
        std::vector<std::shared_ptr<SymbolicNode>> ops;
        ops.push_back(node.clone());
        ops.push_back(std::make_shared<VariableNode>(var));
        result = std::make_shared<MultiplyNode>(ops);
    }

    void visit(VariableNode& node) override {
        // int(x) dx = x^2/2
        if (node.name == var) {
            auto two = std::make_shared<NumberNode>(BigInt(2));
            auto x_pow_2 = std::make_shared<PowerNode>(node.clone(), two);
            // Use Rational for 1/2
            auto half = std::make_shared<NumberNode>(Rational(1, 2)); 
            std::vector<std::shared_ptr<SymbolicNode>> ops = {half, x_pow_2};
            result = std::make_shared<MultiplyNode>(ops);
        } else {
            // int(a) dx = a*x (constant variable)
            std::vector<std::shared_ptr<SymbolicNode>> ops;
            ops.push_back(node.clone());
            ops.push_back(std::make_shared<VariableNode>(var));
            result = std::make_shared<MultiplyNode>(ops);
        }
    }

    void visit(AddNode& node) override {
        // int(f + g) = int(f) + int(g)
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for (const auto& op : node.operands) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = std::make_shared<AddNode>(new_ops);
    }

    void visit(MultiplyNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> c_ops;
        std::vector<std::shared_ptr<SymbolicNode>> x_ops;
        
        for (const auto& op : node.operands) {
            if (std::dynamic_pointer_cast<NumberNode>(op)) {
                c_ops.push_back(op);
            } else {
                x_ops.push_back(op);
            }
        }
        
        if (x_ops.empty()) {
             // Constant
             std::vector<std::shared_ptr<SymbolicNode>> final_ops = node.operands;
             final_ops.push_back(std::make_shared<VariableNode>(var));
             result = std::make_shared<MultiplyNode>(final_ops);
             return;
        }

        if (x_ops.size() == 1) {
            // c * f(x)
            x_ops[0]->accept(*this);
            auto int_f = result;
            
            std::vector<std::shared_ptr<SymbolicNode>> res_ops = c_ops;
            res_ops.push_back(int_f);
            result = std::make_shared<MultiplyNode>(res_ops);
            return;
        }

        // Falback: return integral(expr)
        result = std::make_shared<NumberNode>(BigInt(0)); // Stub
    }
    
    void visit(PowerNode& node) override {
        // int(x^n)
        bool base_is_x = false;
        if (auto v = std::dynamic_pointer_cast<VariableNode>(node.base)) {
            if (v->name == var) base_is_x = true;
        }
        
        if (base_is_x) {
            // Check exponent for -1
            if (auto num = std::dynamic_pointer_cast<NumberNode>(node.exponent)) {
                 if (std::holds_alternative<double>(num->value) && std::get<double>(num->value) == -1.0) {
                     // ln(x)
                     std::vector<std::shared_ptr<SymbolicNode>> args = {node.base};
                     result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args);
                     return;
                 }
                 if (std::holds_alternative<BigInt>(num->value) && std::get<BigInt>(num->value) == BigInt(-1)) {
                     std::vector<std::shared_ptr<SymbolicNode>> args = {node.base};
                     result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args);
                     return;
                 }
                 // x^(n+1)/(n+1)
                 auto one = std::make_shared<NumberNode>(BigInt(1));
                 std::vector<std::shared_ptr<SymbolicNode>> add_ops = {node.exponent, one};
                 auto n_plus_1 = std::make_shared<AddNode>(add_ops);
                 
                 // Normalize n+1 for display
                 NormalizationVisitor norm_exp;
                 n_plus_1->accept(norm_exp);
                 auto n_plus_1_sched = norm_exp.get_result();
                 
                 auto new_pow = std::make_shared<PowerNode>(node.base, n_plus_1_sched);
                 
                 // 1/(n+1) -> (n+1)^-1
                 auto minus_one = std::make_shared<NumberNode>(BigInt(-1));
                 auto denom = std::make_shared<PowerNode>(n_plus_1_sched, minus_one);
                 
                 std::vector<std::shared_ptr<SymbolicNode>> mul_ops = {denom, new_pow};
                 result = std::make_shared<MultiplyNode>(mul_ops);
                 return;
            }
        }
        
        // Fallback for non-x base or complex exponent
        std::vector<std::shared_ptr<SymbolicNode>> ops;
        ops.push_back(node.clone());
        ops.push_back(std::make_shared<VariableNode>(var));
        result = std::make_shared<MultiplyNode>(ops);
    }
    
    void visit(FunctionNode& node) override {
        // Special cases: sin, cos, exp
        bool arg_is_x = false;
        if (node.arguments.size() == 1) {
            if (auto v = std::dynamic_pointer_cast<VariableNode>(node.arguments[0])) {
                if (v->name == var) arg_is_x = true;
            }
        }

        if (arg_is_x) {
            std::vector<std::shared_ptr<SymbolicNode>> args = node.arguments;
            if (node.type == FunctionNode::FuncType::Sin) {
                // -cos(x)
                auto cos_node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Cos, args);
                auto minus_one = std::make_shared<NumberNode>(BigInt(-1));
                std::vector<std::shared_ptr<SymbolicNode>> ops = {minus_one, cos_node};
                result = std::make_shared<MultiplyNode>(ops);
                return;
            } else if (node.type == FunctionNode::FuncType::Cos) {
                // sin(x)
                result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Sin, args);
                return;
            } else if (node.type == FunctionNode::FuncType::Exp) {
                result = node.clone();
                return;
            }
        }
        
        result = node.clone(); 
    }

    void visit(MatrixNode& node) override {
        // Element-wise integration
        // ... (Similar to SubstituteVisitor)
        result = node.clone(); // Stub
    }
};


// =========================================================
// SymbolicExpr Implementation
// =========================================================

// Comparison
int SymbolicExpr::compare(const std::shared_ptr<SymbolicExpr>& other) const {
    if (!root || !other->root) return 0;
    // SymbolicNode::compare is public and handles type priority and same-type comparison
    return root->compare(*other->root);
}

// Substitution
std::shared_ptr<SymbolicExpr> SymbolicExpr::substitute(const std::string& var_name, const std::shared_ptr<SymbolicExpr>& value) const {
    if (!root) return nullptr;
    SubstituteVisitor v(var_name, value->root);
    root->accept(v);
    
    // Auto-simplify
    NormalizationVisitor norm;
    v.get_result()->accept(norm);
    
    return std::make_shared<SymbolicExpr>(norm.get_result());
}

// Expansion
std::shared_ptr<SymbolicExpr> SymbolicExpr::expand() const {
    if (!root) return nullptr;
    // TODO: Implement proper ExpandVisitor
    return std::make_shared<SymbolicExpr>(root->clone());
}

// String Representation
std::string SymbolicExpr::to_string() const {
    PrintVisitor printer;
    if (root) {
        root->accept(printer);
        return printer.get_result();
    }
    return "null";
}

// Double Conversion
double SymbolicExpr::to_double() const {
    if (!root) return 0.0;
    // Simple evaluation for single NumberNodes
    if (auto num = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<double>(num->value)) return std::get<double>(num->value);
        if (std::holds_alternative<::BigInt>(num->value)) return std::get<::BigInt>(num->value).to_double();
        if (std::holds_alternative<::Rational>(num->value)) return std::get<::Rational>(num->value).to_double();
    }
    return 0.0; 
}

bool SymbolicExpr::is_zero() const {
    if (!root) return false;
    if (auto num = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<double>(num->value)) return std::get<double>(num->value) == 0.0;
        if (std::holds_alternative<::BigInt>(num->value)) return std::get<::BigInt>(num->value).to_int() == 0;
        if (std::holds_alternative<::Rational>(num->value)) return std::get<::Rational>(num->value).get_numerator().to_int() == 0;
    }
    return false;
}

bool SymbolicExpr::is_one() const {
    if (!root) return false;
    if (auto num = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<double>(num->value)) return std::get<double>(num->value) == 1.0;
        if (std::holds_alternative<::BigInt>(num->value)) return std::get<::BigInt>(num->value).to_int() == 1;
        if (std::holds_alternative<::Rational>(num->value)) return std::get<::Rational>(num->value).to_double() == 1.0;
    }
    return false;
}

// Simplification
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify() const {
    if (!root) return nullptr;
    NormalizationVisitor v;
    root->accept(v);
    return std::make_shared<SymbolicExpr>(v.get_result());
}

// Differentiation
std::shared_ptr<SymbolicExpr> SymbolicExpr::differentiate(const std::string& var_name) const {
    if (!root) return nullptr;
    DifferentiationVisitor v(var_name);
    root->accept(v);
    
    // Auto-simplify result to avoid 0+1, 1*x etc.
    NormalizationVisitor norm;
    v.get_result()->accept(norm);
    
    return std::make_shared<SymbolicExpr>(norm.get_result());
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::differentiate_legacy(const std::string& var_name) const {
    return differentiate(var_name);
}

// Solving moved to solve_visitor.cpp
/* 
std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::solve(std::shared_ptr<SymbolicExpr> eq, const std::string& var_name) {
   // Stub removed
   return {};
} 
*/


// Factorization (Stub)
std::shared_ptr<SymbolicExpr> SymbolicExpr::factor() const {
    return simplify();
}

// Calculus Implementation
std::shared_ptr<SymbolicExpr> SymbolicExpr::limit(const std::string& var, const std::shared_ptr<SymbolicExpr>& point) const {
    if (!root) return nullptr;
    
    // Check for fraction form: N * D^-1
    // This is a naive check targeting the L'Hopital test case
    std::shared_ptr<SymbolicExpr> num_expr = nullptr;
    std::shared_ptr<SymbolicExpr> den_expr = nullptr;
    
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(root)) {
        std::vector<std::shared_ptr<SymbolicNode>> num_nodes;
        std::vector<std::shared_ptr<SymbolicNode>> den_nodes;
        
        for(const auto& op : mul->operands) {
             bool is_den = false;
             if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                 if (auto sn = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                      // Check for -1
                      if (std::holds_alternative<double>(sn->value) && std::get<double>(sn->value) == -1.0) is_den = true;
                      else if (std::holds_alternative<BigInt>(sn->value) && std::get<BigInt>(sn->value) == BigInt(-1)) is_den = true;
                      else if (std::holds_alternative<Rational>(sn->value) && std::get<Rational>(sn->value) == Rational(-1)) is_den = true;
                      
                      if(is_den) {
                          den_nodes.push_back(pow->base);
                      }
                 }
             }
             if (!is_den) {
                 num_nodes.push_back(op);
             }
        }
        
        if (!den_nodes.empty()) {
            std::shared_ptr<SymbolicNode> n_node, d_node;
            if (num_nodes.empty()) n_node = std::make_shared<NumberNode>(BigInt(1));
            else if (num_nodes.size() == 1) n_node = num_nodes[0];
            else n_node = std::make_shared<MultiplyNode>(num_nodes);
            
            if (den_nodes.size() == 1) d_node = den_nodes[0];
            else d_node = std::make_shared<MultiplyNode>(den_nodes);
            
            num_expr = std::make_shared<SymbolicExpr>(n_node);
            den_expr = std::make_shared<SymbolicExpr>(d_node);
        }
    }
    
    if (num_expr && den_expr) {
        // Evaluate at point
        auto val_n = num_expr->substitute(var, point);
        auto val_d = den_expr->substitute(var, point);
        
        // Check if 0/0
        // Using to_double for quick check (dangerous but valid for simple cases)
        double vn = val_n->to_double();
        double vd = val_d->to_double();
        
        if (std::abs(vn) < 1e-9 && std::abs(vd) < 1e-9) {
             // L'Hopital
             auto dN = num_expr->differentiate(var);
             auto dD = den_expr->differentiate(var);
             
             // New expression: dN/dD
             // ratio = dN * dD^-1
             std::vector<std::shared_ptr<SymbolicNode>> ops;
             ops.push_back(dN->root);
             
             auto minus_one = std::make_shared<NumberNode>(BigInt(-1));
             auto dD_inv = std::make_shared<PowerNode>(dD->root, minus_one);
             ops.push_back(dD_inv);
             
             auto ratio_node = std::make_shared<MultiplyNode>(ops);
             auto ratio = std::make_shared<SymbolicExpr>(ratio_node);
             
             return ratio->limit(var, point);
        }
    }

    // Direct substitution
    SubstituteVisitor v(var, point->root);
    root->accept(v);
    auto res = v.get_result();
    
    NormalizationVisitor norm;
    res->accept(norm);
    
    return std::make_shared<SymbolicExpr>(norm.get_result());
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::integrate(const std::string& var) const {
    if (!root) return nullptr;
    IntegrationVisitor v(var);
    root->accept(v);
    
    NormalizationVisitor norm;
    v.get_result()->accept(norm);
    
    return std::make_shared<SymbolicExpr>(norm.get_result());
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::series(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, int order) const {
    if (!root) return nullptr;
    return std::make_shared<SymbolicExpr>(root->clone());
}

// Static Helpers Stubs - ONLY those declared but not defined inline in header
std::shared_ptr<SymbolicExpr> SymbolicExpr::divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den) {
    return SymbolicExpr::multiply(num, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::make_integral(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    return expr;
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::make_limit(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var, const std::shared_ptr<SymbolicExpr>& point) {
    return expr;
}

/*
// Matrix Stubs moved to solve_visitor.cpp (or implemented there)
std::shared_ptr<SymbolicExpr> SymbolicExpr::determinant(const std::shared_ptr<SymbolicExpr>& mat) { 
    return SymbolicExpr::number(0); 
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::transpose(const std::shared_ptr<SymbolicExpr>& mat) { 
    return mat; 
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::inverse(const std::shared_ptr<SymbolicExpr>& mat) { 
    return mat; 
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::rref(const std::shared_ptr<SymbolicExpr>& mat) { 
    return mat; 
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::charpoly(const std::shared_ptr<SymbolicExpr>& mat, const std::string& lambda) { 
    return SymbolicExpr::number(0); 
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::eigenvalues(const std::shared_ptr<SymbolicExpr>& mat) { 
    return SymbolicExpr::number(0); 
}
std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::vector<std::shared_ptr<SymbolicExpr>>>> SymbolicExpr::eigenvectors(const std::shared_ptr<SymbolicExpr>& mat) { 
    return {}; 
}
*/
// solve_system moved to solve_visitor.cpp

// Polynomial Stubs
std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_gcd(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    return SymbolicExpr::number(1);
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_resultant(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b, const std::string& var) {
    return SymbolicExpr::number(0);
}

// Internal Simplify Stubs (Legacy)
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_sqrt() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_multiply() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_add() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_power() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_sin() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_cos() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_tan() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_ln() const { return simplify(); }

// Deprecated Accessors Implementation
SymbolicExpr::Type SymbolicExpr::get_type() const {
    if (!root) return Type::Number;
    
    // Check dynamic types
    if (std::dynamic_pointer_cast<NumberNode>(root)) return Type::Number;
    if (std::dynamic_pointer_cast<VariableNode>(root)) return Type::Variable;
    if (std::dynamic_pointer_cast<AddNode>(root)) return Type::Add;
    if (std::dynamic_pointer_cast<MultiplyNode>(root)) return Type::Multiply;
    if (std::dynamic_pointer_cast<PowerNode>(root)) return Type::Power;
    
    // Use FunctionNode::type instead of func_type
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(root)) {
        switch (func->type) {
            case FunctionNode::FuncType::Sin: return Type::Sin;
            case FunctionNode::FuncType::Cos: return Type::Cos;
            case FunctionNode::FuncType::Tan: return Type::Tan;
            case FunctionNode::FuncType::Ln: return Type::Ln;
            case FunctionNode::FuncType::Log: return Type::Log;
            case FunctionNode::FuncType::Exp: return Type::Power; 
            case FunctionNode::FuncType::Abs: return Type::Abs;
            case FunctionNode::FuncType::Sqrt: return Type::Sqrt;
            default: return Type::Variable; 
        }
    }
    
    if (auto mat = std::dynamic_pointer_cast<MatrixNode>(root)) {
        if (mat->rows == 1 || mat->cols == 1) return Type::Vector;
        return Type::Matrix;
    }

    return Type::Number;
}

std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::get_operands() const {
    std::vector<std::shared_ptr<SymbolicExpr>> ops;
    if (!root) return ops;

    if (auto add = std::dynamic_pointer_cast<AddNode>(root)) {
        for (const auto& op : add->operands) ops.push_back(std::make_shared<SymbolicExpr>(op));
    } else if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(root)) {
        for (const auto& op : mul->operands) ops.push_back(std::make_shared<SymbolicExpr>(op));
    } else if (auto pow = std::dynamic_pointer_cast<PowerNode>(root)) {
        ops.push_back(std::make_shared<SymbolicExpr>(pow->base));
        ops.push_back(std::make_shared<SymbolicExpr>(pow->exponent));
    } else if (auto func = std::dynamic_pointer_cast<FunctionNode>(root)) {
        for (const auto& arg : func->arguments) ops.push_back(std::make_shared<SymbolicExpr>(arg));
    }
    return ops;
}

// Header declares: std::variant<int, ::BigInt, ::Rational> get_number_value() const;
std::variant<int, ::BigInt, ::Rational> SymbolicExpr::get_number_value() const { 
    if (auto node = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<BigInt>(node->value)) return std::get<BigInt>(node->value);
        if (std::holds_alternative<Rational>(node->value)) return std::get<Rational>(node->value);
        if (std::holds_alternative<double>(node->value)) return (int)std::get<double>(node->value); // Approx or cast
    }
    return 0; // Returns variant with int=0
}

std::string SymbolicExpr::get_identifier() const {
    if (auto v = std::dynamic_pointer_cast<VariableNode>(root)) return v->name;
    return "";
}

