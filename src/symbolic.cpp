#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include "../include/symbolic.hpp"
#include "../include/integration.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include "../include/visitors/differentiation_visitor.hpp"
#include "../include/visitors/expand_visitor.hpp"
#include "../include/visitors/limit_visitor.hpp"
#include "../include/poly_utils.hpp"
#include "../include/matcher.hpp"
#include "../include/integration.hpp"





class VariablesVisitor : public SymbolicVisitor {
public:
    std::set<std::string> vars;

    void visit(NumberNode& node) override {} 
    void visit(VariableNode& node) override {
        vars.insert(node.name);
    }
    void visit(AddNode& node) override {
        for(auto& op : node.operands) op->accept(*this);
    }
    void visit(MultiplyNode& node) override {
        for(auto& op : node.operands) op->accept(*this);
    }
    void visit(PowerNode& node) override {
        node.base->accept(*this);
        node.exponent->accept(*this);
    }
    void visit(FunctionNode& node) override {
        for(auto& arg : node.arguments) arg->accept(*this);
    }
    void visit(MatrixNode& node) override {
        
    }
};

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
        
        std::vector<std::shared_ptr<SymbolicNode>> ops;
        ops.push_back(node.clone());
        ops.push_back(std::make_shared<VariableNode>(var));
        result = std::make_shared<MultiplyNode>(ops);
    }

    void visit(VariableNode& node) override {
        if (node.name == var) {
            auto two = std::make_shared<NumberNode>(BigInt(2));
            auto x_pow_2 = std::make_shared<PowerNode>(node.clone(), two);
            auto half = std::make_shared<NumberNode>(Rational(1, 2));
            result = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{half, x_pow_2});
        } else {
            result = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{node.clone(), std::make_shared<VariableNode>(var)});
        }
    }

    void visit(AddNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for (auto& op : node.operands) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = std::make_shared<AddNode>(new_ops);
    }

    void visit(MultiplyNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> x_ops, c_ops;
        for (auto& op : node.operands) {
            if (lamina::depends_on_var(op, var)) {
                x_ops.push_back(op);
            } else {
                c_ops.push_back(op);
            }
        }

        if (x_ops.empty()) {
            result = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{node.clone(), std::make_shared<VariableNode>(var)});
            return;
        }

        if (x_ops.size() == 1) {
            x_ops[0]->accept(*this);
            auto int_f = result;
            if (c_ops.empty()) {
                result = int_f;
            } else {
                std::vector<std::shared_ptr<SymbolicNode>> res_ops = c_ops;
                res_ops.push_back(int_f);
                result = std::make_shared<MultiplyNode>(res_ops);
            }
            return;
        }
        
        result = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{node.clone(), std::make_shared<VariableNode>(var)});
    }
    
    void visit(PowerNode& node) override {
        
        bool base_is_x = false;
        if (auto v = std::dynamic_pointer_cast<VariableNode>(node.base)) {
            if (v->name == var) base_is_x = true;
        }
        
        if (base_is_x) {
            
            if (auto num = std::dynamic_pointer_cast<NumberNode>(node.exponent)) {
                 if (std::holds_alternative<double>(num->value) && std::get<double>(num->value) == -1.0) {
                     
                     std::vector<std::shared_ptr<SymbolicNode>> args = {node.base};
                     result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args);
                     return;
                 }
                 if (std::holds_alternative<BigInt>(num->value) && std::get<BigInt>(num->value) == BigInt(-1)) {
                     std::vector<std::shared_ptr<SymbolicNode>> args = {node.base};
                     result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args);
                     return;
                 }
                 
                 auto one = std::make_shared<NumberNode>(BigInt(1));
                 std::vector<std::shared_ptr<SymbolicNode>> add_ops = {node.exponent, one};
                 auto n_plus_1 = std::make_shared<AddNode>(add_ops);
                 
                 
                 NormalizationVisitor norm_exp;
                 n_plus_1->accept(norm_exp);
                 auto n_plus_1_sched = norm_exp.get_result();
                 
                 auto new_pow = std::make_shared<PowerNode>(node.base, n_plus_1_sched);
                 
                 
                 auto minus_one = std::make_shared<NumberNode>(BigInt(-1));
                 auto denom = std::make_shared<PowerNode>(n_plus_1_sched, minus_one);
                 
                 std::vector<std::shared_ptr<SymbolicNode>> mul_ops = {denom, new_pow};
                 result = std::make_shared<MultiplyNode>(mul_ops);
                 return;
            }
        }
        
        
        std::vector<std::shared_ptr<SymbolicNode>> ops;
        ops.push_back(node.clone());
        ops.push_back(std::make_shared<VariableNode>(var));
        result = std::make_shared<MultiplyNode>(ops);
    }
    
    void visit(FunctionNode& node) override {
        
        bool arg_is_x = false;
        if (node.arguments.size() == 1) {
            if (auto v = std::dynamic_pointer_cast<VariableNode>(node.arguments[0])) {
                if (v->name == var) arg_is_x = true;
            }
        }

        if (arg_is_x) {
            std::vector<std::shared_ptr<SymbolicNode>> args = node.arguments;
            if (node.type == FunctionNode::FuncType::Sin) {
                
                auto cos_node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Cos, args);
                auto minus_one = std::make_shared<NumberNode>(BigInt(-1));
                std::vector<std::shared_ptr<SymbolicNode>> ops = {minus_one, cos_node};
                result = std::make_shared<MultiplyNode>(ops);
                return;
            } else if (node.type == FunctionNode::FuncType::Cos) {
                
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
        
        
        result = node.clone(); 
    }
};







int SymbolicExpr::compare(const std::shared_ptr<SymbolicExpr>& other) const {
    if (!root || !other->root) return 0;
    
    return root->compare(*other->root);
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::   substitute(const std::string& var_name, const std::shared_ptr<SymbolicExpr>& value) const {
    if (!root) return nullptr;
    SubstituteVisitor v(var_name, value->root);
    root->accept(v);
    
    
    NormalizationVisitor norm;
    v.get_result()->accept(norm);
    
    return std::make_shared<SymbolicExpr>(norm.get_result());
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::expand() const {
    if (!root) return nullptr;
    
    ExpandVisitor v;
    root->accept(v);
    
    auto result_node = v.get_result();
    if (!result_node) return nullptr;
    
    // Normalize after expansion to combine like terms like x^2 + x^2 + x^2
    return std::make_shared<SymbolicExpr>(result_node)->simplify();
}


std::string SymbolicExpr::to_string() const {
    PrintVisitor printer;
    if (root) {
        root->accept(printer);
        return printer.get_result();
    }
    return "null";
}


double SymbolicExpr::to_double() const {
    if (!root) return 0.0;
    
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


std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify() const {
    if (!root) return nullptr;
    NormalizationVisitor v;
    root->accept(v);
    return std::make_shared<SymbolicExpr>(v.get_result());
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_trig() const {
    auto res = simplify();
    if (!res->root) return nullptr;

    static lamina::RewriteEngine engine;
    static bool init = false;
    
    if (!init) {
        init = true;
        using namespace lamina;
        auto x_val = wildcard("x");
        auto x = std::make_shared<SymbolicExpr>(x_val);
        
        // Rule: sin(x)^2 + cos(x)^2 -> 1
        auto sinx = SymbolicExpr::sin(x);
        auto cosx = SymbolicExpr::cos(x);
        auto n2 = SymbolicExpr::number(2);
        auto sin2 = SymbolicExpr::power(sinx, n2);
        auto cos2 = SymbolicExpr::power(cosx, n2);
        
        // Match sin(x)^2 + cos(x)^2 -> 1
        auto pat1 = SymbolicExpr::add(sin2, cos2);
        engine.add_rule(Rule(*pat1, *SymbolicExpr::number(1), {"x"}));
        
        // Match cos(x)^2 + sin(x)^2 -> 1
        auto pat2 = SymbolicExpr::add(cos2, sin2);
        engine.add_rule(Rule(*pat2, *SymbolicExpr::number(1), {"x"}));

        // Rule: sin(2*x) -> 2*sin(x)*cos(x)
        auto two_x = SymbolicExpr::multiply(n2, x);
        auto sin2x = SymbolicExpr::sin(two_x);
        auto two_sin_cos = SymbolicExpr::multiply(n2, 
            SymbolicExpr::multiply(sinx, cosx));
        engine.add_rule(Rule(*sin2x, *two_sin_cos, {"x"}));
        
        // Rule: cos(2*x) -> cos(x)^2 - sin(x)^2
        auto cos2x = SymbolicExpr::cos(two_x);
        auto cos2_sub_sin2 = SymbolicExpr::add(cos2, 
            SymbolicExpr::multiply(sin2, SymbolicExpr::number(-1)));
        engine.add_rule(Rule(*cos2x, *cos2_sub_sin2, {"x"}));
    }
    
    auto simplified = engine.apply(*res);
    auto result_ptr = std::make_shared<SymbolicExpr>(simplified);
    return result_ptr->simplify();
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::differentiate(const std::string& var_name) const {
    if (!root) return nullptr;
    DifferentiationVisitor v(var_name);
    root->accept(v);
    
    
    NormalizationVisitor norm;
    v.get_result()->accept(norm);
    
    return std::make_shared<SymbolicExpr>(norm.get_result());
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::differentiate_legacy(const std::string& var_name) const {
    return differentiate(var_name);
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::factor() const {
    auto simp = simplify();
    if (!simp || !simp->root) return simp;

    if (auto add_node = std::dynamic_pointer_cast<AddNode>(simp->root)) {
        // Try common factor
        std::shared_ptr<SymbolicExpr> common = nullptr;
        for (const auto& op : add_node->operands) {
             auto expr_op = std::make_shared<SymbolicExpr>(op);
             if (!common) common = expr_op;
             else common = poly_gcd(common, expr_op);
        }
        
        if (common && !common->is_one() && !common->is_zero()) {
             std::vector<std::shared_ptr<SymbolicNode>> new_ops;
             for (const auto& op : add_node->operands) {
                  auto term = std::make_shared<SymbolicExpr>(op);
                  // Manually divide: term / common
                  auto inv_common = power(common, number(-1));
                  auto quot = multiply(term, inv_common);
                  quot = quot->simplify();
                  new_ops.push_back(quot->root);
             }
             auto new_sum = std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(new_ops));
             
             // Construct Multiply(common, sum) WITHOUT calling simplify() on result (to avoid re-expansion)
             std::vector<std::shared_ptr<SymbolicNode>> final_ops = {common->root, new_sum->root};
             return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(final_ops));
        }
        
        // Try quadratic
        VariablesVisitor vv;
        simp->root->accept(vv);
        if (vv.vars.size() == 1) {
             std::string var = *vv.vars.begin();
             try {
                 auto poly = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(simp, var);
                 if (poly.degree() == 2) {
                      auto solutions = solve(simp, var);
                      if (solutions.size() == 2) {
                           auto leading = poly.coeffs[2].val; 
                           if (!leading) leading = number(1);
                           leading = leading->simplify();
                           
                           auto x_node = std::make_shared<SymbolicExpr>(std::make_shared<VariableNode>(var));
                           
                           auto t1 = SymbolicExpr::add(x_node, multiply(solutions[0], number(-1)))->simplify();
                           auto t2 = SymbolicExpr::add(x_node, multiply(solutions[1], number(-1)))->simplify();
                           
                           std::vector<std::shared_ptr<SymbolicNode>> factors;
                           if (!leading->is_one()) factors.push_back(leading->root);
                           factors.push_back(t1->root);
                           factors.push_back(t2->root);
                           
                           return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(factors));
                      }
                 }
             } catch (...) {}
        }
    }
    
    return simp;
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::limit(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, const std::string& direction) const {
    if (!root) return nullptr;
    
    LimitVisitor v(var, point->root, direction);
    root->accept(v);
    
    if (v.get_result()) {
        return std::make_shared<SymbolicExpr>(v.get_result());
    }
    
    return nullptr;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::integrate(const std::string& var) const {
    if (!root) return nullptr;
    
    // Use the advanced Integrator class instead of the basic IntegrationVisitor
    lamina::Integrator integrator;
    SymbolicExpr result = integrator.integrate(*this, var);
    
    // The Integrator returns a value, we return a shared_ptr
    auto res_ptr = std::make_shared<SymbolicExpr>(result);
    return res_ptr->simplify(); 
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::series(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, int order) const {
    if (!root) return nullptr;
    return std::make_shared<SymbolicExpr>(root->clone());
}


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



std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_gcd(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    if (!a || !b) return nullptr;
    
    try {
        struct VarVisitor : public SymbolicVisitor {
            std::set<std::string> vars;
            void visit(NumberNode&) override {}
            void visit(VariableNode& n) override { vars.insert(n.name); }
            void visit(AddNode& n) override { for(auto& op : n.operands) op->accept(*this); }
            void visit(MultiplyNode& n) override { for(auto& op : n.operands) op->accept(*this); }
            void visit(PowerNode& n) override { n.base->accept(*this); n.exponent->accept(*this); }
            void visit(FunctionNode& n) override { for(auto& arg : n.arguments) arg->accept(*this); }
            void visit(MatrixNode& n) override {}
            void visit(RelationalNode& n) override { n.left->accept(*this); n.right->accept(*this); }
        } vv;
        if (a->root) a->root->accept(vv);
        if (b->root) b->root->accept(vv);

        if (vv.vars.empty()) return SymbolicExpr::number(1);
        std::string var = *vv.vars.begin();
        
        auto pa = lamina::symbolic_to_poly<BigInt>(a, var);
        auto pb = lamina::symbolic_to_poly<BigInt>(b, var);
        auto g = lamina::Polynomial<BigInt>::gcd(pa, pb);
        return lamina::poly_to_symbolic(g);
    } catch (...) {
        return SymbolicExpr::number(1);
    }
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_resultant(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b, const std::string& var) {
    try {
        auto pb = lamina::symbolic_to_poly<Rational>(b, var);
        auto pa = lamina::symbolic_to_poly<Rational>(a, var);
        
        if (pb.degree() == 1) {
            auto m = pb.coeffs[1];
            auto c = pb.coeffs[0];
            
            auto root_val = SymbolicExpr::number(-c/m);
            auto substitution = a->substitute(var, root_val);
            auto factor = SymbolicExpr::power(SymbolicExpr::number(m), SymbolicExpr::number((int)pa.degree()));
            
            return SymbolicExpr::multiply(factor, substitution)->simplify();
        }
    } catch (...) {}
    
    return SymbolicExpr::number(0);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_sqrt() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_multiply() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_add() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_power() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_sin() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_cos() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_tan() const { return simplify(); }
std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_ln() const { return simplify(); }


SymbolicExpr::Type SymbolicExpr::get_type() const {
    if (!root) return Type::Number;
    
    if (std::dynamic_pointer_cast<NumberNode>(root)) return Type::Number;
    if (std::dynamic_pointer_cast<VariableNode>(root)) return Type::Variable;
    if (std::dynamic_pointer_cast<AddNode>(root)) return Type::Add;
    if (std::dynamic_pointer_cast<MultiplyNode>(root)) return Type::Multiply;
    if (std::dynamic_pointer_cast<PowerNode>(root)) return Type::Power;
    
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
            case FunctionNode::FuncType::Atan2: return Type::Atan2;
            case FunctionNode::FuncType::ArcSin: return Type::ArcSin;
            case FunctionNode::FuncType::ArcCos: return Type::ArcCos;
            case FunctionNode::FuncType::ArcTan: return Type::ArcTan;
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


std::variant<int, ::BigInt, ::Rational> SymbolicExpr::get_number_value() const { 
    if (auto node = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<BigInt>(node->value)) return std::get<BigInt>(node->value);
        if (std::holds_alternative<Rational>(node->value)) return std::get<Rational>(node->value);
        if (std::holds_alternative<double>(node->value)) return (int)std::get<double>(node->value); 
    }
    return 0; 
}

std::string SymbolicExpr::get_identifier() const {
    if (auto v = std::dynamic_pointer_cast<VariableNode>(root)) return v->name;
    return "";
}


