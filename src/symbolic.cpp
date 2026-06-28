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
#include "../include/assumption_context.hpp"
#include "../include/interval.hpp"
#include "../include/transcendental_factor.hpp"
#include "../include/solve_polynomial.hpp"
#include "../include/multivariate_factor.hpp"
#include "../include/calculus_utils.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"

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
    void visit(RelationalNode& node) override {
        node.left->accept(*this);
        node.right->accept(*this);
    }
    void visit(LogicalNode& node) override {
        node.left->accept(*this);
        node.right->accept(*this);
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

    void visit(RelationalNode& node) override {
        node.left->accept(*this);
        auto new_left = result;
        node.right->accept(*this);
        auto new_right = result;
        result = std::make_shared<RelationalNode>(new_left, new_right, node.op);
    }

    void visit(LogicalNode& node) override {
        node.left->accept(*this);
        auto new_left = result;
        node.right->accept(*this);
        auto new_right = result;
        result = std::make_shared<LogicalNode>(new_left, new_right, node.op);
    }

    void visit(SummationNode& node) override {
        // 指标变量是绑定变量：仅当替换变量不是指标时才替换通项。
        std::shared_ptr<SymbolicNode> new_body;
        if (node.index_var == var_name) {
            new_body = node.body->clone();
        } else {
            node.body->accept(*this);
            new_body = result;
        }
        node.lower_bound->accept(*this);
        auto new_lo = result;
        node.upper_bound->accept(*this);
        auto new_hi = result;
        result = std::make_shared<SummationNode>(new_body, node.index_var, new_lo, new_hi);
    }

    void visit(ProductNode_Op& node) override {
        std::shared_ptr<SymbolicNode> new_body;
        if (node.index_var == var_name) {
            new_body = node.body->clone();
        } else {
            node.body->accept(*this);
            new_body = result;
        }
        node.lower_bound->accept(*this);
        auto new_lo = result;
        node.upper_bound->accept(*this);
        auto new_hi = result;
        result = std::make_shared<ProductNode_Op>(new_body, node.index_var, new_lo, new_hi);
    }

    void visit(PiecewiseNode& node) override {
        std::vector<PiecewiseNode::Branch> new_branches;
        for (const auto& b : node.branches) {
            PiecewiseNode::Branch nb;
            b.expression->accept(*this);
            nb.expression = result;
            b.condition->accept(*this);
            nb.condition = result;
            new_branches.push_back(nb);
        }
        std::shared_ptr<SymbolicNode> new_def = nullptr;
        if (node.default_expr) {
            node.default_expr->accept(*this);
            new_def = result;
        }
        result = std::make_shared<PiecewiseNode>(std::move(new_branches), new_def);
    }

    void visit(ComplexNode& node) override {
        node.real->accept(*this);
        auto new_real = result;
        node.imag->accept(*this);
        auto new_imag = result;
        result = std::make_shared<ComplexNode>(new_real, new_imag);
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
            std::vector<std::shared_ptr<SymbolicNode>> args = {node.clone(), std::make_shared<VariableNode>(var)};
            result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args);
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
            std::vector<std::shared_ptr<SymbolicNode>> args = {node.clone(), std::make_shared<VariableNode>(var)};
            result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args);
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

        std::vector<std::shared_ptr<SymbolicNode>> args = {node.clone(), std::make_shared<VariableNode>(var)};
        result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args);
    }

    void visit(PowerNode& node) override {

        bool base_is_x = false;
        if (auto v = std::dynamic_pointer_cast<VariableNode>(node.base)) {
            if (v->name == var) base_is_x = true;
        }

        if (base_is_x) {

            if (auto num = std::dynamic_pointer_cast<NumberNode>(node.exponent)) {
                  if (std::holds_alternative<lmmc_real_t>(num->value)) {
                      int eq;
                      lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(num->value), -1.0, 1e-9, 1e-9, &eq);
                      if (eq) {
                          std::vector<std::shared_ptr<SymbolicNode>> args = {node.base};
                          result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, args);
                          return;
                      }
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

        std::vector<std::shared_ptr<SymbolicNode>> fallback_args = {node.clone(), std::make_shared<VariableNode>(var)};
        result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, fallback_args);
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

        std::vector<std::shared_ptr<SymbolicNode>> fallback_args = {node.clone(), std::make_shared<VariableNode>(var)};
        result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, fallback_args);
    }

    void visit(MatrixNode& node) override {

        std::vector<std::shared_ptr<SymbolicNode>> fallback_args = {node.clone(), std::make_shared<VariableNode>(var)};
        result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, fallback_args);
    }
};

int SymbolicExpr::compare(const std::shared_ptr<SymbolicExpr>& other) const {
    if (!root || !other->root) return 0;

    return root->compare(*other->root);
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::substitute(const std::string& var_name, const std::shared_ptr<SymbolicExpr>& value) const {
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

lmmc_real_t SymbolicExpr::to_numeric() const {
    if (!root) return 0.0;

    if (auto num = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value)) return std::get<lmmc_real_t>(num->value);
        if (std::holds_alternative<::BigInt>(num->value)) return (lmmc_real_t)std::get<::BigInt>(num->value).to_double();
        if (std::holds_alternative<::Rational>(num->value)) return (lmmc_real_t)std::get<::Rational>(num->value).to_double();
    }

    if (auto var = std::dynamic_pointer_cast<VariableNode>(root)) {
        if (var->name == "pi") return LMMC_CONST_PI;
        if (var->name == "e") return std::exp((lmmc_real_t)1.0);
        return 0.0;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(root)) {
        lmmc_real_t sum = 0.0;
        for (const auto& op : add->operands) {
            sum += SymbolicExpr(op).to_numeric();
        }
        return sum;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(root)) {
        lmmc_real_t prod = 1.0;
        for (const auto& op : mul->operands) {
            prod *= SymbolicExpr(op).to_numeric();
        }
        return prod;
    }

    if (auto pow = std::dynamic_pointer_cast<PowerNode>(root)) {
        lmmc_real_t base_val = SymbolicExpr(pow->base).to_numeric();
        lmmc_real_t exp_val = SymbolicExpr(pow->exponent).to_numeric();
        return std::pow(base_val, exp_val);
    }

    if (auto func = std::dynamic_pointer_cast<FunctionNode>(root)) {
        if (func->arguments.size() == 1) {
            SymbolicExpr arg_expr(func->arguments[0]);
            lmmc_real_t arg_val = arg_expr.to_numeric();
            lmmc_real_t result;

            switch (func->type) {
                case FunctionNode::FuncType::Sin:
                    result = std::sin(arg_val);
                    return result;
                case FunctionNode::FuncType::Cos:
                    result = std::cos(arg_val);
                    return result;
                case FunctionNode::FuncType::Tan:
                    result = std::tan(arg_val);
                    return result;
                case FunctionNode::FuncType::Sec:
                    return 1.0 / std::cos(arg_val);
                case FunctionNode::FuncType::Csc:
                    return 1.0 / std::sin(arg_val);
                case FunctionNode::FuncType::Cot:
                    return 1.0 / std::tan(arg_val);
                case FunctionNode::FuncType::Exp:
                    result = std::exp(arg_val);
                    return result;
                case FunctionNode::FuncType::Ln:
                    result = std::log(arg_val);
                    return result;
                case FunctionNode::FuncType::Sqrt:
                    result = std::sqrt(arg_val);
                    return result;
                case FunctionNode::FuncType::Abs:
                    result = std::abs(arg_val);
                    return result;
                case FunctionNode::FuncType::ArcSin:
                    return std::asin(arg_val);
                case FunctionNode::FuncType::ArcCos:
                    return std::acos(arg_val);
                case FunctionNode::FuncType::ArcTan:
                    return std::atan(arg_val);
                case FunctionNode::FuncType::Sinh:
                    return std::sinh(arg_val);
                case FunctionNode::FuncType::Cosh:
                    return std::cosh(arg_val);
                case FunctionNode::FuncType::Tanh:
                    return std::tanh(arg_val);
                case FunctionNode::FuncType::LambertW:
                    if (lmmc_lambertw(arg_val, &result) == LMMC_STATUS_OK) {
                        return result;
                    }
                    break;
                default:
                    break;
            }
        }

        if (func->arguments.size() == 2 && func->type == FunctionNode::FuncType::Atan2) {
            SymbolicExpr y_expr(func->arguments[0]);
            SymbolicExpr x_expr(func->arguments[1]);
            lmmc_real_t y_val = y_expr.to_numeric();
            lmmc_real_t x_val = x_expr.to_numeric();
            return std::atan2(y_val, x_val);
        }
    }

    return 0.0;
}

bool SymbolicExpr::is_zero() const {
    if (!root) return false;
    if (auto num = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            int eq;
            lmmc_double_nearly_equal(v, 0.0, &eq);
            return eq != 0;
        }
        if (std::holds_alternative<::BigInt>(num->value)) return std::get<::BigInt>(num->value).to_int() == 0;
        if (std::holds_alternative<::Rational>(num->value)) return std::get<::Rational>(num->value).get_numerator().to_int() == 0;
    }
    return false;
}

bool SymbolicExpr::is_one() const {
    if (!root) return false;
    if (auto num = std::dynamic_pointer_cast<NumberNode>(root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            int eq;
            lmmc_double_nearly_equal(v, 1.0, &eq);
            return eq != 0;
        }
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

        auto sinx = SymbolicExpr::sin(x);
        auto cosx = SymbolicExpr::cos(x);
        auto n2 = SymbolicExpr::number(2);
        auto sin2 = SymbolicExpr::power(sinx, n2);
        auto cos2 = SymbolicExpr::power(cosx, n2);

        auto pat1 = SymbolicExpr::add(sin2, cos2);
        engine.add_rule(Rule(*pat1, *SymbolicExpr::number(1), {"x"}));

        auto pat2 = SymbolicExpr::add(cos2, sin2);
        engine.add_rule(Rule(*pat2, *SymbolicExpr::number(1), {"x"}));

        auto two_x = SymbolicExpr::multiply(n2, x);
        auto sin2x = SymbolicExpr::sin(two_x);
        auto two_sin_cos = SymbolicExpr::multiply(n2,
            SymbolicExpr::multiply(sinx, cosx));
        engine.add_rule(Rule(*sin2x, *two_sin_cos, {"x"}));

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

// ---------------------------------------------------------------------------
// 多元因式分解辅助函数
// ---------------------------------------------------------------------------

/**
 * @internal
 * @brief 判断符号节点是否为多项式表达式（不含超越函数、负指数等）
 */
static bool is_poly_expr_node(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return false;
    if (std::dynamic_pointer_cast<NumberNode>(node)) return true;
    if (std::dynamic_pointer_cast<VariableNode>(node)) return true;

    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (const auto& op : add->operands) {
            if (!is_poly_expr_node(op)) return false;
        }
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (const auto& op : mul->operands) {
            if (!is_poly_expr_node(op)) return false;
        }
        return true;
    }

    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
            long long exp_val = 0;
            if (std::holds_alternative<BigInt>(exp_num->value)) {
                auto bi = std::get<BigInt>(exp_num->value);
                if (bi.IsNegative()) return false;
                exp_val = bi.to_int();
            } else if (std::holds_alternative<Rational>(exp_num->value)) {
                auto r = std::get<Rational>(exp_num->value);
                if (!r.is_integer()) return false;
                auto bi = r.to_BigInt();
                if (bi.IsNegative()) return false;
                exp_val = bi.to_int();
            } else if (std::holds_alternative<lmmc_real_t>(exp_num->value)) {
                lmmc_real_t d = std::get<lmmc_real_t>(exp_num->value);
                if (!std::isfinite(d) || d < 0 || d != std::floor(d)) return false;
                exp_val = static_cast<long long>(d);
            } else {
                return false;
            }
            if (exp_val < 0 || exp_val > 100) return false;
            return is_poly_expr_node(pow->base);
        }
        return false;
    }

    return false;
}

/**
 * @internal
 * @brief 递归将符号节点转换为 MultiPoly
 *
 * 假设节点已通过 is_poly_expr_node 验证为多项式。
 * @param[in] node 符号节点
 * @param[in] vars 变量名列表（确定单项式各分量的含义）
 * @return 对应的 MultiPoly
 */
static lamina::MultiPoly symbolic_node_to_multipoly(
    const std::shared_ptr<SymbolicNode>& node,
    const std::vector<std::string>& vars)
{
    if (!node) return lamina::MultiPoly(Rational(0), vars);

    if (auto num = std::dynamic_pointer_cast<NumberNode>(node)) {
        Rational coeff(0);
        if (std::holds_alternative<BigInt>(num->value)) {
            coeff = Rational(std::get<BigInt>(num->value));
        } else if (std::holds_alternative<Rational>(num->value)) {
            coeff = std::get<Rational>(num->value);
        } else if (std::holds_alternative<lmmc_real_t>(num->value)) {
            coeff = Rational::from_double(std::get<lmmc_real_t>(num->value));
        }
        return lamina::MultiPoly(coeff, vars);
    }

    if (auto var_node = std::dynamic_pointer_cast<VariableNode>(node)) {
        lamina::Monomial mono(vars.size(), 0);
        for (size_t i = 0; i < vars.size(); ++i) {
            if (vars[i] == var_node->name) {
                mono[i] = 1;
                break;
            }
        }
        std::vector<lamina::MultiPoly::Term> terms;
        terms.push_back({mono, Rational(1)});
        return lamina::MultiPoly(std::move(terms), vars);
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        lamina::MultiPoly result(Rational(0), vars);
        for (const auto& op : add->operands) {
            result = result + symbolic_node_to_multipoly(op, vars);
        }
        return result;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        lamina::MultiPoly result(Rational(1), vars);
        for (const auto& op : mul->operands) {
            result = result * symbolic_node_to_multipoly(op, vars);
        }
        return result;
    }

    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        auto base_poly = symbolic_node_to_multipoly(pow->base, vars);
        int exp_val = 0;
        if (auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
            if (std::holds_alternative<BigInt>(exp_num->value)) {
                exp_val = std::get<BigInt>(exp_num->value).to_int();
            } else if (std::holds_alternative<Rational>(exp_num->value)) {
                exp_val = std::get<Rational>(exp_num->value).to_BigInt().to_int();
            } else if (std::holds_alternative<lmmc_real_t>(exp_num->value)) {
                exp_val = static_cast<int>(std::get<lmmc_real_t>(exp_num->value));
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

        // 系数部分
        if (!(coeff == Rational(1)) || lamina::total_degree(mono) == 0) {
            if (coeff == Rational(-1) && lamina::total_degree(mono) > 0) {
                factors.push_back(SymbolicExpr::number(-1));
            } else {
                factors.push_back(SymbolicExpr::number(coeff));
            }
        }

        // 变量部分
        for (size_t i = 0; i < vars.size() && i < mono.size(); ++i) {
            if (mono[i] == 0) continue;
            auto var_expr = std::make_shared<SymbolicExpr>(
                std::make_shared<VariableNode>(vars[i]));
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

std::shared_ptr<SymbolicExpr> SymbolicExpr::factor() const {
    auto simp = simplify();
    if (!simp || !simp->root) return simp;

    if (auto add_node = std::dynamic_pointer_cast<AddNode>(simp->root)) {

        // 步骤 1：提取各项的公因式（GCD）
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

                  auto inv_common = power(common, number(-1));
                  auto quot = multiply(term, inv_common);
                  quot = quot->simplify();
                  new_ops.push_back(quot->root);
             }
             auto new_sum = std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(new_ops));

             // 递归分解余下的和式
             auto factored_sum = new_sum->factor();

             std::vector<std::shared_ptr<SymbolicNode>> final_ops = {common->root, factored_sum->root};
             return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(final_ops));
        }

        // 步骤 2：一元多项式分解（支持任意次数）
        VariablesVisitor vv;
        simp->root->accept(vv);
        if (vv.vars.size() == 1) {
             std::string var = *vv.vars.begin();
             try {
                 auto poly = lamina::symbolic_to_poly<Rational>(simp, var);
                 int deg = poly.degree();

                 if (deg >= 2) {
                      // 使用有理根定理逐步分解
                      auto roots = lamina::find_rational_roots(poly);

                      if (!roots.empty()) {
                           auto x_node = std::make_shared<SymbolicExpr>(std::make_shared<VariableNode>(var));
                           auto leading = poly.lead_coeff();

                           std::vector<std::shared_ptr<SymbolicNode>> factors;

                           // 首项系数
                           if (!(leading == Rational(1))) {
                               factors.push_back(number(leading)->root);
                           }

                           // 从根构建线性因子 (x - r)
                           for (const auto& r : roots) {
                               std::shared_ptr<SymbolicExpr> linear_factor;
                               if (r == Rational(0)) {
                                   linear_factor = x_node;
                               } else {
                                   auto neg_r = number(Rational(0) - r);
                                   linear_factor = SymbolicExpr::add(x_node, neg_r)->simplify();
                               }
                               factors.push_back(linear_factor->root);
                           }

                           // 计算余下的不可约因子：原多项式 / 已分解因子的乘积
                           lamina::Polynomial<Rational> factored_product({leading}, var);
                           for (const auto& r : roots) {
                               lamina::Polynomial<Rational> lin({Rational(0) - r, Rational(1)}, var);
                               factored_product = factored_product * lin;
                           }

                           auto [quotient, remainder] = poly.div_mod(factored_product);

                           if (remainder.is_zero() && quotient.degree() >= 1) {
                               // 余下的不可约因子
                               if (quotient.degree() == 1 && quotient.lead_coeff() == Rational(1)) {
                                   // 一次因子直接加入
                                   auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                   factors.push_back(q_expr->root);
                               } else if (quotient.degree() >= 2) {
                                   // 尝试递归分解余下部分
                                   auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                   auto q_factored = q_expr->factor();
                                   factors.push_back(q_factored->root);
                               } else {
                                   auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                   factors.push_back(q_expr->root);
                               }
                           } else if (!remainder.is_zero()) {
                               // 除法有余数 → 回退到原始方法
                               goto try_solve_quadratic;
                           }

                           if (factors.size() > 1) {
                               return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(factors));
                           } else if (factors.size() == 1) {
                               return std::make_shared<SymbolicExpr>(factors[0]);
                           }
                      }
                 }
             } catch (...) {}

             try_solve_quadratic:
             // 后备：二次多项式通过求解方程分解
             try {
                 auto poly = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(simp, *vv.vars.begin());
                 if (poly.degree() == 2) {
                      std::string var = *vv.vars.begin();
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

        // 步骤 3：多元多项式 — 先尝试 factor_multivariate，再回退逐变量分解
        if (vv.vars.size() > 1) {
            // 3a: 检测是否为多项式表达式，若是则使用 MultiPoly 路径
            if (is_poly_expr_node(simp->root)) {
                try {
                    std::vector<std::string> var_list(vv.vars.begin(), vv.vars.end());
                    auto mpoly = symbolic_node_to_multipoly(simp->root, var_list);

                    if (!mpoly.is_zero() && !mpoly.is_constant()) {
                        auto result = lamina::factor_multivariate(mpoly);

                        // 检查分解是否产生了多个因子
                        if (!result.factors.empty()) {
                            std::vector<std::shared_ptr<SymbolicNode>> factor_nodes;

                            // 常数因子
                            if (!(result.constant == Rational(1))) {
                                factor_nodes.push_back(
                                    SymbolicExpr::number(result.constant)->root);
                            }

                            // 各不可约因子
                            for (size_t i = 0; i < result.factors.size(); ++i) {
                                auto factor_expr = multipoly_to_symbolic(result.factors[i]);
                                if (!factor_expr || factor_expr->is_zero()) continue;

                                int mult = (i < result.multiplicities.size())
                                    ? result.multiplicities[i] : 1;

                                if (mult == 1) {
                                    factor_nodes.push_back(factor_expr->root);
                                } else {
                                    auto pow_expr = SymbolicExpr::power(
                                        factor_expr, SymbolicExpr::number(mult));
                                    factor_nodes.push_back(pow_expr->root);
                                }
                            }

                            if (factor_nodes.size() > 1) {
                                return std::make_shared<SymbolicExpr>(
                                    std::make_shared<MultiplyNode>(std::move(factor_nodes)));
                            } else if (factor_nodes.size() == 1) {
                                return std::make_shared<SymbolicExpr>(factor_nodes[0]);
                            }
                        }
                    }
                } catch (...) {
                    // MultiPoly 路径失败，回退到逐变量分解
                }
            }

            // 3b: 回退 — 逐变量尝试有理根分解
            for (const auto& var : vv.vars) {
                try {
                    auto poly = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(simp, var);
                    if (poly.degree() >= 2) {
                        // 尝试用 Rational 系数做分解
                        auto poly_r = lamina::symbolic_to_poly<Rational>(simp, var);
                        auto roots = lamina::find_rational_roots(poly_r);
                        if (!roots.empty()) {
                            auto x_node = std::make_shared<SymbolicExpr>(std::make_shared<VariableNode>(var));
                            std::vector<std::shared_ptr<SymbolicNode>> factors;

                            auto leading = poly_r.lead_coeff();
                            if (!(leading == Rational(1))) {
                                factors.push_back(number(leading)->root);
                            }

                            for (const auto& r : roots) {
                                std::shared_ptr<SymbolicExpr> linear_factor;
                                if (r == Rational(0)) {
                                    linear_factor = x_node;
                                } else {
                                    auto neg_r = number(Rational(0) - r);
                                    linear_factor = SymbolicExpr::add(x_node, neg_r)->simplify();
                                }
                                factors.push_back(linear_factor->root);
                            }

                            // 计算余下因子
                            lamina::Polynomial<Rational> factored_product({leading}, var);
                            for (const auto& r : roots) {
                                lamina::Polynomial<Rational> lin({Rational(0) - r, Rational(1)}, var);
                                factored_product = factored_product * lin;
                            }
                            auto [quotient, remainder] = poly_r.div_mod(factored_product);
                            if (remainder.is_zero() && quotient.degree() >= 1) {
                                auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                factors.push_back(q_expr->root);
                            } else if (remainder.is_zero() && !quotient.is_zero() && quotient.degree() == 0) {
                                // 常数商 → 已完全分解
                                if (!(quotient.coeffs[0] == Rational(1))) {
                                    factors.push_back(number(quotient.coeffs[0])->root);
                                }
                            }

                            if (factors.size() > 1) {
                                return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(factors));
                            }
                        }
                    }
                } catch (...) {
                    continue;
                }
            }
        }
    }

    // 超越因式分解后备路径：当标准多项式分解无法处理时，
    // 检测表达式是否含超越函数，若是则尝试超越因式分解。
    {
        VariablesVisitor tv;
        simp->root->accept(tv);
        std::string target_var;
        if (!tv.vars.empty()) {
            target_var = *tv.vars.begin();
        }

        if (!target_var.empty()) {
            try {
                auto trans_factors = lamina::factor_transcendental(simp, target_var);
                if (trans_factors.size() > 1) {
                    // 超越分解成功：组装乘积表达式
                    std::vector<std::shared_ptr<SymbolicNode>> factor_nodes;
                    factor_nodes.reserve(trans_factors.size());
                    for (const auto& f : trans_factors) {
                        if (f && f->root) {
                            factor_nodes.push_back(f->root);
                        }
                    }
                    if (factor_nodes.size() > 1) {
                        return std::make_shared<SymbolicExpr>(
                            std::make_shared<MultiplyNode>(std::move(factor_nodes)));
                    }
                }
            } catch (...) {
                // 超越分解失败，返回原表达式
            }
        }
    }

    return simp;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::cancel() const {
    auto simp = simplify();
    if (!simp || !simp->root) return simp;

    // 辅助 lambda：从乘积节点中分离分子因子和分母因子。
    // 分母因子 = 含负指数的 PowerNode。
    auto separate_num_den = [](const std::shared_ptr<SymbolicNode>& node,
                               std::vector<std::shared_ptr<SymbolicNode>>& num_out,
                               std::vector<std::shared_ptr<SymbolicNode>>& den_out) {
        auto classify = [&](const std::shared_ptr<SymbolicNode>& factor) {
            if (auto pow = std::dynamic_pointer_cast<PowerNode>(factor)) {
                if (auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                    bool is_negative = false;
                    if (std::holds_alternative<BigInt>(exp_num->value)) {
                        is_negative = std::get<BigInt>(exp_num->value).IsNegative();
                    } else if (std::holds_alternative<lmmc_real_t>(exp_num->value)) {
                        is_negative = std::get<lmmc_real_t>(exp_num->value) < 0;
                    } else if (std::holds_alternative<Rational>(exp_num->value)) {
                        is_negative = std::get<Rational>(exp_num->value).get_numerator().IsNegative();
                    }
                    if (is_negative) {
                        // 取绝对值指数
                        std::shared_ptr<SymbolicNode> pos_exp;
                        if (std::holds_alternative<BigInt>(exp_num->value)) {
                            pos_exp = std::make_shared<NumberNode>(BigInt(0) - std::get<BigInt>(exp_num->value));
                        } else if (std::holds_alternative<lmmc_real_t>(exp_num->value)) {
                            pos_exp = std::make_shared<NumberNode>(-std::get<lmmc_real_t>(exp_num->value));
                        } else {
                            auto r = std::get<Rational>(exp_num->value);
                            pos_exp = std::make_shared<NumberNode>(Rational(BigInt(0) - r.get_numerator(), r.get_denominator()));
                        }
                        auto pos_exp_num = std::dynamic_pointer_cast<NumberNode>(pos_exp);
                        if (pos_exp_num && pos_exp_num->is_one()) {
                            den_out.push_back(pow->base);
                        } else {
                            den_out.push_back(std::make_shared<PowerNode>(pow->base, pos_exp));
                        }
                        return;
                    }
                }
            }
            num_out.push_back(factor);
        };

        if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
            for (const auto& op : mul->operands) {
                classify(op);
            }
        } else {
            classify(node);
        }
    };

    // 辅助 lambda：从因子列表构建乘积表达式
    auto build_product = [](const std::vector<std::shared_ptr<SymbolicNode>>& factors) -> std::shared_ptr<SymbolicExpr> {
        if (factors.empty()) return SymbolicExpr::number(1);
        if (factors.size() == 1) return std::make_shared<SymbolicExpr>(factors[0]);
        return std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(factors));
    };

    // 策略 1：表达式是 AddNode，各项可能含公共分母因子。
    // 例如 simplify 后 (x²-1)/(x-1) 变为 -1*(x-1)^-1 + x²*(x-1)^-1
    // 需要提取公共分母，重组为 (分子之和)/分母 再做 GCD 约分。
    if (auto add_node = std::dynamic_pointer_cast<AddNode>(simp->root)) {
        // 对每个加法项分离分子/分母
        struct TermInfo {
            std::shared_ptr<SymbolicExpr> numerator;
            std::shared_ptr<SymbolicExpr> denominator;
        };
        std::vector<TermInfo> terms;
        bool has_denominator = false;

        for (const auto& term : add_node->operands) {
            std::vector<std::shared_ptr<SymbolicNode>> t_num, t_den;
            separate_num_den(term, t_num, t_den);
            auto num_expr = build_product(t_num)->simplify();
            auto den_expr = build_product(t_den)->simplify();
            if (!den_expr->is_one()) has_denominator = true;
            terms.push_back({num_expr, den_expr});
        }

        if (has_denominator) {
            // 计算公共分母（所有项分母的 LCM，简化处理：乘积）
            // 对于常见情况（所有项分母相同），直接提取。
            // 检查是否所有分母相同
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
                // 所有项分母相同：分子直接相加
                combined_den = terms[0].denominator;
                std::vector<std::shared_ptr<SymbolicNode>> num_ops;
                for (const auto& t : terms) {
                    if (t.numerator && t.numerator->root) {
                        num_ops.push_back(t.numerator->root);
                    }
                }
                if (num_ops.empty()) return SymbolicExpr::number(0);
                if (num_ops.size() == 1) combined_num = std::make_shared<SymbolicExpr>(num_ops[0]);
                else combined_num = std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(num_ops));
                combined_num = combined_num->simplify();
            } else {
                // 分母不同：通分（乘以其他项的分母）
                // 计算总分母 = 所有不同分母的乘积
                std::shared_ptr<SymbolicExpr> total_den = SymbolicExpr::number(1);
                // 收集不同的分母
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

                // 每项分子乘以 (总分母 / 该项分母)
                std::vector<std::shared_ptr<SymbolicNode>> num_ops;
                for (const auto& t : terms) {
                    auto factor = divide(total_den, t.denominator)->simplify();
                    auto adjusted_num = SymbolicExpr::multiply(t.numerator, factor)->simplify();
                    if (adjusted_num && adjusted_num->root) {
                        num_ops.push_back(adjusted_num->root);
                    }
                }
                if (num_ops.empty()) return SymbolicExpr::number(0);
                if (num_ops.size() == 1) combined_num = std::make_shared<SymbolicExpr>(num_ops[0]);
                else combined_num = std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(num_ops));
                combined_num = combined_num->expand()->simplify();
            }

            // 对 combined_num / combined_den 做多项式 GCD 约分
            combined_den = combined_den->expand()->simplify();

            VariablesVisitor vv;
            if (combined_num->root) combined_num->root->accept(vv);
            if (combined_den->root) combined_den->root->accept(vv);

            auto cur_num = combined_num;
            auto cur_den = combined_den;

            for (const auto& var : vv.vars) {
                try {
                    auto num_expanded = cur_num->expand()->simplify();
                    auto den_expanded = cur_den->expand()->simplify();

                    // 使用 SymbolicPolyCoeff 保留其他变量作为符号系数
                    auto poly_num = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(num_expanded, var);
                    auto poly_den = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(den_expanded, var);

                    if (poly_num.is_zero()) return SymbolicExpr::number(0);
                    if (poly_den.is_zero()) return simp;

                    // 尝试用 Rational 系数做 GCD（纯数值系数时更可靠）
                    auto poly_num_r = lamina::symbolic_to_poly<Rational>(num_expanded, var);
                    auto poly_den_r = lamina::symbolic_to_poly<Rational>(den_expanded, var);

                    // 检查 Rational 转换是否丢失了信息（多元情况）
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
                        // 多元情况：使用 SymbolicPolyCoeff 做 GCD
                        auto g = lamina::Polynomial<lamina::SymbolicPolyCoeff>::gcd(poly_num, poly_den);
                        if (g.degree() >= 1) {
                            auto [q_num, r_num] = poly_num.div_mod(g);
                            auto [q_den, r_den] = poly_den.div_mod(g);
                            cur_num = lamina::poly_to_symbolic(q_num)->simplify();
                            cur_den = lamina::poly_to_symbolic(q_den)->simplify();
                        }
                    }
                } catch (...) {
                    continue;
                }
            }

            if (cur_den->is_one()) return cur_num;

            // 分母为 -1 时取负
            if (cur_den->root) {
                auto diff = SymbolicExpr::add(std::make_shared<SymbolicExpr>(cur_den->root), SymbolicExpr::number(-1))->simplify();
                if (diff->is_zero()) {
                    return SymbolicExpr::multiply(SymbolicExpr::number(-1), cur_num)->simplify();
                }
            }

            return divide(cur_num, cur_den)->simplify();
        }
    }

    // 策略 2：表达式是 MultiplyNode，直接分离分子/分母。
    std::vector<std::shared_ptr<SymbolicNode>> num_factors;
    std::vector<std::shared_ptr<SymbolicNode>> den_factors;
    separate_num_den(simp->root, num_factors, den_factors);

    if (den_factors.empty()) return simp;

    auto numerator = build_product(num_factors)->simplify();
    auto denominator = build_product(den_factors)->simplify();

    if (denominator->is_one()) return numerator;

    VariablesVisitor vv;
    if (numerator->root) numerator->root->accept(vv);
    if (denominator->root) denominator->root->accept(vv);

    if (vv.vars.empty()) {
        return divide(numerator, denominator)->simplify();
    }

    auto cur_num = numerator;
    auto cur_den = denominator;

    for (const auto& var : vv.vars) {
        try {
            auto num_expanded = cur_num->expand()->simplify();
            auto den_expanded = cur_den->expand()->simplify();

            auto poly_num_r = lamina::symbolic_to_poly<Rational>(num_expanded, var);
            auto poly_den_r = lamina::symbolic_to_poly<Rational>(den_expanded, var);

            // 检查 Rational 转换是否丢失了信息
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
        } catch (...) {
            continue;
        }
    }

    if (cur_den->is_one()) return cur_num;

    if (cur_den->root) {
        auto diff = SymbolicExpr::add(std::make_shared<SymbolicExpr>(cur_den->root), SymbolicExpr::number(-1))->simplify();
        if (diff->is_zero()) {
            return SymbolicExpr::multiply(SymbolicExpr::number(-1), cur_num)->simplify();
        }
    }

    return divide(cur_num, cur_den)->simplify();
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

    lamina::Integrator integrator;
    SymbolicExpr result = integrator.integrate(*this, var);

    auto res_ptr = std::make_shared<SymbolicExpr>(result);
    return res_ptr->simplify();
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::series(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, int order, const lamina::AssumptionContext* ctx) const {
    if (!root || !point) return nullptr;
    if (order < 0) return SymbolicExpr::number(0);
    if (ctx) {
    }

    auto x = SymbolicExpr::variable(var);
    auto neg_point = SymbolicExpr::multiply(SymbolicExpr::number(-1), point);
    auto delta = SymbolicExpr::add(x, neg_point);

    std::vector<std::shared_ptr<SymbolicNode>> terms;
    auto deriv = std::make_shared<SymbolicExpr>(root->clone());

    for (int n = 0; n <= order; ++n) {
        if (n > 0) {
            deriv = deriv->differentiate(var);
            if (!deriv) break;
            deriv = deriv->simplify();
        }

        auto coeff = deriv->substitute(var, point);
        if (!coeff) break;
        coeff = coeff->simplify();

        auto term = coeff;
        if (n > 0) {
            term = SymbolicExpr::multiply(term, SymbolicExpr::power(delta, SymbolicExpr::number(n)));
        }

        if (n > 1) {
            BigInt fact = BigInt::factorial(static_cast<unsigned int>(n));
            auto inv_fact = SymbolicExpr::number(Rational(BigInt(1), fact));
            term = SymbolicExpr::multiply(term, inv_fact);
        }

        terms.push_back(term->root);
    }

    if (terms.empty()) return SymbolicExpr::number(0);
    auto sum = std::make_shared<AddNode>(terms);
    return std::make_shared<SymbolicExpr>(sum)->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den) {
    return SymbolicExpr::multiply(num, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::make_integral(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr) return nullptr;
    auto var_node = std::make_shared<VariableNode>(var);
    std::vector<std::shared_ptr<SymbolicNode>> args = {expr->root, var_node};
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args));
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::make_limit(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var, const std::shared_ptr<SymbolicExpr>& point) {
    if (!expr || !point) return nullptr;
    auto var_node = std::make_shared<VariableNode>(var);
    std::vector<std::shared_ptr<SymbolicNode>> args = {expr->root, var_node, point->root};
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(FunctionNode::FuncType::Limit, args));
}

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
            void visit(LogicalNode& n) override { n.left->accept(*this); n.right->accept(*this); }
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
        auto pa = lamina::symbolic_to_poly<Rational>(a, var);
        auto pb = lamina::symbolic_to_poly<Rational>(b, var);

        int m = pa.degree();
        int n = pb.degree();

        // 退化情形：任一为零多项式 → 结式为 0
        if (m < 0 || n < 0) return SymbolicExpr::number(0);
        // 任一为非零常数：Res(a, const c) = c^deg(other)（两者都常数时为 1）
        if (m == 0 && n == 0) return SymbolicExpr::number(1);
        if (n == 0) {
            // Res(a, c) = c^m
            Rational c = pb.coeffs[0];
            Rational r(1);
            for (int i = 0; i < m; ++i) r = r * c;
            return SymbolicExpr::number(r);
        }
        if (m == 0) {
            Rational c = pa.coeffs[0];
            Rational r(1);
            for (int i = 0; i < n; ++i) r = r * c;
            return SymbolicExpr::number(r);
        }

        // 构造 (m+n)×(m+n) Sylvester 矩阵。
        // 前 n 行由 a 的系数移位排列，后 m 行由 b 的系数移位排列。
        // 行内按降幂排列：列 0 对应 x^(m+n-1) 系数。
        int N = m + n;
        std::vector<std::vector<Rational>> S(N, std::vector<Rational>(N, Rational(0)));

        // a 的降幂系数：a_m, a_{m-1}, ..., a_0  （pa.coeffs[m..0]）
        for (int row = 0; row < n; ++row) {
            for (int k = 0; k <= m; ++k) {
                // a 的 x^(m-k) 系数 = pa.coeffs[m-k]
                S[row][row + k] = pa.coeffs[m - k];
            }
        }
        // b 的降幂系数
        for (int row = 0; row < m; ++row) {
            for (int k = 0; k <= n; ++k) {
                S[n + row][row + k] = pb.coeffs[n - k];
            }
        }

        // Bareiss 无分式高斯消元求行列式（精确 Rational 运算）。
        Rational det(1);
        Rational prev(1);
        int sign = 1;
        for (int i = 0; i < N; ++i) {
            // 主元为 0 时寻找下方非零行交换
            if (S[i][i] == Rational(0)) {
                int swap_row = -1;
                for (int r = i + 1; r < N; ++r) {
                    if (!(S[r][i] == Rational(0))) { swap_row = r; break; }
                }
                if (swap_row < 0) return SymbolicExpr::number(0); // 奇异 → 结式 0
                std::swap(S[i], S[swap_row]);
                sign = -sign;
            }
            for (int r = i + 1; r < N; ++r) {
                for (int c = i + 1; c < N; ++c) {
                    Rational num = S[r][c] * S[i][i] - S[r][i] * S[i][c];
                    S[r][c] = num / prev;
                }
                S[r][i] = Rational(0);
            }
            prev = S[i][i];
        }
        det = S[N - 1][N - 1];
        if (sign < 0) det = Rational(0) - det;

        return SymbolicExpr::number(det);
    } catch (...) {}

    return SymbolicExpr::number(0);
}



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
            case FunctionNode::FuncType::LambertW: return Type::Variable;
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
        if (std::holds_alternative<lmmc_real_t>(node->value)) return (int)std::get<lmmc_real_t>(node->value);
    }
    return 0;
}

std::string SymbolicExpr::get_identifier() const {
    if (auto v = std::dynamic_pointer_cast<VariableNode>(root)) return v->name;
    return "";
}
