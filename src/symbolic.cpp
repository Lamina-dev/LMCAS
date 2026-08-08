#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include "../include/symbolic.hpp"
#include "../include/symbolic_ast.hpp"
#include "../include/integration.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include "../include/visitors/differentiation_visitor.hpp"
#include "../include/visitors/expand_visitor.hpp"
#include "../include/visitors/limit_visitor.hpp"
#include "../include/poly_utils.hpp"
#include "poly_utils_internal.hpp"
#include "../include/matcher.hpp"
#include "../include/integration.hpp"
#include "../include/assumption_context.hpp"
#include "../include/interval.hpp"
#include "../include/transcendental_factor.hpp"
#include "../include/solve_polynomial.hpp"
#include "../include/multivariate_factor.hpp"
#include "../include/calculus_utils.hpp"
#include "../include/numeric_evaluation.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"

std::shared_ptr<SymbolicExpr> SymbolicExpr::number(int n) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(BigInt(n)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::number(long long n) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(BigInt(n)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::number(double n) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(n)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::number(const BigInt& value) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(value));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::number(const Rational& value) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(value));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::infinity(int sign) {
    auto node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<const SymbolicNode>>{});
    auto expression = lamina::detail::make_expression_ptr(node);
    return sign < 0 ? multiply(number(-1), expression) : expression;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::sqrt(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("sqrt operand cannot be null");
    }
    auto half = lamina::detail::make_node<NumberNode>(Rational(1, 2));
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<PowerNode>(operand->impl_->root, half));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::multiply(
    std::shared_ptr<SymbolicExpr> left,
    std::shared_ptr<SymbolicExpr> right) {
    if (!left || !left->impl_->root || !right || !right->impl_->root) {
        throw std::invalid_argument("multiply operands cannot be null");
    }
    std::vector<std::shared_ptr<const SymbolicNode>> operands{
        left->impl_->root, right->impl_->root};
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<MultiplyNode>(std::move(operands)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::add(
    std::shared_ptr<SymbolicExpr> left,
    std::shared_ptr<SymbolicExpr> right) {
    if (!left || !left->impl_->root || !right || !right->impl_->root) {
        throw std::invalid_argument("add operands cannot be null");
    }
    std::vector<std::shared_ptr<const SymbolicNode>> operands{
        left->impl_->root, right->impl_->root};
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<AddNode>(std::move(operands)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::power(
    std::shared_ptr<SymbolicExpr> base,
    std::shared_ptr<SymbolicExpr> exponent) {
    if (!base || !base->impl_->root || !exponent || !exponent->impl_->root) {
        throw std::invalid_argument("power operands cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<PowerNode>(base->impl_->root, exponent->impl_->root));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::sin(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("sin operand cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::cos(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("cos operand cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Cos,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::tan(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("tan operand cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Tan,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::ln(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("ln operand cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Ln,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::exp(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("exp operand cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::lambertw(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("lambertw operand cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::LambertW,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::log(
    std::shared_ptr<SymbolicExpr> value,
    std::shared_ptr<SymbolicExpr> base) {
    if (!value || !value->impl_->root || !base || !base->impl_->root) {
        throw std::invalid_argument("log operands cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Log,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                value->impl_->root, base->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::atan2(
    std::shared_ptr<SymbolicExpr> y,
    std::shared_ptr<SymbolicExpr> x) {
    if (!y || !y->impl_->root || !x || !x->impl_->root) {
        throw std::invalid_argument("atan2 operands cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Atan2,
            std::vector<std::shared_ptr<const SymbolicNode>>{y->impl_->root, x->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::root_of(
    std::shared_ptr<SymbolicExpr> polynomial,
    const std::string& variable_name,
    int index) {
    if (!polynomial || !polynomial->impl_->root) {
        throw std::invalid_argument("root_of polynomial cannot be null");
    }
    auto variable_expression = variable(variable_name);
    auto index_expression = number(index);
    std::vector<std::shared_ptr<const SymbolicNode>> arguments{
        polynomial->impl_->root, variable_expression->impl_->root, index_expression->impl_->root};
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::RootOf, std::move(arguments)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::eq(
    std::shared_ptr<SymbolicExpr> left,
    std::shared_ptr<SymbolicExpr> right) {
    if (!left || !left->impl_->root || !right || !right->impl_->root) {
        throw std::invalid_argument("eq operands cannot be null");
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<RelationalNode>(
            left->impl_->root, right->impl_->root, RelationalNode::Op::EQ));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::integral(
    std::shared_ptr<SymbolicExpr> operand,
    const std::string& variable_name) {
    return operand ? operand->integrate(variable_name) : nullptr;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::limit_func(
    std::shared_ptr<SymbolicExpr> operand,
    const std::string& variable_name,
    std::shared_ptr<SymbolicExpr> target) {
    return operand ? operand->limit(variable_name, std::move(target)) : nullptr;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::matrix(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& elements) {
    if (elements.empty()) {
        throw std::invalid_argument("matrix requires at least one row");
    }
    const std::size_t column_count = elements.front().size();
    if (column_count == 0) {
        throw std::invalid_argument("matrix requires at least one column");
    }

    std::vector<std::vector<std::shared_ptr<const SymbolicNode>>> node_elements;
    node_elements.reserve(elements.size());
    for (const auto& row : elements) {
        if (row.size() != column_count) {
            throw std::invalid_argument(
                "matrix rows must have the same number of columns");
        }
        std::vector<std::shared_ptr<const SymbolicNode>> node_row;
        node_row.reserve(row.size());
        for (const auto& element : row) {
            if (!element || !element->impl_->root) {
                throw std::invalid_argument("matrix elements cannot be null");
            }
            node_row.push_back(element->impl_->root);
        }
        node_elements.push_back(std::move(node_row));
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<MatrixNode>(std::move(node_elements)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::variable(
    const std::string& name) {
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<VariableNode>(name));
}

bool SymbolicExpr::is_number() const {
    return impl_->root->is_number();
}

bool SymbolicExpr::get_number_value_is_zero() const {
    return impl_->root->is_zero();
}

bool SymbolicExpr::is_big_int() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    return number_node &&
        std::holds_alternative<BigInt>(number_node->value());
}

bool SymbolicExpr::is_rational() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    return number_node &&
        std::holds_alternative<Rational>(number_node->value());
}

bool SymbolicExpr::is_int() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (!number_node) return false;
    if (std::holds_alternative<BigInt>(number_node->value())) return true;
    if (std::holds_alternative<Rational>(number_node->value())) {
        return std::get<Rational>(number_node->value()).is_integer();
    }

    const lmmc_real_t value = std::get<lmmc_real_t>(number_node->value());
    const lmmc_real_t rounded = std::round(value);
    int equal = 0;
    lmmc_double_nearly_equal_tol(value, rounded, 1e-12, 1e-12, &equal);
    return equal != 0;
}

std::variant<int, BigInt, Rational> SymbolicExpr::get_number() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (!number_node) {
        throw std::runtime_error("Expression is not a number");
    }
    if (std::holds_alternative<BigInt>(number_node->value())) {
        return std::get<BigInt>(number_node->value());
    }
    if (std::holds_alternative<Rational>(number_node->value())) {
        return std::get<Rational>(number_node->value());
    }
    return static_cast<int>(std::get<lmmc_real_t>(number_node->value()));
}

int SymbolicExpr::get_int() const {
    if (!is_int()) {
        throw std::runtime_error("Expression is not an integer");
    }
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (std::holds_alternative<BigInt>(number_node->value())) {
        return std::get<BigInt>(number_node->value()).to_int();
    }
    if (std::holds_alternative<Rational>(number_node->value())) {
        return std::get<Rational>(number_node->value()).to_BigInt().to_int();
    }
    return static_cast<int>(std::get<lmmc_real_t>(number_node->value()));
}

BigInt SymbolicExpr::get_big_int() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (!number_node ||
        !std::holds_alternative<BigInt>(number_node->value())) {
        throw std::runtime_error("Expression is not a BigInt");
    }
    return std::get<BigInt>(number_node->value());
}

Rational SymbolicExpr::get_rational() const {
    auto number_node = std::dynamic_pointer_cast<const NumberNode>(impl_->root);
    if (!number_node ||
        !std::holds_alternative<Rational>(number_node->value())) {
        throw std::runtime_error("Expression is not a Rational");
    }
    return std::get<Rational>(number_node->value());
}

Rational SymbolicExpr::convert_rational() const {
    if (is_rational()) return get_rational();
    if (is_big_int()) return Rational(get_big_int());
    return Rational(0);
}

class VariablesVisitor : public lamina::detail::SymbolicVisitor {
    std::set<std::string> bound_vars;

    bool is_bound(const std::string& name) const {
        return bound_vars.find(name) != bound_vars.end();
    }

public:
    std::set<std::string> vars;

    void visit(const NumberNode&) override {}
    void visit(const VariableNode& node) override {
        if (!is_bound(node.name())) {
            vars.insert(node.name());
        }
    }
    void visit(const AddNode& node) override {
        for(auto& op : node.operands()) op->accept(*this);
    }
    void visit(const MultiplyNode& node) override {
        for(auto& op : node.operands()) op->accept(*this);
    }
    void visit(const PowerNode& node) override {
        node.base()->accept(*this);
        node.exponent()->accept(*this);
    }
    void visit(const FunctionNode& node) override {
        for(auto& arg : node.arguments()) arg->accept(*this);
    }
    void visit(const MatrixNode&) override {

    }
    void visit(const RelationalNode& node) override {
        node.left()->accept(*this);
        node.right()->accept(*this);
    }
    void visit(const LogicalNode& node) override {
        node.left()->accept(*this);
        if (node.right()) node.right()->accept(*this);
    }
    void visit(const PiecewiseNode& node) override {
        for (const auto& branch : node.branches()) {
            branch.expression->accept(*this);
            branch.condition->accept(*this);
        }
        if (node.default_expr()) node.default_expr()->accept(*this);
    }
    void visit(const SummationNode& node) override {
        node.lower_bound()->accept(*this);
        node.upper_bound()->accept(*this);
        bool inserted = bound_vars.insert(node.index_var()).second;
        node.body()->accept(*this);
        if (inserted) bound_vars.erase(node.index_var());
    }
    void visit(const ProductNode_Op& node) override {
        node.lower_bound()->accept(*this);
        node.upper_bound()->accept(*this);
        bool inserted = bound_vars.insert(node.index_var()).second;
        node.body()->accept(*this);
        if (inserted) bound_vars.erase(node.index_var());
    }
    void visit(const TransformNode& node) override {
        bool inserted = bound_vars.insert(node.source_var()).second;
        node.body()->accept(*this);
        if (inserted) bound_vars.erase(node.source_var());
        if (!is_bound(node.target_var())) vars.insert(node.target_var());
    }
    void visit(const QuantifierNode& node) override {
        node.domain()->accept(*this);
        bool inserted = bound_vars.insert(node.bound_var()).second;
        node.predicate()->accept(*this);
        if (inserted) bound_vars.erase(node.bound_var());
    }
    void visit(const SetBuilderNode& node) override {
        node.domain()->accept(*this);
        bool inserted = bound_vars.insert(node.element_var()).second;
        node.predicate()->accept(*this);
        if (inserted) bound_vars.erase(node.element_var());
    }
    void visit(const ComplexNode& node) override {
        node.real()->accept(*this);
        node.imag()->accept(*this);
    }
};

class SubstituteVisitor : public lamina::detail::SymbolicVisitor {
    std::string var_name;
    std::shared_ptr<const SymbolicNode> new_val;
public:
    std::shared_ptr<const SymbolicNode> result;

    SubstituteVisitor(std::string v, std::shared_ptr<const SymbolicNode> val)
        : var_name(std::move(v)), new_val(std::move(val)) {}

    std::shared_ptr<const SymbolicNode> get_result() const { return result; }

    void visit(const NumberNode& node) override {
        result = node.clone();
    }

    void visit(const VariableNode& node) override {
        if (node.name() == var_name) {
            result = new_val->clone();
        } else {
            result = node.clone();
        }
    }

    void visit(const AddNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (const auto& op : node.operands()) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = lamina::detail::make_node<AddNode>(new_ops);
    }

    void visit(const MultiplyNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (const auto& op : node.operands()) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = lamina::detail::make_node<MultiplyNode>(new_ops);
    }

    void visit(const PowerNode& node) override {
        node.base()->accept(*this);
        auto new_base = result;
        node.exponent()->accept(*this);
        auto new_exp = result;
        result = lamina::detail::make_node<PowerNode>(new_base, new_exp);
    }

    void visit(const FunctionNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (const auto& arg : node.arguments()) {
            arg->accept(*this);
            new_args.push_back(result);
        }
        result = lamina::detail::make_node<FunctionNode>(node.type(), new_args);
    }

    void visit(const MatrixNode& node) override {

        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
            const auto& dense = std::get<MatrixNode::DenseStorage>(node.storage());
            MatrixNode::DenseStorage new_dense;
            for(const auto& e : dense) {
                e->accept(*this);
                new_dense.push_back(result);
            }
            result = lamina::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_dense);
        } else {
            const auto& sparse = std::get<MatrixNode::SparseStorage>(node.storage());
            MatrixNode::SparseStorage new_sparse;
            for(const auto& [idx, val] : sparse) {
                val->accept(*this);
                new_sparse[idx] = result;
            }
            result = lamina::detail::make_node<MatrixNode>(node.rows(), node.cols(), new_sparse);
        }
    }

    void visit(const RelationalNode& node) override {
        node.left()->accept(*this);
        auto new_left = result;
        node.right()->accept(*this);
        auto new_right = result;
        result = lamina::detail::make_node<RelationalNode>(new_left, new_right, node.op());
    }

    void visit(const LogicalNode& node) override {
        node.left()->accept(*this);
        auto new_left = result;
        std::shared_ptr<const SymbolicNode> new_right = nullptr;
        if (node.right()) {
            node.right()->accept(*this);
            new_right = result;
        }
        result = lamina::detail::make_node<LogicalNode>(new_left, new_right, node.op());
    }

    void visit(const SummationNode& node) override {
        /// 指标变量是绑定变量：仅当替换变量不是指标时才替换通项。
        std::shared_ptr<const SymbolicNode> new_body;
        if (node.index_var() == var_name) {
            new_body = node.body()->clone();
        } else {
            node.body()->accept(*this);
            new_body = result;
        }
        node.lower_bound()->accept(*this);
        auto new_lo = result;
        node.upper_bound()->accept(*this);
        auto new_hi = result;
        result = lamina::detail::make_node<SummationNode>(new_body, node.index_var(), new_lo, new_hi);
    }

    void visit(const ProductNode_Op& node) override {
        std::shared_ptr<const SymbolicNode> new_body;
        if (node.index_var() == var_name) {
            new_body = node.body()->clone();
        } else {
            node.body()->accept(*this);
            new_body = result;
        }
        node.lower_bound()->accept(*this);
        auto new_lo = result;
        node.upper_bound()->accept(*this);
        auto new_hi = result;
        result = lamina::detail::make_node<ProductNode_Op>(new_body, node.index_var(), new_lo, new_hi);
    }

    void visit(const PiecewiseNode& node) override {
        std::vector<PiecewiseNode::Branch> new_branches;
        for (const auto& b : node.branches()) {
            PiecewiseNode::Branch nb;
            b.expression->accept(*this);
            nb.expression = result;
            b.condition->accept(*this);
            nb.condition = result;
            new_branches.push_back(nb);
        }
        std::shared_ptr<const SymbolicNode> new_def = nullptr;
        if (node.default_expr()) {
            node.default_expr()->accept(*this);
            new_def = result;
        }
        result = lamina::detail::make_node<PiecewiseNode>(std::move(new_branches), new_def);
    }

    void visit(const ComplexNode& node) override {
        node.real()->accept(*this);
        auto new_real = result;
        node.imag()->accept(*this);
        auto new_imag = result;
        result = lamina::detail::make_node<ComplexNode>(new_real, new_imag);
    }

    void visit(const TransformNode& node) override {
        std::shared_ptr<const SymbolicNode> new_body;
        if (node.source_var() == var_name) {
            new_body = node.body()->clone();
        } else {
            node.body()->accept(*this);
            new_body = result;
        }
        result = lamina::detail::make_node<TransformNode>(
            node.transform_type(), new_body, node.source_var(), node.target_var());
    }

    void visit(const QuantifierNode& node) override {
        node.domain()->accept(*this);
        auto new_domain = result;

        std::shared_ptr<const SymbolicNode> new_predicate;
        if (node.bound_var() == var_name) {
            new_predicate = node.predicate()->clone();
        } else {
            node.predicate()->accept(*this);
            new_predicate = result;
        }

        result = lamina::detail::make_node<QuantifierNode>(
            node.quantifier_type(), node.bound_var(), new_domain, new_predicate);
    }

    void visit(const SetBuilderNode& node) override {
        node.domain()->accept(*this);
        auto new_domain = result;

        std::shared_ptr<const SymbolicNode> new_predicate;
        if (node.element_var() == var_name) {
            new_predicate = node.predicate()->clone();
        } else {
            node.predicate()->accept(*this);
            new_predicate = result;
        }

        result = lamina::detail::make_node<SetBuilderNode>(
            node.element_var(), new_domain, new_predicate);
    }
};

class IntegrationVisitor : public lamina::detail::SymbolicVisitor {
    std::string var;
    std::shared_ptr<const SymbolicNode> unevaluated_integral(const std::shared_ptr<const SymbolicNode>& node) const {
        std::vector<std::shared_ptr<const SymbolicNode>> args = {node->clone(), lamina::detail::make_node<VariableNode>(var)};
        return lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args);
    }
public:
    std::shared_ptr<const SymbolicNode> result;

    IntegrationVisitor(std::string v) : var(std::move(v)) {}

    std::shared_ptr<const SymbolicNode> get_result() const { return result; }

    void visit(const NumberNode& node) override {

        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        ops.push_back(node.clone());
        ops.push_back(lamina::detail::make_node<VariableNode>(var));
        result = lamina::detail::make_node<MultiplyNode>(ops);
    }

    void visit(const VariableNode& node) override {
        if (node.name() == var) {
            auto two = lamina::detail::make_node<NumberNode>(BigInt(2));
            auto x_pow_2 = lamina::detail::make_node<PowerNode>(node.clone(), two);
            auto half = lamina::detail::make_node<NumberNode>(Rational(1, 2));
            result = lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{half, x_pow_2});
        } else {
            std::vector<std::shared_ptr<const SymbolicNode>> args = {node.clone(), lamina::detail::make_node<VariableNode>(var)};
            result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args);
        }
    }

    void visit(const AddNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (auto& op : node.operands()) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = lamina::detail::make_node<AddNode>(new_ops);
    }

    void visit(const MultiplyNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> x_ops, c_ops;
        for (auto& op : node.operands()) {
            if (lamina::depends_on_var(op, var)) {
                x_ops.push_back(op);
            } else {
                c_ops.push_back(op);
            }
        }

        if (x_ops.empty()) {
            std::vector<std::shared_ptr<const SymbolicNode>> args = {node.clone(), lamina::detail::make_node<VariableNode>(var)};
            result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args);
            return;
        }

        if (x_ops.size() == 1) {
            x_ops[0]->accept(*this);
            auto int_f = result;
            if (c_ops.empty()) {
                result = int_f;
            } else {
                std::vector<std::shared_ptr<const SymbolicNode>> res_ops = c_ops;
                res_ops.push_back(int_f);
                result = lamina::detail::make_node<MultiplyNode>(res_ops);
            }
            return;
        }

        std::vector<std::shared_ptr<const SymbolicNode>> args = {node.clone(), lamina::detail::make_node<VariableNode>(var)};
        result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args);
    }

    void visit(const PowerNode& node) override {

        bool base_is_x = false;
        if (auto v = std::dynamic_pointer_cast<const VariableNode>(node.base())) {
            if (v->name() == var) base_is_x = true;
        }

        if (base_is_x) {

            if (auto num = std::dynamic_pointer_cast<const NumberNode>(node.exponent())) {
                  if (std::holds_alternative<lmmc_real_t>(num->value())) {
                      int eq;
                      lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(num->value()), -1.0, 1e-9, 1e-9, &eq);
                      if (eq) {
                          std::vector<std::shared_ptr<const SymbolicNode>> args = {node.base()};
                          result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, args);
                          return;
                      }
                  }
                  if (std::holds_alternative<BigInt>(num->value()) && std::get<BigInt>(num->value()) == BigInt(-1)) {
                     std::vector<std::shared_ptr<const SymbolicNode>> args = {node.base()};
                     result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, args);
                     return;
                 }

                 auto one = lamina::detail::make_node<NumberNode>(BigInt(1));
                 std::vector<std::shared_ptr<const SymbolicNode>> add_ops = {node.exponent(), one};
                 auto n_plus_1 = lamina::detail::make_node<AddNode>(add_ops);

                 NormalizationVisitor norm_exp;
                 n_plus_1->accept(norm_exp);
                 auto n_plus_1_sched = norm_exp.get_result();

                 auto new_pow = lamina::detail::make_node<PowerNode>(node.base(), n_plus_1_sched);

                 auto minus_one = lamina::detail::make_node<NumberNode>(BigInt(-1));
                 auto denom = lamina::detail::make_node<PowerNode>(n_plus_1_sched, minus_one);

                 std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {denom, new_pow};
                 result = lamina::detail::make_node<MultiplyNode>(mul_ops);
                 return;
            }
        }

        std::vector<std::shared_ptr<const SymbolicNode>> fallback_args = {node.clone(), lamina::detail::make_node<VariableNode>(var)};
        result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, fallback_args);
    }

    void visit(const FunctionNode& node) override {

        bool arg_is_x = false;
        if (node.arguments().size() == 1) {
            if (auto v = std::dynamic_pointer_cast<const VariableNode>(node.arguments()[0])) {
                if (v->name() == var) arg_is_x = true;
            }
        }

        if (arg_is_x) {
            std::vector<std::shared_ptr<const SymbolicNode>> args = node.arguments();
            if (node.type() == FunctionNode::FuncType::Sin) {

                auto cos_node = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cos, args);
                auto minus_one = lamina::detail::make_node<NumberNode>(BigInt(-1));
                std::vector<std::shared_ptr<const SymbolicNode>> ops = {minus_one, cos_node};
                result = lamina::detail::make_node<MultiplyNode>(ops);
                return;
            } else if (node.type() == FunctionNode::FuncType::Cos) {

                result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sin, args);
                return;
            } else if (node.type() == FunctionNode::FuncType::Exp) {
                result = node.clone();
                return;
            }
        }

        std::vector<std::shared_ptr<const SymbolicNode>> fallback_args = {node.clone(), lamina::detail::make_node<VariableNode>(var)};
        result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, fallback_args);
    }

    void visit(const MatrixNode& node) override {

        result = unevaluated_integral(node.clone());
    }

    void visit(const RelationalNode& node) override {
        result = unevaluated_integral(node.clone());
    }

    void visit(const LogicalNode& node) override {
        result = unevaluated_integral(node.clone());
    }

    void visit(const PiecewiseNode& node) override {
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(node.branches().size());
        for (const auto& branch : node.branches()) {
            branch.expression->accept(*this);
            branches.push_back({result, branch.condition->clone()});
        }
        std::shared_ptr<const SymbolicNode> default_expr = nullptr;
        if (node.default_expr()) {
            node.default_expr()->accept(*this);
            default_expr = result;
        }
        result = lamina::detail::make_node<PiecewiseNode>(std::move(branches), default_expr);
    }

    void visit(const SummationNode& node) override {
        result = unevaluated_integral(node.clone());
    }

    void visit(const ProductNode_Op& node) override {
        result = unevaluated_integral(node.clone());
    }

    void visit(const TransformNode& node) override {
        result = unevaluated_integral(node.clone());
    }

    void visit(const QuantifierNode& node) override {
        result = unevaluated_integral(node.clone());
    }

    void visit(const SetBuilderNode& node) override {
        result = unevaluated_integral(node.clone());
    }

    void visit(const ComplexNode& node) override {
        node.real()->accept(*this);
        auto real_int = result;
        node.imag()->accept(*this);
        auto imag_int = result;
        result = lamina::detail::make_node<ComplexNode>(real_int, imag_int);
    }
};

int SymbolicExpr::compare(const std::shared_ptr<SymbolicExpr>& other) const {
    if (!impl_->root || !lamina::detail::node(other)) return 0;

    return impl_->root->compare(*lamina::detail::node(other));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::substitute(const std::string& var_name, const std::shared_ptr<SymbolicExpr>& value) const {
    if (!impl_->root) return nullptr;
    SubstituteVisitor v(var_name, lamina::detail::node(value));
    impl_->root->accept(v);

    NormalizationVisitor norm;
    v.get_result()->accept(norm);

    return lamina::detail::make_expression_ptr(norm.get_result());
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::expand() const {
    if (!impl_->root) return nullptr;

    ExpandVisitor v;
    impl_->root->accept(v);

    auto result_node = v.get_result();
    if (!result_node) return nullptr;

    return lamina::detail::make_expression_ptr(result_node)->simplify();
}

std::string SymbolicExpr::to_string() const {
    PrintVisitor printer;
    if (impl_->root) {
        impl_->root->accept(printer);
        return printer.get_result();
    }
    return "null";
}

lmmc_real_t SymbolicExpr::to_numeric() const {
    auto evaluated = lamina::evaluate_numeric(*this);
    if (!evaluated) {
        throw std::runtime_error("numeric evaluation failed: " + evaluated.error().message);
    }
    return static_cast<lmmc_real_t>(evaluated.value().value);
}

bool SymbolicExpr::is_zero() const {
    if (!impl_->root) return false;
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(impl_->root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value())) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value());
            int eq;
            lmmc_double_nearly_equal(v, 0.0, &eq);
            return eq != 0;
        }
        if (std::holds_alternative<::BigInt>(num->value())) return std::get<::BigInt>(num->value()).to_int() == 0;
        if (std::holds_alternative<::Rational>(num->value())) return std::get<::Rational>(num->value()).get_numerator().to_int() == 0;
    }
    return false;
}

bool SymbolicExpr::is_one() const {
    if (!impl_->root) return false;
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(impl_->root)) {
        if (std::holds_alternative<lmmc_real_t>(num->value())) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value());
            int eq;
            lmmc_double_nearly_equal(v, 1.0, &eq);
            return eq != 0;
        }
        if (std::holds_alternative<::BigInt>(num->value())) return std::get<::BigInt>(num->value()).to_int() == 1;
        if (std::holds_alternative<::Rational>(num->value())) return std::get<::Rational>(num->value()).to_double() == 1.0;
    }
    return false;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify() const {
    if (!impl_->root) return nullptr;
    NormalizationVisitor v;
    impl_->root->accept(v);
    return lamina::detail::make_expression_ptr(v.get_result());
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::simplify_trig() const {
    auto res = simplify();
    if (!lamina::detail::node(res)) return nullptr;

    static lamina::RewriteEngine engine;
    static bool init = false;

    if (!init) {
        init = true;
        using namespace lamina;
        auto x_val = wildcard("x");
        auto x = lamina::detail::make_expression_ptr(x_val);

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
    auto result_ptr = lamina::detail::make_expression_ptr(simplified);
    return result_ptr->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::differentiate(const std::string& var_name) const {
    if (!impl_->root) return nullptr;
    DifferentiationVisitor v(var_name);
    impl_->root->accept(v);

    NormalizationVisitor norm;
    v.get_result()->accept(norm);

    return lamina::detail::make_expression_ptr(norm.get_result());
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::differentiate_legacy(const std::string& var_name) const {
    return differentiate(var_name);
}

/// 多元因式分解辅助函数

/**
 * @internal
 * @brief 判断符号节点是否为多项式表达式（不含超越函数、负指数等）
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

std::shared_ptr<SymbolicExpr> SymbolicExpr::factor() const {
    auto simp = simplify();
    if (!simp || !lamina::detail::node(simp)) return simp;

    if (auto add_node = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(simp))) {

        /// 步骤 1：提取各项的公因式（GCD）
        std::shared_ptr<SymbolicExpr> common = nullptr;
        for (const auto& op : add_node->operands()) {
             auto expr_op = lamina::detail::make_expression_ptr(op);
             if (!common) common = expr_op;
             else common = poly_gcd(common, expr_op);
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
             auto factored_sum = new_sum->factor();

             std::vector<std::shared_ptr<const SymbolicNode>> final_ops = {lamina::detail::node(common), lamina::detail::node(factored_sum)};
             return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(final_ops));
        }

        /// 步骤 2：一元多项式分解（支持任意次数）
        VariablesVisitor vv;
        lamina::detail::node(simp)->accept(vv);
        if (vv.vars.size() == 1) {
             std::string var = *vv.vars.begin();
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

                           /// 计算余下的不可约因子：原多项式 / 已分解因子的乘积
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
                                   auto q_factored = q_expr->factor();
                                   factors.push_back(lamina::detail::node(q_factored));
                               } else {
                                   auto q_expr = lamina::poly_to_symbolic(quotient)->simplify();
                                   factors.push_back(lamina::detail::node(q_expr));
                               }
                           } else if (!remainder.is_zero()) {
                               /// 除法有余数 → 回退到原始方法
                               goto try_solve_quadratic;
                           }

                           if (factors.size() > 1) {
                               return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(factors));
                           } else if (factors.size() == 1) {
                               return lamina::detail::make_expression_ptr(factors[0]);
                           }
                      }
                 }
             } catch (...) {}

             try_solve_quadratic:
             /// 后备：二次多项式通过求解方程分解
             try {
                 auto poly = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(simp, *vv.vars.begin());
                 if (poly.degree() == 2) {
                      std::string var = *vv.vars.begin();
                      auto solutions = solve(simp, var);
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
             } catch (...) {}
        }

        /// 步骤 3：多元多项式 — 先尝试 factor_multivariate，再回退逐变量分解
        if (vv.vars.size() > 1) {
            /// 3a: 检测是否为多项式表达式，若是则使用 MultiPoly 路径
            if (is_poly_expr_node(lamina::detail::node(simp))) {
                try {
                    std::vector<std::string> var_list(vv.vars.begin(), vv.vars.end());
                    auto mpoly = symbolic_node_to_multipoly(lamina::detail::node(simp), var_list);

                    if (!mpoly.is_zero() && !mpoly.is_constant()) {
                        auto result = lamina::factor_multivariate(mpoly);

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
                } catch (...) {
                    /// MultiPoly 路径失败，回退到逐变量分解
                }
            }

            /// 3b: 回退 — 逐变量尝试有理根分解
            for (const auto& var : vv.vars) {
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
                                /// 常数商 → 已完全分解
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
                } catch (...) {
                    continue;
                }
            }
        }
    }

    /// 超越因式分解后备路径：当标准多项式分解无法处理时，
    /// 检测表达式是否含超越函数，若是则尝试超越因式分解。
    {
        VariablesVisitor tv;
        lamina::detail::node(simp)->accept(tv);
        std::string target_var;
        if (!tv.vars.empty()) {
            target_var = *tv.vars.begin();
        }

        if (!target_var.empty()) {
            try {
                auto trans_factors = lamina::factor_transcendental(simp, target_var);
                if (trans_factors.size() > 1) {
                    /// 超越分解成功：组装乘积表达式
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
            } catch (...) {
                /// 超越分解失败，返回原表达式
            }
        }
    }

    return simp;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::cancel() const {
    auto simp = simplify();
    if (!simp || !lamina::detail::node(simp)) return simp;

    /// 辅助 lambda：从乘积节点中分离分子因子和分母因子。
    /// 分母因子 = 含负指数的 PowerNode。
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

    /// 辅助 lambda：从因子列表构建乘积表达式
    auto build_product = [](const std::vector<std::shared_ptr<const SymbolicNode>>& factors) -> std::shared_ptr<SymbolicExpr> {
        if (factors.empty()) return SymbolicExpr::number(1);
        if (factors.size() == 1) return lamina::detail::make_expression_ptr(factors[0]);
        return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(factors));
    };

    /// 策略 1：表达式是 AddNode，各项可能含公共分母因子。
    /// 例如 simplify 后 (x²-1)/(x-1) 变为 -1*(x-1)^-1 + x²*(x-1)^-1
    /// 需要提取公共分母，重组为 (分子之和)/分母 再做 GCD 约分。
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
            /// 计算公共分母（所有项分母的 LCM，简化处理：乘积）
            /// 对于常见情况（所有项分母相同），直接提取。
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
                /// 所有项分母相同：分子直接相加
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
                /// 分母不同：通分（乘以其他项的分母）
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

            VariablesVisitor vv;
            if (lamina::detail::node(combined_num)) lamina::detail::node(combined_num)->accept(vv);
            if (lamina::detail::node(combined_den)) lamina::detail::node(combined_den)->accept(vv);

            auto cur_num = combined_num;
            auto cur_den = combined_den;

            for (const auto& var : vv.vars) {
                try {
                    auto num_expanded = cur_num->expand()->simplify();
                    auto den_expanded = cur_den->expand()->simplify();

                    /// 使用 SymbolicPolyCoeff 保留其他变量作为符号系数
                    auto poly_num = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(num_expanded, var);
                    auto poly_den = lamina::symbolic_to_poly<lamina::SymbolicPolyCoeff>(den_expanded, var);

                    if (poly_num.is_zero()) return SymbolicExpr::number(0);
                    if (poly_den.is_zero()) return simp;

                    /// 尝试用 Rational 系数做 GCD（纯数值系数时更可靠）
                    auto poly_num_r = lamina::symbolic_to_poly<Rational>(num_expanded, var);
                    auto poly_den_r = lamina::symbolic_to_poly<Rational>(den_expanded, var);

                    /// 检查 Rational 转换是否丢失了信息（多元情况）
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
                        /// 多元情况：使用 SymbolicPolyCoeff 做 GCD
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

    /// 策略 2：表达式是 MultiplyNode，直接分离分子/分母。
    std::vector<std::shared_ptr<const SymbolicNode>> num_factors;
    std::vector<std::shared_ptr<const SymbolicNode>> den_factors;
    separate_num_den(lamina::detail::node(simp), num_factors, den_factors);

    if (den_factors.empty()) return simp;

    auto numerator = build_product(num_factors)->simplify();
    auto denominator = build_product(den_factors)->simplify();

    if (denominator->is_one()) return numerator;

    VariablesVisitor vv;
    if (lamina::detail::node(numerator)) lamina::detail::node(numerator)->accept(vv);
    if (lamina::detail::node(denominator)) lamina::detail::node(denominator)->accept(vv);

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
        } catch (...) {
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

std::shared_ptr<SymbolicExpr> SymbolicExpr::limit(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, const std::string& direction) const {
    if (!impl_->root) return nullptr;

    LimitVisitor v(var, lamina::detail::node(point), direction);
    impl_->root->accept(v);

    if (v.get_result()) {
        return lamina::detail::make_expression_ptr(v.get_result());
    }

    return nullptr;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::integrate(const std::string& var) const {
    if (!impl_->root) return nullptr;

    lamina::Integrator integrator;
    SymbolicExpr result = integrator.integrate(*this, var);

    auto res_ptr = lamina::detail::make_expression_ptr(result);
    return res_ptr->simplify();
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::series(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, int order, const lamina::AssumptionContext* ctx) const {
    if (!impl_->root || !point) return nullptr;
    if (order < 0) return SymbolicExpr::number(0);
    if (ctx) {
    }

    auto x = SymbolicExpr::variable(var);
    auto neg_point = SymbolicExpr::multiply(SymbolicExpr::number(-1), point);
    auto delta = SymbolicExpr::add(x, neg_point);

    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    auto deriv = lamina::detail::make_expression_ptr(impl_->root->clone());

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

        terms.push_back(lamina::detail::node(term));
    }

    if (terms.empty()) return SymbolicExpr::number(0);
    auto sum = lamina::detail::make_node<AddNode>(terms);
    return lamina::detail::make_expression_ptr(sum)->simplify();
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den) {
    return SymbolicExpr::multiply(num, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::make_integral(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr) return nullptr;
    auto var_node = lamina::detail::make_node<VariableNode>(var);
    std::vector<std::shared_ptr<const SymbolicNode>> args = {lamina::detail::node(expr), var_node};
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args));
}
std::shared_ptr<SymbolicExpr> SymbolicExpr::make_limit(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var, const std::shared_ptr<SymbolicExpr>& point) {
    if (!expr || !point) return nullptr;
    auto var_node = lamina::detail::make_node<VariableNode>(var);
    std::vector<std::shared_ptr<const SymbolicNode>> args = {lamina::detail::node(expr), var_node, lamina::detail::node(point)};
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Limit, args));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::poly_gcd(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    if (!a || !b) return nullptr;

    try {
        struct VarVisitor : public lamina::detail::SymbolicVisitor {
            std::set<std::string> vars;
            void visit(const NumberNode&) override {}
            void visit(const VariableNode& n) override { vars.insert(n.name()); }
            void visit(const AddNode& n) override { for(auto& op : n.operands()) op->accept(*this); }
            void visit(const MultiplyNode& n) override { for(auto& op : n.operands()) op->accept(*this); }
            void visit(const PowerNode& n) override { n.base()->accept(*this); n.exponent()->accept(*this); }
            void visit(const FunctionNode& n) override { for(auto& arg : n.arguments()) arg->accept(*this); }
            void visit(const MatrixNode&) override {}
            void visit(const RelationalNode& n) override { n.left()->accept(*this); n.right()->accept(*this); }
            void visit(const LogicalNode& n) override { n.left()->accept(*this); if (n.right()) n.right()->accept(*this); }
            void visit(const PiecewiseNode& n) override {
                for (const auto& branch : n.branches()) {
                    branch.expression->accept(*this);
                    branch.condition->accept(*this);
                }
                if (n.default_expr()) n.default_expr()->accept(*this);
            }
            void visit(const SummationNode& n) override {
                n.lower_bound()->accept(*this);
                n.upper_bound()->accept(*this);
                n.body()->accept(*this);
                vars.erase(n.index_var());
            }
            void visit(const ProductNode_Op& n) override {
                n.lower_bound()->accept(*this);
                n.upper_bound()->accept(*this);
                n.body()->accept(*this);
                vars.erase(n.index_var());
            }
            void visit(const TransformNode& n) override {
                n.body()->accept(*this);
                vars.erase(n.source_var());
                vars.insert(n.target_var());
            }
            void visit(const QuantifierNode& n) override {
                n.domain()->accept(*this);
                n.predicate()->accept(*this);
                vars.erase(n.bound_var());
            }
            void visit(const SetBuilderNode& n) override {
                n.domain()->accept(*this);
                n.predicate()->accept(*this);
                vars.erase(n.element_var());
            }
            void visit(const ComplexNode& n) override {
                n.real()->accept(*this);
                n.imag()->accept(*this);
            }
        } vv;
        if (lamina::detail::node(a)) lamina::detail::node(a)->accept(vv);
        if (lamina::detail::node(b)) lamina::detail::node(b)->accept(vv);

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

        /// 退化情形：任一为零多项式 → 结式为 0
        if (m < 0 || n < 0) return SymbolicExpr::number(0);
        /// 任一为非零常数：Res(a, const c) = c^deg(other)（两者都常数时为 1）
        if (m == 0 && n == 0) return SymbolicExpr::number(1);
        if (n == 0) {
            /// Res(a, c) = c^m
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

        /// 构造 (m+n)×(m+n) Sylvester 矩阵。
        /// 前 n 行由 a 的系数移位排列，后 m 行由 b 的系数移位排列。
        /// 行内按降幂排列：列 0 对应 x^(m+n-1) 系数。
        int N = m + n;
        std::vector<std::vector<Rational>> S(N, std::vector<Rational>(N, Rational(0)));

        /// a 的降幂系数：a_m, a_{m-1}, ..., a_0  （pa.coeffs[m..0]）
        for (int row = 0; row < n; ++row) {
            for (int k = 0; k <= m; ++k) {
                /// a 的 x^(m-k) 系数 = pa.coeffs[m-k]
                S[row][row + k] = pa.coeffs[m - k];
            }
        }
        /// b 的降幂系数
        for (int row = 0; row < m; ++row) {
            for (int k = 0; k <= n; ++k) {
                S[n + row][row + k] = pb.coeffs[n - k];
            }
        }

        /// Bareiss 无分式高斯消元求行列式（精确 Rational 运算）。
        Rational det(1);
        Rational prev(1);
        int sign = 1;
        for (int i = 0; i < N; ++i) {
            /// 主元为 0 时寻找下方非零行交换
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
    if (!impl_->root) return Type::Number;

    if (std::dynamic_pointer_cast<const NumberNode>(impl_->root)) return Type::Number;
    if (std::dynamic_pointer_cast<const VariableNode>(impl_->root)) return Type::Variable;
    if (std::dynamic_pointer_cast<const AddNode>(impl_->root)) return Type::Add;
    if (std::dynamic_pointer_cast<const MultiplyNode>(impl_->root)) return Type::Multiply;
    if (std::dynamic_pointer_cast<const PowerNode>(impl_->root)) return Type::Power;

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(impl_->root)) {
        switch (func->type()) {
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

    if (auto mat = std::dynamic_pointer_cast<const MatrixNode>(impl_->root)) {
        if (mat->rows() == 1 || mat->cols() == 1) return Type::Vector;
        return Type::Matrix;
    }

    return Type::Number;
}

std::vector<std::shared_ptr<SymbolicExpr>> SymbolicExpr::get_operands() const {
    std::vector<std::shared_ptr<SymbolicExpr>> ops;
    if (!impl_->root) return ops;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(impl_->root)) {
        for (const auto& op : add->operands()) ops.push_back(lamina::detail::make_expression_ptr(op));
    } else if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(impl_->root)) {
        for (const auto& op : mul->operands()) ops.push_back(lamina::detail::make_expression_ptr(op));
    } else if (auto pow = std::dynamic_pointer_cast<const PowerNode>(impl_->root)) {
        ops.push_back(lamina::detail::make_expression_ptr(pow->base()));
        ops.push_back(lamina::detail::make_expression_ptr(pow->exponent()));
    } else if (auto func = std::dynamic_pointer_cast<const FunctionNode>(impl_->root)) {
        for (const auto& arg : func->arguments()) ops.push_back(lamina::detail::make_expression_ptr(arg));
    }
    return ops;
}

std::variant<int, ::BigInt, ::Rational> SymbolicExpr::get_number_value() const {
    if (auto node = std::dynamic_pointer_cast<const NumberNode>(impl_->root)) {
        if (std::holds_alternative<BigInt>(node->value())) return std::get<BigInt>(node->value());
        if (std::holds_alternative<Rational>(node->value())) return std::get<Rational>(node->value());
        if (std::holds_alternative<lmmc_real_t>(node->value())) return (int)std::get<lmmc_real_t>(node->value());
    }
    return 0;
}

std::string SymbolicExpr::get_identifier() const {
    if (auto v = std::dynamic_pointer_cast<const VariableNode>(impl_->root)) return v->name();
    return "";
}
