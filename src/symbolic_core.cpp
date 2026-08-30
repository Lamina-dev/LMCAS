#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include "../include/symbolic.hpp"
#include "../include/symbolic_ast.hpp"
#include "../include/internal/expression_analysis.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include "../include/visitors/differentiation_visitor.hpp"
#include "../include/visitors/expand_visitor.hpp"
#include "../include/matcher.hpp"
#include "../include/numeric_evaluation.hpp"
#include "../include/root_of_utils.hpp"
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
    if (index < 0) throw std::invalid_argument("root_of index cannot be negative");
    auto result = lamina::make_rootof_checked(
        polynomial, variable_name, static_cast<std::size_t>(index));
    if (!result) throw std::invalid_argument(result.error().message);
    return std::move(result.value());
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


int SymbolicExpr::compare(const std::shared_ptr<SymbolicExpr>& other) const {
    if (!impl_->root || !lamina::detail::node(other)) return 0;

    return impl_->root->compare(*lamina::detail::node(other));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::substitute(
    const std::string& var_name,
    const std::shared_ptr<SymbolicExpr>& value) const {
    if (!value) return nullptr;
    auto substituted = lamina::substitute_free(
        impl_->root, var_name, lamina::detail::node(value));
    NormalizationVisitor normalization;
    substituted->accept(normalization);
    return lamina::detail::make_expression_ptr(normalization.get_result());
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

    static const lamina::RewriteEngine engine = [] {
        lamina::RewriteEngine configured;
        using namespace lamina;
        auto x_val = wildcard("x");
        auto x = lamina::detail::make_expression_ptr(x_val);

        auto sinx = SymbolicExpr::sin(x);
        auto cosx = SymbolicExpr::cos(x);
        auto n2 = SymbolicExpr::number(2);
        auto sin2 = SymbolicExpr::power(sinx, n2);
        auto cos2 = SymbolicExpr::power(cosx, n2);

        auto pat1 = SymbolicExpr::add(sin2, cos2);
        configured.add_rule(Rule(*pat1, *SymbolicExpr::number(1), {"x"}));

        auto pat2 = SymbolicExpr::add(cos2, sin2);
        configured.add_rule(Rule(*pat2, *SymbolicExpr::number(1), {"x"}));

        auto two_x = SymbolicExpr::multiply(n2, x);
        auto sin2x = SymbolicExpr::sin(two_x);
        auto two_sin_cos = SymbolicExpr::multiply(n2,
            SymbolicExpr::multiply(sinx, cosx));
        configured.add_rule(Rule(*sin2x, *two_sin_cos, {"x"}));

        auto cos2x = SymbolicExpr::cos(two_x);
        auto cos2_sub_sin2 = SymbolicExpr::add(cos2,
            SymbolicExpr::multiply(sin2, SymbolicExpr::number(-1)));
        configured.add_rule(Rule(*cos2x, *cos2_sub_sin2, {"x"}));
        return configured;
    }();

    lamina::ComputationContext context;
    auto rewritten = engine.apply_checked(*res, context);
    if (!rewritten) return res;
    auto simplified = std::move(rewritten.value());
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


std::shared_ptr<SymbolicExpr> SymbolicExpr::divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den) {
    return SymbolicExpr::multiply(num, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
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
