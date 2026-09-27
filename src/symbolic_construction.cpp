#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_construction.hpp"
#include "root_of_identity.hpp"
#include <exception>
#include <new>
#include <utility>

namespace LMCAS {

ExpressionResult detail::constant_expression(const char* name) {
    constexpr const char* operation = "LMCAS.constant";
    try {
        auto expression = SymbolicExpr::variable(name);
        if (!expression || !LMCAS::detail::node(expression)) {
            return ExpressionResult::failure(CasErrc::InternalInvariant,
                                             "constant factory returned null",
                                             operation);
        }
        return ExpressionResult::success(std::move(expression));
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                         "constant allocation failed",
                                         operation);
    } catch (const std::exception& error) {
        return ExpressionResult::failure(CasErrc::InvalidArgument, error.what(),
                                         operation);
    }
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::number(int n) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(BigInt(n)));
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::number(double n) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(n)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::number(const BigInt& value) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(value));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::number(const Rational& value) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(value));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::infinity(int sign) {
    auto node = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Infinity,
        std::vector<std::shared_ptr<const SymbolicNode>>{});
    auto expression = LMCAS::detail::make_expression_ptr(node);
    return sign < 0 ? multiply(number(-1), expression) : expression;
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::sqrt(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("sqrt operand cannot be null");
    }
    auto half = LMCAS::detail::make_node<NumberNode>(Rational(1, 2));
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<PowerNode>(operand->impl_->root, half));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::multiply(
    std::shared_ptr<SymbolicExpr> left,
    std::shared_ptr<SymbolicExpr> right) {
    if (!left || !left->impl_->root || !right || !right->impl_->root) {
        throw std::invalid_argument("multiply operands cannot be null");
    }
    std::vector<std::shared_ptr<const SymbolicNode>> operands{
        left->impl_->root, right->impl_->root};
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<MultiplyNode>(std::move(operands)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::add(
    std::shared_ptr<SymbolicExpr> left,
    std::shared_ptr<SymbolicExpr> right) {
    if (!left || !left->impl_->root || !right || !right->impl_->root) {
        throw std::invalid_argument("add operands cannot be null");
    }
    std::vector<std::shared_ptr<const SymbolicNode>> operands{
        left->impl_->root, right->impl_->root};
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<AddNode>(std::move(operands)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::power(
    std::shared_ptr<SymbolicExpr> base,
    std::shared_ptr<SymbolicExpr> exponent) {
    if (!base || !base->impl_->root || !exponent || !exponent->impl_->root) {
        throw std::invalid_argument("power operands cannot be null");
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<PowerNode>(base->impl_->root, exponent->impl_->root));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::sin(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("sin operand cannot be null");
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::cos(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("cos operand cannot be null");
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Cos,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::tan(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("tan operand cannot be null");
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Tan,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::ln(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("ln operand cannot be null");
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Ln,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::exp(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("exp operand cannot be null");
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::lambertw(
    std::shared_ptr<SymbolicExpr> operand) {
    if (!operand || !operand->impl_->root) {
        throw std::invalid_argument("lambertw operand cannot be null");
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::LambertW,
            std::vector<std::shared_ptr<const SymbolicNode>>{operand->impl_->root}));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::log(
    std::shared_ptr<SymbolicExpr> value,
    std::shared_ptr<SymbolicExpr> base) {
    if (!value || !value->impl_->root || !base || !base->impl_->root) {
        throw std::invalid_argument("log operands cannot be null");
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
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
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
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
    auto result = LMCAS::make_rootof_checked(
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
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<RelationalNode>(
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
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<MatrixNode>(std::move(node_elements)));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::variable(
    const std::string& name) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<VariableNode>(name));
}

std::shared_ptr<SymbolicExpr> SymbolicExpr::divide(const std::shared_ptr<SymbolicExpr>& num, const std::shared_ptr<SymbolicExpr>& den) {
    return SymbolicExpr::multiply(num, SymbolicExpr::power(den, SymbolicExpr::number(-1)));
}

}
