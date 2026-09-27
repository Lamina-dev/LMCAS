#include "internal/expr_parser.hpp"
#include "internal/expr_common.hpp"
#include "root_of_identity.hpp"

#include <exception>
#include <utility>

namespace LMCAS::expr_detail {
namespace {

using UnaryFunction = ExprResult (*)(const ExprPtr&);

UnaryFunction find_unary_function(const std::string& name) {
    static const struct {
        const char* name;
        UnaryFunction function;
    } functions[] = {
        {"sin", LMCAS::sin},
        {"cos", LMCAS::cos},
        {"tan", LMCAS::tan},
        {"asin", LMCAS::asin},
        {"acos", LMCAS::acos},
        {"atan", LMCAS::atan},
        {"sqrt", LMCAS::sqrt},
        {"exp", LMCAS::exp},
        {"ln", LMCAS::log},
        {"log", LMCAS::log},
        {"log10", LMCAS::log10},
        {"floor", LMCAS::floor},
        {"ceil", LMCAS::ceil},
        {"round", LMCAS::round},
        {"abs", LMCAS::abs},
        {"real", LMCAS::real},
        {"imag", LMCAS::imag},
        {"conj", LMCAS::conj}
    };
    for (const auto& entry : functions) {
        if (name == entry.name) {
            return entry.function;
        }
    }
    return nullptr;
}

}

ExprResult ExprParser::apply_log(const std::vector<ExprPtr>& arguments) {
    try {
        auto result = SymbolicExpr::log(arguments[0], arguments[1]);
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expr_detail::expr_common::expression_failure(CasErrc::ResourceLimit,
                                  "log expression allocation failed",
                                  kParseOperation);
    } catch (const std::exception& error) {
        return expr_detail::expr_common::expression_failure(CasErrc::InternalInvariant, error.what(),
                                  kParseOperation);
    }
}

ExprResult ExprParser::apply_atan2(const std::vector<ExprPtr>& arguments) {
    try {
        auto result = SymbolicExpr::atan2(arguments[0], arguments[1]);
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expr_detail::expr_common::expression_failure(CasErrc::ResourceLimit,
                                  "atan2 expression allocation failed",
                                  kParseOperation);
    } catch (const std::exception& error) {
        return expr_detail::expr_common::expression_failure(CasErrc::InternalInvariant, error.what(),
                                  kParseOperation);
    }
}

ExprResult ExprParser::apply_integral(const std::vector<ExprPtr>& arguments) {
    auto variable = std::dynamic_pointer_cast<const VariableNode>(
        LMCAS::detail::node(arguments[1]));
    if (!variable) {
        return fail("Integral variable must be a symbol");
    }
    try {
        auto node = LMCAS::detail::make_node<IntegralNode>(
            LMCAS::detail::node(arguments[0]), variable->name(),
            arguments.size() == 4
                ? LMCAS::detail::node(arguments[2])
                : nullptr,
            arguments.size() == 4
                ? LMCAS::detail::node(arguments[3])
                : nullptr);
        return ExprResult::success(
            LMCAS::detail::make_expression_ptr(std::move(node)));
    } catch (const std::bad_alloc&) {
        return expr_detail::expr_common::expression_failure(CasErrc::ResourceLimit,
        "Integral expression allocation failed",
        kParseOperation);
    } catch (const std::exception& error) {
        return expr_detail::expr_common::expression_failure(CasErrc::InternalInvariant, error.what(), kParseOperation);
    }
}

ExprResult ExprParser::apply_limit(const std::vector<ExprPtr>& arguments) {
    auto variable = std::dynamic_pointer_cast<const VariableNode>(
        LMCAS::detail::node(arguments[1]));
    auto direction = std::dynamic_pointer_cast<const VariableNode>(
        LMCAS::detail::node(arguments[3]));
    if (!variable || !direction) {
        return fail("limit variable and direction must be symbols");
    }
    LimitDirection parsed_direction;
    if (direction->name() == "both") {
        parsed_direction = LimitDirection::Both;
    } else if (direction->name() == "below" ||
               direction->name() == "left" ||
               direction->name() == "from_below") {
        parsed_direction = LimitDirection::FromBelow;
    } else if (direction->name() == "above" ||
               direction->name() == "right" ||
               direction->name() == "from_above") {
        parsed_direction = LimitDirection::FromAbove;
    } else {
        return fail("invalid limit direction");
    }
    return ExprResult::success(LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<LimitNode>(
            LMCAS::detail::node(arguments[0]), variable->name(),
            LMCAS::detail::node(arguments[2]), parsed_direction)));
}

ExprResult ExprParser::apply_rootof(const std::vector<ExprPtr>& arguments) {
    auto variable = std::dynamic_pointer_cast<const VariableNode>(
        LMCAS::detail::node(arguments[1]));
    auto index = std::dynamic_pointer_cast<const NumberNode>(
        LMCAS::detail::node(arguments[2]));
    if (!variable || !index ||
        std::holds_alternative<lmmc_real_t>(index->value())) {
        return fail("RootOf variable and index are invalid");
    }
    BigInt index_value;
    if (std::holds_alternative<BigInt>(index->value())) {
        index_value = std::get<BigInt>(index->value());
    } else {
        const auto& rational = std::get<Rational>(index->value());
        if (!rational.is_integer()) {
            return fail("RootOf index must be an integer");
        }
        index_value = rational.to_bigint();
    }
    const auto converted = index_value.try_to_int64();
    if (!converted || *converted < 0) {
        return fail("RootOf index is negative or too large");
    }
    auto root = make_rootof_checked(
        arguments[0], variable->name(),
        static_cast<std::size_t>(*converted));
    if (!root) {
        return ExprResult::failure(root.error());
    }
    return ExprResult::success(std::move(root.value()));
}

ExprResult ExprParser::apply_clamp(const std::vector<ExprPtr>& arguments) {
    return clamp(arguments[0], arguments[1], arguments[2]);
}

ExprResult ExprParser::apply_function(const std::string& name,
                                    const std::vector<ExprPtr>& arguments) {
    if (auto unary = find_unary_function(name)) {
        if (arguments.size() == 1) {
            return unary(arguments[0]);
        }
        if (name == "log" && arguments.size() == 2) {
            return apply_log(arguments);
        }
        return fail("invalid argument count for function '" + name + "'");
    }
    using Constructor = ExprResult (ExprParser::*)(const std::vector<ExprPtr>&);
    static const struct {
        const char* name;
        Constructor construct;
        std::size_t arity;
        std::size_t alternative_arity;
    } functions[] = {
        {"atan2", &ExprParser::apply_atan2, 2, 2},
        {"clamp", &ExprParser::apply_clamp, 3, 3},
        {"Integral", &ExprParser::apply_integral, 2, 4},
        {"integral", &ExprParser::apply_integral, 2, 4},
        {"Limit", &ExprParser::apply_limit, 4, 4},
        {"limit", &ExprParser::apply_limit, 4, 4},
        {"RootOf", &ExprParser::apply_rootof, 3, 3},
        {"rootof", &ExprParser::apply_rootof, 3, 3}
    };
    for (const auto& entry : functions) {
        if (name == entry.name) {
            if (arguments.size() != entry.arity &&
                arguments.size() != entry.alternative_arity) {
                return fail("invalid argument count for function '" + name + "'");
            }
            return (this->*entry.construct)(arguments);
        }
    }
    if (!arguments.empty() && (name == "max" || name == "min")) {
        return function_node(name, arguments,
            name == "max" ? FunctionNode::FuncType::Max : FunctionNode::FuncType::Min);
    }
    return uninterpreted_function(name, arguments);
}

ExprResult ExprParser::uninterpreted_function(
    const std::string& name, const std::vector<ExprPtr>& arguments) {
    try {
        std::vector<std::shared_ptr<const SymbolicNode>> nodes;
        nodes.reserve(arguments.size());
        for (const auto& argument : arguments) {
            if (!argument || !LMCAS::detail::node(argument)) {
                return expr_detail::expr_common::expression_failure(CasErrc::InvalidArgument,
                                          "function arguments cannot be null",
                                          kParseOperation);
            }
            nodes.push_back(LMCAS::detail::node(argument));
        }
        auto node = LMCAS::detail::make_node<UninterpretedFunctionNode>(
            name, std::move(nodes));
        return ExprResult::success(LMCAS::detail::make_expression_ptr(std::move(node)));
    } catch (const std::bad_alloc&) {
        return expr_detail::expr_common::expression_failure(CasErrc::ResourceLimit,
                                  "function expression allocation failed",
                                  kParseOperation);
    } catch (const std::exception& error) {
        return expr_detail::expr_common::expression_failure(CasErrc::InternalInvariant, error.what(),
                                  kParseOperation);
    }
}

}
