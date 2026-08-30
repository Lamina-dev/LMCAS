#include "lsr_expr.hpp"

#include <exception>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>

#include "assumption_context.hpp"
#include "lsr_expr_internal.hpp"
#include "poly_utils.hpp"
#include "symbolic_ast.hpp"

namespace lamina::lsr {
namespace {

Result<EqvProfile> eqv_profile_failure(std::string message) {
    return Result<EqvProfile>::failure(
        CasErrc::UnsupportedExpression,
        std::move(message),
        kEquivalentProfileOperation);
}

Result<void> validate_eqv_options(const EqvOptions& options) {
    if (options.budget.max_rewrite_steps == 0 ||
        options.budget.max_rewrite_depth == 0 ||
        options.budget.max_node_growth_factor == 0) {
        return Result<void>::failure(
            CasErrc::ResourceLimit,
            "equivalence rewrite budget exhausted before normalization",
            kEquivalentOperation);
    }
    return Result<void>::success();
}

bool exact_integer_node(const std::shared_ptr<const SymbolicNode>& node,
                        int expected) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return false;
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()) == BigInt(expected);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()) == Rational(expected);
    }
    return false;
}

bool trig_square_argument(const std::shared_ptr<const SymbolicNode>& node,
                          FunctionNode::FuncType type,
                          std::shared_ptr<const SymbolicNode>& argument) {
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power || !exact_integer_node(power->exponent(), 2)) return false;
    auto function = std::dynamic_pointer_cast<const FunctionNode>(power->base());
    if (!function || function->type() != type ||
        function->arguments().size() != 1) {
        return false;
    }
    argument = function->arguments()[0];
    return true;
}

std::optional<ExprPtr> unwrap_trig_negated_argument(
    const std::shared_ptr<const SymbolicNode>& node) {
    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) return std::nullopt;

    bool found_negative_one = false;
    std::vector<std::shared_ptr<const SymbolicNode>> remaining;
    remaining.reserve(multiply->operands().size());
    for (const auto& operand : multiply->operands()) {
        if (!found_negative_one && exact_integer_node(operand, -1)) {
            found_negative_one = true;
            continue;
        }
        remaining.push_back(operand);
    }
    if (!found_negative_one || remaining.empty()) return std::nullopt;
    if (remaining.size() == 1) {
        return lamina::detail::make_expression_ptr(remaining.front());
    }
    return lamina::detail::make_expression_ptr(
        SymbolicFactory::create_multiply(std::move(remaining)))->simplify();
}

ExprPtr rewrite_trig_basic_identity(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return nullptr;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> rewritten;
        rewritten.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            auto child = rewrite_trig_basic_identity(operand);
            rewritten.push_back(lamina::detail::node(child));
        }

        std::vector<bool> used(rewritten.size(), false);
        std::vector<std::shared_ptr<const SymbolicNode>> result_nodes;
        for (std::size_t i = 0; i < rewritten.size(); ++i) {
            if (used[i]) continue;
            std::shared_ptr<const SymbolicNode> sin_arg;
            std::shared_ptr<const SymbolicNode> cos_arg;
            const bool is_sin_square = trig_square_argument(
                rewritten[i], FunctionNode::FuncType::Sin, sin_arg);
            const bool is_cos_square = trig_square_argument(
                rewritten[i], FunctionNode::FuncType::Cos, cos_arg);
            bool matched = false;
            for (std::size_t j = i + 1; j < rewritten.size(); ++j) {
                if (used[j]) continue;
                std::shared_ptr<const SymbolicNode> other_arg;
                if (is_sin_square &&
                    trig_square_argument(rewritten[j],
                                         FunctionNode::FuncType::Cos,
                                         other_arg) &&
                    sin_arg->equals(*other_arg)) {
                    matched = true;
                } else if (is_cos_square &&
                           trig_square_argument(rewritten[j],
                                                FunctionNode::FuncType::Sin,
                                                other_arg) &&
                           cos_arg->equals(*other_arg)) {
                    matched = true;
                }
                if (matched) {
                    used[i] = true;
                    used[j] = true;
                    result_nodes.push_back(
                        lamina::detail::node(SymbolicExpr::number(1)));
                    break;
                }
            }
            if (!matched && !used[i]) {
                result_nodes.push_back(rewritten[i]);
            }
        }

        if (result_nodes.empty()) {
            return SymbolicExpr::number(0);
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_add(std::move(result_nodes)))->simplify();
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(multiply->operands().size());
        for (const auto& operand : multiply->operands()) {
            auto child = rewrite_trig_basic_identity(operand);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_multiply(std::move(operands)))->simplify();
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = rewrite_trig_basic_identity(power->base());
        auto exponent = rewrite_trig_basic_identity(power->exponent());
        return SymbolicExpr::power(base, exponent)->simplify();
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(function->arguments().size());
        for (const auto& argument : function->arguments()) {
            auto child = rewrite_trig_basic_identity(argument);
            args.push_back(lamina::detail::node(child));
        }

        if (args.size() == 1) {
            auto positive_arg = unwrap_trig_negated_argument(args[0]);
            if (positive_arg && *positive_arg) {
                if (function->type() == FunctionNode::FuncType::Sin) {
                    return SymbolicExpr::multiply(
                        SymbolicExpr::number(-1),
                        SymbolicExpr::sin(*positive_arg))->simplify();
                }
                if (function->type() == FunctionNode::FuncType::Cos) {
                    return SymbolicExpr::cos(*positive_arg)->simplify();
                }
            }
        }

        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(function->type(),
                                                    std::move(args)))->simplify();
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real_part = rewrite_trig_basic_identity(complex_node->real());
        auto imag_part = rewrite_trig_basic_identity(complex_node->imag());
        auto value = complex(real_part, imag_part);
        if (!value) throw std::runtime_error(value.error().message);
        return value.value()->simplify();
    }

    return lamina::detail::make_expression_ptr(node);
}

ExprPtr rewrite_exp_log_basic_identity(
    const std::shared_ptr<const SymbolicNode>& node,
    const AssumptionContext* assumptions) {
    if (!node) return nullptr;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            auto child = rewrite_exp_log_basic_identity(operand, assumptions);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_add(std::move(operands)))->simplify();
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(multiply->operands().size());
        for (const auto& operand : multiply->operands()) {
            auto child = rewrite_exp_log_basic_identity(operand, assumptions);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_multiply(std::move(operands)))->simplify();
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = rewrite_exp_log_basic_identity(power->base(), assumptions);
        auto exponent = rewrite_exp_log_basic_identity(power->exponent(), assumptions);
        return SymbolicExpr::power(base, exponent)->simplify();
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(function->arguments().size());
        for (const auto& argument : function->arguments()) {
            auto child = rewrite_exp_log_basic_identity(argument, assumptions);
            args.push_back(lamina::detail::node(child));
        }

        if (function->type() == FunctionNode::FuncType::Exp &&
            args.size() == 1 && exact_integer_node(args[0], 0)) {
            return SymbolicExpr::number(1);
        }
        if (function->type() == FunctionNode::FuncType::Ln &&
            args.size() == 1 && exact_integer_node(args[0], 1)) {
            return SymbolicExpr::number(0);
        }
        if (function->type() == FunctionNode::FuncType::Exp &&
            args.size() == 1) {
            auto inner_ln = std::dynamic_pointer_cast<const FunctionNode>(args[0]);
            if (inner_ln && inner_ln->type() == FunctionNode::FuncType::Ln &&
                inner_ln->arguments().size() == 1) {
                auto ln_arg = lamina::detail::make_expression_ptr(
                    inner_ln->arguments()[0]);
                const bool known_positive = inner_ln->arguments()[0]->is_positive() ||
                    (assumptions &&
                     detail::propagate_result(assumptions->is_positive(*ln_arg)) ==
                         Tribool::True);
                if (known_positive) {
                    return ln_arg->simplify();
                }
            }
        }

        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(function->type(),
                                                    std::move(args)))->simplify();
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real_part = rewrite_exp_log_basic_identity(complex_node->real(),
                                                        assumptions);
        auto imag_part = rewrite_exp_log_basic_identity(complex_node->imag(),
                                                        assumptions);
        auto value = complex(real_part, imag_part);
        if (!value) throw std::runtime_error(value.error().message);
        return value.value()->simplify();
    }

    return lamina::detail::make_expression_ptr(node);
}

void collect_variable_names(
    const std::shared_ptr<const SymbolicNode>& node,
    std::set<std::string>& variables) {
    if (!node) return;
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        variables.insert(variable->name());
        return;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            collect_variable_names(operand, variables);
        }
        return;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            collect_variable_names(operand, variables);
        }
        return;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        collect_variable_names(power->base(), variables);
        collect_variable_names(power->exponent(), variables);
        return;
    }
    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        collect_variable_names(complex_node->real(), variables);
        collect_variable_names(complex_node->imag(), variables);
    }
}

Result<std::optional<bool>> prove_rational_polynomial_equivalence(
    const ExprPtr& difference,
    ComputationContext& context,
    const EqvOptions& options) {
    if (!difference) {
        return Result<std::optional<bool>>::failure(
            CasErrc::InternalInvariant,
            "equivalence difference is null",
            kEquivalentOperation);
    }

    if (options.budget.max_rewrite_steps < 4) {
        return Result<std::optional<bool>>::failure(
            CasErrc::ResourceLimit,
            "equivalence rewrite budget exhausted before polynomial normalization",
            kEquivalentOperation);
    }

    auto step = context.consume_steps(4, kEquivalentOperation);
    if (!step) return Result<std::optional<bool>>::failure(step.error());

    auto expanded = difference->expand();
    if (!expanded || !lamina::detail::node(expanded)) {
        return Result<std::optional<bool>>::failure(
            CasErrc::InternalInvariant,
            "equivalence expansion returned null",
            kEquivalentOperation);
    }

    std::set<std::string> variables;
    collect_variable_names(lamina::detail::node(expanded), variables);
    if (variables.size() > 1) {
        return Result<std::optional<bool>>::success(std::nullopt);
    }
    const std::string variable = variables.empty() ? "x" : *variables.begin();

    auto recognized = recognize_rational_polynomial(*expanded, variable, context);
    if (!recognized) {
        return Result<std::optional<bool>>::failure(recognized.error());
    }
    if (!recognized.value()) {
        return Result<std::optional<bool>>::success(std::nullopt);
    }
    return Result<std::optional<bool>>::success(recognized.value()->is_zero());
}

ExprPtr canonicalize_lsr_complex_product(const SymbolicExpr& expression) {
    const auto& node = lamina::detail::node(expression);

    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (is_imaginary_unit_name(variable->name())) {
            auto unit = imaginary_unit();
            if (!unit) throw std::runtime_error(unit.error().message);
            return unit.value();
        }
        return lamina::detail::make_expression_ptr(node);
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            auto canonical_operand = canonicalize_lsr_complex_product(
                *lamina::detail::make_expression_ptr(operand));
            operands.push_back(lamina::detail::node(canonical_operand));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_add(std::move(operands)))->simplify();
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = exact_small_integer_node(power->exponent(), 0, 16);
        if (exponent) {
            auto canonical_base = canonicalize_lsr_complex_product(
                *lamina::detail::make_expression_ptr(power->base()));
            if (std::dynamic_pointer_cast<const ComplexNode>(
                    lamina::detail::node(canonical_base))) {
                auto result = SymbolicExpr::number(1);
                for (int i = 0; i < *exponent; ++i) {
                    result = canonicalize_lsr_complex_product(
                        *SymbolicExpr::multiply(result, canonical_base));
                }
                return result->simplify();
            }
        }
        return lamina::detail::make_expression_ptr(node);
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real_part = canonicalize_lsr_complex_product(
            *lamina::detail::make_expression_ptr(complex_node->real()));
        auto imag_part = canonicalize_lsr_complex_product(
            *lamina::detail::make_expression_ptr(complex_node->imag()));
        auto value = complex(real_part, imag_part);
        if (!value) throw std::runtime_error(value.error().message);
        return value.value()->simplify();
    }

    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) {
        return lamina::detail::make_expression_ptr(node);
    }

    bool saw_complex = false;
    auto real = SymbolicExpr::number(1);
    auto imag = SymbolicExpr::number(0);

    for (const auto& operand : multiply->operands()) {
        auto canonical_operand = canonicalize_lsr_complex_product(
            *lamina::detail::make_expression_ptr(operand));
        const auto& operand_node = lamina::detail::node(canonical_operand);
        ExprPtr factor_real;
        ExprPtr factor_imag;
        if (auto complex_operand = std::dynamic_pointer_cast<const ComplexNode>(
                operand_node)) {
            saw_complex = true;
            factor_real = lamina::detail::make_expression_ptr(complex_operand->real());
            factor_imag = lamina::detail::make_expression_ptr(complex_operand->imag());
        } else {
            factor_real = canonical_operand;
            factor_imag = SymbolicExpr::number(0);
        }

        auto ac = SymbolicExpr::multiply(real, factor_real);
        auto bd = SymbolicExpr::multiply(imag, factor_imag);
        auto ad = SymbolicExpr::multiply(real, factor_imag);
        auto bc = SymbolicExpr::multiply(imag, factor_real);
        auto next_real = SymbolicExpr::add(
            ac, SymbolicExpr::multiply(SymbolicExpr::number(-1), bd))->simplify();
        auto next_imag = SymbolicExpr::add(ad, bc)->simplify();
        real = std::move(next_real);
        imag = std::move(next_imag);
    }

    if (!saw_complex) {
        return lamina::detail::make_expression_ptr(node);
    }
    auto result = complex(real, imag);
    if (!result) {
        throw std::runtime_error(result.error().message);
    }
    return result.value()->simplify();
}

} // namespace

Result<EqvProfile> eqv_profile_from_name(const std::string& name) {
    if (name == "Core") {
        return Result<EqvProfile>::success(EqvProfile::Core);
    }
    if (name == "Trig-Basic") {
        return Result<EqvProfile>::success(EqvProfile::TrigBasic);
    }
    if (name == "ExpLog-Basic") {
        return Result<EqvProfile>::success(EqvProfile::ExpLogBasic);
    }
    return eqv_profile_failure("unsupported equivalence profile: " + name);
}

Result<void> set_eqv_profile(EqvOptions& options,
                             const std::string& name) {
    auto profile = eqv_profile_from_name(name);
    if (!profile) {
        return Result<void>::failure(profile.error());
    }
    options.profile = profile.value();
    return Result<void>::success();
}

Result<void> set_eqv_budget(EqvOptions& options,
                            std::size_t steps,
                            std::size_t depth,
                            std::size_t growth) {
    EqvOptions candidate = options;
    candidate.budget.max_rewrite_steps = steps;
    candidate.budget.max_rewrite_depth = depth;
    candidate.budget.max_node_growth_factor = growth;
    auto valid = validate_eqv_options(candidate);
    if (!valid) {
        return valid;
    }
    options = candidate;
    return Result<void>::success();
}

bool structurally_equal(const SymbolicExpr& lhs, const SymbolicExpr& rhs) {
    const auto& left = lamina::detail::node(lhs);
    const auto& right = lamina::detail::node(rhs);
    if (!left || !right) return left == right;
    return left->equals(*right);
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context,
                             const EqvOptions& options) {
    auto options_valid = validate_eqv_options(options);
    if (!options_valid) return Result<bool>::failure(options_valid.error());

    auto step = context.consume_steps(1, kEquivalentOperation);
    if (!step) return Result<bool>::failure(step.error());
    try {
        auto lhs_dimension = dimension_of(lhs);
        if (!lhs_dimension) return Result<bool>::failure(lhs_dimension.error());
        auto rhs_dimension = dimension_of(rhs);
        if (!rhs_dimension) return Result<bool>::failure(rhs_dimension.error());
        if (lhs_dimension.value() != rhs_dimension.value()) {
            return Result<bool>::success(false);
        }
        auto lhs_ptr = std::make_shared<SymbolicExpr>(lhs);
        auto rhs_ptr = std::make_shared<SymbolicExpr>(rhs);
        if (std::dynamic_pointer_cast<const QuantityNode>(lamina::detail::node(lhs))) {
            auto stripped = strip_unit(lhs_ptr, UnitStripMode::BaseValue, context);
            if (!stripped) return Result<bool>::failure(stripped.error());
            lhs_ptr = stripped.value();
        }
        if (std::dynamic_pointer_cast<const QuantityNode>(lamina::detail::node(rhs))) {
            auto stripped = strip_unit(rhs_ptr, UnitStripMode::BaseValue, context);
            if (!stripped) return Result<bool>::failure(stripped.error());
            rhs_ptr = stripped.value();
        }
        auto canonical_lhs = canonicalize_lsr_complex_product(*lhs_ptr);
        auto canonical_rhs = canonicalize_lsr_complex_product(*rhs_ptr);
        if (structurally_equal(*canonical_lhs, *canonical_rhs)) {
            return Result<bool>::success(true);
        }
        auto difference = SymbolicExpr::add(
            canonical_lhs,
            SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                   canonical_rhs));
        if (!difference) {
            return Result<bool>::failure(CasErrc::InternalInvariant,
                                         "equivalence difference construction failed",
                                         kEquivalentOperation);
        }
        if (difference->simplify()->is_zero()) {
            return Result<bool>::success(true);
        }

        auto polynomial_proof = prove_rational_polynomial_equivalence(
            difference, context, options);
        if (!polynomial_proof) {
            return Result<bool>::failure(polynomial_proof.error());
        }
        if (polynomial_proof.value()) {
            return Result<bool>::success(*polynomial_proof.value());
        }
        if (options.profile == EqvProfile::TrigBasic) {
            if (options.budget.max_rewrite_steps < 8) {
                return Result<bool>::failure(
                    CasErrc::ResourceLimit,
                    "equivalence rewrite budget exhausted before Trig-Basic normalization",
                    kEquivalentOperation);
            }
            auto trig_step = context.consume_steps(8, kEquivalentOperation);
            if (!trig_step) return Result<bool>::failure(trig_step.error());
            auto trig_lhs = rewrite_trig_basic_identity(lamina::detail::node(lhs));
            auto trig_rhs = rewrite_trig_basic_identity(lamina::detail::node(rhs));
            EqvOptions core_options = options;
            core_options.profile = EqvProfile::Core;
            return equivalent_core(*trig_lhs, *trig_rhs, context,
                                   core_options);
        }
        if (options.profile == EqvProfile::ExpLogBasic) {
            if (options.budget.max_rewrite_steps < 8) {
                return Result<bool>::failure(
                    CasErrc::ResourceLimit,
                    "equivalence rewrite budget exhausted before ExpLog-Basic normalization",
                    kEquivalentOperation);
            }
            auto exp_log_step = context.consume_steps(8, kEquivalentOperation);
            if (!exp_log_step) return Result<bool>::failure(exp_log_step.error());
            auto exp_log_lhs = rewrite_exp_log_basic_identity(
                lamina::detail::node(lhs), context.assumptions().get());
            auto exp_log_rhs = rewrite_exp_log_basic_identity(
                lamina::detail::node(rhs), context.assumptions().get());
            EqvOptions core_options = options;
            core_options.profile = EqvProfile::Core;
            return equivalent_core(*exp_log_lhs, *exp_log_rhs, context,
                                   core_options);
        }
        return Result<bool>::success(false);
    } catch (const std::bad_alloc&) {
        return Result<bool>::failure(CasErrc::ResourceLimit,
                                     "equivalence check allocation failed",
                                     kEquivalentOperation);
    } catch (const std::exception& error) {
        return Result<bool>::failure(CasErrc::Inconclusive, error.what(),
                                     kEquivalentOperation);
    }
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context) {
    return equivalent_core(lhs, rhs, context, EqvOptions{});
}

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context,
                        const EqvOptions& options) {
    auto checked = equivalent_core(lhs, rhs, context, options);
    if (checked) return checked;

    const auto code = checked.error().code;
    if (code == CasErrc::ResourceLimit ||
        code == CasErrc::Inconclusive ||
        code == CasErrc::UnsupportedExpression) {
        (void)context.add_diagnostic(
            Diagnostic{DiagnosticSeverity::Warning,
                       kEquivalentOperation,
                       checked.error().message});
        return Result<bool>::success(false);
    }
    return Result<bool>::failure(checked.error());
}

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context) {
    return equivalent(lhs, rhs, context, EqvOptions{});
}

} // namespace lamina::lsr
