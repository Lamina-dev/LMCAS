#include "residual_verification.hpp"

#include "internal/equivalence_engine.hpp"
#include "internal/exact_matrix.hpp"
#include "internal/exact_matrix_support.hpp"
#include "internal/expression_analysis.hpp"
#include "poly_utils.hpp"
#include "irrational.hpp"
#include "multivariate_poly.hpp"
#include "internal/normalization_utils.hpp"

#include <algorithm>
#include <exception>
#include <new>
#include <limits>
#include <optional>
#include <utility>


namespace LMCAS {
namespace {

using RootReplacementResult =
    Result<std::optional<detail::SymbolicNodePtr>>;

class ResidualRecursionScope {
public:
    explicit ResidualRecursionScope(ComputationContext& context)
        : context_(context) {}
    ~ResidualRecursionScope() { context_.leave_recursion(); }
private:
    ComputationContext& context_;
};

/**
 * @brief 在原公共定义域上证明形式有理恒等式。
 * 将不透明标量子式视为不定元，保留其点值可能为零的情况；
 * 交叉相乘免去因式分解及符号系数 GCD。
 */
class FormalResidual {
    struct Fraction { MultiPoly numerator; MultiPoly denominator; };
    using Node = detail::SymbolicNodePtr;
    ComputationContext& context_;
    std::vector<Node> atoms_;
    std::vector<std::string> variables_;
    std::vector<std::optional<Rational>> radical_squares_;
    std::vector<Node> radical_radicands_;
    std::vector<std::pair<Node, Node>> canonical_radicals_;
    std::optional<CasError> failure_;

    bool checked(Result<void> result) {
        if (result) return true;
        failure_ = result.error();
        return false;
    }
    std::optional<std::int64_t> integer_power(const PowerNode& power) {
        BigInt exponent;
        if (!try_get_integer_value(
                std::dynamic_pointer_cast<const NumberNode>(power.exponent()), exponent))
            return std::nullopt;
        return exponent.try_to_int64();
    }
    bool collect_operands(
        const std::vector<Node>& operands, bool canonicalize_radicals) {
        for (const auto& child : operands) {
            if (!collect(child, canonicalize_radicals)) {
                return false;
            }
        }
        return true;
    }
    std::optional<bool> collect_integer_power(
        const PowerNode& power, bool canonicalize_radicals) {
        const auto exponent = integer_power(power);
        if (!exponent) {
            return std::nullopt;
        }
        if (*exponent < -1000 || *exponent > 1000) {
            return false;
        }
        return collect(power.base(), canonicalize_radicals);
    }
    bool collect(const Node& node, bool canonicalize_radicals = true) {
        if (!checked(context_.consume_steps(1, "residual.formal"))) {
            return false;
        }
        if (!checked(context_.enter_recursion("residual.formal"))) {
            return false;
        }
        ResidualRecursionScope scope(context_);
        for (const auto& atom : atoms_) {
            if (atom->equals(*node)) {
                return true;
            }
        }
        if (auto number = dynamic_cast<const NumberNode*>(node.get())) {
            return !std::holds_alternative<lmmc_real_t>(number->value());
        }
        if (auto sum = dynamic_cast<const AddNode*>(node.get())) {
            return collect_operands(
                sum->operands(), canonicalize_radicals);
        }
        if (auto product = dynamic_cast<const MultiplyNode*>(node.get())) {
            return collect_operands(
                product->operands(), canonicalize_radicals);
        }
        if (auto power = dynamic_cast<const PowerNode*>(node.get())) {
            auto collected =
                collect_integer_power(*power, canonicalize_radicals);
            if (collected) {
                return *collected;
            }
        }
        return register_atom(node, canonicalize_radicals);
    }
    Node radical_radicand(const Node& node) {
        if (auto function = dynamic_cast<const FunctionNode*>(node.get())) {
            if (function->type() == FunctionNode::FuncType::Sqrt &&
                function->arguments().size() == 1) {
                return function->arguments()[0];
            }
        } else if (auto power = dynamic_cast<const PowerNode*>(node.get())) {
            auto exponent = detail::matrix_kernel::exact_rational(
                detail::make_expression_ptr(power->exponent()));
            if (exponent && *exponent == Rational(1, 2)) {
                return power->base();
            }
        }
        return nullptr;
    }
    std::optional<Rational> radical_square(const Node& node) {
        const auto radicand = radical_radicand(node);
        return radicand
            ? detail::matrix_kernel::exact_rational(
                  detail::make_expression_ptr(radicand))
            : std::nullopt;
    }

    bool register_atom(const Node& node, bool canonicalize_radicals) {
        if (!dynamic_cast<const VariableNode*>(node.get()) &&
            !dynamic_cast<const FunctionNode*>(node.get()) &&
            !dynamic_cast<const PowerNode*>(node.get()) &&
            !dynamic_cast<const RootOfNode*>(node.get())) { return false; }
        for (const auto& atom : atoms_) { if (atom->equals(*node)) { return true; } }
        for (const auto& radical : canonical_radicals_)
            { if (radical.first->equals(*node)) { return true; } }
        return register_new_atom(node, canonicalize_radicals);
    }
    std::optional<bool> register_canonical_radical(
        const Node& node, const Rational& square) {
        const auto& numerator = square.get_numerator();
        const auto& denominator = square.get_denominator();
        if (!checked(context_.require_integer_bits(
                numerator.bit_length() + denominator.bit_length() + 2,
                "residual.formal"))) {
            return false;
        }
        auto root = Irrational::sqrt_checked(
            numerator * denominator, 1.0, context_);
        if (!root) {
            failure_ = root.error();
            return false;
        }
        if (!checked(context_.reserve_nodes(5, "residual.formal"))) {
            return false;
        }
        auto canonical = detail::node(denominator == BigInt(1)
            ? root.value().to_symbolic()
            : (root.value() * Rational(BigInt(1), denominator)).to_symbolic());
        if (canonical->equals(*node)) {
            return std::nullopt;
        }
        if (!collect(canonical, false) ||
            !checked(context_.require_expansion_terms(
                atoms_.size() + canonical_radicals_.size() + 1,
                "residual.formal"))) {
            return false;
        }
        canonical_radicals_.emplace_back(node, std::move(canonical));
        return true;
    }

    bool register_regular_atom(
        const Node& node, std::optional<Rational> square) {
        const auto radicand = radical_radicand(node);
        if (radicand && !square && !collect(radicand)) {
            return false;
        }
        if (!checked(context_.require_expansion_terms(
                atoms_.size() + canonical_radicals_.size() + 1,
                "residual.formal"))) {
            return false;
        }
        atoms_.push_back(node);
        variables_.push_back(
            "_formal" + std::to_string(variables_.size()));
        radical_squares_.push_back(std::move(square));
        radical_radicands_.push_back(
            radical_squares_.back() ? nullptr : radicand);
        return true;
    }

    bool register_new_atom(const Node& node, bool canonicalize_radicals) {
        auto square = radical_square(node);
        if (canonicalize_radicals && square &&
            !square->get_numerator().is_negative()) {
            auto registered =
                register_canonical_radical(node, *square);
            if (registered) {
                return *registered;
            }
        }
        return register_regular_atom(node, std::move(square));
    }
    MultiPoly constant(const Rational& value) const { return MultiPoly(value, variables_); }
    std::optional<MultiPoly> multiply(const MultiPoly& a, const MultiPoly& b) {
        const auto count = static_cast<std::size_t>(a.num_terms()) *
            static_cast<std::size_t>(b.num_terms());
        if (!checked(context_.require_expansion_terms(count, "residual.formal")) ||
            !checked(context_.consume_steps(count + 1, "residual.formal"))) return std::nullopt;
        std::size_t a_bits = 0, b_bits = 0;
        std::int64_t a_degree = 0, b_degree = 0;
        for (const auto& term : a.terms()) {
            a_bits = std::max(a_bits, term.second.get_numerator().bit_length() +
                term.second.get_denominator().bit_length());
            std::int64_t degree = 0;
            for (int exponent : term.first) degree += exponent;
            a_degree = std::max(a_degree, degree);
        }
        for (const auto& term : b.terms()) {
            b_bits = std::max(b_bits, term.second.get_numerator().bit_length() +
                term.second.get_denominator().bit_length());
            std::int64_t degree = 0;
            for (int exponent : term.first) degree += exponent;
            b_degree = std::max(b_degree, degree);
        }
        if (a_degree + b_degree > std::numeric_limits<int>::max()) {
            checked(Result<void>::failure(CasErrc::ResourceLimit,
                "formal polynomial degree exceeds its representation", "residual.formal"));
            return std::nullopt;
        }
        if (!checked(context_.require_integer_bits(
                (a_bits + b_bits + 1) * (count + 1), "residual.formal"))) return std::nullopt;
        return a * b;
    }
    std::optional<Fraction> combine(Fraction a, Fraction b, bool sum) {
        auto denominator = multiply(a.denominator, b.denominator);
        if (!denominator) return std::nullopt;
        if (!sum) {
            auto numerator = multiply(a.numerator, b.numerator);
            if (!numerator) return std::nullopt;
            return Fraction{std::move(*numerator), std::move(*denominator)};
        }
        auto left = multiply(a.numerator, b.denominator);
        auto right = multiply(b.numerator, a.denominator);
        if (!left || !right) return std::nullopt;
        if (!checked(context_.require_expansion_terms(
                left->num_terms() + right->num_terms(), "residual.formal"))) return std::nullopt;
        std::size_t left_bits = 0, right_bits = 0;
        for (const auto& term : left->terms())
            left_bits = std::max(left_bits, term.second.get_numerator().bit_length() +
                term.second.get_denominator().bit_length());
        for (const auto& term : right->terms())
            right_bits = std::max(right_bits, term.second.get_numerator().bit_length() +
                term.second.get_denominator().bit_length());
        if (!checked(context_.require_integer_bits(
                left_bits + right_bits + 1, "residual.formal"))) return std::nullopt;
        return Fraction{*left + *right, std::move(*denominator)};
    }
    std::optional<Fraction> convert_integer_power(const PowerNode& power, std::int64_t n) {
        auto base = convert(power.base());
        if (!base) { return std::nullopt; }
        if (n < 0) { std::swap(base->numerator, base->denominator); n = -n; }
        std::optional<Fraction> result = Fraction{constant(Rational(1)), constant(Rational(1))};
        while (n) {
            if (n & 1) {
                result = combine(std::move(*result), *base, false);
                if (!result) { return std::nullopt; }
            }
            n >>= 1;
            if (n) {
                base = combine(*base, *base, false);
                if (!base) { return std::nullopt; }
            }
        }
        return result;
    }
    std::optional<Fraction> convert_number(const Node& node) {
        auto value = detail::matrix_kernel::exact_rational(detail::make_expression_ptr(node));
        if (!value) { return std::nullopt; }
        if (!checked(context_.require_integer_bits(
                value->get_numerator().bit_length() + value->get_denominator().bit_length(),
                "residual.formal"))) { return std::nullopt; }
        return Fraction{constant(*value), constant(Rational(1))};
    }
    std::optional<Fraction> convert_operands(const std::vector<Node>& operands, bool sum) {
        std::optional<Fraction> result = Fraction{
            constant(Rational(sum ? 0 : 1)), constant(Rational(1))};
        for (const auto& child : operands) {
            auto operand = convert(child);
            if (!operand) { return std::nullopt; }
            result = combine(std::move(*result), std::move(*operand), sum);
            if (!result) { return std::nullopt; }
        }
        return result;
    }
    std::optional<Fraction> convert(const Node& node) {
        if (!checked(context_.enter_recursion("residual.formal"))) { return std::nullopt; }
        ResidualRecursionScope scope(context_);
        if (dynamic_cast<const NumberNode*>(node.get())) {
            return convert_number(node);
        }
        for (const auto& radical : canonical_radicals_) {
            if (radical.first->equals(*node)) {
                return convert(radical.second);
            }
        }
        for (std::size_t i = 0; i < atoms_.size(); ++i) {
            if (!atoms_[i]->equals(*node)) { continue; }
            Monomial powers(atoms_.size(), 0);
            powers[i] = 1;
            return Fraction{
                MultiPoly({{std::move(powers), Rational(1)}}, variables_),
                constant(Rational(1))};
        }
        auto sum = dynamic_cast<const AddNode*>(node.get());
        auto product = dynamic_cast<const MultiplyNode*>(node.get());
        if (sum || product) {
            return convert_operands(
                sum ? sum->operands() : product->operands(), sum != nullptr);
        }
        if (auto power = dynamic_cast<const PowerNode*>(node.get())) {
            if (auto exponent = integer_power(*power)) {
                return convert_integer_power(*power, *exponent);
            }
        }
        return std::nullopt;
    }
    std::optional<MultiPoly> reduce_radicals(const MultiPoly& polynomial) {
        MultiPoly reduced(Rational(0), variables_);
        for (auto term : polynomial.terms()) {
            std::vector<std::pair<Node, int>> symbolic_squares;
            for (std::size_t i = 0; i < radical_squares_.size(); ++i) {
                if (term.first[i] < 2) { continue; }
                const auto exponent = term.first[i] / 2;
                if (radical_squares_[i]) {
                    const auto& square = *radical_squares_[i];
                    auto bits = term.second.get_numerator().bit_length() +
                        term.second.get_denominator().bit_length() +
                        static_cast<std::size_t>(exponent) *
                        (square.get_numerator().bit_length() +
                         square.get_denominator().bit_length());
                    if (!checked(context_.require_integer_bits(
                            bits, "residual.formal"))) {
                        return std::nullopt;
                    }
                    term.second =
                        term.second * square.power(BigInt(exponent));
                } else if (radical_radicands_[i]) {
                    symbolic_squares.emplace_back(
                        radical_radicands_[i], exponent);
                } else {
                    continue;
                }
                if (!checked(context_.consume_steps(
                        static_cast<std::size_t>(exponent) + 1,
                        "residual.formal"))) {
                    return std::nullopt;
                }
                term.first[i] %= 2;
            }
            MultiPoly expanded({term}, variables_);
            for (const auto& [radicand, exponent] : symbolic_squares) {
                auto square = convert(radicand);
                if (!square ||
                    square->denominator != constant(Rational(1))) {
                    return std::nullopt;
                }
                for (int power = 0; power < exponent; ++power) {
                    auto product = multiply(expanded, square->numerator);
                    if (!product) { return std::nullopt; }
                    expanded = std::move(*product);
                }
            }
            if (!checked(context_.require_expansion_terms(
                    static_cast<std::size_t>(reduced.num_terms()) +
                    static_cast<std::size_t>(expanded.num_terms()),
                    "residual.formal"))) {
                return std::nullopt;
            }
            reduced = reduced + expanded;
        }
        return reduced;
    }
public:
    explicit FormalResidual(ComputationContext& context) : context_(context) {}
    Result<bool> prove(const Node& expression) {
        if (!collect(expression))
            { return failure_ ? Result<bool>::failure(*failure_) : Result<bool>::success(false); }
        auto fraction = convert(expression);
        if (failure_) { return Result<bool>::failure(*failure_); }
        if (!fraction) { return false; }
        const bool has_radicals = std::any_of(
            radical_squares_.begin(), radical_squares_.end(),
            [](const auto& value) { return value.has_value(); }) ||
            std::any_of(
                radical_radicands_.begin(), radical_radicands_.end(),
                [](const auto& value) { return value != nullptr; });
        if (!has_radicals) {
            return fraction->numerator.is_zero() &&
                   !fraction->denominator.is_zero();
        }
        auto numerator = reduce_radicals(fraction->numerator);
        auto denominator = reduce_radicals(fraction->denominator);
        if (failure_) { return Result<bool>::failure(*failure_); }
        return numerator && denominator && numerator->is_zero() && !denominator->is_zero();
    }
};


RootReplacementResult replace_exact_root(
    const detail::SymbolicNodePtr& node, const detail::ExactRootId*& identity,
    ComputationContext& context, const std::string& operation);

RootReplacementResult replace_root_operands(const detail::SymbolicNodePtr& node,
    const std::vector<detail::SymbolicNodePtr>& operands, bool sum,
    const detail::ExactRootId*& identity, ComputationContext& context, const std::string& operation) {
    bool changed = false;
    std::vector<detail::SymbolicNodePtr> replaced;
    for (auto operand = operands.begin(); operand != operands.end(); ++operand) {
        auto child = replace_exact_root(*operand, identity, context, operation);
        if (!child || !child.value()) { return child; }
        if (!changed && *child.value() != *operand) {
            auto reserved = context.reserve_nodes(1, operation);
            if (!reserved) { return RootReplacementResult::failure(reserved.error()); }
            replaced.reserve(operands.size());
            replaced.insert(replaced.end(), operands.begin(), operand);
            changed = true;
        }
        if (changed) { replaced.push_back(std::move(*child.value())); }
    }
    if (!changed) { return RootReplacementResult::success(node); }
    return RootReplacementResult::success(sum
        ? SymbolicFactory::create_add(std::move(replaced))
        : SymbolicFactory::create_multiply(std::move(replaced)));
}

RootReplacementResult replace_root_power(const detail::SymbolicNodePtr& node,
    const PowerNode& power, const detail::ExactRootId*& identity,
    ComputationContext& context, const std::string& operation) {
    auto base = replace_exact_root(power.base(), identity, context, operation);
    if (!base || !base.value()) { return base; }
    auto exponent = replace_exact_root(power.exponent(), identity, context, operation);
    if (!exponent || !exponent.value()) { return exponent; }
    if (*base.value() == power.base() && *exponent.value() == power.exponent()) {
        return RootReplacementResult::success(node);
    }
    auto reserved = context.reserve_nodes(1, operation);
    if (!reserved) { return RootReplacementResult::failure(reserved.error()); }
    return RootReplacementResult::success(
        SymbolicFactory::create_power(*base.value(), *exponent.value()));
}

RootReplacementResult replace_exact_root(
    const detail::SymbolicNodePtr& node,
    const detail::ExactRootId*& identity,
    ComputationContext& context,
    const std::string& operation) {
    auto entered = context.enter_recursion(operation);
    if (!entered) { return RootReplacementResult::failure(entered.error()); }
    ResidualRecursionScope scope(context);
    if (const auto* root = dynamic_cast<const RootOfNode*>(node.get())) {
        if (identity && identity->compare(root->exact_id()) != 0) {
            return RootReplacementResult::success(std::nullopt);
        }
        identity = &root->exact_id();
        auto reserved = context.reserve_nodes(1, operation);
        if (!reserved) { return RootReplacementResult::failure(reserved.error()); }
        return RootReplacementResult::success(
            SymbolicFactory::create_variable("_residual_root"));
    }
    if (dynamic_cast<const NumberNode*>(node.get())) {
        return RootReplacementResult::success(node);
    }
    const auto* sum = dynamic_cast<const AddNode*>(node.get());
    const auto* product = dynamic_cast<const MultiplyNode*>(node.get());
    if (sum || product) {
        const auto& operands = sum ? sum->operands() : product->operands();
        return replace_root_operands(node, operands, sum != nullptr, identity, context, operation);
    }
    if (const auto* power = dynamic_cast<const PowerNode*>(node.get())) {
        return replace_root_power(node, *power, identity, context, operation);
    }
    return RootReplacementResult::success(std::nullopt);
}

Result<void> check_residual_coefficient(
    const Rational& value, ComputationContext& context,
    const std::string& operation) {
    return context.require_integer_bits(std::max(
        value.get_numerator().bit_length(),
        value.get_denominator().bit_length()), operation);
}

Result<bool> has_zero_root_remainder(
    Polynomial<Rational> remainder,
    const Polynomial<Rational>& defining,
    ComputationContext& context,
    const std::string& operation) {
    auto terms = context.require_expansion_terms(
        std::max(remainder.coeffs.size(), defining.coeffs.size()), operation);
    if (!terms) return Result<bool>::failure(terms.error());
    for (const auto& coefficient : remainder.coeffs) {
        auto checked = check_residual_coefficient(coefficient, context, operation);
        if (!checked) return Result<bool>::failure(checked.error());
    }
    for (const auto& coefficient : defining.coeffs) {
        auto checked = check_residual_coefficient(coefficient, context, operation);
        if (!checked) return Result<bool>::failure(checked.error());
    }
    while (remainder.degree() >= defining.degree()) {
        const auto offset = remainder.coeffs.size() - defining.coeffs.size();
        const auto& factor = remainder.coeffs.back();
        for (std::size_t index = 0; index + 1 < defining.coeffs.size(); ++index) {
            auto step = context.consume_steps(1, operation);
            if (!step) return Result<bool>::failure(step.error());
            auto product = factor * defining.coeffs[index];
            auto checked = check_residual_coefficient(product, context, operation);
            if (!checked) return Result<bool>::failure(checked.error());
            auto& coefficient = remainder.coeffs[offset + index];
            coefficient = coefficient - product;
            checked = check_residual_coefficient(coefficient, context, operation);
            if (!checked) return Result<bool>::failure(checked.error());
        }
        remainder.coeffs.pop_back();
        remainder.trim();
    }
    return Result<bool>::success(remainder.is_zero());
}

ResidualCheckResult prove_exact_root_residual(
    const ExprPtr& expression, ComputationContext& context,
    const std::string& operation) {
    const detail::ExactRootId* identity = nullptr;
    auto replaced = replace_exact_root(
        detail::node(expression), identity, context, operation);
    if (!replaced) return ResidualCheckResult::failure(replaced.error());
    if (!replaced.value() || !identity) {
        return ResidualCheckResult::success(UnprovedResidual{expression});
    }
    auto polynomial = recognize_rational_polynomial(
        *detail::make_expression_ptr(*replaced.value()),
        "_residual_root", context);
    if (!polynomial) return ResidualCheckResult::failure(polynomial.error());
    if (polynomial.value()) {
        auto zero = has_zero_root_remainder(
            std::move(*polynomial.value()), identity->polynomial, context, operation);
        if (!zero) return ResidualCheckResult::failure(zero.error());
        if (zero.value()) {
            return ResidualCheckResult::success(ProvedZeroResidual{
                ExactResidualProof{expression}});
        }
    }
    return ResidualCheckResult::success(UnprovedResidual{expression});
}

ExprPtr normalized(ExprPtr expression) {
    if (!expression) { return expression; }
    auto simplified = expression->simplify();
    return simplified ? simplified : expression;
}

ResidualCheckResult classify_residual(
    const ExprPtr& expression,
    ComputationContext& context,
    const std::string& operation) {
    auto root_proof = prove_exact_root_residual(expression, context, operation);
    if (!root_proof ||
        std::holds_alternative<ProvedZeroResidual>(root_proof.value())) {
        return root_proof;
    }
    auto proof = detail::classify_exact_zero(
        expression, context, operation);
    if (!proof) { return ResidualCheckResult::failure(proof.error()); }
    if (proof.value() == detail::ZeroProof::Zero) {
        return ResidualCheckResult::success(ProvedZeroResidual{
            ExactNormalizationProof{expression}});
    }
    const auto variables = free_variables(detail::node(expression));
    if (variables.size() == 1) {
        auto polynomial = recognize_rational_polynomial(
            *expression, *variables.begin(), context);
        if (!polynomial) {
            return ResidualCheckResult::failure(polynomial.error());
        }
        if (polynomial.value() && !polynomial.value()->is_zero()) {
            return ResidualCheckResult::success(
                ProvedNonIdentityResidual{expression});
        }
    }
    if (proof.value() == detail::ZeroProof::NonZero) {
        return ResidualCheckResult::success(
            ProvedNonIdentityResidual{expression});
    }
    return ResidualCheckResult::success(UnprovedResidual{expression});
}

bool is_decided(const ResidualCheck& result) {
    return !std::holds_alternative<UnprovedResidual>(result);
}

struct ResidualProfiles {
    bool trig = false;
    bool exp_log = false;
};

void collect_residual_profiles(
    const detail::SymbolicNodePtr& node, ResidualProfiles& profiles);

void collect_residual_profile_children(
    const std::vector<detail::SymbolicNodePtr>& children,
    ResidualProfiles& profiles) {
    for (const auto& child : children) {
        collect_residual_profiles(child, profiles);
    }
}

void collect_function_profile(
    const FunctionNode& function, ResidualProfiles& profiles) {
    profiles.trig = profiles.trig ||
        function.type() == FunctionNode::FuncType::Sin ||
        function.type() == FunctionNode::FuncType::Cos;
    profiles.exp_log = profiles.exp_log ||
        function.type() == FunctionNode::FuncType::Exp ||
        function.type() == FunctionNode::FuncType::Ln;
    collect_residual_profile_children(function.arguments(), profiles);
}

void collect_residual_profiles(
    const detail::SymbolicNodePtr& node, ResidualProfiles& profiles) {
    if (!node || (profiles.trig && profiles.exp_log)) {
        return;
    }
    if (const auto function =
            std::dynamic_pointer_cast<const FunctionNode>(node)) {
        collect_function_profile(*function, profiles);
        return;
    }
    if (const auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        collect_residual_profile_children(add->operands(), profiles);
        return;
    }
    if (const auto multiply =
            std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        collect_residual_profile_children(multiply->operands(), profiles);
        return;
    }
    if (const auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        collect_residual_profiles(power->base(), profiles);
        collect_residual_profiles(power->exponent(), profiles);
        return;
    }
    if (const auto complex =
            std::dynamic_pointer_cast<const ComplexNode>(node)) {
        collect_residual_profiles(complex->real(), profiles);
        collect_residual_profiles(complex->imag(), profiles);
    }
}

ResidualCheckResult prove_with_profiles(
    const ExprPtr& current, ComputationContext& context,
    const EqvOptions& options) {
    ResidualProfiles required;
    collect_residual_profiles(detail::node(current), required);
    const auto zero = SymbolicExpr::number(0);
    for (const auto profile : {
             EqvProfile::Core, EqvProfile::TrigBasic, EqvProfile::ExpLogBasic}) {
        if ((profile == EqvProfile::TrigBasic && !required.trig) ||
            (profile == EqvProfile::ExpLogBasic && !required.exp_log)) {
            continue;
        }
        auto profile_options = options;
        profile_options.profile = profile;
        auto equivalent = equivalent_core(*current, *zero, context, profile_options);
        if (!equivalent) {
            if (equivalent.error().code == CasErrc::UnsupportedExpression ||
                equivalent.error().code == CasErrc::Inconclusive) {
                continue;
            }
            return ResidualCheckResult::failure(equivalent.error());
        }
        if (equivalent.value()) {
            return ResidualCheckResult::success(ProvedZeroResidual{
                ExactResidualProof{current}});
        }
    }
    return ResidualCheckResult::success(UnprovedResidual{current});
}

ResidualCheckResult prove_normalized_residual(const ExprPtr& residual,
    ComputationContext& context, const EqvOptions& options, const std::string& operation) {
    auto current = normalized(residual);
    auto first = classify_residual(current, context, operation);
    if (!first || is_decided(first.value())) { return first; }
    current = normalized(current->expand());
    auto expanded = classify_residual(current, context, operation);
    if (!expanded || is_decided(expanded.value())) { return expanded; }
    /**
     * @brief 通分证明受计算预算约束，替代符号系数的 Euclid 及因子搜索。
     * 支持绑定整数名同时出现在参数分母中的解族。
     */
    auto formal = FormalResidual(context).prove(detail::node(current));
    if (!formal) { return ResidualCheckResult::failure(formal.error()); }
    if (formal.value()) {
        return ResidualCheckResult::success(ProvedZeroResidual{ExactResidualProof{current}});
    }
    return prove_with_profiles(current, context, options);
}

}

ResidualCheckResult check_zero_residual(
    const ExprPtr& residual,
    ComputationContext& context,
    const LMCAS::EqvOptions& options) {
    constexpr const char* operation = "residual.check_zero";
    if (!residual || !detail::node(residual)) {
        return ResidualCheckResult::failure(
            CasErrc::InvalidArgument,
            "residual expression cannot be null", operation);
    }
    auto access = context.consume_steps(1, operation);
    if (!access) { return ResidualCheckResult::failure(access.error()); }

    try {
        auto root_proof = prove_exact_root_residual(residual, context, operation);
        if (!root_proof ||
            std::holds_alternative<ProvedZeroResidual>(root_proof.value())) {
            return root_proof;
        }
        auto formal = FormalResidual(context).prove(detail::node(residual));
        if (!formal) { return ResidualCheckResult::failure(formal.error()); }
        if (formal.value())
            { return ResidualCheckResult::success(ProvedZeroResidual{ExactResidualProof{residual}}); }
        return prove_normalized_residual(residual, context, options, operation);
    } catch (const CasError& error) {
        return ResidualCheckResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ResidualCheckResult::failure(
            CasErrc::ResourceLimit,
            "residual verification allocation failed", operation);
    } catch (const std::exception& error) {
        return ResidualCheckResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}


}
