#include "residual_verification.hpp"

#include "internal/equivalence_engine.hpp"
#include "internal/exact_matrix.hpp"
#include "internal/exact_matrix_support.hpp"
#include "internal/expression_analysis.hpp"
#include "poly_utils.hpp"
#include "irrational.hpp"
#include "multivariate_poly.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/expression_transform.hpp"
#include "internal/rewrite_budget.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/facts_query.hpp"

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
    detail::RewriteBudget& budget_;
    std::size_t canonical_nodes_ = 0;
    std::vector<Node> atoms_;
    std::vector<std::string> variables_;
    std::vector<std::optional<Rational>> radical_squares_;
    std::vector<Node> radical_radicands_;
    std::vector<std::pair<Node, Node>> canonical_atoms_;
    std::optional<CasError> failure_;
    bool inexact_ = false;
    struct Angle {
        Node base;
        Rational scale;
        bool hyperbolic;
    };
    std::vector<Angle> angles_;

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
            auto exponent = detail::exact_rational_value(power->exponent());
            if (exponent && *exponent == Rational(1, 2)) {
                return power->base();
            }
        }
        return nullptr;
    }
    std::optional<Rational> radical_square(const Node& node) {
        const auto radicand = radical_radicand(node);
        return radicand ? detail::exact_rational_value(radicand) : std::nullopt;
    }
    Node number(int value) {
        budget_.require_nodes(1);
        return SymbolicFactory::create_number(BigInt(value));
    }
    Node power(const Node& base, int exponent) {
        auto degree = number(exponent);
        normalization_check_children(&budget_, 1, base, degree);
        return detail::make_node<PowerNode>(base, std::move(degree));
    }
    Node function(FunctionNode::FuncType type, const Node& argument) {
        normalization_check_children(&budget_, 1, argument);
        return detail::make_node<FunctionNode>(type, std::vector<Node>{argument});
    }
    Node product(std::vector<Node> operands) {
        return SymbolicFactory::create_multiply(std::move(operands), &budget_);
    }
    Node sum(std::vector<Node> operands) {
        return SymbolicFactory::create_add(std::move(operands), &budget_);
    }
    Node square_relation(const Node& node) {
        if (auto radicand = radical_radicand(node)) { return radicand; }
        const auto* fn = dynamic_cast<const FunctionNode*>(node.get());
        if (!fn || fn->arguments().size() != 1) { return nullptr; }
        using FT = FunctionNode::FuncType;
        const auto& argument = fn->arguments()[0];
        if (fn->type() == FT::Cos) {
            return sum({number(1), product({
                number(-1), power(function(FT::Sin, argument), 2)})});
        }
        if (fn->type() == FT::Cosh) {
            return sum({number(1), power(function(FT::Sinh, argument), 2)});
        }
        if (fn->type() == FT::Abs) {
            std::optional<detail::AssumptionFacts> assumed;
            if (context_.assumptions()) { assumed.emplace(*context_.assumptions()); }
            const auto& facts = assumed
                ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
            auto real = detail::query_real_value(argument, facts, context_);
            if (!real) { failure_ = real.error(); return nullptr; }
            if (real.value() == Tribool::True) { return power(argument, 2); }
        }
        return nullptr;
    }
    bool register_replacement(const Node& node, Node replacement) {
        budget_.measure(replacement);
        if (contains_inexact_number(replacement) || !collect(replacement)) { return false; }
        if (!checked(context_.require_expansion_terms(
                atoms_.size() + canonical_atoms_.size() + 1, "residual.formal"))) {
            return false;
        }
        canonical_nodes_ = budget_.append_size(canonical_nodes_, budget_.measure(replacement));
        canonical_atoms_.emplace_back(node, std::move(replacement));
        return true;
    }
    std::pair<Rational, Node> coefficient(const Node& argument) {
        auto value = Rational(1);
        const auto* mul = dynamic_cast<const MultiplyNode*>(argument.get());
        if (!mul) { return {value, argument}; }
        std::vector<Node> factors;
        for (const auto& child : mul->operands()) {
            if (auto exact = detail::exact_rational_value(child)) { value = value * *exact; }
            else { factors.push_back(child); }
        }
        return {value, product(std::move(factors))};
    }
    void remember_angle(const FunctionNode& fn) {
        using FT = FunctionNode::FuncType;
        const bool hyperbolic = fn.type() == FT::Sinh || fn.type() == FT::Cosh ||
            fn.type() == FT::Tanh;
        if (fn.arguments().size() != 1 ||
            (!hyperbolic && fn.type() != FT::Sin && fn.type() != FT::Cos &&
             fn.type() != FT::Tan && fn.type() != FT::Cot &&
             fn.type() != FT::Sec && fn.type() != FT::Csc)) { return; }
        auto argument = detail::normalize_node(fn.arguments()[0], budget_);
        if (contains_inexact_number(argument)) { inexact_ = true; return; }
        auto [scale, base] = coefficient(argument);
        if (scale < Rational(0)) { scale = -scale; }
        if (scale.is_zero()) { return; }
        for (const auto& angle : angles_) {
            budget_.consume();
            if (angle.hyperbolic == hyperbolic && angle.scale == scale &&
                angle.base->equals(*base)) { return; }
        }
        if (!checked(context_.require_expansion_terms(angles_.size() + 1,
                                                     "residual.formal"))) { return; }
        canonical_nodes_ = budget_.append_size(canonical_nodes_, budget_.measure(base));
        angles_.push_back({std::move(base), std::move(scale), hyperbolic});
    }
    class AngleCollector final : public detail::RecursiveSymbolicVisitor {
        FormalResidual& owner_;
    public:
        explicit AngleCollector(FormalResidual& owner)
            : RecursiveSymbolicVisitor(&owner.budget_), owner_(owner) {}
        void visit(const FunctionNode& fn) override {
            owner_.remember_angle(fn);
            RecursiveSymbolicVisitor::visit(fn);
        }
    };
    std::int64_t double_angle_multiple(
        const Rational& scale, const Node& base, bool hyperbolic) {
        std::int64_t multiple = 1;
        for (const auto& angle : angles_) {
            budget_.consume();
            if (angle.hyperbolic != hyperbolic || !angle.base->equals(*base)) { continue; }
            auto ratio = scale / angle.scale;
            if (ratio.get_denominator() != BigInt(1)) { continue; }
            auto integer = ratio.get_numerator().try_to_int64();
            if (integer && *integer > multiple && *integer <= 1000 &&
                (*integer & (*integer - 1)) == 0) { multiple = *integer; }
        }
        return multiple;
    }
    Node rewrite_radical(const Node& node) {
        auto radicand = radical_radicand(node);
        if (!radicand || detail::exact_rational_value(radicand)) { return nullptr; }
        auto canonical = detail::expand_node(radicand, budget_);
        if (contains_inexact_number(canonical)) { inexact_ = true; return nullptr; }
        Node term = canonical;
        if (const auto* add = dynamic_cast<const AddNode*>(canonical.get())) {
            for (const auto& child : add->operands()) {
                if (!dynamic_cast<const NumberNode*>(child.get())) { term = child; break; }
            }
        }
        auto scale = coefficient(term).first;
        if (scale < Rational(0)) { scale = -scale; }
        if (scale.is_zero()) { return nullptr; }
        if (scale != Rational(1)) {
            canonical = detail::expand_node(product({
                SymbolicFactory::create_number(Rational(1) / scale), canonical}), budget_);
        }
        if (scale == Rational(1) && canonical->equals(*radicand)) { return nullptr; }
        auto half = SymbolicFactory::create_number(Rational(1, 2));
        normalization_check_children(&budget_, 1, canonical, half);
        auto radical = detail::make_node<PowerNode>(canonical, half);
        if (scale == Rational(1)) { return radical; }
        auto scalar = SymbolicFactory::create_number(scale);
        normalization_check_children(&budget_, 1, scalar, half);
        return product({detail::make_node<PowerNode>(scalar, half), radical});
    }
    Node rewrite_function(const FunctionNode& fn) {
        if (fn.arguments().size() != 1) { return nullptr; }
        using FT = FunctionNode::FuncType;
        const auto& argument = fn.arguments()[0];
        if (fn.type() == FT::Sqrt) {
            auto half = SymbolicFactory::create_number(Rational(1, 2));
            normalization_check_children(&budget_, 1, argument, half);
            return detail::make_node<PowerNode>(argument, std::move(half));
        }
        if (fn.type() == FT::Tan || fn.type() == FT::Cot ||
            fn.type() == FT::Sec || fn.type() == FT::Csc) {
            auto sine = function(FT::Sin, argument);
            auto cosine = function(FT::Cos, argument);
            if (fn.type() == FT::Tan) { return product({sine, power(cosine, -1)}); }
            if (fn.type() == FT::Cot) { return product({cosine, power(sine, -1)}); }
            return power(fn.type() == FT::Sec ? cosine : sine, -1);
        }
        if (fn.type() == FT::Tanh) {
            return product({function(FT::Sinh, argument),
                            power(function(FT::Cosh, argument), -1)});
        }
        if (fn.type() != FT::Sin && fn.type() != FT::Cos &&
            fn.type() != FT::Sinh && fn.type() != FT::Cosh &&
            fn.type() != FT::Exp) { return nullptr; }
        auto canonical = detail::normalize_node(argument, budget_);
        if (contains_inexact_number(canonical)) { inexact_ = true; return nullptr; }
        if (!canonical->equals(*argument)) { return function(fn.type(), canonical); }
        if (fn.type() == FT::Exp) {
            if (auto exact = detail::exact_rational_value(argument);
                exact && exact->is_zero()) { return number(1); }
            if (const auto* add = dynamic_cast<const AddNode*>(argument.get())) {
                std::vector<Node> factors;
                for (const auto& child : add->operands()) {
                    factors.push_back(function(FT::Exp, child));
                }
                return product(std::move(factors));
            }
        }
        auto [scale, base] = coefficient(argument);
        auto integer = scale.get_denominator() == BigInt(1)
            ? scale.get_numerator().try_to_int64() : std::optional<std::int64_t>{};
        if (fn.type() == FT::Exp) {
            if (integer && *integer != 1 && *integer >= -1000 && *integer <= 1000) {
                return power(function(FT::Exp, base), static_cast<int>(*integer));
            }
            if (scale < Rational(0)) {
                auto positive = product({SymbolicFactory::create_number(-scale), base});
                return power(function(FT::Exp, positive), -1);
            }
            return nullptr;
        }
        const bool odd = fn.type() == FT::Sin || fn.type() == FT::Sinh;
        if (scale < Rational(0)) {
            auto positive = product({SymbolicFactory::create_number(-scale), base});
            auto value = function(fn.type(), positive);
            return odd ? product({number(-1), value}) : value;
        }
        const bool hyperbolic = fn.type() == FT::Sinh || fn.type() == FT::Cosh;
        if (double_angle_multiple(scale, base, hyperbolic) == 1) { return nullptr; }
        auto half = product({SymbolicFactory::create_number(scale / Rational(2)), base});
        auto sine = function(hyperbolic ? FT::Sinh : FT::Sin, half);
        auto cosine = function(hyperbolic ? FT::Cosh : FT::Cos, half);
        if (odd) { return product({number(2), sine, cosine}); }
        return hyperbolic
            ? sum({power(cosine, 2), power(sine, 2)})
            : sum({power(cosine, 2), product({number(-1), power(sine, 2)})});
    }

    bool register_atom(const Node& node, bool canonicalize_radicals) {
        if (!dynamic_cast<const VariableNode*>(node.get()) &&
            !dynamic_cast<const FunctionNode*>(node.get()) &&
            !dynamic_cast<const PowerNode*>(node.get()) &&
            !dynamic_cast<const RootOfNode*>(node.get())) { return false; }
        for (const auto& atom : atoms_) { if (atom->equals(*node)) { return true; } }
        for (const auto& atom : canonical_atoms_)
            { if (atom.first->equals(*node)) { return true; } }
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
        auto canonical = detail::node(denominator == BigInt(1)
            ? root.value().to_symbolic()
            : (root.value() * Rational(BigInt(1), denominator)).to_symbolic());
        if (canonical->equals(*node)) {
            return std::nullopt;
        }
        if (!collect(canonical, false) ||
            !checked(context_.require_expansion_terms(
                atoms_.size() + canonical_atoms_.size() + 1,
                "residual.formal"))) {
            return false;
        }
        canonical_nodes_ = budget_.append_size(canonical_nodes_, budget_.measure(canonical));
        canonical_atoms_.emplace_back(node, std::move(canonical));
        return true;
    }

    bool register_regular_atom(
        const Node& node, std::optional<Rational> square) {
        const auto radicand = square_relation(node);
        if (radicand && !square && !collect(radicand)) {
            return false;
        }
        if (failure_) { return false; }
        if (radicand && !square) {
            canonical_nodes_ = budget_.append_size(canonical_nodes_, budget_.measure(radicand));
        }
        if (!checked(context_.require_expansion_terms(
                atoms_.size() + canonical_atoms_.size() + 1,
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
        if (const auto* fn = dynamic_cast<const FunctionNode*>(node.get())) {
            if (auto replacement = rewrite_function(*fn)) {
                return register_replacement(node, std::move(replacement));
            }
            if (inexact_) { return false; }
        }
        if (const auto* exponentiation = dynamic_cast<const PowerNode*>(node.get())) {
            auto exponent = detail::exact_rational_value(exponentiation->exponent());
            if (exponent && exponent->get_denominator() == BigInt(2) &&
                *exponent != Rational(1, 2)) {
                auto degree = exponent->get_numerator().try_to_int64();
                if (degree && *degree >= -1000 && *degree <= 1000) {
                    auto half = SymbolicFactory::create_number(Rational(1, 2));
                    normalization_check_children(&budget_, 1, exponentiation->base(), half);
                    auto radical = detail::make_node<PowerNode>(exponentiation->base(), half);
                    return register_replacement(node, power(radical, static_cast<int>(*degree)));
                }
            }
        }
        if (auto replacement = rewrite_radical(node)) {
            return register_replacement(node, std::move(replacement));
        }
        if (inexact_) { return false; }
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
        auto value = detail::exact_rational_value(node);
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
        for (const auto& atom : canonical_atoms_) {
            if (atom.first->equals(*node)) {
                return convert(atom.second);
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
        if (!checked(context_.enter_recursion("residual.formal"))) { return std::nullopt; }
        ResidualRecursionScope scope(context_);
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
    explicit FormalResidual(detail::RewriteBudget& budget)
        : context_(budget.context()), budget_(budget) {}
    Result<bool> prove(const Node& expression) {
        AngleCollector angles(*this);
        expression->accept(angles);
        if (failure_) { return Result<bool>::failure(*failure_); }
        if (inexact_) { return false; }
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
    detail::RewriteBudget& budget, const std::string& operation);

RootReplacementResult replace_root_operands(const detail::SymbolicNodePtr& node,
    const std::vector<detail::SymbolicNodePtr>& operands, bool sum,
    const detail::ExactRootId*& identity, detail::RewriteBudget& budget, const std::string& operation) {
    bool changed = false;
    std::vector<detail::SymbolicNodePtr> replaced;
    std::size_t nodes = 1;
    for (auto operand = operands.begin(); operand != operands.end(); ++operand) {
        auto child = replace_exact_root(*operand, identity, budget, operation);
        if (!child || !child.value()) { return child; }
        nodes = budget.append_size(nodes, budget.measure(*child.value()));
        if (!changed && *child.value() != *operand) {
            replaced.reserve(operands.size());
            replaced.insert(replaced.end(), operands.begin(), operand);
            changed = true;
        }
        if (changed) { replaced.push_back(std::move(*child.value())); }
    }
    if (!changed) { return RootReplacementResult::success(node); }
    return RootReplacementResult::success(sum
        ? SymbolicFactory::create_add(std::move(replaced), &budget)
        : SymbolicFactory::create_multiply(std::move(replaced), &budget));
}

RootReplacementResult replace_root_power(const detail::SymbolicNodePtr& node,
    const PowerNode& power, const detail::ExactRootId*& identity,
    detail::RewriteBudget& budget, const std::string& operation) {
    auto base = replace_exact_root(power.base(), identity, budget, operation);
    if (!base || !base.value()) { return base; }
    auto exponent = replace_exact_root(power.exponent(), identity, budget, operation);
    if (!exponent || !exponent.value()) { return exponent; }
    if (*base.value() == power.base() && *exponent.value() == power.exponent()) {
        return RootReplacementResult::success(node);
    }
    budget.append_size(budget.append_size(1, budget.measure(*base.value())),
                       budget.measure(*exponent.value()));
    return RootReplacementResult::success(
        SymbolicFactory::create_power(*base.value(), *exponent.value()));
}

RootReplacementResult replace_exact_root(
    const detail::SymbolicNodePtr& node,
    const detail::ExactRootId*& identity,
    detail::RewriteBudget& budget,
    const std::string& operation) {
    auto& context = budget.context();
    auto entered = context.enter_recursion(operation);
    if (!entered) { return RootReplacementResult::failure(entered.error()); }
    ResidualRecursionScope scope(context);
    if (const auto* root = dynamic_cast<const RootOfNode*>(node.get())) {
        if (identity && identity->compare(root->exact_id()) != 0) {
            return RootReplacementResult::success(std::nullopt);
        }
        identity = &root->exact_id();
        budget.require_nodes(1);
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
        return replace_root_operands(node, operands, sum != nullptr, identity, budget, operation);
    }
    if (const auto* power = dynamic_cast<const PowerNode*>(node.get())) {
        return replace_root_power(node, *power, identity, budget, operation);
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
    const ExprPtr& expression, detail::RewriteBudget& budget,
    const std::string& operation) {
    auto& context = budget.context();
    const detail::ExactRootId* identity = nullptr;
    auto replaced = replace_exact_root(
        detail::node(expression), identity, budget, operation);
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


ResidualCheckResult classify_residual(
    const ExprPtr& expression,
    detail::RewriteBudget& budget,
    const std::string& operation) {
    auto& context = budget.context();
    auto root_proof = prove_exact_root_residual(expression, budget, operation);
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
    const auto variables = free_variables(detail::node(expression), &budget);
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
        auto equivalent = detail::equivalent_impl(*current, *zero, context, profile_options);
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
    detail::RewriteBudget& budget, const EqvOptions& options, const std::string& operation) {
    auto& context = budget.context();
    std::optional<detail::AssumptionFacts> assumed;
    if (context.assumptions()) { assumed.emplace(*context.assumptions()); }
    const auto& facts = assumed
        ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
    auto current = detail::make_expression_ptr(
        detail::normalize_node(detail::node(residual), budget, facts));
    if (contains_inexact_number(detail::node(current))) {
        return ResidualCheckResult::success(UnprovedResidual{current});
    }
    auto first = classify_residual(current, budget, operation);
    if (!first || is_decided(first.value())) { return first; }
    current = detail::make_expression_ptr(
        detail::expand_node(detail::node(current), budget, facts));
    if (contains_inexact_number(detail::node(current))) {
        return ResidualCheckResult::success(UnprovedResidual{current});
    }
    auto expanded = classify_residual(current, budget, operation);
    if (!expanded || is_decided(expanded.value())) { return expanded; }
    /**
     * @brief 通分证明受计算预算约束，替代符号系数的 Euclid 及因子搜索。
     * 支持绑定整数名同时出现在参数分母中的解族。
     */
    auto formal = FormalResidual(budget).prove(detail::node(current));
    if (!formal) { return ResidualCheckResult::failure(formal.error()); }
    if (formal.value()) {
        return ResidualCheckResult::success(ProvedZeroResidual{ExactResidualProof{current}});
    }
    return prove_with_profiles(current, context, options);
}

const PiecewiseNode* residual_piecewise(
    const detail::SymbolicNodePtr& node, detail::RewriteBudget& budget) {
    detail::RewriteScope scope(&budget);
    if (const auto* piecewise = dynamic_cast<const PiecewiseNode*>(node.get())) {
        return piecewise;
    }
    const auto* add = dynamic_cast<const AddNode*>(node.get());
    const auto* mul = dynamic_cast<const MultiplyNode*>(node.get());
    if (add || mul) {
        for (const auto& child : add ? add->operands() : mul->operands()) {
            if (auto piecewise = residual_piecewise(child, budget)) { return piecewise; }
        }
    }
    if (const auto* power = dynamic_cast<const PowerNode*>(node.get())) {
        return residual_piecewise(power->base(), budget);
    }
    return nullptr;
}

class ResidualBranch final : public detail::SymbolicRewriter {
    const PiecewiseNode& selected_;
    std::size_t branch_;
public:
    ResidualBranch(const PiecewiseNode& selected, std::size_t branch,
                   detail::RewriteBudget& budget)
        : SymbolicRewriter(&budget), selected_(selected), branch_(branch) {}
    void visit(const PiecewiseNode& node) override {
        if (node.branches().size() != selected_.branches().size() ||
            static_cast<bool>(node.default_expr()) !=
                static_cast<bool>(selected_.default_expr())) {
            set_result(current());
            return;
        }
        for (std::size_t index = 0; index < node.branches().size(); ++index) {
            if (!node.branches()[index].condition->equals(
                    *selected_.branches()[index].condition)) {
                set_result(current());
                return;
            }
        }
        set_result(branch_ == node.branches().size()
            ? node.default_expr() : node.branches()[branch_].expression);
    }
};

std::optional<ResidualCheckResult> prove_piecewise_residual(
    const ExprPtr& residual, detail::RewriteBudget& budget,
    const EqvOptions& options) {
    const auto* piecewise = residual_piecewise(detail::node(residual), budget);
    if (!piecewise) { return std::nullopt; }
    auto& context = budget.context();
    auto entered = context.enter_recursion("residual.piecewise");
    if (!entered) { return ResidualCheckResult::failure(entered.error()); }
    ResidualRecursionScope scope(context);
    const auto count = piecewise->branches().size() +
        static_cast<std::size_t>(piecewise->default_expr() != nullptr);
    for (std::size_t index = 0; index < count; ++index) {
        ResidualBranch select(*piecewise, index, budget);
        auto branch = select.rewrite(detail::node(residual));
        budget.measure(branch);
        auto proof = check_zero_residual(detail::make_expression_ptr(branch), context, options);
        if (!proof) { return proof; }
        if (!std::holds_alternative<ProvedZeroResidual>(proof.value())) {
            return ResidualCheckResult::success(UnprovedResidual{residual});
        }
    }
    return ResidualCheckResult::success(
        ProvedZeroResidual{ExactResidualProof{residual}});
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
        detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
                                     context.limits().max_ast_nodes, operation);
        budget.measure(detail::node(residual));
        if (contains_inexact_number(detail::node(residual))) {
            return ResidualCheckResult::success(UnprovedResidual{residual});
        }
        if (auto piecewise = prove_piecewise_residual(residual, budget, options)) {
            return std::move(*piecewise);
        }
        auto root_proof = prove_exact_root_residual(residual, budget, operation);
        if (!root_proof ||
            std::holds_alternative<ProvedZeroResidual>(root_proof.value())) {
            return root_proof;
        }
        auto formal = FormalResidual(budget).prove(detail::node(residual));
        if (!formal) { return ResidualCheckResult::failure(formal.error()); }
        if (formal.value())
            { return ResidualCheckResult::success(ProvedZeroResidual{ExactResidualProof{residual}}); }
        return prove_normalized_residual(residual, budget, options, operation);
    } catch (const CasError& error) {
        return ResidualCheckResult::failure(error);
    } catch (const std::length_error& error) {
        return ResidualCheckResult::failure(CasErrc::ResourceLimit, error.what(), operation);
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
