#include "internal/integration_definite_support.hpp"
#include "internal/exact_sturm.hpp"
#include "polynomial_conversion.hpp"
#include "root_of_identity.hpp"

#include <limits>

namespace LMCAS::detail::definite {
static int sign(const Rational& value) {
    return value < Rational(0) ? -1 : value > Rational(0) ? 1 : 0;
}

std::optional<Rational> exact_number(const Node& node, bool allow_binary) {
    const auto* number = dynamic_cast<const NumberNode*>(node.get());
    if (!number) return std::nullopt;
    if (const auto* n = std::get_if<BigInt>(&number->value())) return Rational(*n);
    if (const auto* n = std::get_if<Rational>(&number->value())) return *n;
    if (allow_binary) {
        const auto n = std::get<lmmc_real_t>(number->value());
        if (std::isfinite(n)) return Rational::from_double(n);
    }
    return std::nullopt;
}

/**
 * @brief 识别绝对值节点及公共接口的平方开根表示。
 * @pre 调用方须先证明参数为实数，再在开区间上使用 |u| = ±u。
 */
Node absolute_argument(const SymbolicNode& node) {
    Node radicand;
    if (const auto* function = dynamic_cast<const FunctionNode*>(&node)) {
        if (function->arguments().size() != 1) return nullptr;
        if (function->type() == FunctionNode::FuncType::Abs) return function->arguments()[0];
        if (function->type() != FunctionNode::FuncType::Sqrt) return nullptr;
        radicand = function->arguments()[0];
    } else if (const auto* power = dynamic_cast<const PowerNode*>(&node)) {
        auto exponent = exact_number(power->exponent(), true);
        if (!exponent || *exponent != Rational(1, 2)) return nullptr;
        radicand = power->base();
    } else return nullptr;
    const auto* square = dynamic_cast<const PowerNode*>(radicand.get());
    if (!square) return nullptr;
    auto exponent = exact_number(square->exponent(), true);
    return exponent && *exponent == Rational(2) ? square->base() : nullptr;
}

int infinity_sign(const Node& node) {
    if (const auto* f = dynamic_cast<const FunctionNode*>(node.get())) {
        if (f->type() == FunctionNode::FuncType::Infinity) return 1;
    }
    if (const auto* m = dynamic_cast<const MultiplyNode*>(node.get())) {
        int result = 1;
        bool infinity = false;
        for (const auto& factor : m->operands()) {
            if (infinity_sign(factor)) {
                if (infinity) return 0;
                infinity = true;
                result *= infinity_sign(factor);
            } else {
                auto value = exact_number(factor, true);
                if (!value || sign(*value) == 0) return 0;
                result *= sign(*value);
            }
        }
        return infinity ? result : 0;
    }
    return 0;
}

Result<Boundary> boundary_checked(const SymbolicExpr& expression,
    ComputationContext& context) {
    if (int infinity = infinity_sign(detail::node(expression)))
        return Boundary{infinity, {}};
    if (auto rational = exact_number(detail::node(expression), true))
        return Boundary{0, {Poly({-*rational, Rational(1)}), *rational, *rational, 0, 1}};
    if (dynamic_cast<const NumberNode*>(detail::node(expression).get()))
        return Result<Boundary>::failure(CasErrc::DomainError,
            "nonfinite numeric endpoint requires a typed infinity", operation);
    if (const auto* root = dynamic_cast<const RootOfNode*>(detail::node(expression).get())) {
        auto value = detail::make_exact_real_algebraic(
            root->exact_id().polynomial, root->index(), 1, context);
        if (!value) return Result<Boundary>::failure(value.error());
        return Boundary{0, std::move(value.value())};
    }
    return undecided<Boundary>();
}

Result<int> compare_boundaries(const Boundary& a, const Boundary& b,
    ComputationContext& context) {
    if (a.infinity || b.infinity)
        { return a.infinity < b.infinity ? -1 : a.infinity > b.infinity ? 1 : 0; }
    auto left = a.finite;
    auto right = b.finite;
    right.polynomial.variable_name = left.polynomial.variable_name;
    if (!(left.upper < right.lower) && !(right.upper < left.lower)) {
        auto equal = detail::equal_exact_real_algebraic(left, right, context);
        if (!equal) { return Result<int>::failure(equal.error()); }
        if (equal.value()) { return 0; }
        while (!(left.upper < right.lower) && !(right.upper < left.lower)) {
            auto refined = detail::refine_exact_real_algebraic(left, context, operation);
            if (!refined) { return Result<int>::failure(refined.error()); }
            refined = detail::refine_exact_real_algebraic(right, context, operation);
            if (!refined) { return Result<int>::failure(refined.error()); }
        }
    }
    return detail::compare_exact_real_algebraic(std::move(left), std::move(right), context);
}

static Result<void> polynomial_budget(const Poly& polynomial, ComputationContext& context) {
    auto terms = context.require_expansion_terms(polynomial.coeffs.size(), operation);
    if (!terms) return terms;
    for (const auto& coefficient : polynomial.coeffs) {
        auto checked = context.require_integer_bits(std::max(
            coefficient.get_numerator().bit_length(),
            coefficient.get_denominator().bit_length()), operation);
        if (!checked) return checked;
    }
    return context.consume_steps(polynomial.coeffs.size() + 1, operation);
}

Result<Poly> product_checked(const Poly& a, const Poly& b, ComputationContext& context) {
    const auto size = a.coeffs.size() + b.coeffs.size();
    auto budget = context.require_expansion_terms(size, operation);
    if (!budget) return Result<Poly>::failure(budget.error());
    if (size > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        return Result<Poly>::failure(CasErrc::ResourceLimit,
            "polynomial degree exceeds representation", operation);
    auto step = context.consume_steps(a.coeffs.size() * b.coeffs.size() + 1, operation);
    if (!step) return Result<Poly>::failure(step.error());
    auto value = a * b;
    budget = polynomial_budget(value, context);
    if (!budget) return Result<Poly>::failure(budget.error());
    return value;
}

static Result<Poly> gcd_checked(Poly a, Poly b, ComputationContext& context) {
    while (!b.is_zero()) {
        auto step = context.consume_steps(a.coeffs.size() * b.coeffs.size() + 1, operation);
        if (!step) return Result<Poly>::failure(step.error());
        auto remainder = a.div_mod(b).second;
        auto checked = polynomial_budget(remainder, context);
        if (!checked) return Result<Poly>::failure(checked.error());
        a = std::move(b);
        b = std::move(remainder);
    }
    return a.make_monic();
}

static RationalResult rational_operands_checked(const std::vector<Node>& operands,
    bool sum, const std::string& var, ComputationContext& context) {
    RationalFunction result{Poly(Rational(sum ? 0 : 1), var), Poly(Rational(1), var)};
    for (const auto& operand : operands) {
        auto child = rational_checked(operand, var, context);
        if (!child || !child.value()) { return child; }
        auto n = product_checked(result.numerator,
            sum ? child.value()->denominator : child.value()->numerator, context);
        if (!n) { return RationalResult::failure(n.error()); }
        if (sum) {
            auto other = product_checked(child.value()->numerator, result.denominator, context);
            if (!other) { return RationalResult::failure(other.error()); }
            n.value() = n.value() + other.value();
        }
        auto d = product_checked(result.denominator, child.value()->denominator, context);
        if (!d) { return RationalResult::failure(d.error()); }
        result = {std::move(n.value()), std::move(d.value())};
    }
    return std::optional<RationalFunction>{std::move(result)};
}

static RationalResult rational_power_checked(const PowerNode& power,
    const std::string& var, ComputationContext& context) {
    auto exponent = exact_number(power.exponent());
    if (!exponent || !exponent->is_integer()) { return std::optional<RationalFunction>{}; }
    auto base = rational_checked(power.base(), var, context);
    if (!base || !base.value()) { return base; }
    BigInt count = exponent->to_bigint();
    if (count < BigInt(0)) {
        std::swap(base.value()->numerator, base.value()->denominator);
        count = -count;
    }
    RationalFunction result{Poly(Rational(1), var), Poly(Rational(1), var)};
    while (count != BigInt(0)) {
        if (count % BigInt(2) != BigInt(0)) {
            auto n = product_checked(result.numerator, base.value()->numerator, context);
            auto d = product_checked(result.denominator, base.value()->denominator, context);
            if (!n) { return RationalResult::failure(n.error()); }
            if (!d) { return RationalResult::failure(d.error()); }
            result = {std::move(n.value()), std::move(d.value())};
        }
        count = count / BigInt(2);
        if (count == BigInt(0)) { break; }
        auto n = product_checked(base.value()->numerator, base.value()->numerator, context);
        auto d = product_checked(base.value()->denominator, base.value()->denominator, context);
        if (!n) { return RationalResult::failure(n.error()); }
        if (!d) { return RationalResult::failure(d.error()); }
        base.value() = RationalFunction{std::move(n.value()), std::move(d.value())};
    }
    return std::optional<RationalFunction>{std::move(result)};
}
RationalResult rational_checked(const Node& node, const std::string& var,
    ComputationContext& context) {
    auto entered = context.enter_recursion(operation);
    if (!entered) { return RationalResult::failure(entered.error()); }
    struct Leave { ComputationContext& context; ~Leave() { context.leave_recursion(); } } leave{context};
    auto polynomial = recognize_rational_polynomial(detail::expression_from_node(node), var, context);
    if (!polynomial) { return RationalResult::failure(polynomial.error()); }
    if (polynomial.value()) { return std::optional<RationalFunction>{
        RationalFunction{std::move(*polynomial.value()), Poly(Rational(1), var)}}; }
    const auto* sum = dynamic_cast<const AddNode*>(node.get());
    const auto* product = dynamic_cast<const MultiplyNode*>(node.get());
    if (sum || product) {
        return rational_operands_checked(
            sum ? sum->operands() : product->operands(), sum != nullptr, var, context);
    }
    if (const auto* power = dynamic_cast<const PowerNode*>(node.get())) {
        return rational_power_checked(*power, var, context);
    }
    return std::optional<RationalFunction>{};
}

static Result<void> recover_rational_cut(Algebraic& root, const BigInt& denominator_bound,
    ComputationContext& context) {
    if (root.is_rational()) { return Result<void>::success(); }
    auto bits = context.require_integer_bits(2 * denominator_bound.bit_length() + 2, operation);
    if (!bits) { return bits; }
    const Rational separation(BigInt(1), BigInt(2) * denominator_bound * denominator_bound);
    while (root.upper - root.lower >= separation) {
        auto refined = detail::refine_exact_real_algebraic(root, context, operation);
        if (!refined) { return refined; }
        if (root.is_rational()) { return Result<void>::success(); }
    }
    /**
     * @brief 有理根定理给出分母 q 的上界；当前包围宽度保证根 p/q 位于中点连分数的渐近分数中。
     * 通过连分数恢复有理分割点，免去整数因子枚举。
     */
    const Rational midpoint = (root.lower + root.upper) / Rational(2);
    BigInt numerator = midpoint.get_numerator(), denominator = midpoint.get_denominator();
    BigInt old_p(0), p(1), old_q(1), q(0);
    while (denominator != BigInt(0)) {
        auto step = context.consume_steps(1, operation);
        if (!step) { return step; }
        BigInt quotient = numerator / denominator;
        BigInt remainder = numerator % denominator;
        if (remainder < BigInt(0)) { quotient = quotient - BigInt(1); remainder = remainder + denominator; }
        BigInt next_p = quotient * p + old_p;
        BigInt next_q = quotient * q + old_q;
        if (next_q > denominator_bound) { break; }
        Rational candidate(next_p, next_q);
        if (candidate >= root.lower && candidate <= root.upper &&
            root.polynomial.eval(candidate) == Rational(0)) {
            root.lower = root.upper = candidate;
            break;
        }
        old_p = std::move(p); p = std::move(next_p);
        old_q = std::move(q); q = std::move(next_q);
        numerator = std::move(denominator); denominator = std::move(remainder);
    }
    return Result<void>::success();
}

Result<std::vector<Algebraic>> roots_checked(const Poly& polynomial, ComputationContext& context) {
    if (polynomial.degree() <= 0) return std::vector<Algebraic>{};
    auto repeated = gcd_checked(polynomial, polynomial.differentiate(), context);
    if (!repeated) return Result<std::vector<Algebraic>>::failure(repeated.error());
    auto square_free = polynomial.div_mod(repeated.value()).first.make_monic();
    BigInt denominator_bound(1);
    for (const auto& coefficient : square_free.coeffs) {
        const auto& denominator = coefficient.get_denominator();
        denominator_bound = denominator_bound / BigInt::gcd(denominator_bound, denominator) * denominator;
        auto checked = context.require_integer_bits(denominator_bound.bit_length(), operation);
        if (!checked) return Result<std::vector<Algebraic>>::failure(checked.error());
    }
    auto isolated = detail::isolate_real_roots_exact(square_free, context, operation);
    if (!isolated) return Result<std::vector<Algebraic>>::failure(isolated.error());
    std::vector<Algebraic> roots;
    for (std::size_t i = 0; i < isolated.value().size(); ++i) {
        Algebraic root{square_free, isolated.value()[i].first, isolated.value()[i].second, i, 1};
        auto recovered = recover_rational_cut(root, denominator_bound, context);
        if (!recovered) return Result<std::vector<Algebraic>>::failure(recovered.error());
        roots.push_back(std::move(root));
    }
    return roots;
}

Result<std::optional<int>> polynomial_sign_checked(const Poly& polynomial,
    const Boundary& a, const Boundary& b, ComputationContext& context) {
    if (polynomial.is_zero()) return std::optional<int>{0};
    if (polynomial.degree() == 0) return std::optional<int>{sign(polynomial.lead_coeff())};
    auto roots = roots_checked(polynomial, context);
    if (!roots) return Result<std::optional<int>>::failure(roots.error());
    for (const auto& root : roots.value()) {
        auto left = compare_boundaries(a, Boundary{0, root}, context);
        auto right = compare_boundaries(Boundary{0, root}, b, context);
        if (!left) return Result<std::optional<int>>::failure(left.error());
        if (!right) return Result<std::optional<int>>::failure(right.error());
        if (left.value() < 0 && right.value() < 0) return std::optional<int>{};
    }
    if (a.infinity < 0)
        return std::optional<int>{sign(polynomial.lead_coeff()) *
            (polynomial.degree() % 2 ? -1 : 1)};
    if (b.infinity > 0) return std::optional<int>{sign(polynomial.lead_coeff())};
    auto left = a.finite;
    auto right = b.finite;
    while (!(left.upper < right.lower)) {
        auto l = detail::refine_exact_real_algebraic(left, context, operation);
        if (!l) return Result<std::optional<int>>::failure(l.error());
        auto r = detail::refine_exact_real_algebraic(right, context, operation);
        if (!r) return Result<std::optional<int>>::failure(r.error());
    }
    return std::optional<int>{sign(polynomial.eval((left.upper + right.lower) / Rational(2)))};
}

static Result<int> polynomial_at_root_sign(Poly polynomial, Algebraic root,
    ComputationContext& context) {
    if (polynomial.is_zero()) { return 0; }
    if (root.is_rational()) { return sign(polynomial.eval(root.lower)); }
    auto common = gcd_checked(polynomial, root.polynomial, context);
    if (!common) { return Result<int>::failure(common.error()); }
    if (common.value().degree() > 0) {
        auto zeros = detail::count_real_roots_exact(common.value(), root.lower, root.upper, context, operation);
        if (!zeros) { return Result<int>::failure(zeros.error()); }
        if (zeros.value()) { return 0; }
    }
    for (;;) {
        auto zeros = detail::count_real_roots_exact(polynomial, root.lower, root.upper, context, operation);
        if (!zeros) { return Result<int>::failure(zeros.error()); }
        if (!zeros.value()) { return sign(polynomial.eval(root.lower)); }
        auto refined = detail::refine_exact_real_algebraic(root, context, operation);
        if (!refined) { return Result<int>::failure(refined.error()); }
    }
}

static bool include_divergence(int& divergence, int contribution) {
    if (divergence && divergence != contribution) { return false; }
    divergence = contribution;
    return true;
}

static Result<void> rational_pole_divergence(const RationalFunction& rational,
    const Algebraic& root, const Boundary& a, const Boundary& b,
    int& divergence, ComputationContext& context) {
    auto left = compare_boundaries(a, Boundary{0, root}, context);
    auto right = compare_boundaries(Boundary{0, root}, b, context);
    if (!left) { return Result<void>::failure(left.error()); }
    if (!right) { return Result<void>::failure(right.error()); }
    if (left.value() > 0 || right.value() > 0) { return Result<void>::success(); }
    Poly derivative = rational.denominator;
    std::size_t multiplicity = 0;
    int denominator_sign = 0;
    do {
        derivative = derivative.differentiate();
        ++multiplicity;
        auto value = polynomial_at_root_sign(derivative, root, context);
        if (!value) { return Result<void>::failure(value.error()); }
        denominator_sign = value.value();
    } while (!denominator_sign);
    auto numerator_sign = polynomial_at_root_sign(rational.numerator, root, context);
    if (!numerator_sign) { return Result<void>::failure(numerator_sign.error()); }
    int coefficient_sign = numerator_sign.value() * denominator_sign;
    if (!coefficient_sign) { return Result<void>::failure(CasErrc::InternalInvariant,
        "reduced rational pole has a zero numerator", operation); }
    if (left.value() < 0 && !include_divergence(divergence, coefficient_sign * (multiplicity % 2 ? -1 : 1)))
        { return divergent<void>(); }
    if (right.value() < 0 && !include_divergence(divergence, coefficient_sign)) { return divergent<void>(); }
    return Result<void>::success();
}

Result<int> rational_divergence_checked(RationalFunction& rational,
    const Boundary& a, const Boundary& b, ComputationContext& context) {
    if (rational.denominator.is_zero()) { return divergent<int>(); }
    auto gcd = gcd_checked(rational.numerator, rational.denominator, context);
    if (!gcd) { return Result<int>::failure(gcd.error()); }
    rational.numerator = rational.numerator.div_mod(gcd.value()).first;
    rational.denominator = rational.denominator.div_mod(gcd.value()).first;
    if (rational.numerator.is_zero()) { return 0; }
    auto roots = roots_checked(rational.denominator, context);
    if (!roots) { return Result<int>::failure(roots.error()); }
    int divergence = 0;
    for (const auto& root : roots.value()) {
        auto pole = rational_pole_divergence(rational, root, a, b, divergence, context);
        if (!pole) { return Result<int>::failure(pole.error()); }
    }
    const int degree = rational.numerator.degree() - rational.denominator.degree();
    if (degree >= -1) {
        int coefficient = sign(rational.numerator.lead_coeff()) * sign(rational.denominator.lead_coeff());
        if (a.infinity < 0 && !include_divergence(divergence, coefficient * (degree % 2 ? -1 : 1))) { return divergent<int>(); }
        if (b.infinity > 0 && !include_divergence(divergence, coefficient)) { return divergent<int>(); }
    }
    return divergence;
}
}
