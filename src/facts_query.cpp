#include "internal/facts_query.hpp"
#include "computation_context.hpp"
#include "internal/exact_constant_bounds.hpp"
#include "internal/exact_root.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/symbolic_ast.hpp"

#include <cmath>
#include <new>
#include <utility>

namespace LMCAS::detail {
namespace {
using Node = std::shared_ptr<const SymbolicNode>;
using T = Tribool;
using FT = FunctionNode::FuncType;

// Fold closed rational arithmetic without normalization or fact queries.
std::optional<Rational> closed_rational_value(
    const Node& node, Domain domain, ComputationContext& context);

std::optional<Rational> closed_rational_power(
    const PowerNode& power, Domain domain, ComputationContext& context,
    ExactBoundArithmetic& arithmetic) {
    auto base = closed_rational_value(power.base(), domain, context);
    if (!base) { return std::nullopt; }
    auto exponent = closed_rational_value(power.exponent(), domain, context);
    if (!exponent || !exponent->is_integer()) { return std::nullopt; }
    BigInt count = exponent->to_bigint();
    if (base->is_zero() && (count.is_negative() ||
        (count.is_zero() && domain == Domain::Real))) { return std::nullopt; }
    if (count.is_negative()) {
        *base = arithmetic.div(Rational(1), *base);
        count = -count;
    }
    Rational value(1);
    while (!count.is_zero()) {
        arithmetic.step();
        if (count.is_odd()) { value = arithmetic.mul(value, *base); }
        count = count >> 1;
        if (!count.is_zero()) { *base = arithmetic.mul(*base, *base); }
    }
    return value;
}

std::optional<Rational> closed_rational_value(
    const Node& node, Domain domain, ComputationContext& context) {
    ExactBoundArithmetic arithmetic(context, "facts.constant_exponent");
    arithmetic.require(context.enter_recursion("facts.constant_exponent"));
    struct Leave {
        ComputationContext& context;
        ~Leave() { context.leave_recursion(); }
    } leave{context};
    if (const auto* number = dynamic_cast<const NumberNode*>(node.get())) {
        if (const auto* value = std::get_if<double>(&number->value());
            value && !std::isfinite(*value)) { return std::nullopt; }
        auto value = exact_number_as_rational(*number);
        arithmetic.bits(std::max(value.get_numerator().bit_length(),
                                 value.get_denominator().bit_length()));
        return value;
    }
    const auto* sum = dynamic_cast<const AddNode*>(node.get());
    const auto* product = dynamic_cast<const MultiplyNode*>(node.get());
    if (sum || product) {
        Rational value(product ? 1 : 0);
        for (const auto& operand : sum ? sum->operands() : product->operands()) {
            auto child = closed_rational_value(operand, domain, context);
            if (!child) { return std::nullopt; }
            value = product ? arithmetic.mul(value, *child) : arithmetic.add(value, *child);
        }
        return value;
    }
    const auto* power = dynamic_cast<const PowerNode*>(node.get());
    if (!power) { return std::nullopt; }
    return closed_rational_power(*power, domain, context, arithmetic);
}

T conjunction(T a, T b) {
    if (a == T::False || b == T::False) return T::False;
    return a == T::True && b == T::True ? T::True : T::Unknown;
}

T disjunction(T a, T b) {
    if (a == T::True || b == T::True) return T::True;
    return a == T::False && b == T::False ? T::False : T::Unknown;
}

struct ValueFacts {
    T defined = T::Unknown;
    T real = T::Unknown;
    T nonzero = T::Unknown;
    T positive = T::Unknown;
    T nonnegative = T::Unknown;
};

// Traverse once without normalization, floating evaluation, or callbacks.
class DomainTraversal {
public:
    DomainTraversal(const FactsQuery& facts, Domain domain,
                    ComputationContext& context, bool project = false)
        : facts_(facts), domain_(domain), context_(context), project_(project) {}

    std::optional<CasError> error;
    bool expressible = true;
    std::vector<std::shared_ptr<SymbolicExpr>> conditions;

    ValueFacts walk(const Node& node, bool values = false) {
        if (error) return {};
        if (!node) {
            error = CasError{CasErrc::InvalidArgument, "null expression", "facts.domain"};
            return {};
        }
        auto entered = context_.enter_recursion("facts.domain");
        if (!entered) {
            error = entered.error();
            return {};
        }
        struct Leave {
            ComputationContext& context;
            ~Leave() { context.leave_recursion(); }
        } leave{context_};
        auto result = structural(node, values);
        // Value queries need definedness; projections retain conditions such as x != 0.
        if (result.defined == T::False || (result.defined != T::True && !project_)) {
            result.real = result.nonzero = result.positive = result.nonnegative = T::Unknown;
        } else if (result.defined == T::True && values) {
            infer(node, result);
        }
        return result;
    }

    Node false_condition() { return relation(number(0), number(1), RelationOp::EQ); }

private:
    const FactsQuery& facts_;
    Domain domain_;
    ComputationContext& context_;
    bool project_;

    T checked(Result<T> result) {
        if (!result) {
            if (!error) error = result.error();
            return T::Unknown;
        }
        return result.value();
    }

    template<class N, class... Args>
    Node make(Args&&... args) {
        if (error) return {};
        auto reserved = context_.reserve_nodes(1, "facts.domain_condition");
        if (!reserved) {
            error = reserved.error();
            return {};
        }
        return make_node<N>(std::forward<Args>(args)...);
    }

    Node number(int n) { return make<NumberNode>(BigInt(n)); }
    Node relation(Node a, Node b, RelationOp op) {
        return make<RelationalNode>(std::move(a), std::move(b), op);
    }
    Node logical(Node a, Node b, LogicalNode::Op op) {
        return make<LogicalNode>(std::move(a), std::move(b), op);
    }
    Node function(FT type, const Node& argument) {
        return make<FunctionNode>(type, std::vector<Node>{argument});
    }
    Node multiply(Node a, Node b) {
        std::vector<Node> operands;
        operands.reserve(2);
        operands.push_back(std::move(a));
        operands.push_back(std::move(b));
        return make<MultiplyNode>(std::move(operands));
    }

    template<class Predicate>
    void require(ValueFacts& result, T truth, Predicate&& predicate) {
        result.defined = conjunction(result.defined, truth);
        if (project_ && truth == T::Unknown && !error) {
            auto condition = predicate();
            if (!error && condition) conditions.push_back(make_expression_ptr(std::move(condition)));
        }
    }

    static void implications(ValueFacts& result) {
        if (result.nonnegative == T::True && result.real != T::False) {
            result.real = T::True;
            if (result.nonzero == T::True) result.positive = T::True;
        }
        if (result.real == T::False) {
            result.positive = result.nonnegative = T::False;
        } else if (result.positive == T::True) {
            result.real = result.nonzero = result.nonnegative = T::True;
        } else if (result.nonzero == T::False) {
            result.real = result.nonnegative = T::True;
            result.positive = T::False;
        } else if (result.real == T::True && result.nonnegative == T::False) {
            result.nonzero = T::True;
            result.positive = T::False;
        }
    }

    void infer(const Node& node, ValueFacts& result) {
        implications(result);
        if (result.positive == T::Unknown) { result.positive = checked(facts_.is_positive(node, context_)); }
        if (error) { return; }
        implications(result);
        if (result.real == T::Unknown) { result.real = checked(facts_.is_real(node, context_)); }
        if (error) { return; }
        if (result.nonzero == T::Unknown) { result.nonzero = checked(facts_.is_nonzero(node, context_)); }
        if (error) { return; }
        if (result.nonnegative == T::Unknown) { result.nonnegative = checked(facts_.is_nonnegative(node, context_)); }
        if (error) { return; }
        if (result.positive == T::Unknown || result.nonnegative == T::Unknown) {
            auto negative = checked(facts_.is_negative(node, context_));
            if (error) { return; }
            if (negative == T::True) {
                result.real = result.nonzero = T::True;
                result.positive = result.nonnegative = T::False;
            }
        }
        implications(result);
    }

    struct ArithmeticReduction {
        T all_real = T::True, all_nonnegative = T::True, all_nonzero = T::True;
        bool any_positive = false, any_zero = false, all_zero = true, signs_known = true;
        unsigned negative_parity = 0;
        std::size_t nonreal_count = 0;
        bool classified_reality = true, other_factors_nonzero = true;

        ValueFacts finish(ValueFacts result, bool product) const {
            if (all_real == T::True) { result.real = T::True; }
            if (nonreal_count == 1 && classified_reality && (!product || other_factors_nonzero))
                { result.real = T::False; }
            if (product) {
                if (any_zero) { result.nonzero = T::False; }
                else if (all_nonzero == T::True) { result.nonzero = T::True; }
                if (signs_known) {
                    result.positive = result.nonnegative = negative_parity ? T::False : T::True;
                }
            } else if (all_nonnegative == T::True) {
                result.nonnegative = T::True;
                if (any_positive) { result.positive = T::True; }
                if (all_zero) { result.nonzero = T::False; }
            }
            return result;
        }
    };

    ValueFacts arithmetic(const std::vector<Node>& operands, bool product, bool values) {
        ValueFacts result;
        result.defined = T::True;
        ArithmeticReduction reduction;
        for (const auto& operand : operands) {
            auto child = walk(operand, values);
            if (error) { return {}; }
            result.defined = conjunction(result.defined, child.defined);
            if (!values) { continue; }
            reduction.all_real = conjunction(reduction.all_real, child.real);
            reduction.all_nonnegative = conjunction(reduction.all_nonnegative, child.nonnegative);
            reduction.all_nonzero = conjunction(reduction.all_nonzero, child.nonzero);
            reduction.any_positive |= child.positive == T::True;
            reduction.any_zero |= child.nonzero == T::False;
            reduction.all_zero &= child.nonzero == T::False;
            if (child.real == T::False) { ++reduction.nonreal_count; }
            else if (child.real == T::Unknown) { reduction.classified_reality = false; }
            else if (child.nonzero != T::True) { reduction.other_factors_nonzero = false; }
            if (child.positive == T::True) {
                continue;
            } else if (child.real == T::True && child.nonnegative == T::False) {
                reduction.negative_parity ^= 1;
            } else {
                reduction.signs_known = false;
            }
        }
        if (!values) { return result; }
        return reduction.finish(result, product);
    }

    ValueFacts power(const PowerNode& node, bool values) {
        const auto exponent_number = std::dynamic_pointer_cast<const NumberNode>(node.exponent());
        BigInt integer;
        const bool integral = exponent_number && try_get_integer_value(exponent_number, integer);
        bool needs_base_value = values || !integral;
        if (!needs_base_value) {
            needs_base_value = integer.is_negative() ||
                (integer.is_zero() && domain_ == Domain::Real);
        }
        auto base = walk(node.base(), needs_base_value);
        if (error) { return {}; }
        auto exponent = walk(node.exponent(), !integral);
        if (error) { return {}; }
        ValueFacts result;
        result.defined = conjunction(base.defined, exponent.defined);
        if (integral) {
            return integral_power(node, base, integer, result);
        }
        if (exponent_number) {
            return nonintegral_power(node, base, exponent, result, true);
        }
        if (exponent.nonzero == T::False) {
            return integral_power(node, base, BigInt(0), result);
        }
        const auto exact_exponent = closed_rational_value(node.exponent(), domain_, context_);
        if (exact_exponent && exact_exponent->is_integer()) {
            return integral_power(node, base, exact_exponent->to_bigint(), result);
        }
        return nonintegral_power(node, base, exponent, result, exact_exponent.has_value());
    }

    ValueFacts integral_power(const PowerNode& node, const ValueFacts& base, const BigInt& integer, ValueFacts result) {
        const bool zero = integer == BigInt(0);
        if (integer.is_negative() || (zero && domain_ == Domain::Real)) {
            require(result, base.nonzero, [&] {
                return relation(node.base(), number(0), RelationOp::NEQ);
            });
        }
        if (zero) {
            result.real = result.nonzero = result.positive = result.nonnegative = T::True;
        } else {
            if (base.real == T::True) { result.real = T::True; }
            const T nonzero = project_ && integer.is_negative() ? T::True : base.nonzero;
            result.nonzero = nonzero;
            if (base.real == T::True) {
                if (integer % BigInt(2) == BigInt(0)) {
                    result.nonnegative = T::True;
                    result.positive = nonzero;
                } else {
                    result.positive = base.positive;
                    result.nonnegative = base.nonnegative;
                }
            }
        }
        return result;
    }

    T nonintegral_domain(const ValueFacts& base, const ValueFacts& exponent, bool known_noninteger, bool& unresolved_integer_domain) {
        T allowed;
        unresolved_integer_domain = false;
        if (domain_ == Domain::Real) {
            allowed = disjunction(conjunction(base.positive, exponent.real),
                                  conjunction(base.nonnegative, exponent.positive));
            unresolved_integer_domain = !known_noninteger && base.nonnegative != T::True;
            if (unresolved_integer_domain && allowed != T::True) { allowed = T::Unknown; }
            if (known_noninteger && base.real == T::True && base.nonnegative == T::False)
                { allowed = T::False; }
            if (base.nonzero == T::False && exponent.positive == T::False &&
                exponent.real == T::True) { allowed = T::False; }
        } else {
            const T zero_nonnegative = conjunction(
                base.nonzero == T::False ? T::True :
                    base.nonzero == T::True ? T::False : T::Unknown,
                exponent.nonnegative);
            allowed = disjunction(base.nonzero, zero_nonnegative);
        }
        return allowed;
    }

    ValueFacts nonintegral_power(const PowerNode& node, const ValueFacts& base, const ValueFacts& exponent, ValueFacts result, bool known_noninteger) {
        bool unresolved_integer_domain;
        const T allowed = nonintegral_domain(base, exponent, known_noninteger, unresolved_integer_domain);
        if (allowed == T::Unknown &&
            (exponent.real != T::True || unresolved_integer_domain)) {
            expressible = false;
            result.defined = conjunction(result.defined, T::Unknown);
        } else {
            require(result, allowed, [&] {
                return logical(
                    relation(node.base(), number(0),
                        domain_ == Domain::Real ? RelationOp::GT : RelationOp::NEQ),
                    logical(relation(node.base(), number(0), RelationOp::EQ),
                            relation(node.exponent(), number(0),
                                domain_ == Domain::Real ? RelationOp::GT : RelationOp::GEQ),
                            LogicalNode::Op::And),
                    LogicalNode::Op::Or);
            });
        }
        return nonintegral_value(base, exponent, result, known_noninteger);
    }

    ValueFacts nonintegral_value(const ValueFacts& base, const ValueFacts& exponent,
                                ValueFacts result, bool known_noninteger) const {
        result.nonzero = base.nonzero;
        if (base.nonzero == T::False && domain_ == Domain::Complex &&
            exponent.positive != T::True) {
            return zero_base_power_value(exponent, result);
        }
        if (base.positive == T::True && exponent.real == T::True) {
            result.real = result.positive = result.nonnegative = T::True;
        } else if (base.nonnegative == T::True && exponent.positive == T::True) {
            result.real = result.nonnegative = T::True;
            result.positive = base.nonzero;
        } else if (known_noninteger && base.real == T::True && base.nonnegative == T::False) {
            result.real = T::False;
        }
        return result;
    }

    ValueFacts zero_base_power_value(const ValueFacts& exponent, ValueFacts result) const {
        result.nonzero = T::Unknown;
        if (exponent.real != T::True ||
            (exponent.nonnegative != T::True && !project_)) { return result; }
        result.real = result.nonnegative = T::True;
        if (project_ && exponent.positive == T::False) {
            result.nonzero = result.positive = T::True;
        }
        return result;
    }

    ValueFacts elementary(const FunctionNode& node, bool values) {
        const auto& arguments = node.arguments();
        if (node.type() == FT::Infinity) { return {T::False}; }
        const bool logarithm_base = node.type() == FT::Log && arguments.size() == 2;
        if (arguments.size() != 1 && !logarithm_base) {
            expressible = false;
            return {};
        }
        const auto& argument = arguments[0];
        const bool restricted = restricted_argument(node.type());
        auto child = walk(argument, values || restricted);
        if (error) { return {}; }
        ValueFacts result;
        result.defined = child.defined;
        return elementary_value(node, child, result);
    }

    bool restricted_argument(FT type) const {
        if (type == FT::Ln || type == FT::Log || type == FT::Tan ||
            type == FT::Sec || type == FT::Cot || type == FT::Csc) {
            return true;
        }
        if (type == FT::ArcTan || type == FT::Sgn) {
            return domain_ == Domain::Complex;
        }
        if (type == FT::Sqrt || type == FT::ArcSin || type == FT::ArcCos) {
            return domain_ == Domain::Real;
        }
        return false;
    }

    ValueFacts elementary_value(const FunctionNode& node, const ValueFacts& child, ValueFacts result) {
        switch (node.type()) {
        case FT::Exp: {
            result.nonzero = T::True;
            if (child.real == T::True)
                { result.real = result.positive = result.nonnegative = T::True; }
            return result;
        }
        case FT::Sin: { break; }
        case FT::Cos: { break; }
        case FT::Sinh: { break; }
        case FT::Cosh: { break; }
        case FT::Erf: { break; }
        case FT::Si: { break; }
        default: { return restricted_function(node, child, result); }
        }
        if (child.real == T::True) { result.real = T::True; }
        if (child.nonzero == T::False) {
            result.nonzero = node.type() == FT::Cos || node.type() == FT::Cosh
                ? T::True : T::False;
            if (result.nonzero == T::True) { result.positive = T::True; }
        }
        return result;
    }

    ValueFacts restricted_function(const FunctionNode& node, const ValueFacts& child, ValueFacts result) {
        switch (node.type()) {
        case FT::Ln: { return logarithm(node.arguments(), child, result); }
        case FT::Log: { return logarithm(node.arguments(), child, result); }
        case FT::Sqrt: { return real_valued_function(node, child, result); }
        case FT::Abs: { return real_valued_function(node, child, result); }
        case FT::ComplexAbs: { return real_valued_function(node, child, result); }
        case FT::Tan: { return trigonometric_domain(node, child, result); }
        case FT::Sec: { return trigonometric_domain(node, child, result); }
        case FT::Cot: { return trigonometric_domain(node, child, result); }
        case FT::Csc: { return trigonometric_domain(node, child, result); }
        case FT::ArcSin: { return trigonometric_domain(node, child, result); }
        case FT::ArcCos: { return trigonometric_domain(node, child, result); }
        case FT::Sgn: { return real_valued_function(node, child, result); }
        case FT::ArcTan: { return trigonometric_domain(node, child, result); }
        case FT::LambertW: { return lambert_w_domain(node, child, result); }
        default: {
            // Unsupported singularities and conditional branches remain undecidable.
            result.defined = conjunction(result.defined, T::Unknown);
            expressible = false;
            break;
        }
        }
        return result;
    }

    ValueFacts real_valued_function(const FunctionNode& node, const ValueFacts& child, ValueFacts result) {
        if (node.type() == FT::Sqrt) {
            if (domain_ == Domain::Real) {
                require(result, child.nonnegative, [&] {
                    return relation(node.arguments()[0], number(0), RelationOp::GEQ);
                });
            }
            result.nonzero = child.nonzero;
            if (child.nonnegative == T::True || (project_ && domain_ == Domain::Real)) {
                result.real = result.nonnegative = T::True;
                result.positive = child.nonzero;
            } else if (child.real == T::True && child.nonnegative == T::False) {
                result.real = T::False;
            }
            return result;
        }
        if (node.type() == FT::Abs || node.type() == FT::ComplexAbs) {
            result.real = result.nonnegative = T::True;
            result.nonzero = result.positive = child.nonzero;
            return result;
        }
        if (domain_ == Domain::Real || child.real == T::True) {
            if (child.defined == T::True) { result.real = T::True; }
            result.nonzero = child.nonzero;
            result.positive = child.positive;
            result.nonnegative = child.nonnegative;
        } else {
            result.defined = conjunction(result.defined, T::Unknown);
            expressible = false;
        }
        return result;
    }
    ValueFacts lambert_w_domain(const FunctionNode& node, const ValueFacts& child,
                                ValueFacts result) {
        if (child.real == T::True && child.nonnegative == T::True) {
            result.real = result.nonnegative = T::True;
            result.positive = child.positive;
            result.nonzero = child.nonzero;
            return result;
        }
        if (domain_ == Domain::Real) {
            const T in_principal_domain = child.nonnegative == T::True ? T::True : T::Unknown;
            require(result, in_principal_domain, [&] {
                auto lower = multiply(number(-1), function(FT::Exp, number(-1)));
                return relation(node.arguments()[0], std::move(lower), RelationOp::GEQ);
            });
            if (in_principal_domain == T::True || project_) {
                result.real = T::True;
                result.nonnegative = child.nonnegative;
                result.positive = child.positive;
                result.nonzero = child.nonzero;
            }
            return result;
        }
        result.defined = conjunction(result.defined, T::Unknown);
        expressible = false;
        return result;
    }

    ValueFacts logarithm(const std::vector<Node>& arguments, const ValueFacts& child, ValueFacts result) {
        const auto& argument = arguments[0];
        require(result, domain_ == Domain::Real ? child.positive : child.nonzero, [&] {
            return relation(argument, number(0),
                domain_ == Domain::Real ? RelationOp::GT : RelationOp::NEQ);
        });
        if (child.positive == T::True || (project_ && domain_ == Domain::Real))
            { result.real = T::True; }
        if (const auto* numeric = dynamic_cast<const NumberNode*>(argument.get())) {
            if (numeric->is_one()) { result.nonzero = T::False; }
            else if (numeric->is_positive()) {
                result.nonzero = T::True;
                result.positive = result.nonnegative =
                    exact_number_as_rational(*numeric) > Rational(1) ? T::True : T::False;
            } else if (!numeric->is_zero()) {
                result.real = T::False;
                result.nonzero = T::True;
            }
        }
        if (arguments.size() == 2) { return logarithm_base(arguments[1], child, result); }
        return result;
    }
    ValueFacts logarithm_base(const Node& argument, const ValueFacts& child, ValueFacts result) {
        auto base = walk(argument, true);
        if (error) { return {}; }
        result.defined = conjunction(result.defined, base.defined);
        require(result, domain_ == Domain::Real ? base.positive : base.nonzero, [&] {
            return relation(argument, number(0),
                domain_ == Domain::Real ? RelationOp::GT : RelationOp::NEQ);
        });
        T different_from_one = T::Unknown;
        if (const auto* numeric = dynamic_cast<const NumberNode*>(argument.get()))
            { different_from_one = numeric->is_one() ? T::False : T::True; }
        require(result, different_from_one, [&] {
            return relation(argument, number(1), RelationOp::NEQ);
        });
        result.real = (project_ && domain_ == Domain::Real) ||
            (child.positive == T::True && base.positive == T::True)
            ? T::True : T::Unknown;
        result.positive = result.nonnegative = T::Unknown;
        return result;
    }

    ValueFacts trigonometric_domain(const FunctionNode& node, const ValueFacts& child, ValueFacts result) {
        const auto& argument = node.arguments()[0];
        if (node.type() == FT::ArcSin || node.type() == FT::ArcCos) {
            if (domain_ == Domain::Real) {
                T inside = T::Unknown;
                if (const auto* numeric = dynamic_cast<const NumberNode*>(argument.get())) {
                    const auto value = exact_number_as_rational(*numeric);
                    inside = value >= Rational(-1) && value <= Rational(1) ? T::True : T::False;
                }
                require(result, inside, [&] {
                    return logical(relation(argument, number(-1), RelationOp::GEQ),
                                   relation(argument, number(1), RelationOp::LEQ),
                                   LogicalNode::Op::And);
                });
            }
            return result;
        }
        if (node.type() == FT::ArcTan) {
            if (child.real == T::True) {
                result.real = T::True;
                result.nonzero = child.nonzero;
                result.positive = child.positive;
                result.nonnegative = child.nonnegative;
            } else if (domain_ == Domain::Complex) {
                result.defined = conjunction(result.defined, T::Unknown);
                expressible = false;
            }
            return result;
        }
        const bool cosine = node.type() == FT::Tan || node.type() == FT::Sec;
        T nonzero = T::Unknown;
        if (child.nonzero == T::False) { nonzero = cosine ? T::True : T::False; }
        require(result, nonzero, [&] {
            return relation(function(cosine ? FT::Cos : FT::Sin, argument),
                            number(0), RelationOp::NEQ);
        });
        if (child.real == T::True) { result.real = T::True; }
        return result;
    }

    ValueFacts structural(const Node& node, bool values) {
        if (const auto* number = dynamic_cast<const NumberNode*>(node.get())) {
            const T nz = number->is_zero() ? T::False : T::True;
            const T positive = number->is_positive() ? T::True : T::False;
            return {T::True, T::True, nz, positive,
                    nz == T::False || positive == T::True ? T::True : T::False};
        }
        if (const auto* variable = dynamic_cast<const VariableNode*>(node.get())) {
            return variable_facts(*variable);
        }
        if (const auto* add = dynamic_cast<const AddNode*>(node.get()))
            { return arithmetic(add->operands(), false, values); }
        if (const auto* multiply = dynamic_cast<const MultiplyNode*>(node.get()))
            { return arithmetic(multiply->operands(), true, values); }
        if (const auto* pow = dynamic_cast<const PowerNode*>(node.get())) { return power(*pow, values); }
        if (const auto* fn = dynamic_cast<const FunctionNode*>(node.get())) { return elementary(*fn, values); }
        if (const auto* complex = dynamic_cast<const ComplexNode*>(node.get())) {
            return complex_facts(*complex);
        }
        if (const auto* root = dynamic_cast<const RootOfNode*>(node.get())) {
            return root_facts(*root, values);
        }
        if (const auto* relation = dynamic_cast<const RelationalNode*>(node.get())) {
            auto left = walk(relation->left());
            if (error) { return {}; }
            auto right = walk(relation->right());
            return {conjunction(left.defined, right.defined)};
        }
        expressible = false;
        return {};
    }

    ValueFacts variable_facts(const VariableNode& variable) {
        if (variable.name() == "pi" || variable.name() == "π" ||
            variable.name() == "e" || variable.name() == "phi")
            { return {T::True, T::True, T::True, T::True, T::True}; }
        if (is_imaginary_unit_name(variable.name()))
            { return {domain_ == Domain::Real ? T::False : T::True,
                    T::False, T::True, T::False, T::False}; }
        return {T::True};
    }

    ValueFacts complex_facts(const ComplexNode& complex) {
        auto real = walk(complex.real(), true);
        if (error) { return {}; }
        auto imag = walk(complex.imag(), true);
        if (error) { return {}; }
        ValueFacts result;
        result.defined = conjunction(conjunction(real.defined, imag.defined),
                                     conjunction(real.real, imag.real));
        if (real.real == T::Unknown || imag.real == T::Unknown) { expressible = false; }
        if (domain_ == Domain::Real) {
            const T zero = imag.nonzero == T::False ? T::True :
                imag.nonzero == T::True ? T::False : T::Unknown;
            require(result, zero, [&] {
                return relation(complex.imag(), number(0), RelationOp::EQ);
            });
        }
        result.nonzero = disjunction(real.nonzero, imag.nonzero);
        result.real = imag.nonzero == T::False ? T::True :
            imag.nonzero == T::True ? T::False : T::Unknown;
        if (result.real == T::True) {
            result.positive = real.positive;
            result.nonnegative = real.nonnegative;
        }
        return result;
    }

    ValueFacts root_facts(const RootOfNode& root, bool values) {
        ValueFacts result;
        result.defined = T::True; /**< ExactRootId 已在构造时校验。 */
        if (values || domain_ == Domain::Real) {
            auto isolation = isolate_exact_root(root.exact_id(), context_, "facts.rootof");
            if (!isolation) {
                error = isolation.error();
                return {};
            }
            if (const auto* real = std::get_if<RealIsolation>(&isolation.value())) {
                result.real = T::True;
                if (real->value.lower > Rational(0)) { result.positive = T::True; }
                else if (real->value.upper < Rational(0)) { result.nonnegative = T::False; }
                else if (real->value.lower == Rational(0) && real->value.upper == Rational(0))
                    { result.nonzero = T::False; }
            } else {
                result.real = T::False;
                result.nonzero = T::True;
                if (domain_ == Domain::Real) { result.defined = T::False; }
            }
        }
        return result;
    }
};

Result<Tribool> query_value(const Node& node, const FactsQuery& facts, Domain domain,
                           T ValueFacts::*member, bool values, ComputationContext& context) {
    if (domain != Domain::Real && domain != Domain::Complex) {
        return Result<T>::failure(CasErrc::InvalidArgument,
            "value domain must be Real or Complex", "facts.domain");
    }
    try {
        DomainTraversal traversal(facts, domain, context);
        auto result = traversal.walk(node, values);
        if (traversal.error) return Result<T>::failure(*traversal.error);
        return result.*member;
    } catch (const CasError& error) {
        return Result<T>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<T>::failure(CasErrc::ResourceLimit,
                                  "value facts allocation failed", "facts.domain");
    }
}

}

Result<Tribool> query_definedness(const Node& node, const FactsQuery& facts, Domain domain, ComputationContext& context) {
    return query_value(node, facts, domain, &ValueFacts::defined, false, context);
}
Result<Tribool> query_nonzero_value(const Node& node, const FactsQuery& facts, Domain domain, ComputationContext& context) {
    return query_value(node, facts, domain, &ValueFacts::nonzero, true, context);
}
Result<Tribool> query_real_value(const Node& node, const FactsQuery& facts, ComputationContext& context) {
    return query_value(node, facts, Domain::Complex, &ValueFacts::real, true, context);
}
Result<Tribool> query_positive_value(const Node& node, const FactsQuery& facts, ComputationContext& context) {
    return query_value(node, facts, Domain::Complex, &ValueFacts::positive, true, context);
}
Result<Tribool> query_nonnegative_value(const Node& node, const FactsQuery& facts, ComputationContext& context) {
    return query_value(node, facts, Domain::Complex, &ValueFacts::nonnegative, true, context);
}

Result<std::optional<std::vector<std::shared_ptr<SymbolicExpr>>>> domain_constraints(
    const Node& node, const FactsQuery& facts, Domain domain, ComputationContext& context) {
    using Output = std::optional<std::vector<std::shared_ptr<SymbolicExpr>>>;
    if (domain != Domain::Real && domain != Domain::Complex) {
        return Result<Output>::failure(CasErrc::InvalidArgument,
            "value domain must be Real or Complex", "facts.domain");
    }
    try {
        DomainTraversal traversal(facts, domain, context, true);
        auto result = traversal.walk(node);
        if (traversal.error) return Result<Output>::failure(*traversal.error);
        if (result.defined == T::False) {
            auto impossible = traversal.false_condition();
            if (traversal.error) return Result<Output>::failure(*traversal.error);
            traversal.conditions.clear();
            traversal.conditions.push_back(make_expression_ptr(std::move(impossible)));
        } else if (!traversal.expressible) {
            return Output{};
        }
        return Output{std::move(traversal.conditions)};
    } catch (const CasError& error) {
        return Result<Output>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<Output>::failure(CasErrc::ResourceLimit,
                                      "domain projection allocation failed", "facts.domain");
    }
}

}
