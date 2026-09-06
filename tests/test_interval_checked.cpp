#include "interval.hpp"
#include "test_common.hpp"

int main() {
    using namespace LMCAS;

    TEST_CASE("Checked interval normalization preserves exact large endpoints");
    const BigInt two_to_53("9007199254740992");
    const BigInt next_integer = two_to_53 + BigInt(1);
    std::vector<Interval> points{
        Interval::point(SymbolicExpr::number(next_integer)),
        Interval::point(SymbolicExpr::number(two_to_53))
    };
    ComputationContext exact_context;
    auto normalized = normalize_intervals_checked(std::move(points), exact_context);
    EXPECT_TRUE(normalized && normalized.value().size() == 2,
                "distinct integers beyond double precision do not merge");
    if (normalized && normalized.value().size() == 2) {
        EXPECT_TRUE(normalized.value()[0].lower.value->to_string() == two_to_53.to_string(),
                    "large exact endpoints sort by mathematical value");
        EXPECT_TRUE(normalized.value()[1].lower.value->to_string() == next_integer.to_string(),
                    "the adjacent large integer remains distinct");
    }

    TEST_CASE("Checked interval comparison uses exact IEEE values");
    auto exact_third = SymbolicExpr::number(Rational(1, 3));
    auto approximate_third = SymbolicExpr::number(1.0 / 3.0);
    Interval forward{Endpoint::closed(approximate_third), Endpoint::closed(exact_third)};
    Interval reverse{Endpoint::closed(exact_third), Endpoint::closed(approximate_third)};
    ComputationContext forward_context;
    auto forward_empty = interval_is_empty_checked(forward, forward_context);
    EXPECT_TRUE(forward_empty && !forward_empty.value(),
                "IEEE one-third is ordered below exact one-third");
    ComputationContext reverse_context;
    auto reverse_empty = interval_is_empty_checked(reverse, reverse_context);
    EXPECT_TRUE(reverse_empty && reverse_empty.value(),
                "reversing exact and approximate one-third is empty");

    TEST_CASE("Checked interval comparison proves simple quadratic surd order");
    Interval exact_surd{
        Endpoint::closed(SymbolicExpr::sqrt(SymbolicExpr::number(2))),
        Endpoint::closed(SymbolicExpr::number(2))
    };
    ComputationContext surd_context;
    auto surd_empty = interval_is_empty_checked(exact_surd, surd_context);
    EXPECT_TRUE(surd_empty && !surd_empty.value(),
                "exact sqrt(2) is proven smaller than 2 without numeric evaluation");

    Interval distinct_extensions{
        Endpoint::closed(SymbolicExpr::sqrt(SymbolicExpr::number(2))),
        Endpoint::closed(SymbolicExpr::sqrt(SymbolicExpr::number(3)))
    };
    ComputationContext distinct_context;
    auto distinct_empty = interval_is_empty_checked(
        distinct_extensions, distinct_context);
    EXPECT_TRUE(distinct_empty && !distinct_empty.value(),
                "sqrt(2) is certified smaller than sqrt(3)");

    Interval exact_numeric_expression{
        Endpoint::closed(SymbolicExpr::add(
            SymbolicExpr::number(Rational(1, 3)),
            SymbolicExpr::number(Rational(1, 6)))),
        Endpoint::closed(SymbolicExpr::number(Rational(1, 2)))
    };
    ComputationContext simplified_context;
    auto simplified_empty = interval_is_empty_checked(
        exact_numeric_expression, simplified_context);
    EXPECT_TRUE(simplified_empty && !simplified_empty.value(),
                "exact endpoint arithmetic is simplified before comparison");

    TEST_CASE("Square-root normalization preserves exactness domains");
    auto exact_irrational = SymbolicExpr::sqrt(SymbolicExpr::number(2))->simplify();
    auto exact_irrational_number = exact_irrational
        ? std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(exact_irrational))
        : nullptr;
    EXPECT_TRUE(exact_irrational &&
                    (!exact_irrational_number ||
                     !std::holds_alternative<lmmc_real_t>(
                         exact_irrational_number->value())),
                "sqrt of an exact non-square never becomes ApproxReal");

    auto exact_rational_root = SymbolicExpr::sqrt(
        SymbolicExpr::number(Rational(4, 9)))->simplify();
    auto exact_rational_number = exact_rational_root
        ? std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(exact_rational_root))
        : nullptr;
    EXPECT_TRUE(exact_rational_number &&
                    std::holds_alternative<Rational>(exact_rational_number->value()) &&
                    std::get<Rational>(exact_rational_number->value()) == Rational(2, 3),
                "sqrt of an exact rational square remains exact");

    auto approximate_root = SymbolicExpr::sqrt(
        SymbolicExpr::number(2.0))->simplify();
    auto approximate_number = approximate_root
        ? std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(approximate_root))
        : nullptr;
    EXPECT_TRUE(approximate_number &&
                    std::holds_alternative<lmmc_real_t>(approximate_number->value()) &&
                    std::abs(std::get<lmmc_real_t>(approximate_number->value()) -
                             std::sqrt(2.0)) < 1e-12,
                "sqrt of an explicit ApproxReal remains an approximate operation");

    TEST_CASE("Elementary-function normalization keeps exact inputs symbolic");
    auto tiny_exact = SymbolicExpr::number(Rational(
        BigInt(1), BigInt("10000000000000")));
    const std::vector<std::pair<FunctionNode::FuncType,
                                std::shared_ptr<SymbolicExpr>>> exact_functions{
        {FunctionNode::FuncType::Sin, SymbolicExpr::sin(tiny_exact)},
        {FunctionNode::FuncType::Cos, SymbolicExpr::cos(tiny_exact)},
        {FunctionNode::FuncType::Tan, SymbolicExpr::tan(tiny_exact)},
        {FunctionNode::FuncType::Exp, SymbolicExpr::exp(tiny_exact)},
        {FunctionNode::FuncType::Ln, SymbolicExpr::ln(tiny_exact)}
    };
    for (const auto& [type, expression] : exact_functions) {
        auto normalized = expression->simplify();
        auto function = normalized
            ? std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(normalized))
            : nullptr;
        EXPECT_TRUE(function && function->type() == type,
                    "tiny exact elementary-function argument is not tolerance-folded");
    }

    auto exact_abs = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                LMCAS::detail::make_node<NumberNode>(Rational(-7, 3))}))->simplify();
    auto exact_abs_number = exact_abs
        ? std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(exact_abs))
        : nullptr;
    EXPECT_TRUE(exact_abs_number &&
                    std::holds_alternative<Rational>(exact_abs_number->value()) &&
                    std::get<Rational>(exact_abs_number->value()) == Rational(7, 3),
                "abs of an exact Rational remains exact");

    auto exact_lambert = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::LambertW,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                LMCAS::detail::make_node<NumberNode>(Rational(1))}))->simplify();
    auto exact_lambert_function = exact_lambert
        ? std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(exact_lambert))
        : nullptr;
    EXPECT_TRUE(exact_lambert_function &&
                    exact_lambert_function->type() == FunctionNode::FuncType::LambertW,
                "LambertW of an exact non-identity input stays symbolic");

    auto approximate_lambert = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::LambertW,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                LMCAS::detail::make_node<NumberNode>(1.0)}))->simplify();
    auto approximate_lambert_number = approximate_lambert
        ? std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(approximate_lambert))
        : nullptr;
    EXPECT_TRUE(approximate_lambert_number &&
                    std::holds_alternative<lmmc_real_t>(
                        approximate_lambert_number->value()),
                "LambertW evaluates only for an explicit ApproxReal input");

    TEST_CASE("Normalization does not narrow oversized exact exponents");
    const BigInt oversized_exponent("9223372036854775808");
    auto oversized_power = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<PowerNode>(
            LMCAS::detail::make_node<NumberNode>(BigInt(2)),
            LMCAS::detail::make_node<NumberNode>(oversized_exponent)))->simplify();
    auto preserved_power = oversized_power
        ? std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(oversized_power))
        : nullptr;
    auto preserved_exponent = preserved_power
        ? std::dynamic_pointer_cast<const NumberNode>(preserved_power->exponent())
        : nullptr;
    EXPECT_TRUE(preserved_exponent &&
                    std::holds_alternative<BigInt>(preserved_exponent->value()) &&
                    std::get<BigInt>(preserved_exponent->value()) == oversized_exponent,
                "oversized BigInt exponent remains exact and unevaluated");

    TEST_CASE("Checked interval errors are explicit");
    Interval symbolic{
        Endpoint::closed(SymbolicExpr::variable("a")),
        Endpoint::closed(SymbolicExpr::number(1))
    };
    ComputationContext symbolic_context;
    auto symbolic_result = interval_is_empty_checked(symbolic, symbolic_context);
    EXPECT_TRUE(!symbolic_result &&
                    symbolic_result.error().code == CasErrc::UnboundSymbol,
                "unbound symbolic endpoints return UnboundSymbol");

    Interval invalid_domain{
        Endpoint::closed(SymbolicExpr::ln(SymbolicExpr::number(-1))),
        Endpoint::closed(SymbolicExpr::number(1))
    };
    ComputationContext domain_context;
    auto domain_result = interval_is_empty_checked(invalid_domain, domain_context);
    EXPECT_TRUE(!domain_result && domain_result.error().code == CasErrc::DomainError,
                "invalid endpoint domains return DomainError");

    Endpoint malformed_infinity = Endpoint::neg_inf();
    malformed_infinity.is_open = false;
    Interval malformed{malformed_infinity, Endpoint::pos_inf()};
    ComputationContext malformed_context;
    auto malformed_result = interval_is_empty_checked(malformed, malformed_context);
    EXPECT_TRUE(!malformed_result &&
                    malformed_result.error().code == CasErrc::InvalidArgument,
                "closed infinity endpoints are rejected");

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context({}, cancellation);
    auto cancelled = interval_contains_checked(
        Interval::point(SymbolicExpr::number(0)), 0.0, cancelled_context);
    EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                "interval membership observes cancellation");

    auto default_symbolic_contains = interval_contains_checked(symbolic, 0.0);
    EXPECT_TRUE(!default_symbolic_contains &&
                    default_symbolic_contains.error().code == CasErrc::UnboundSymbol,
                "default-context interval membership preserves UnboundSymbol");

    auto default_symbolic_empty = interval_is_empty_checked(symbolic);
    EXPECT_TRUE(!default_symbolic_empty &&
                    default_symbolic_empty.error().code == CasErrc::UnboundSymbol,
                "default-context interval emptiness preserves UnboundSymbol");

    EXPECT_TRUE(!symbolic.contains(0.0),
                "legacy interval membership unwraps endpoint errors to false");
    EXPECT_TRUE(!symbolic.is_empty(),
                "legacy interval emptiness unwraps endpoint errors conservatively");

    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext limited_context(limits);
    auto limited = interval_is_empty_checked(
        Interval::point(SymbolicExpr::number(0)), limited_context);
    EXPECT_TRUE(!limited && limited.error().code == CasErrc::ResourceLimit,
                "endpoint comparisons consume the shared step budget");

    TEST_CASE("Checked IntervalUnion union preserves exact adjacent large points");
    IntervalUnion large_left = IntervalUnion::from_single(
        Interval::point(SymbolicExpr::number(two_to_53)));
    IntervalUnion large_right = IntervalUnion::from_single(
        Interval::point(SymbolicExpr::number(next_integer)));
    ComputationContext union_context;
    auto large_union = large_left.unite_checked(large_right, union_context);
    EXPECT_TRUE(large_union && large_union.value().intervals().size() == 2,
                "checked union does not merge integers that only collide as doubles");
    if (large_union && large_union.value().intervals().size() == 2) {
        EXPECT_TRUE(large_union.value().intervals()[0].lower.value->to_string() ==
                        two_to_53.to_string(),
                    "checked union keeps the smaller exact integer first");
        EXPECT_TRUE(large_union.value().intervals()[1].lower.value->to_string() ==
                        next_integer.to_string(),
                    "checked union keeps the adjacent exact integer separate");
    }

    auto default_large_union = large_left.unite_checked(large_right);
    EXPECT_TRUE(default_large_union && default_large_union.value().intervals().size() == 2,
                "default-context checked union preserves exact adjacent large points");

    auto legacy_large_union = large_left.unite(large_right);
    EXPECT_TRUE(legacy_large_union.intervals().size() == 2,
                "legacy union now unwraps the checked exact endpoint path");

    IntervalUnion constructed_large_union({
        Interval::point(SymbolicExpr::number(next_integer)),
        Interval::point(SymbolicExpr::number(two_to_53))
    });
    EXPECT_TRUE(constructed_large_union.intervals().size() == 2,
                "legacy constructor normalization preserves exact adjacent large points");
    if (constructed_large_union.intervals().size() == 2) {
        EXPECT_TRUE(constructed_large_union.intervals()[0].lower.value->to_string() ==
                        two_to_53.to_string(),
                    "constructor normalization sorts exact endpoints by mathematical value");
        EXPECT_TRUE(constructed_large_union.intervals()[1].lower.value->to_string() ==
                        next_integer.to_string(),
                    "constructor normalization keeps adjacent exact endpoints separate");
    }

    TEST_CASE("Checked IntervalUnion intersection and complement use exact endpoints");
    IntervalUnion closed_open = IntervalUnion::from_single(Interval{
        Endpoint::closed(SymbolicExpr::number(1)),
        Endpoint::open(SymbolicExpr::number(2))
    });
    IntervalUnion open_closed = IntervalUnion::from_single(Interval{
        Endpoint::open(SymbolicExpr::number(2)),
        Endpoint::closed(SymbolicExpr::number(3))
    });
    ComputationContext intersection_context;
    auto intersection = closed_open.intersect_checked(open_closed, intersection_context);
    EXPECT_TRUE(intersection && intersection.value().is_empty(),
                "disjoint intervals sharing only an open boundary have empty intersection");

    IntervalUnion punctured = closed_open.unite(open_closed);
    ComputationContext complement_context;
    auto complement = punctured.complement_checked(complement_context);
    EXPECT_TRUE(complement && complement.value().intervals().size() == 3,
                "complement preserves the missing boundary point");
    if (complement && complement.value().intervals().size() == 3) {
        const auto& middle = complement.value().intervals()[1];
        EXPECT_TRUE(!middle.lower.is_open && !middle.upper.is_open,
                    "the gap at 2 is closed on both sides");
        EXPECT_TRUE(middle.lower.value->to_string() == "2" &&
                        middle.upper.value->to_string() == "2",
                    "the middle complement interval is the singleton {2}");
    }

    TEST_CASE("Checked IntervalUnion operations propagate endpoint errors");
    ComputationContext symbolic_factory_context;
    auto symbolic_factory = IntervalUnion::from_intervals_checked(
        {symbolic}, symbolic_factory_context);
    EXPECT_TRUE(!symbolic_factory &&
                    symbolic_factory.error().code == CasErrc::UnboundSymbol,
                "checked IntervalUnion factory propagates endpoint errors");

    IntervalUnion symbolic_union = IntervalUnion::from_single(symbolic);
    ComputationContext symbolic_union_context;
    auto symbolic_union_result = symbolic_union.unite_checked(
        IntervalUnion::empty(), symbolic_union_context);
    EXPECT_TRUE(!symbolic_union_result &&
                    symbolic_union_result.error().code == CasErrc::UnboundSymbol,
                "checked union propagates symbolic endpoint errors");

    auto default_symbolic_union = symbolic_union.unite_checked(IntervalUnion::empty());
    EXPECT_TRUE(!default_symbolic_union &&
                    default_symbolic_union.error().code == CasErrc::UnboundSymbol,
                "default-context checked union propagates symbolic endpoint errors");

    auto legacy_symbolic_union = symbolic_union.unite(IntervalUnion::empty());
    EXPECT_TRUE(legacy_symbolic_union.is_empty(),
                "legacy union unwraps endpoint errors to an empty compatibility result");

    auto legacy_symbolic_complement = symbolic_union.complement();
    EXPECT_TRUE(legacy_symbolic_complement.is_empty(),
                "legacy complement unwraps endpoint errors to an empty compatibility result");

    CancellationToken union_cancellation;
    union_cancellation.cancel();
    ComputationContext cancelled_union_context({}, union_cancellation);
    auto cancelled_factory = IntervalUnion::from_intervals_checked(
        {Interval::point(SymbolicExpr::number(0))}, cancelled_union_context);
    EXPECT_TRUE(!cancelled_factory &&
                    cancelled_factory.error().code == CasErrc::Cancelled,
                "checked IntervalUnion factory observes cancellation");
    auto cancelled_union = large_left.intersect_checked(
        large_right, cancelled_union_context);
    EXPECT_TRUE(!cancelled_union &&
                    cancelled_union.error().code == CasErrc::Cancelled,
                "checked IntervalUnion operations observe cancellation");

    return TEST_REPORT();
}
