#include "expr.hpp"
#include "internal/symbolic_ast.hpp"
#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

using namespace LMCAS;

namespace {

ExprPtr round_trip(const ExprPtr& expression) {
    auto encoded = serialize_expr(expression);
    if (!encoded) {
        ADD_FAILURE() << encoded.error().message;
        return nullptr;
    }
    auto decoded = parse_serialized_expr(encoded.value());
    if (!decoded) {
        ADD_FAILURE() << decoded.error().message;
        return nullptr;
    }
    return decoded.value();
}

void expect_same_value(const ExprPtr& expression) {
    auto decoded = round_trip(expression);
    ASSERT_NE(decoded, nullptr);
    EXPECT_TRUE(structurally_equal(*decoded, *expression));
}

TEST(ExprSerialization, ExactNumbers) {
    const BigInt huge("-12345678901234567890123456789012345678901234567890");
    auto integer_value = round_trip(SymbolicExpr::number(huge));
    ASSERT_NE(integer_value, nullptr);
    ASSERT_TRUE(integer_value->is_big_int());
    EXPECT_EQ(integer_value->get_big_int(), huge);

    const Rational fraction(
        BigInt("123456789012345678901234567890123456789"),
        BigInt("700000000000000000000000000000000000003"));
    auto rational_value = round_trip(SymbolicExpr::number(fraction));
    ASSERT_NE(rational_value, nullptr);
    ASSERT_TRUE(rational_value->is_rational());
    EXPECT_EQ(rational_value->get_rational(), fraction);
}

TEST(ExprSerialization, Binary64NegativeZero) {
    auto negative = approx_real(-0.0);
    ASSERT_TRUE(negative);
    auto encoded = serialize_expr(negative.value());
    ASSERT_TRUE(encoded);
    EXPECT_NE(encoded.value().find("8000000000000000"), std::string::npos);

    auto positive = approx_real(0.0);
    ASSERT_TRUE(positive);
    auto positive_bytes = serialize_expr(positive.value());
    ASSERT_TRUE(positive_bytes);
    EXPECT_NE(encoded.value(), positive_bytes.value());

    auto decoded = parse_serialized_expr(encoded.value());
    ASSERT_TRUE(decoded);
    auto value = evalf(*decoded.value());
    ASSERT_TRUE(value);
    EXPECT_EQ(value.value().value, 0.0);
    EXPECT_TRUE(std::signbit(value.value().value));
}

TEST(ExprSerialization, SymbolNamesAndConstants) {
    auto constant = pi();
    ASSERT_TRUE(constant);
    auto variable = SymbolicExpr::variable("pi");
    auto parsed_constant = round_trip(constant.value());
    auto parsed_variable = round_trip(variable);
    ASSERT_NE(parsed_constant, nullptr);
    ASSERT_NE(parsed_variable, nullptr);
    EXPECT_FALSE(symbol_name(parsed_constant).has_value());
    EXPECT_EQ(symbol_name(parsed_variable), std::optional<std::string_view>("pi"));
    EXPECT_FALSE(structurally_equal(*parsed_constant, *parsed_variable));
    auto constant_bytes = serialize_expr(constant.value());
    auto variable_bytes = serialize_expr(variable);
    ASSERT_TRUE(constant_bytes);
    ASSERT_TRUE(variable_bytes);
    EXPECT_NE(constant_bytes.value(), variable_bytes.value());

    auto unbound = evalf(*parsed_variable);
    ASSERT_FALSE(unbound);
    EXPECT_EQ(unbound.error().code, CasErrc::UnboundSymbol);
    const NumericBindings bindings{{"pi", 7.0}};
    auto constant_value = evalf(*parsed_constant, bindings);
    auto variable_value = evalf(*parsed_variable, bindings);
    ASSERT_TRUE(constant_value);
    ASSERT_TRUE(variable_value);
    EXPECT_NEAR(constant_value.value().value, LMMC_CONST_PI, 1e-15);
    EXPECT_EQ(variable_value.value().value, 7.0);

    const std::string binary_name("a\0b", 3);
    auto binary_symbol = round_trip(SymbolicExpr::variable(binary_name));
    ASSERT_NE(binary_symbol, nullptr);
    auto name = symbol_name(binary_symbol);
    ASSERT_TRUE(name.has_value());
    EXPECT_EQ(*name, binary_name);
}

TEST(ExprSerialization, FunctionsRelationsAndFiniteSets) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto call = function("f", {x, SymbolicExpr::number(Rational(3, 7))});
    ASSERT_TRUE(call);
    auto condition = logical_and(
        gt(call.value(), SymbolicExpr::number(0)).value(),
        logical_not(eq(y, SymbolicExpr::number(2)).value()).value());
    ASSERT_TRUE(condition);
    expect_same_value(condition.value());
    expect_same_value(SymbolicExpr::log(x, SymbolicExpr::number(10)));

    auto ordered = relation(x, y, RelationOp::LEQ);
    ASSERT_TRUE(ordered);
    auto parsed_relation = round_trip(ordered.value());
    ASSERT_NE(parsed_relation, nullptr);
    EXPECT_EQ(relation_op(parsed_relation), RelationOp::LEQ);
    EXPECT_TRUE(structurally_equal(*parsed_relation, *ordered.value()));

    auto set = finite_set({SymbolicExpr::number(2), x, SymbolicExpr::number(-1)});
    ASSERT_TRUE(set);
    expect_same_value(set.value());
}

TEST(ExprSerialization, ArithmeticComplexAndInterval) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::power(
        SymbolicExpr::add(x, SymbolicExpr::number(Rational(1, 3))),
        SymbolicExpr::number(-2));
    expect_same_value(expression);

    auto imaginary = complex(SymbolicExpr::number(BigInt("98765432109876543210")),
                             SymbolicExpr::number(Rational(-2, 5)));
    ASSERT_TRUE(imaginary);
    expect_same_value(imaginary.value());

    auto range = interval(SymbolicExpr::number(-2),
                          SymbolicExpr::number(5), false, true);
    ASSERT_TRUE(range);
    expect_same_value(range.value());
    auto in_range = membership(x, range.value());
    ASSERT_TRUE(in_range);
    expect_same_value(in_range.value());
}

TEST(ExprSerialization, BoundNamesAndScopes) {
    auto k = SymbolicFactory::create_variable("k");
    auto x = SymbolicFactory::create_variable("x");
    auto zero = SymbolicFactory::create_number(BigInt(0));
    auto n = SymbolicFactory::create_variable("n");
    auto body = SymbolicFactory::create_add({k, x});
    auto sum = detail::make_expression_ptr(
        detail::make_node<SummationNode>(body, "k", zero, n));
    auto product = detail::make_expression_ptr(
        detail::make_node<ProductNode>(body, "k", zero, n));
    expect_same_value(sum);
    expect_same_value(product);

    auto decoded = round_trip(sum);
    ASSERT_NE(decoded, nullptr);
    auto changed_free = decoded->substitute("x", SymbolicExpr::number(5));
    auto expected = sum->substitute("x", SymbolicExpr::number(5));
    ASSERT_NE(changed_free, nullptr);
    ASSERT_NE(expected, nullptr);
    EXPECT_TRUE(structurally_equal(*changed_free, *expected));
    auto changed_bound = decoded->substitute("k", SymbolicExpr::number(5));
    ASSERT_NE(changed_bound, nullptr);
    EXPECT_TRUE(structurally_equal(*changed_bound, *decoded));

    auto predicate = detail::make_node<RelationalNode>(k, x, RelationOp::GT);
    auto quantified = detail::make_expression_ptr(
        detail::make_node<QuantifierNode>(QuantifierNode::Type::ForAll,
                                          "k", n, predicate));
    auto builder = detail::make_expression_ptr(
        detail::make_node<SetBuilderNode>("k", n, predicate));
    expect_same_value(quantified);
    expect_same_value(builder);
}

TEST(ExprSerialization, MatrixLayoutIsNotPartOfTheValue) {
    auto zero = SymbolicFactory::create_number(BigInt(0));
    auto five = SymbolicFactory::create_number(BigInt(5));
    MatrixNode::DenseStorage entries{zero, zero, zero, zero, five, zero};
    MatrixNode::SparseStorage nonzero{{4, five}};
    auto dense = detail::make_expression_ptr(
        detail::make_node<MatrixNode>(2, 3, std::move(entries)));
    auto sparse = detail::make_expression_ptr(
        detail::make_node<MatrixNode>(2, 3, std::move(nonzero)));
    auto dense_bytes = serialize_expr(dense);
    auto sparse_bytes = serialize_expr(sparse);
    ASSERT_TRUE(dense_bytes);
    ASSERT_TRUE(sparse_bytes);
    EXPECT_EQ(dense_bytes.value(), sparse_bytes.value());

    auto decoded = parse_serialized_expr(sparse_bytes.value());
    ASSERT_TRUE(decoded);
    auto expected = SymbolicExpr::matrix({
        {SymbolicExpr::number(0), SymbolicExpr::number(0), SymbolicExpr::number(0)},
        {SymbolicExpr::number(0), SymbolicExpr::number(5), SymbolicExpr::number(0)}});
    EXPECT_TRUE(structurally_equal(*decoded.value(), *expected));
}

TEST(ExprSerialization, QuantitiesEncodePhysicalValueAndDimension) {
    const auto length = DimensionSignature::base("length");
    ComputationContext context;
    auto scaled = with_unit_definition(
        SymbolicExpr::number(Rational(3, 2)), "large-length",
        UnitDefinition{length, Rational(100)}, context);
    ASSERT_TRUE(scaled);
    auto base = with_unit_definition(
        SymbolicExpr::number(150), "base-length",
        UnitDefinition{length, Rational(1)}, context);
    ASSERT_TRUE(base);
    auto scaled_bytes = serialize_expr(scaled.value());
    auto base_bytes = serialize_expr(base.value());
    ASSERT_TRUE(scaled_bytes);
    ASSERT_TRUE(base_bytes);
    EXPECT_EQ(scaled_bytes.value(), base_bytes.value());

    auto decoded = parse_serialized_expr(scaled_bytes.value());
    ASSERT_TRUE(decoded);
    auto dimension = dimension_of(*decoded.value());
    ASSERT_TRUE(dimension);
    EXPECT_EQ(dimension.value(), length);
    auto physical_value = strip_to_base_value(decoded.value(), context);
    ASSERT_TRUE(physical_value);
    EXPECT_EQ(physical_value.value()->simplify()->convert_rational(), Rational(150));

    auto different_dimension = with_unit_definition(
        SymbolicExpr::number(150), "base-time",
        UnitDefinition{DimensionSignature::base("time"), Rational(1)}, context);
    ASSERT_TRUE(different_dimension);
    auto different_bytes = serialize_expr(different_dimension.value());
    ASSERT_TRUE(different_bytes);
    EXPECT_NE(scaled_bytes.value(), different_bytes.value());
}

TEST(ExprSerialization, RootIntegralAndDirectionalLimit) {
    auto x = SymbolicExpr::variable("x");
    auto polynomial = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::number(-2));
    auto root = SymbolicExpr::root_of(polynomial, "x", 1);
    auto decoded_root = round_trip(root);
    ASSERT_NE(decoded_root, nullptr);
    EXPECT_TRUE(structurally_equal(*decoded_root, *root));
    EXPECT_EQ(decoded_root->to_string(), root->to_string());

    auto body = SymbolicExpr::add(x, SymbolicExpr::variable("y"));
    auto integral = SymbolicExpr::make_integral(body, "x");
    expect_same_value(integral);
    auto bounded = detail::make_expression_ptr(detail::make_node<IntegralNode>(
        detail::node(body), "x", detail::node(SymbolicExpr::number(0)),
        detail::node(SymbolicExpr::number(4))));
    expect_same_value(bounded);

    auto limit = detail::make_expression_ptr(detail::make_node<LimitNode>(
        detail::node(body), "x", detail::node(SymbolicExpr::number(0)),
        LimitDirection::FromAbove));
    expect_same_value(limit);
}

TEST(ExprSerialization, MalformedAndUnknownInput) {
    auto encoded = serialize_expr(SymbolicExpr::variable("x"));
    ASSERT_TRUE(encoded);
    const std::string& valid = encoded.value();
    auto version = parse_serialized_expr(
        std::string("LMCAS_EXPR/2\n") + valid.substr(std::string("LMCAS_EXPR/1\n").size()));
    ASSERT_FALSE(version);
    EXPECT_EQ(version.error().code, CasErrc::ParseError);

    auto truncated = parse_serialized_expr(valid.substr(0, valid.size() - 1));
    ASSERT_FALSE(truncated);
    EXPECT_EQ(truncated.error().code, CasErrc::ParseError);
    auto trailing = parse_serialized_expr(valid + "extra");
    ASSERT_FALSE(trailing);
    EXPECT_EQ(trailing.error().code, CasErrc::ParseError);
    auto bad_length = parse_serialized_expr("LMCAS_EXPR/1\n99:x,");
    ASSERT_FALSE(bad_length);
    EXPECT_EQ(bad_length.error().code, CasErrc::ParseError);
    auto unknown = parse_serialized_expr("LMCAS_EXPR/1\n7:unknown,");
    ASSERT_FALSE(unknown);
    EXPECT_EQ(unknown.error().code, CasErrc::UnsupportedExpression);
}

TEST(ExprSerialization, ContextErrorsOnEncodingAndDecoding) {
    auto source = SymbolicExpr::add(
        SymbolicExpr::variable("x"), SymbolicExpr::number(1));
    auto encoded = serialize_expr(source);
    ASSERT_TRUE(encoded);

    ResourceLimits steps;
    steps.max_steps = 0;
    ComputationContext no_encode_steps(steps);
    auto encode_limited = serialize_expr(source, no_encode_steps);
    ASSERT_FALSE(encode_limited);
    EXPECT_EQ(encode_limited.error().code, CasErrc::ResourceLimit);
    ComputationContext no_decode_steps(steps);
    auto decode_limited = parse_serialized_expr(encoded.value(), no_decode_steps);
    ASSERT_FALSE(decode_limited);
    EXPECT_EQ(decode_limited.error().code, CasErrc::ResourceLimit);

    ResourceLimits bytes;
    bytes.max_input_bytes = encoded.value().size() - 1;
    ComputationContext short_input(bytes);
    auto input_limited = parse_serialized_expr(encoded.value(), short_input);
    ASSERT_FALSE(input_limited);
    EXPECT_EQ(input_limited.error().code, CasErrc::ResourceLimit);

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_encode({}, cancellation);
    auto encode_cancelled = serialize_expr(source, cancelled_encode);
    ASSERT_FALSE(encode_cancelled);
    EXPECT_EQ(encode_cancelled.error().code, CasErrc::Cancelled);
    ComputationContext cancelled_decode({}, cancellation);
    auto decode_cancelled = parse_serialized_expr(encoded.value(), cancelled_decode);
    ASSERT_FALSE(decode_cancelled);
    EXPECT_EQ(decode_cancelled.error().code, CasErrc::Cancelled);
}

} // namespace
