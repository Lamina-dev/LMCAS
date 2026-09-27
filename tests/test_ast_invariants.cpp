#include "test_common.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "expr.hpp"

#include <limits>
#include <stdexcept>
#include <thread>
#include <type_traits>

using namespace LMCAS;

namespace {

bool node_depends_on(
    const std::shared_ptr<const SymbolicNode> &node,
    const std::string &variable) {
    return LMCAS::contains(
        LMCAS::detail::expression_from_node(node), variable);
}

template <typename Fn>
void expect_invalid(Fn &&fn, const std::string &label) {
    bool rejected = false;
    try {
        fn();
    } catch (const std::invalid_argument &) {
        rejected = true;
    } catch (const std::length_error &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << label;
}

TEST(AstInvariants, NullRoot) {
    expect_invalid([]() {
        (void)LMCAS::detail::expression_from_node(std::shared_ptr<const SymbolicNode>{});
    },
                   "SymbolicExpr rejects null root construction");
}

TEST(AstInvariants, NullChildren) {
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<AddNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{one, nullptr});
    },
                   "AddNode rejects null operands");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{one, nullptr});
    },
                   "MultiplyNode rejects null operands");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<PowerNode>(one, nullptr);
    },
                   "PowerNode rejects null exponent");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin,
            std::vector<std::shared_ptr<const SymbolicNode>>{nullptr});
    },
                   "FunctionNode rejects null arguments");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<ComplexNode>(one, nullptr);
    },
                   "ComplexNode rejects null part");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<RelationalNode>(
            one, nullptr, RelationalNode::Op::EQ);
    },
                   "RelationalNode rejects null operand");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<LogicalNode>(
            one, nullptr, LogicalNode::Op::And);
    },
                   "LogicalNode rejects null binary operand");
}

TEST(AstInvariants, InvalidPiecewiseAndIterationChildren) {
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{});
    },
                   "PiecewiseNode rejects empty branch list");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{{one, nullptr}});
    },
                   "PiecewiseNode rejects null condition");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<SummationNode>(nullptr, "k", zero, one);
    },
                   "SummationNode rejects null body");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<SummationNode>(one, "", zero, one);
    },
                   "SummationNode rejects empty index variable");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<ProductNode>(one, "k", zero, nullptr);
    },
                   "ProductNode rejects null bound");
}

TEST(AstInvariants, InvalidTransformAndSetChildren) {
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    auto condition = LMCAS::detail::make_node<RelationalNode>(one, zero, RelationalNode::Op::GT);
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<TransformNode>(
            TransformNode::TransformType::Laplace, nullptr, "t",
            SymbolicFactory::create_variable("s"));
    },
                   "TransformNode rejects null body");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<TransformNode>(
            TransformNode::TransformType::Laplace, one, "",
            SymbolicFactory::create_variable("s"));
    },
                   "TransformNode rejects empty source variable");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<QuantifierNode>(
            QuantifierNode::Type::ForAll, "x", nullptr, condition);
    },
                   "QuantifierNode rejects null domain");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<SetBuilderNode>("", one, condition);
    },
                   "SetBuilderNode rejects empty element variable");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<FiniteSetNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{nullptr});
    },
                   "FiniteSetNode rejects null elements");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<IntervalNode>(nullptr, one, true, true);
    },
                   "IntervalNode rejects null endpoints");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MembershipNode>(one, nullptr);
    },
                   "MembershipNode rejects null sets");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<QuantityNode>(
            one, LMCAS::DimensionSignature::base("m"), Rational(0), "m");
    },
                   "QuantityNode rejects zero scale");
}

TEST(AstInvariants, NullFactoryChildren) {
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    expect_invalid([&]() {
        (void)SymbolicFactory::create_add({one, nullptr});
    },
                   "create_add rejects null operands");
    expect_invalid([&]() {
        (void)SymbolicFactory::create_multiply({one, nullptr});
    },
                   "create_multiply rejects null operands");
    expect_invalid([&]() {
        (void)SymbolicFactory::create_power(one, nullptr);
    },
                   "create_power rejects null operands");
    expect_invalid([&]() {
        (void)SymbolicFactory::create_complex(one, nullptr);
    },
                   "create_complex rejects null operands");
}

TEST(AstInvariants, NonfiniteNumbers) {
    expect_invalid([&]() {
        (void)SymbolicExpr::number(
            std::numeric_limits<double>::quiet_NaN());
    },
                   "SymbolicExpr rejects NaN");
    expect_invalid([&]() {
        (void)SymbolicExpr::number(
            std::numeric_limits<double>::infinity());
    },
                   "SymbolicExpr rejects positive infinity");
    expect_invalid([&]() {
        (void)SymbolicExpr::number(
            -std::numeric_limits<double>::infinity());
    },
                   "SymbolicExpr rejects negative infinity");
    EXPECT_TRUE((SymbolicExpr::infinity() != nullptr)) << "explicit symbolic infinity remains valid";
}

TEST(AstInvariants, MatrixShapeErrors) {
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{});
    },
                   "MatrixNode rejects empty grid");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{{one}, {one, one}});
    },
                   "MatrixNode rejects ragged grid");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{{one, nullptr}});
    },
                   "MatrixNode rejects null grid elements");
    expect_invalid([&]() {
        MatrixNode::DenseStorage dense = {one};
        (void)LMCAS::detail::make_node<MatrixNode>(1, 2, std::move(dense));
    },
                   "MatrixNode rejects dense storage size mismatch");
    expect_invalid([&]() {
        MatrixNode::SparseStorage sparse;
        sparse[2] = one;
        (void)LMCAS::detail::make_node<MatrixNode>(1, 2, std::move(sparse));
    },
                   "MatrixNode rejects sparse index out of bounds");
    expect_invalid([&]() {
        MatrixNode::DenseStorage dense;
        (void)LMCAS::detail::make_node<MatrixNode>(
            std::numeric_limits<size_t>::max(), 2, std::move(dense));
    },
                   "MatrixNode rejects dimension overflow");
}

TEST(AstInvariants, NullExpressions) {
    auto x = SymbolicExpr::variable("x");
    std::shared_ptr<SymbolicExpr> null_expr;
    expect_invalid([&]() {
        (void)SymbolicExpr::add(x, null_expr);
    },
                   "SymbolicExpr::add rejects null expression");
    expect_invalid([&]() {
        (void)SymbolicExpr::power(null_expr, x);
    },
                   "SymbolicExpr::power rejects null expression");
    expect_invalid([&]() {
        (void)SymbolicExpr::matrix({{x, null_expr}});
    },
                   "SymbolicExpr::matrix rejects null expression");
}

TEST(AstInvariants, SummationBinding) {
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto k_node = LMCAS::detail::make_node<VariableNode>("k");
    auto n_node = LMCAS::detail::make_node<VariableNode>("n");
    auto sum_body = SymbolicFactory::create_add({k_node, x_node});
    auto sum = LMCAS::detail::make_node<SummationNode>(sum_body, "k", zero, n_node);
    EXPECT_TRUE((!node_depends_on(sum, "k"))) << "Summation bound variable is not free";
    EXPECT_TRUE((node_depends_on(sum, "x"))) << "Summation body free variable is detected";
    EXPECT_TRUE((node_depends_on(sum, "n"))) << "Summation bound free variable is detected";

    auto sum_expr = LMCAS::detail::make_expression_ptr(sum);
    auto substituted_bound = sum_expr->substitute("k", SymbolicExpr::number(5));
    auto substituted_sum = std::dynamic_pointer_cast<const SummationNode>(LMCAS::detail::node(substituted_bound));
    EXPECT_TRUE((substituted_sum != nullptr)) << "Substitution preserves SummationNode";
    EXPECT_TRUE((node_depends_on(substituted_sum->body(), "k"))) << "Substitution does not replace bound summation variable in body";

    auto substituted_free = sum_expr->substitute("x", SymbolicExpr::number(5));
    auto substituted_free_sum = std::dynamic_pointer_cast<const SummationNode>(LMCAS::detail::node(substituted_free));
    EXPECT_TRUE((substituted_free_sum != nullptr)) << "Free substitution preserves SummationNode";
    EXPECT_TRUE((!node_depends_on(substituted_free_sum->body(), "x"))) << "Substitution replaces free variable inside summation body";
}

TEST(AstInvariants, PiecewiseExpansion) {
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto y_node = LMCAS::detail::make_node<VariableNode>("y");
    auto x_plus_one = SymbolicFactory::create_add({x_node, one});
    auto y_plus_one = SymbolicFactory::create_add({y_node, one});
    auto product = SymbolicFactory::create_multiply({x_plus_one, y_plus_one});
    auto positive_x = LMCAS::detail::make_node<RelationalNode>(x_node, zero, RelationalNode::Op::GT);
    auto piecewise = LMCAS::detail::make_node<PiecewiseNode>(
        std::vector<PiecewiseNode::Branch>{{product, positive_x}});
    auto expanded = LMCAS::detail::make_expression_ptr(piecewise)->expand();
    auto expanded_piecewise = std::dynamic_pointer_cast<const PiecewiseNode>(LMCAS::detail::node(expanded));
    EXPECT_TRUE((expanded_piecewise != nullptr)) << "Expand preserves PiecewiseNode";
    EXPECT_TRUE((std::dynamic_pointer_cast<const AddNode>(expanded_piecewise->branches()[0].expression) != nullptr)) << "Expand traverses and expands PiecewiseNode branch expression";
}

TEST(AstInvariants, QuantifierBinding) {
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto domain = LMCAS::detail::make_node<VariableNode>("R");
    auto y_node_for_predicate = LMCAS::detail::make_node<VariableNode>("y");
    auto predicate = LMCAS::detail::make_node<RelationalNode>(x_node, y_node_for_predicate, RelationalNode::Op::GT);
    auto quantified = LMCAS::detail::make_node<QuantifierNode>(
        QuantifierNode::Type::ForAll, "x", domain, predicate);
    auto quantified_expr = LMCAS::detail::make_expression_ptr(quantified);
    auto q_bound_sub = quantified_expr->substitute("x", SymbolicExpr::number(3));
    auto q_bound = std::dynamic_pointer_cast<const QuantifierNode>(LMCAS::detail::node(q_bound_sub));
    EXPECT_TRUE((q_bound != nullptr)) << "Substitution preserves QuantifierNode";
    EXPECT_TRUE((node_depends_on(q_bound->predicate(), "x"))) << "Substitution does not replace bound quantifier variable";
    auto q_free_sub = quantified_expr->substitute("y", SymbolicExpr::number(3));
    auto q_free = std::dynamic_pointer_cast<const QuantifierNode>(LMCAS::detail::node(q_free_sub));
    EXPECT_TRUE((q_free != nullptr)) << "Free substitution preserves QuantifierNode";
    EXPECT_TRUE((!node_depends_on(q_free->predicate(), "y"))) << "Substitution replaces free variable in quantifier predicate";
}

TEST(AstInvariants, SumProductCaptureAvoidance) {
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto k_node = LMCAS::detail::make_node<VariableNode>("k");
    auto n_node = LMCAS::detail::make_node<VariableNode>("n");
    auto replacement_k = SymbolicFactory::create_variable("k");
    auto capture_source = LMCAS::detail::make_node<SummationNode>(
        SymbolicFactory::create_add({k_node, x_node}), "k", x_node, n_node);
    auto capture_result = std::dynamic_pointer_cast<const SummationNode>(
        LMCAS::substitute_free(capture_source, "x", replacement_k));
    EXPECT_TRUE((capture_result && capture_result->index_var() != "k")) << "Summation binder alpha-renames to avoid capture";
    EXPECT_TRUE((capture_result && node_depends_on(capture_result->body(), "k"))) << "Replacement variable remains free after alpha-renaming";
    EXPECT_TRUE((capture_result && node_depends_on(capture_result->lower_bound(), "k"))) << "Summation bounds remain outside binder scope";
    auto product_source = LMCAS::detail::make_node<ProductNode>(
        SymbolicFactory::create_add({k_node, x_node}), "k", x_node, n_node);
    auto product_result = std::dynamic_pointer_cast<const ProductNode>(
        LMCAS::substitute_free(product_source, "x", replacement_k));
    EXPECT_TRUE((product_result && product_result->index_var() != "k" &&
                 node_depends_on(product_result->lower_bound(), "k")))
        << "Product alpha-renames its body binder but rewrites bounds";
}

TEST(AstInvariants, IntegralTransformCaptureAvoidance) {
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto k_node = LMCAS::detail::make_node<VariableNode>("k");
    auto n_node = LMCAS::detail::make_node<VariableNode>("n");
    auto replacement_k = SymbolicFactory::create_variable("k");
    auto integral = LMCAS::detail::make_node<IntegralNode>(
        SymbolicFactory::create_add({x_node, k_node}), "k", x_node, n_node);
    auto substituted_integral = std::dynamic_pointer_cast<const IntegralNode>(
        LMCAS::substitute_free(integral, "x", replacement_k));
    EXPECT_TRUE((substituted_integral && substituted_integral->variable() != "k")) << "Integral binder alpha-renames to avoid capture";
    EXPECT_TRUE((substituted_integral &&
                 node_depends_on(substituted_integral->lower(), "k")))
        << "Integral bounds remain outside binder scope";

    auto transform = LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::Laplace,
        SymbolicFactory::create_add({x_node, k_node}), "k", x_node);
    auto substituted_transform = std::dynamic_pointer_cast<const TransformNode>(
        LMCAS::substitute_free(transform, "x", replacement_k));
    EXPECT_TRUE((substituted_transform &&
                 substituted_transform->source_var() != "k"))
        << "Transform source binder alpha-renames to avoid capture";
    EXPECT_TRUE((substituted_transform &&
                 node_depends_on(substituted_transform->target(), "k")))
        << "Transform target is a free structural child";
}

TEST(AstInvariants, SetQuantifierCaptureAvoidance) {
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto k_node = LMCAS::detail::make_node<VariableNode>("k");
    auto replacement_k = SymbolicFactory::create_variable("k");
    auto set_builder = LMCAS::detail::make_node<SetBuilderNode>(
        "k", x_node,
        LMCAS::detail::make_node<RelationalNode>(
            x_node, k_node, RelationalNode::Op::GT));
    auto substituted_set = std::dynamic_pointer_cast<const SetBuilderNode>(
        LMCAS::substitute_free(set_builder, "x", replacement_k));
    EXPECT_TRUE((substituted_set && substituted_set->element_var() != "k")) << "Set-builder binder alpha-renames to avoid capture";
    EXPECT_TRUE((substituted_set && node_depends_on(substituted_set->domain(), "k"))) << "Set-builder domain remains outside binder scope";
    auto quantified_domain = LMCAS::detail::make_node<QuantifierNode>(
        QuantifierNode::Type::ForAll, "x", x_node,
        LMCAS::detail::make_node<RelationalNode>(
            x_node, zero, RelationalNode::Op::GT));
    auto substituted_quantified_domain =
        std::dynamic_pointer_cast<const QuantifierNode>(
            LMCAS::substitute_free(
                quantified_domain, "x",
                SymbolicFactory::create_number(BigInt(9))));
    EXPECT_TRUE((substituted_quantified_domain &&
                 substituted_quantified_domain->domain()->is_number() &&
                 node_depends_on(
                     substituted_quantified_domain->predicate(), "x")))
        << "Quantifier domain is outside its predicate binder";
}

TEST(AstInvariants, NestedShadowing) {
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto k_node = LMCAS::detail::make_node<VariableNode>("k");
    auto n_node = LMCAS::detail::make_node<VariableNode>("n");
    auto nested_shadow = LMCAS::detail::make_node<SummationNode>(
        LMCAS::detail::make_node<ProductNode>(
            SymbolicFactory::create_add({k_node, x_node}), "k", zero, n_node),
        "k", zero, n_node);
    auto shadow_result = std::dynamic_pointer_cast<const SummationNode>(
        LMCAS::substitute_free(
            nested_shadow, "k", SymbolicFactory::create_number(BigInt(4))));
    auto shadow_product = shadow_result
                              ? std::dynamic_pointer_cast<const ProductNode>(shadow_result->body())
                              : nullptr;
    EXPECT_TRUE((shadow_product &&
                 node_depends_on(shadow_product->body(), "k")))
        << "Nested shadowing prevents substitution in both scoped bodies";
}

TEST(AstInvariants, NestedFreeVariableLookup) {
    auto k = SymbolicFactory::create_variable("k");
    auto x = SymbolicFactory::create_variable("x");
    auto zero = SymbolicFactory::create_number(BigInt(0));
    auto one = SymbolicFactory::create_number(BigInt(1));
    auto inner = detail::make_node<ProductNode>(
        SymbolicFactory::create_add({k, x}), "k", zero, k);
    EXPECT_TRUE((expression_depends_on_variable(inner, "k"))) << "a binder's own upper bound can contain its name freely";
    auto outer = detail::make_node<SummationNode>(inner, "k", zero, one);
    EXPECT_TRUE((!expression_depends_on_variable(outer, "k"))) << "the enclosing binder binds the inner binder's outside children";
    EXPECT_TRUE((expression_depends_on_variable(outer, "x"))) << "unrelated names remain free through repeated shadowing";
    EXPECT_TRUE((free_variables(outer) == std::set<std::string>{"x"})) << "free-variable collection agrees with nested single-name lookup";
    auto free_outer_bound = detail::make_node<SummationNode>(inner, "k", zero, k);
    EXPECT_TRUE((expression_depends_on_variable(free_outer_bound, "k") &&
                 free_variables(free_outer_bound) == (std::set<std::string>{"k", "x"})))
        << "the outermost bound remains outside both lexical bodies";
}

TEST(AstInvariants, LimitCaptureAvoidance) {
    auto k = SymbolicFactory::create_variable("k");
    auto x = SymbolicFactory::create_variable("x");
    auto source = detail::make_node<LimitNode>(
        SymbolicFactory::create_add({k, x}), "k", x, LimitDirection::FromAbove);
    auto replaced = std::dynamic_pointer_cast<const LimitNode>(substitute_free(source, "x", k));
    EXPECT_TRUE((replaced && replaced->variable() != "k" &&
                 expression_depends_on_variable(replaced->body(), "k")))
        << "alpha-renaming protects the replacement from the limit's binder";
    EXPECT_TRUE((replaced && replaced->point()->equals(*k) &&
                 free_variables(replaced) == std::set<std::string>{"k"}))
        << "the point is substituted outside scope without leaking the fresh binder";
    auto bound_replaced = std::dynamic_pointer_cast<const LimitNode>(
        substitute_free(source, "k", SymbolicFactory::create_number(BigInt(7))));
    EXPECT_TRUE((bound_replaced && bound_replaced->body()->equals(*source->body()))) << "substitution of a bound limit name leaves its body unchanged";
}

TEST(AstInvariants, NestedMatrixTraversal) {
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto k_node = LMCAS::detail::make_node<VariableNode>("k");
    auto n_node = LMCAS::detail::make_node<VariableNode>("n");
    auto nested_sum = LMCAS::detail::make_node<SummationNode>(
        SymbolicFactory::create_add({k_node, x_node}), "k", zero, n_node);
    auto uninterpreted = LMCAS::detail::make_node<UninterpretedFunctionNode>(
        "f", std::vector<std::shared_ptr<const SymbolicNode>>{nested_sum});
    auto matrix = LMCAS::detail::make_node<MatrixNode>(
        1, 1, MatrixNode::DenseStorage{uninterpreted});
    const auto matrix_free = LMCAS::free_variables(matrix);
    EXPECT_TRUE((matrix_free.count("x") == 1 && matrix_free.count("n") == 1 &&
                 matrix_free.count("k") == 0))
        << "Free-variable traversal reaches matrix and function arguments";
    auto matrix_substituted = LMCAS::substitute_free(
        matrix, "x", SymbolicFactory::create_number(BigInt(7)));
    EXPECT_TRUE((!LMCAS::expression_depends_on_variable(matrix_substituted, "x"))) << "Substitution reaches matrix and function arguments";
}

TEST(AstInvariants, TraversalDepth) {
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    std::shared_ptr<const SymbolicNode> deep = x_node;
    for (int depth = 0; depth < 205; ++depth) {
        deep = LMCAS::detail::make_node<PowerNode>(deep, one);
    }
    bool depth_rejected = false;
    try {
        LMCAS::detail::RecursiveSymbolicVisitor traversal;
        deep->accept(traversal);
    } catch (const std::runtime_error &) {
        depth_rejected = true;
    }
    EXPECT_TRUE((depth_rejected)) << "Recursive traversal rejects excessive depth";
}

TEST(AstInvariants, LimitRootBindingRoundtrip) {
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto y_node_for_predicate = LMCAS::detail::make_node<VariableNode>("y");
    auto scoped_limit = LMCAS::detail::make_node<LimitNode>(
        SymbolicFactory::create_add({x_node, y_node_for_predicate}),
        "x", x_node, LimitDirection::Both);
    const auto limit_free = LMCAS::free_variables(scoped_limit);
    EXPECT_TRUE((limit_free.count("x") == 1 && limit_free.count("y") == 1)) << "Limit binds its body variable but not its point";
    auto root_polynomial_expression = SymbolicExpr::add(
        SymbolicExpr::power(
            SymbolicExpr::variable("x"), SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    auto root_expression = SymbolicExpr::root_of(
        root_polynomial_expression, "x", 1);
    auto root_node = std::dynamic_pointer_cast<const RootOfNode>(
        LMCAS::detail::node(root_expression));
    const auto root_free = LMCAS::free_variables(root_node);
    EXPECT_TRUE((root_free.empty())) << "canonical RootOf identity has no dummy free variable";

    auto limit_node = LMCAS::detail::make_node<LimitNode>(
        SymbolicFactory::create_add({x_node, one}), "x", zero,
        LimitDirection::FromAbove);
    auto limit_expression = LMCAS::detail::make_expression_ptr(limit_node);
    auto parsed_limit = LMCAS::parse_expr(limit_expression->to_string());
    EXPECT_TRUE((parsed_limit &&
                 LMCAS::detail::node(parsed_limit.value())->equals(*limit_node)))
        << "LimitNode survives print/parse round trip";

    EXPECT_TRUE((root_node != nullptr)) << "checked RootOf construction returns a RootOfNode";
    auto parsed_root = LMCAS::parse_expr(root_expression->to_string());
    EXPECT_TRUE((parsed_root &&
                 LMCAS::detail::node(parsed_root.value())->equals(*root_node)))
        << "RootOfNode survives print/parse round trip";
}

TEST(AstInvariants, ConcurrentImmutableReads) {
    auto shared_expr = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    const auto expected_hash = LMCAS::detail::node(shared_expr)->hash();
    const auto expected_text = shared_expr->to_string();
    constexpr std::size_t worker_count = 8;
    std::vector<std::size_t> hashes(worker_count);
    std::vector<std::string> texts(worker_count);
    std::vector<int> equal(worker_count, 0);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (std::size_t i = 0; i < worker_count; ++i) {
        workers.emplace_back([&, i] {
            hashes[i] = LMCAS::detail::node(shared_expr)->hash();
            texts[i] = shared_expr->to_string();
            equal[i] = LMCAS::detail::node(shared_expr)->equals(*LMCAS::detail::node(shared_expr)) ? 1 : 0;
        });
    }
    for (auto &worker : workers)
        worker.join();
    bool deterministic = true;
    for (std::size_t i = 0; i < worker_count; ++i) {
        deterministic = deterministic && hashes[i] == expected_hash &&
                        texts[i] == expected_text && equal[i];
    }
    EXPECT_TRUE((deterministic)) << "shared immutable expressions have deterministic concurrent reads";
}

} // namespace
