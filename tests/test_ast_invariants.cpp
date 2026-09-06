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
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable) {
    return LMCAS::contains(
        LMCAS::detail::expression_from_node(node), variable);
}

template <typename Fn>
void expect_invalid(Fn&& fn, const std::string& label) {
    bool rejected = false;
    try {
        fn();
    } catch (const std::invalid_argument&) {
        rejected = true;
    } catch (const std::length_error&) {
        rejected = true;
    }
    EXPECT_TRUE(rejected, label);
}

class ExhaustiveProbeVisitor : public LMCAS::detail::SymbolicVisitor {
public:
    void visit(const NumberNode&) override {}
    void visit(const VariableNode&) override {}
    void visit(const AddNode&) override {}
    void visit(const MultiplyNode&) override {}
    void visit(const PowerNode&) override {}
    void visit(const FunctionNode&) override {}
    void visit(const UninterpretedFunctionNode&) override {}
    void visit(const MatrixNode&) override {}
    void visit(const RelationalNode&) override {}
    void visit(const LogicalNode&) override {}
    void visit(const PiecewiseNode&) override {}
    void visit(const SummationNode&) override {}
    void visit(const ProductNode&) override {}
    void visit(const TransformNode&) override {}
    void visit(const QuantifierNode&) override {}
    void visit(const SetBuilderNode&) override {}
    void visit(const ComplexNode&) override {}
    void visit(const FiniteSetNode&) override {}
    void visit(const IntervalNode&) override {}
    void visit(const MembershipNode&) override {}
    void visit(const QuantityNode&) override {}
    void visit(const IntegralNode&) override {}
    void visit(const LimitNode&) override {}
    void visit(const RootOfNode&) override {}
};

static_assert(!std::is_abstract<ExhaustiveProbeVisitor>::value,
              "LMCAS::detail::SymbolicVisitor implementations must explicitly cover every node type");

template <typename T, typename = void>
struct HasPublicRoot : std::false_type {};

template <typename T>
struct HasPublicRoot<T, std::void_t<decltype(std::declval<T>().root)>>
    : std::true_type {};

#define LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(TraitName, MemberName)             \
    template <typename T, typename = void>                                  \
    struct TraitName : std::false_type {};                                  \
    template <typename T>                                                   \
    struct TraitName<T, std::void_t<decltype(std::declval<T>().MemberName)>> \
        : std::true_type {}

LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicRealField, real);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicImagField, imag);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicBranchesField, branches);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicDefaultExprField, default_expr);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicBodyField, body);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicIndexVarField, index_var);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicLowerBoundField, lower_bound);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicUpperBoundField, upper_bound);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicTransformTypeField, transform_type);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicSourceVarField, source_var);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicTargetVarField, target_var);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicQuantifierTypeField, quantifier_type);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicBoundVarField, bound_var);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicElementVarField, element_var);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicDomainField, domain);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicPredicateField, predicate);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicOperandsField, operands);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicBaseField, base);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicExponentField, exponent);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicArgumentsField, arguments);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicTypeField, type);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicRowsField, rows);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicColsField, cols);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicStorageField, storage);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicLeftField, left);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicRightField, right);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicOpField, op);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicValueField, value);
LMCAS_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicNameField, name);

#undef LMCAS_DEFINE_PUBLIC_MEMBER_PROBE

using RootPointer = std::decay_t<decltype(
    LMCAS::detail::node(std::declval<const SymbolicExpr&>()))>;
using RootElement = typename RootPointer::element_type;
using AddOperandPointer = typename std::decay_t<
    decltype(std::declval<const AddNode&>().operands())>::value_type;
static_assert(std::is_const<RootElement>::value,
              "SymbolicExpr must own a const AST root");
static_assert(!HasPublicRoot<SymbolicExpr>::value,
              "SymbolicExpr must not expose a public root handle");
static_assert(!HasPublicRealField<ComplexNode>::value &&
              !HasPublicImagField<ComplexNode>::value,
              "ComplexNode storage must remain private");
static_assert(!HasPublicBranchesField<PiecewiseNode>::value &&
              !HasPublicDefaultExprField<PiecewiseNode>::value,
              "PiecewiseNode storage must remain private");
static_assert(!HasPublicBodyField<SummationNode>::value &&
              !HasPublicIndexVarField<SummationNode>::value &&
              !HasPublicLowerBoundField<SummationNode>::value &&
              !HasPublicUpperBoundField<SummationNode>::value,
              "SummationNode storage must remain private");
static_assert(!HasPublicBodyField<ProductNode>::value &&
              !HasPublicIndexVarField<ProductNode>::value &&
              !HasPublicLowerBoundField<ProductNode>::value &&
              !HasPublicUpperBoundField<ProductNode>::value,
              "ProductNode storage must remain private");
static_assert(!HasPublicTransformTypeField<TransformNode>::value &&
              !HasPublicBodyField<TransformNode>::value &&
              !HasPublicSourceVarField<TransformNode>::value &&
              !HasPublicTargetVarField<TransformNode>::value,
              "TransformNode storage must remain private");
static_assert(!HasPublicQuantifierTypeField<QuantifierNode>::value &&
              !HasPublicBoundVarField<QuantifierNode>::value &&
              !HasPublicDomainField<QuantifierNode>::value &&
              !HasPublicPredicateField<QuantifierNode>::value,
              "QuantifierNode storage must remain private");
static_assert(!HasPublicElementVarField<SetBuilderNode>::value &&
              !HasPublicDomainField<SetBuilderNode>::value &&
              !HasPublicPredicateField<SetBuilderNode>::value,
              "SetBuilderNode storage must remain private");
static_assert(!HasPublicOperandsField<AddNode>::value,
              "AddNode storage must remain private");
static_assert(!HasPublicOperandsField<MultiplyNode>::value,
              "MultiplyNode storage must remain private");
static_assert(!HasPublicBaseField<PowerNode>::value &&
              !HasPublicExponentField<PowerNode>::value,
              "PowerNode storage must remain private");
static_assert(!HasPublicTypeField<FunctionNode>::value &&
              !HasPublicArgumentsField<FunctionNode>::value,
              "FunctionNode storage must remain private");
static_assert(!HasPublicRowsField<MatrixNode>::value &&
              !HasPublicColsField<MatrixNode>::value &&
              !HasPublicStorageField<MatrixNode>::value,
              "MatrixNode storage must remain private");
static_assert(!HasPublicLeftField<RelationalNode>::value &&
              !HasPublicRightField<RelationalNode>::value &&
              !HasPublicOpField<RelationalNode>::value,
              "RelationalNode storage must remain private");
static_assert(!HasPublicLeftField<LogicalNode>::value &&
              !HasPublicRightField<LogicalNode>::value &&
              !HasPublicOpField<LogicalNode>::value,
              "LogicalNode storage must remain private");
static_assert(!HasPublicValueField<NumberNode>::value,
              "NumberNode storage must remain private");
static_assert(!HasPublicNameField<VariableNode>::value,
              "VariableNode storage must remain private");
static_assert(std::is_same<decltype(std::declval<const NumberNode&>().value()),
                           const std::variant<BigInt, Rational, lmmc_real_t>&>::value,
              "NumberNode must expose its value through a const reference");
static_assert(std::is_same<decltype(std::declval<const VariableNode&>().name()),
                           const std::string&>::value,
              "VariableNode must expose its name through a const reference");
static_assert(std::is_same<decltype(std::declval<const ComplexNode&>().real()),
                           const std::shared_ptr<const SymbolicNode>&>::value,
              "AST child accessors must preserve const ownership");
static_assert(std::is_same<decltype(std::declval<const AddNode&>().operands()),
                           const std::vector<std::shared_ptr<const SymbolicNode>>&>::value,
              "variadic AST accessors must preserve const ownership");
static_assert(!std::is_default_constructible<SymbolicExpr>::value,
              "SymbolicExpr must always contain a valid AST root");
static_assert(
    !std::is_constructible<SymbolicExpr,
                           std::shared_ptr<const SymbolicNode>>::value,
    "SymbolicExpr must not expose raw AST construction");
static_assert(!std::is_constructible<SymbolicExpr, SymbolicExpr::Type>::value,
              "legacy type-only construction must not fabricate expressions");
static_assert(!std::is_constructible<NumberNode, BigInt>::value,
              "NumberNode construction must be factory-only");
static_assert(!std::is_constructible<VariableNode, std::string>::value,
              "VariableNode construction must be factory-only");
static_assert(!std::is_constructible<AddNode,
              std::vector<std::shared_ptr<const SymbolicNode>>>::value,
              "AddNode construction must be factory-only");
static_assert(!std::is_constructible<MultiplyNode,
              std::vector<std::shared_ptr<const SymbolicNode>>>::value,
              "MultiplyNode construction must be factory-only");
static_assert(!std::is_constructible<PowerNode,
              std::shared_ptr<const SymbolicNode>,
              std::shared_ptr<const SymbolicNode>>::value,
              "PowerNode construction must be factory-only");
static_assert(!std::is_constructible<FunctionNode, FunctionNode::FuncType,
              std::vector<std::shared_ptr<const SymbolicNode>>>::value,
              "FunctionNode construction must be factory-only");
static_assert(!std::is_constructible<MatrixNode,
              const std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>&>::value,
              "MatrixNode construction must be factory-only");
static_assert(!std::is_constructible<RelationalNode,
              std::shared_ptr<const SymbolicNode>,
              std::shared_ptr<const SymbolicNode>, RelationalNode::Op>::value,
              "RelationalNode construction must be factory-only");
static_assert(!std::is_constructible<LogicalNode,
              std::shared_ptr<const SymbolicNode>,
              std::shared_ptr<const SymbolicNode>, LogicalNode::Op>::value,
              "LogicalNode construction must be factory-only");
static_assert(!std::is_constructible<PiecewiseNode,
              std::vector<PiecewiseNode::Branch>,
              std::shared_ptr<const SymbolicNode>>::value,
              "PiecewiseNode construction must be factory-only");
static_assert(!std::is_constructible<SummationNode,
              std::shared_ptr<const SymbolicNode>, std::string,
              std::shared_ptr<const SymbolicNode>,
              std::shared_ptr<const SymbolicNode>>::value,
              "SummationNode construction must be factory-only");
static_assert(!std::is_constructible<ProductNode,
              std::shared_ptr<const SymbolicNode>, std::string,
              std::shared_ptr<const SymbolicNode>,
              std::shared_ptr<const SymbolicNode>>::value,
              "ProductNode construction must be factory-only");
static_assert(!std::is_constructible<TransformNode,
              TransformNode::TransformType,
              std::shared_ptr<const SymbolicNode>, std::string,
              std::string>::value,
              "TransformNode construction must be factory-only");
static_assert(!std::is_constructible<QuantifierNode, QuantifierNode::Type,
              std::string, std::shared_ptr<const SymbolicNode>,
              std::shared_ptr<const SymbolicNode>>::value,
              "QuantifierNode construction must be factory-only");
static_assert(!std::is_constructible<SetBuilderNode, std::string,
              std::shared_ptr<const SymbolicNode>,
              std::shared_ptr<const SymbolicNode>>::value,
              "SetBuilderNode construction must be factory-only");
static_assert(!std::is_constructible<ComplexNode,
              std::shared_ptr<const SymbolicNode>,
              std::shared_ptr<const SymbolicNode>>::value,
              "ComplexNode construction must be factory-only");
static_assert(std::is_const<typename AddOperandPointer::element_type>::value,
              "AST child edges must point to const nodes");
static_assert(std::is_same<decltype(&SymbolicNode::accept),
                           void (SymbolicNode::*)(LMCAS::detail::SymbolicVisitor&) const>::value,
              "AST traversal must not expose mutable nodes");

} // namespace

int main() {
    TEST_CASE("SymbolicExpr rejects an empty AST root");
    expect_invalid([]() {
        (void)LMCAS::detail::expression_from_node(std::shared_ptr<const SymbolicNode>{});
    }, "SymbolicExpr rejects null root construction");

    TEST_CASE("AST nodes reject null children");
    auto one = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<AddNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{one, nullptr});
    }, "AddNode rejects null operands");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{one, nullptr});
    }, "MultiplyNode rejects null operands");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<PowerNode>(one, nullptr);
    }, "PowerNode rejects null exponent");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin,
            std::vector<std::shared_ptr<const SymbolicNode>>{nullptr});
    }, "FunctionNode rejects null arguments");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<ComplexNode>(one, nullptr);
    }, "ComplexNode rejects null part");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<RelationalNode>(
            one, nullptr, RelationalNode::Op::EQ);
    }, "RelationalNode rejects null operand");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<LogicalNode>(
            one, nullptr, LogicalNode::Op::And);
    }, "LogicalNode rejects null binary operand");

    TEST_CASE("Extended AST nodes reject invalid children");
    auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
    auto condition = LMCAS::detail::make_node<RelationalNode>(one, zero, RelationalNode::Op::GT);
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{});
    }, "PiecewiseNode rejects empty branch list");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{{one, nullptr}});
    }, "PiecewiseNode rejects null condition");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<SummationNode>(nullptr, "k", zero, one);
    }, "SummationNode rejects null body");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<SummationNode>(one, "", zero, one);
    }, "SummationNode rejects empty index variable");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<ProductNode>(one, "k", zero, nullptr);
    }, "ProductNode rejects null bound");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<TransformNode>(
            TransformNode::TransformType::Laplace, nullptr, "t",
            SymbolicFactory::create_variable("s"));
    }, "TransformNode rejects null body");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<TransformNode>(
            TransformNode::TransformType::Laplace, one, "",
            SymbolicFactory::create_variable("s"));
    }, "TransformNode rejects empty source variable");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<QuantifierNode>(
            QuantifierNode::Type::ForAll, "x", nullptr, condition);
    }, "QuantifierNode rejects null domain");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<SetBuilderNode>("", one, condition);
    }, "SetBuilderNode rejects empty element variable");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<FiniteSetNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{nullptr});
    }, "FiniteSetNode rejects null elements");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<IntervalNode>(nullptr, one, true, true);
    }, "IntervalNode rejects null endpoints");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MembershipNode>(one, nullptr);
    }, "MembershipNode rejects null sets");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<QuantityNode>(
            one, LMCAS::DimensionSignature::base("m"), Rational(0), "m");
    }, "QuantityNode rejects zero scale");

    TEST_CASE("SymbolicFactory rejects null children");
    expect_invalid([&]() {
        (void)SymbolicFactory::create_add({one, nullptr});
    }, "create_add rejects null operands");
    expect_invalid([&]() {
        (void)SymbolicFactory::create_multiply({one, nullptr});
    }, "create_multiply rejects null operands");
    expect_invalid([&]() {
        (void)SymbolicFactory::create_power(one, nullptr);
    }, "create_power rejects null operands");
    expect_invalid([&]() {
        (void)SymbolicFactory::create_complex(one, nullptr);
    }, "create_complex rejects null operands");

    TEST_CASE("NumberNode rejects nonfinite approximate values");
    expect_invalid([&]() {
        (void)SymbolicExpr::number(
            std::numeric_limits<double>::quiet_NaN());
    }, "SymbolicExpr rejects NaN");
    expect_invalid([&]() {
        (void)SymbolicExpr::number(
            std::numeric_limits<double>::infinity());
    }, "SymbolicExpr rejects positive infinity");
    expect_invalid([&]() {
        (void)SymbolicExpr::number(
            -std::numeric_limits<double>::infinity());
    }, "SymbolicExpr rejects negative infinity");
    EXPECT_TRUE(SymbolicExpr::infinity() != nullptr,
                "explicit symbolic infinity remains valid");

    TEST_CASE("MatrixNode validates shape and storage");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{});
    }, "MatrixNode rejects empty grid");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{{one}, {one, one}});
    }, "MatrixNode rejects ragged grid");
    expect_invalid([&]() {
        (void)LMCAS::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{{one, nullptr}});
    }, "MatrixNode rejects null grid elements");
    expect_invalid([&]() {
        MatrixNode::DenseStorage dense = {one};
        (void)LMCAS::detail::make_node<MatrixNode>(1, 2, std::move(dense));
    }, "MatrixNode rejects dense storage size mismatch");
    expect_invalid([&]() {
        MatrixNode::SparseStorage sparse;
        sparse[2] = one;
        (void)LMCAS::detail::make_node<MatrixNode>(1, 2, std::move(sparse));
    }, "MatrixNode rejects sparse index out of bounds");
    expect_invalid([&]() {
        MatrixNode::DenseStorage dense;
        (void)LMCAS::detail::make_node<MatrixNode>(
            std::numeric_limits<size_t>::max(), 2, std::move(dense));
    }, "MatrixNode rejects dimension overflow");
    MatrixNode::SparseStorage sparse;
    sparse[0] = one;
    auto sparse_matrix =
        LMCAS::detail::make_node<MatrixNode>(2, 2, std::move(sparse));
    auto first_missing = sparse_matrix->get(0, 1);
    auto second_missing = sparse_matrix->get(1, 0);
    EXPECT_TRUE(first_missing.get() == second_missing.get(),
                "sparse missing coordinates share one immutable zero node");
    EXPECT_TRUE(first_missing.get() == sparse_matrix->get(0, 1).get(),
                "repeated sparse reads return the same zero node");

    TEST_CASE("SymbolicExpr factories reject null expressions");
    auto x = SymbolicExpr::variable("x");
    std::shared_ptr<SymbolicExpr> null_expr;
    expect_invalid([&]() {
        (void)SymbolicExpr::add(x, null_expr);
    }, "SymbolicExpr::add rejects null expression");
    expect_invalid([&]() {
        (void)SymbolicExpr::power(null_expr, x);
    }, "SymbolicExpr::power rejects null expression");
    expect_invalid([&]() {
        (void)SymbolicExpr::matrix({{x, null_expr}});
    }, "SymbolicExpr::matrix rejects null expression");

    TEST_CASE("Visitor-dependent helpers respect bound variables");
    auto x_node = LMCAS::detail::make_node<VariableNode>("x");
    auto k_node = LMCAS::detail::make_node<VariableNode>("k");
    auto n_node = LMCAS::detail::make_node<VariableNode>("n");
    auto sum_body = SymbolicFactory::create_add({k_node, x_node});
    auto sum = LMCAS::detail::make_node<SummationNode>(sum_body, "k", zero, n_node);
    EXPECT_TRUE(!node_depends_on(sum, "k"), "Summation bound variable is not free");
    EXPECT_TRUE(node_depends_on(sum, "x"), "Summation body free variable is detected");
    EXPECT_TRUE(node_depends_on(sum, "n"), "Summation bound free variable is detected");

    auto sum_expr = LMCAS::detail::make_expression_ptr(sum);
    auto substituted_bound = sum_expr->substitute("k", SymbolicExpr::number(5));
    auto substituted_sum = std::dynamic_pointer_cast<const SummationNode>(LMCAS::detail::node(substituted_bound));
    EXPECT_TRUE(substituted_sum != nullptr, "Substitution preserves SummationNode");
    EXPECT_TRUE(node_depends_on(substituted_sum->body(), "k"),
                "Substitution does not replace bound summation variable in body");

    auto substituted_free = sum_expr->substitute("x", SymbolicExpr::number(5));
    auto substituted_free_sum = std::dynamic_pointer_cast<const SummationNode>(LMCAS::detail::node(substituted_free));
    EXPECT_TRUE(substituted_free_sum != nullptr, "Free substitution preserves SummationNode");
    EXPECT_TRUE(!node_depends_on(substituted_free_sum->body(), "x"),
                "Substitution replaces free variable inside summation body");

    TEST_CASE("Extended nodes are traversed by expand");
    auto y_node = LMCAS::detail::make_node<VariableNode>("y");
    auto x_plus_one = SymbolicFactory::create_add({x_node, one});
    auto y_plus_one = SymbolicFactory::create_add({y_node, one});
    auto product = SymbolicFactory::create_multiply({x_plus_one, y_plus_one});
    auto positive_x = LMCAS::detail::make_node<RelationalNode>(x_node, zero, RelationalNode::Op::GT);
    auto piecewise = LMCAS::detail::make_node<PiecewiseNode>(
        std::vector<PiecewiseNode::Branch>{{product, positive_x}});
    auto expanded = LMCAS::detail::make_expression_ptr(piecewise)->expand();
    auto expanded_piecewise = std::dynamic_pointer_cast<const PiecewiseNode>(LMCAS::detail::node(expanded));
    EXPECT_TRUE(expanded_piecewise != nullptr, "Expand preserves PiecewiseNode");
    EXPECT_TRUE(std::dynamic_pointer_cast<const AddNode>(expanded_piecewise->branches()[0].expression) != nullptr,
                "Expand traverses and expands PiecewiseNode branch expression");

    TEST_CASE("Substitution respects quantifier binding");
    auto domain = LMCAS::detail::make_node<VariableNode>("R");
    auto y_node_for_predicate = LMCAS::detail::make_node<VariableNode>("y");
    auto predicate = LMCAS::detail::make_node<RelationalNode>(x_node, y_node_for_predicate, RelationalNode::Op::GT);
    auto quantified = LMCAS::detail::make_node<QuantifierNode>(
        QuantifierNode::Type::ForAll, "x", domain, predicate);
    auto quantified_expr = LMCAS::detail::make_expression_ptr(quantified);
    auto q_bound_sub = quantified_expr->substitute("x", SymbolicExpr::number(3));
    auto q_bound = std::dynamic_pointer_cast<const QuantifierNode>(LMCAS::detail::node(q_bound_sub));
    EXPECT_TRUE(q_bound != nullptr, "Substitution preserves QuantifierNode");
    EXPECT_TRUE(node_depends_on(q_bound->predicate(), "x"),
                "Substitution does not replace bound quantifier variable");
    auto q_free_sub = quantified_expr->substitute("y", SymbolicExpr::number(3));
    auto q_free = std::dynamic_pointer_cast<const QuantifierNode>(LMCAS::detail::node(q_free_sub));
    EXPECT_TRUE(q_free != nullptr, "Free substitution preserves QuantifierNode");
    EXPECT_TRUE(!node_depends_on(q_free->predicate(), "y"),
                "Substitution replaces free variable in quantifier predicate");

    TEST_CASE("Capture-avoiding substitution covers every binder scope");
    auto replacement_k = SymbolicFactory::create_variable("k");
    auto capture_source = LMCAS::detail::make_node<SummationNode>(
        SymbolicFactory::create_add({k_node, x_node}), "k", x_node, n_node);
    auto capture_result = std::dynamic_pointer_cast<const SummationNode>(
        LMCAS::substitute_free(capture_source, "x", replacement_k));
    EXPECT_TRUE(capture_result && capture_result->index_var() != "k",
                "Summation binder alpha-renames to avoid capture");
    EXPECT_TRUE(capture_result && node_depends_on(capture_result->body(), "k"),
                "Replacement variable remains free after alpha-renaming");
    EXPECT_TRUE(capture_result && node_depends_on(capture_result->lower_bound(), "k"),
                "Summation bounds remain outside binder scope");
    auto product_source = LMCAS::detail::make_node<ProductNode>(
        SymbolicFactory::create_add({k_node, x_node}), "k", x_node, n_node);
    auto product_result = std::dynamic_pointer_cast<const ProductNode>(
        LMCAS::substitute_free(product_source, "x", replacement_k));
    EXPECT_TRUE(product_result && product_result->index_var() != "k" &&
                    node_depends_on(product_result->lower_bound(), "k"),
                "Product alpha-renames its body binder but rewrites bounds");


    auto integral = LMCAS::detail::make_node<IntegralNode>(
        SymbolicFactory::create_add({x_node, k_node}), "k", x_node, n_node);
    auto substituted_integral = std::dynamic_pointer_cast<const IntegralNode>(
        LMCAS::substitute_free(integral, "x", replacement_k));
    EXPECT_TRUE(substituted_integral && substituted_integral->variable() != "k",
                "Integral binder alpha-renames to avoid capture");
    EXPECT_TRUE(substituted_integral &&
                    node_depends_on(substituted_integral->lower(), "k"),
                "Integral bounds remain outside binder scope");

    auto transform = LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::Laplace,
        SymbolicFactory::create_add({x_node, k_node}), "k", x_node);
    auto substituted_transform = std::dynamic_pointer_cast<const TransformNode>(
        LMCAS::substitute_free(transform, "x", replacement_k));
    EXPECT_TRUE(substituted_transform &&
                    substituted_transform->source_var() != "k",
                "Transform source binder alpha-renames to avoid capture");
    EXPECT_TRUE(substituted_transform &&
                    node_depends_on(substituted_transform->target(), "k"),
                "Transform target is a free structural child");

    auto set_builder = LMCAS::detail::make_node<SetBuilderNode>(
        "k", x_node,
        LMCAS::detail::make_node<RelationalNode>(
            x_node, k_node, RelationalNode::Op::GT));
    auto substituted_set = std::dynamic_pointer_cast<const SetBuilderNode>(
        LMCAS::substitute_free(set_builder, "x", replacement_k));
    EXPECT_TRUE(substituted_set && substituted_set->element_var() != "k",
                "Set-builder binder alpha-renames to avoid capture");
    EXPECT_TRUE(substituted_set && node_depends_on(substituted_set->domain(), "k"),
                "Set-builder domain remains outside binder scope");
    auto quantified_domain = LMCAS::detail::make_node<QuantifierNode>(
        QuantifierNode::Type::ForAll, "x", x_node,
        LMCAS::detail::make_node<RelationalNode>(
            x_node, zero, RelationalNode::Op::GT));
    auto substituted_quantified_domain =
        std::dynamic_pointer_cast<const QuantifierNode>(
            LMCAS::substitute_free(
                quantified_domain, "x",
                SymbolicFactory::create_number(BigInt(9))));
    EXPECT_TRUE(substituted_quantified_domain &&
                    substituted_quantified_domain->domain()->is_number() &&
                    node_depends_on(
                        substituted_quantified_domain->predicate(), "x"),
                "Quantifier domain is outside its predicate binder");

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
    EXPECT_TRUE(shadow_product &&
                    node_depends_on(shadow_product->body(), "k"),
                "Nested shadowing prevents substitution in both scoped bodies");


    TEST_CASE("Traversal reaches matrices and uninterpreted functions");
    auto nested_sum = LMCAS::detail::make_node<SummationNode>(
        SymbolicFactory::create_add({k_node, x_node}), "k", zero, n_node);
    auto uninterpreted = LMCAS::detail::make_node<UninterpretedFunctionNode>(
        "f", std::vector<std::shared_ptr<const SymbolicNode>>{nested_sum});
    auto matrix = LMCAS::detail::make_node<MatrixNode>(
        1, 1, MatrixNode::DenseStorage{uninterpreted});
    const auto matrix_free = LMCAS::free_variables(matrix);
    EXPECT_TRUE(matrix_free.count("x") == 1 && matrix_free.count("n") == 1 &&
                    matrix_free.count("k") == 0,
                "Free-variable traversal reaches matrix and function arguments");
    auto matrix_substituted = LMCAS::substitute_free(
        matrix, "x", SymbolicFactory::create_number(BigInt(7)));
    EXPECT_TRUE(!LMCAS::expression_depends_on_variable(matrix_substituted, "x"),
                "Substitution reaches matrix and function arguments");

    TEST_CASE("Recursive traversal enforces depth limit");
    std::shared_ptr<const SymbolicNode> deep = x_node;
    for (int depth = 0; depth < 205; ++depth) {
        deep = LMCAS::detail::make_node<PowerNode>(deep, one);
    }
    bool depth_rejected = false;
    try {
        LMCAS::detail::RecursiveSymbolicVisitor traversal;
        deep->accept(traversal);
    } catch (const std::runtime_error&) {
        depth_rejected = true;
    }
    EXPECT_TRUE(depth_rejected, "Recursive traversal rejects excessive depth");
    TEST_CASE("Limit and RootOf expose precise binder scopes");
    auto scoped_limit = LMCAS::detail::make_node<LimitNode>(
        SymbolicFactory::create_add({x_node, y_node_for_predicate}),
        "x", x_node, LimitDirection::Both);
    const auto limit_free = LMCAS::free_variables(scoped_limit);
    EXPECT_TRUE(limit_free.count("x") == 1 && limit_free.count("y") == 1,
                "Limit binds its body variable but not its point");
    auto root_polynomial_expression = SymbolicExpr::add(
        SymbolicExpr::power(
            SymbolicExpr::variable("x"), SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    auto root_expression = SymbolicExpr::root_of(
        root_polynomial_expression, "x", 1);
    auto root_node = std::dynamic_pointer_cast<const RootOfNode>(
        LMCAS::detail::node(root_expression));
    const auto root_free = LMCAS::free_variables(root_node);
    EXPECT_TRUE(root_free.empty(),
                "canonical RootOf identity has no dummy free variable");

    TEST_CASE("Limit and RootOf print parse round trips preserve structure");
    auto limit_node = LMCAS::detail::make_node<LimitNode>(
        SymbolicFactory::create_add({x_node, one}), "x", zero,
        LimitDirection::FromAbove);
    auto limit_expression = LMCAS::detail::make_expression_ptr(limit_node);
    auto parsed_limit = LMCAS::parse_expr(limit_expression->to_string());
    EXPECT_TRUE(parsed_limit &&
                    LMCAS::detail::node(parsed_limit.value())->equals(*limit_node),
                "LimitNode survives print/parse round trip");

    EXPECT_TRUE(root_node != nullptr,
                "checked RootOf construction returns a RootOfNode");
    auto parsed_root = LMCAS::parse_expr(root_expression->to_string());
    EXPECT_TRUE(parsed_root &&
                    LMCAS::detail::node(parsed_root.value())->equals(*root_node),
                "RootOfNode survives print/parse round trip");


    TEST_CASE("Immutable expressions can be shared for concurrent reads");
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
    for (auto& worker : workers) worker.join();
    bool deterministic = true;
    for (std::size_t i = 0; i < worker_count; ++i) {
        deterministic = deterministic && hashes[i] == expected_hash &&
                        texts[i] == expected_text && equal[i];
    }
    EXPECT_TRUE(deterministic,
                "shared immutable expressions have deterministic concurrent reads");

    return TEST_REPORT();
}
