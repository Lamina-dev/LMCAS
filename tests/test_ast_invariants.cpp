#include "test_common.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"

#include <limits>
#include <stdexcept>
#include <thread>
#include <type_traits>

namespace {

bool node_depends_on(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable) {
    return lamina::contains(
        lamina::detail::expression_from_node(node), variable);
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

class ExhaustiveProbeVisitor : public lamina::detail::SymbolicVisitor {
public:
    void visit(const NumberNode&) override {}
    void visit(const VariableNode&) override {}
    void visit(const AddNode&) override {}
    void visit(const MultiplyNode&) override {}
    void visit(const PowerNode&) override {}
    void visit(const FunctionNode&) override {}
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
};

static_assert(!std::is_abstract<ExhaustiveProbeVisitor>::value,
              "lamina::detail::SymbolicVisitor implementations must explicitly cover every node type");

template <typename T, typename = void>
struct HasPublicRoot : std::false_type {};

template <typename T>
struct HasPublicRoot<T, std::void_t<decltype(std::declval<T>().root)>>
    : std::true_type {};

#define LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(TraitName, MemberName)             \
    template <typename T, typename = void>                                  \
    struct TraitName : std::false_type {};                                  \
    template <typename T>                                                   \
    struct TraitName<T, std::void_t<decltype(std::declval<T>().MemberName)>> \
        : std::true_type {}

LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicRealField, real);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicImagField, imag);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicBranchesField, branches);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicDefaultExprField, default_expr);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicBodyField, body);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicIndexVarField, index_var);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicLowerBoundField, lower_bound);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicUpperBoundField, upper_bound);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicTransformTypeField, transform_type);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicSourceVarField, source_var);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicTargetVarField, target_var);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicQuantifierTypeField, quantifier_type);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicBoundVarField, bound_var);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicElementVarField, element_var);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicDomainField, domain);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicPredicateField, predicate);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicOperandsField, operands);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicBaseField, base);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicExponentField, exponent);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicArgumentsField, arguments);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicTypeField, type);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicRowsField, rows);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicColsField, cols);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicStorageField, storage);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicLeftField, left);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicRightField, right);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicOpField, op);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicValueField, value);
LAMINA_DEFINE_PUBLIC_MEMBER_PROBE(HasPublicNameField, name);

#undef LAMINA_DEFINE_PUBLIC_MEMBER_PROBE

using RootPointer = std::decay_t<decltype(
    lamina::detail::node(std::declval<const SymbolicExpr&>()))>;
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
                           void (SymbolicNode::*)(lamina::detail::SymbolicVisitor&) const>::value,
              "AST traversal must not expose mutable nodes");

} // namespace

int main() {
    TEST_CASE("SymbolicExpr rejects an empty AST root");
    expect_invalid([]() {
        (void)lamina::detail::expression_from_node(std::shared_ptr<const SymbolicNode>{});
    }, "SymbolicExpr rejects null root construction");

    TEST_CASE("AST nodes reject null children");
    auto one = lamina::detail::make_node<NumberNode>(BigInt(1));
    expect_invalid([&]() {
        (void)lamina::detail::make_node<AddNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{one, nullptr});
    }, "AddNode rejects null operands");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{one, nullptr});
    }, "MultiplyNode rejects null operands");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<PowerNode>(one, nullptr);
    }, "PowerNode rejects null exponent");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sin,
            std::vector<std::shared_ptr<const SymbolicNode>>{nullptr});
    }, "FunctionNode rejects null arguments");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<ComplexNode>(one, nullptr);
    }, "ComplexNode rejects null part");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<RelationalNode>(
            one, nullptr, RelationalNode::Op::EQ);
    }, "RelationalNode rejects null operand");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<LogicalNode>(
            one, nullptr, LogicalNode::Op::And);
    }, "LogicalNode rejects null binary operand");

    TEST_CASE("Extended AST nodes reject invalid children");
    auto zero = lamina::detail::make_node<NumberNode>(BigInt(0));
    auto condition = lamina::detail::make_node<RelationalNode>(one, zero, RelationalNode::Op::GT);
    expect_invalid([&]() {
        (void)lamina::detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{});
    }, "PiecewiseNode rejects empty branch list");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{{one, nullptr}});
    }, "PiecewiseNode rejects null condition");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<SummationNode>(nullptr, "k", zero, one);
    }, "SummationNode rejects null body");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<SummationNode>(one, "", zero, one);
    }, "SummationNode rejects empty index variable");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<ProductNode>(one, "k", zero, nullptr);
    }, "ProductNode rejects null bound");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<TransformNode>(
            TransformNode::TransformType::Laplace, nullptr, "t", "s");
    }, "TransformNode rejects null body");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<TransformNode>(
            TransformNode::TransformType::Laplace, one, "", "s");
    }, "TransformNode rejects empty source variable");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<QuantifierNode>(
            QuantifierNode::Type::ForAll, "x", nullptr, condition);
    }, "QuantifierNode rejects null domain");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<SetBuilderNode>("", one, condition);
    }, "SetBuilderNode rejects empty element variable");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<FiniteSetNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{nullptr});
    }, "FiniteSetNode rejects null elements");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<IntervalNode>(nullptr, one, true, true);
    }, "IntervalNode rejects null endpoints");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<MembershipNode>(one, nullptr);
    }, "MembershipNode rejects null sets");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<QuantityNode>(
            one, lamina::DimensionSignature::base("m"), Rational(0), "m");
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

    TEST_CASE("MatrixNode validates shape and storage");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{});
    }, "MatrixNode rejects empty grid");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{{one}, {one, one}});
    }, "MatrixNode rejects ragged grid");
    expect_invalid([&]() {
        (void)lamina::detail::make_node<MatrixNode>(
            std::vector<std::vector<std::shared_ptr<const SymbolicNode>>>{{one, nullptr}});
    }, "MatrixNode rejects null grid elements");
    expect_invalid([&]() {
        MatrixNode::DenseStorage dense = {one};
        (void)lamina::detail::make_node<MatrixNode>(1, 2, std::move(dense));
    }, "MatrixNode rejects dense storage size mismatch");
    expect_invalid([&]() {
        MatrixNode::SparseStorage sparse;
        sparse[2] = one;
        (void)lamina::detail::make_node<MatrixNode>(1, 2, std::move(sparse));
    }, "MatrixNode rejects sparse index out of bounds");
    expect_invalid([&]() {
        MatrixNode::DenseStorage dense;
        (void)lamina::detail::make_node<MatrixNode>(
            std::numeric_limits<size_t>::max(), 2, std::move(dense));
    }, "MatrixNode rejects dimension overflow");

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
    auto x_node = lamina::detail::make_node<VariableNode>("x");
    auto k_node = lamina::detail::make_node<VariableNode>("k");
    auto n_node = lamina::detail::make_node<VariableNode>("n");
    auto sum_body = SymbolicFactory::create_add({k_node, x_node});
    auto sum = lamina::detail::make_node<SummationNode>(sum_body, "k", zero, n_node);
    EXPECT_TRUE(!node_depends_on(sum, "k"), "Summation bound variable is not free");
    EXPECT_TRUE(node_depends_on(sum, "x"), "Summation body free variable is detected");
    EXPECT_TRUE(node_depends_on(sum, "n"), "Summation bound free variable is detected");

    auto sum_expr = lamina::detail::make_expression_ptr(sum);
    auto substituted_bound = sum_expr->substitute("k", SymbolicExpr::number(5));
    auto substituted_sum = std::dynamic_pointer_cast<const SummationNode>(lamina::detail::node(substituted_bound));
    EXPECT_TRUE(substituted_sum != nullptr, "Substitution preserves SummationNode");
    EXPECT_TRUE(node_depends_on(substituted_sum->body(), "k"),
                "Substitution does not replace bound summation variable in body");

    auto substituted_free = sum_expr->substitute("x", SymbolicExpr::number(5));
    auto substituted_free_sum = std::dynamic_pointer_cast<const SummationNode>(lamina::detail::node(substituted_free));
    EXPECT_TRUE(substituted_free_sum != nullptr, "Free substitution preserves SummationNode");
    EXPECT_TRUE(!node_depends_on(substituted_free_sum->body(), "x"),
                "Substitution replaces free variable inside summation body");

    TEST_CASE("Extended nodes are traversed by expand");
    auto y_node = lamina::detail::make_node<VariableNode>("y");
    auto x_plus_one = SymbolicFactory::create_add({x_node, one});
    auto y_plus_one = SymbolicFactory::create_add({y_node, one});
    auto product = SymbolicFactory::create_multiply({x_plus_one, y_plus_one});
    auto positive_x = lamina::detail::make_node<RelationalNode>(x_node, zero, RelationalNode::Op::GT);
    auto piecewise = lamina::detail::make_node<PiecewiseNode>(
        std::vector<PiecewiseNode::Branch>{{product, positive_x}});
    auto expanded = lamina::detail::make_expression_ptr(piecewise)->expand();
    auto expanded_piecewise = std::dynamic_pointer_cast<const PiecewiseNode>(lamina::detail::node(expanded));
    EXPECT_TRUE(expanded_piecewise != nullptr, "Expand preserves PiecewiseNode");
    EXPECT_TRUE(std::dynamic_pointer_cast<const AddNode>(expanded_piecewise->branches()[0].expression) != nullptr,
                "Expand traverses and expands PiecewiseNode branch expression");

    TEST_CASE("Substitution respects quantifier binding");
    auto domain = lamina::detail::make_node<VariableNode>("R");
    auto y_node_for_predicate = lamina::detail::make_node<VariableNode>("y");
    auto predicate = lamina::detail::make_node<RelationalNode>(x_node, y_node_for_predicate, RelationalNode::Op::GT);
    auto quantified = lamina::detail::make_node<QuantifierNode>(
        QuantifierNode::Type::ForAll, "x", domain, predicate);
    auto quantified_expr = lamina::detail::make_expression_ptr(quantified);
    auto q_bound_sub = quantified_expr->substitute("x", SymbolicExpr::number(3));
    auto q_bound = std::dynamic_pointer_cast<const QuantifierNode>(lamina::detail::node(q_bound_sub));
    EXPECT_TRUE(q_bound != nullptr, "Substitution preserves QuantifierNode");
    EXPECT_TRUE(node_depends_on(q_bound->predicate(), "x"),
                "Substitution does not replace bound quantifier variable");
    auto q_free_sub = quantified_expr->substitute("y", SymbolicExpr::number(3));
    auto q_free = std::dynamic_pointer_cast<const QuantifierNode>(lamina::detail::node(q_free_sub));
    EXPECT_TRUE(q_free != nullptr, "Free substitution preserves QuantifierNode");
    EXPECT_TRUE(!node_depends_on(q_free->predicate(), "y"),
                "Substitution replaces free variable in quantifier predicate");

    TEST_CASE("Immutable expressions can be shared for concurrent reads");
    auto shared_expr = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    const auto expected_hash = lamina::detail::node(shared_expr)->hash();
    const auto expected_text = shared_expr->to_string();
    constexpr std::size_t worker_count = 8;
    std::vector<std::size_t> hashes(worker_count);
    std::vector<std::string> texts(worker_count);
    std::vector<int> equal(worker_count, 0);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (std::size_t i = 0; i < worker_count; ++i) {
        workers.emplace_back([&, i] {
            hashes[i] = lamina::detail::node(shared_expr)->hash();
            texts[i] = shared_expr->to_string();
            equal[i] = lamina::detail::node(shared_expr)->equals(*lamina::detail::node(shared_expr)) ? 1 : 0;
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
