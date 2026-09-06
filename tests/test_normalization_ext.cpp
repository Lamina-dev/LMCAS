/**
 * @file test_normalization_ext.cpp
 * @brief 测试 NormalizationVisitor 对新节点类型的规范化功能。
 *
 * 覆盖: PiecewiseNode 条件验证、De Morgan 定律、双重否定、
 *       蕴含化简、虚数单位幂次化简、量词化简。
 */
#include "test_common.hpp"
#include "../include/symbolic_ast.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/expr.hpp"
#include "../include/assumption_context.hpp"
#include <cmath>

using namespace LMCAS;

static std::string to_str(const std::shared_ptr<const SymbolicNode>& node) {
    PrintVisitor pv;
    node->accept(pv);
    return pv.get_result();
}

static std::shared_ptr<const SymbolicNode> normalize(const std::shared_ptr<const SymbolicNode>& node) {
    NormalizationVisitor v;
    node->accept(v);
    return v.get_result();
}

int main() {
    TEST_CASE("Approximate square-root functions preserve their numeric domain");
    for (double argument : {1e-30, 1e200, 9.0 + 1e-11, 4.0}) {
        auto root = LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sqrt,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                LMCAS::detail::node(SymbolicExpr::number(argument))});
        auto value = LMCAS::detail::make_expression_ptr(normalize(root));
        EXPECT_TRUE(LMCAS::structurally_equal(
                        *value, *SymbolicExpr::number(std::sqrt(argument))),
                    "approximate function roots retain small, large, and near-integer values");
    }

    TEST_CASE("Power normalization preserves approximate exponents and final range");
    const struct {
        double base;
        double exponent;
    } power_cases[] = {
        {1e200, 0.5000000001},
        {1e110, -3},
        {1e-108, -3},
        {1e200, -2},
        {1e100, -3},
        {2.0, 0x1p63},
        {2.0, -0x1p63}
    };
    for (const auto& test : power_cases) {
        auto power = SymbolicExpr::power(
            SymbolicExpr::number(test.base), SymbolicExpr::number(test.exponent));
        for (auto expression : {power, SymbolicExpr::multiply(
                                   SymbolicExpr::variable("power_factor"), power)}) {
            auto simplified = LMCAS::simplify(expression);
            EXPECT_TRUE(simplified.has_value(), "numeric powers remain valid expressions");
            if (!simplified) continue;
            auto value = LMCAS::evaluate_numeric(*simplified.value(), {{"power_factor", 1.0}});
            const double expected = std::pow(test.base, test.exponent);
            if (std::isfinite(expected)) {
                EXPECT_TRUE(value && value.value().value == expected,
                            "normalization preserves the directly evaluated power");
            } else {
                EXPECT_TRUE(value && value.value().status == LMCAS::NumericStatus::PositiveInfinity,
                            "true overflow remains an explicit numeric result");
            }
        }
    }
    auto exact_reciprocal = SymbolicExpr::power(
        SymbolicExpr::number(10), SymbolicExpr::number(-3))->simplify();
    EXPECT_TRUE(LMCAS::structurally_equal(
                    *exact_reciprocal, *SymbolicExpr::number(Rational(1, 1000))),
                "exact negative powers retain exact rational results");
    auto perturbed_root = SymbolicExpr::power(
        SymbolicExpr::number(-1.0 - 1e-13), SymbolicExpr::number(Rational(1, 2)));
    auto perturbed_square = SymbolicExpr::power(perturbed_root, SymbolicExpr::number(2));
    EXPECT_TRUE(LMCAS::structurally_equal(
                    *perturbed_square->simplify(), *perturbed_square),
                "a nearby negative base is not classified as the imaginary unit");

    TEST_CASE("Nonnegative assumptions preserve nearby non-square exponents");
    LMCAS::AssumptionContext assumptions;
    auto nonnegative = assumptions.assume_sign("root_base", LMCAS::Sign::NonNegative);
    EXPECT_TRUE(nonnegative.has_value(), "nonnegative root-base assumption is accepted");
    const double near_two = 2.0 + 1e-13;
    auto near_square = SymbolicExpr::power(
        SymbolicExpr::variable("root_base"), SymbolicExpr::number(near_two));
    auto root_function = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Sqrt,
        std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(near_square)});
    NormalizationVisitor assumed_normalizer(&assumptions);
    root_function->accept(assumed_normalizer);
    auto assumed_root = LMCAS::detail::make_expression_ptr(assumed_normalizer.get_result());
    auto assumed_value = LMCAS::evaluate_numeric(*assumed_root, {{"root_base", 1e100}});
    EXPECT_TRUE(assumed_value &&
                    assumed_value.value().value == std::sqrt(std::pow(1e100, near_two)),
                "nonnegative assumptions do not round a nearby exponent to two");

    TEST_CASE("PiecewiseNode normalization");
    {
        // piecewise(2+3 if x>0, 1*y if x<0) should normalize expressions
        auto x = LMCAS::detail::make_node<VariableNode>("x");
        auto y = LMCAS::detail::make_node<VariableNode>("y");
        auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));

        auto cond1 = LMCAS::detail::make_node<RelationalNode>(x, zero, RelationalNode::Op::GT);
        auto cond2 = LMCAS::detail::make_node<RelationalNode>(x, zero, RelationalNode::Op::LT);

        // Expression: 2 + 3 (should normalize to 5)
        std::vector<std::shared_ptr<const SymbolicNode>> add_ops = {
            LMCAS::detail::make_node<NumberNode>(BigInt(2)),
            LMCAS::detail::make_node<NumberNode>(BigInt(3))
        };
        auto expr1 = LMCAS::detail::make_node<AddNode>(add_ops);

        // Expression: 1 * y (should normalize to y)
        std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {
            LMCAS::detail::make_node<NumberNode>(BigInt(1)),
            y
        };
        auto expr2 = LMCAS::detail::make_node<MultiplyNode>(mul_ops);

        std::vector<PiecewiseNode::Branch> branches;
        PiecewiseNode::Branch b1; b1.expression = expr1; b1.condition = cond1;
        PiecewiseNode::Branch b2; b2.expression = expr2; b2.condition = cond2;
        branches.push_back(b1);
        branches.push_back(b2);
        auto pw = LMCAS::detail::make_node<PiecewiseNode>(branches);
        auto result = normalize(pw);

        auto pw_result = std::dynamic_pointer_cast<const PiecewiseNode>(result);
        EXPECT_TRUE(pw_result != nullptr, "PiecewiseNode normalization returns PiecewiseNode");
        if (pw_result) {
            // First branch expression should be 5
            auto num = std::dynamic_pointer_cast<const NumberNode>(pw_result->branches()[0].expression);
            EXPECT_TRUE(num != nullptr, "First branch expression normalized to number");
            if (num) {
                EXPECT_TRUE(std::holds_alternative<BigInt>(num->value()) &&
                            std::get<BigInt>(num->value()) == BigInt(5),
                            "First branch expression = 5");
            }
            // Second branch expression should be y
            auto var = std::dynamic_pointer_cast<const VariableNode>(pw_result->branches()[1].expression);
            EXPECT_TRUE(var != nullptr && var->name() == "y", "Second branch expression = y");
        }
    }

    TEST_CASE("De Morgan's law: NOT(A AND B) = NOT(A) OR NOT(B)");
    {
        auto a = LMCAS::detail::make_node<VariableNode>("a");
        auto b = LMCAS::detail::make_node<VariableNode>("b");
        auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));

        auto cond_a = LMCAS::detail::make_node<RelationalNode>(a, zero, RelationalNode::Op::GT);
        auto cond_b = LMCAS::detail::make_node<RelationalNode>(b, zero, RelationalNode::Op::GT);

        // NOT(A AND B)
        auto and_node = LMCAS::detail::make_node<LogicalNode>(cond_a, cond_b, LogicalNode::Op::And);
        auto not_and = LMCAS::detail::make_node<LogicalNode>(and_node, nullptr, LogicalNode::Op::Not);

        auto result = normalize(not_and);
        // Should become (NOT(A) OR NOT(B))
        auto logical = std::dynamic_pointer_cast<const LogicalNode>(result);
        EXPECT_TRUE(logical != nullptr, "De Morgan AND result is LogicalNode");
        if (logical) {
            EXPECT_TRUE(logical->op() == LogicalNode::Op::Or, "De Morgan AND: outer op is Or");
            auto left_not = std::dynamic_pointer_cast<const LogicalNode>(logical->left());
            auto right_not = std::dynamic_pointer_cast<const LogicalNode>(logical->right());
            EXPECT_TRUE(left_not != nullptr && left_not->op() == LogicalNode::Op::Not,
                        "De Morgan AND: left is Not");
            EXPECT_TRUE(right_not != nullptr && right_not->op() == LogicalNode::Op::Not,
                        "De Morgan AND: right is Not");
        }
    }

    TEST_CASE("De Morgan's law: NOT(A OR B) = NOT(A) AND NOT(B)");
    {
        auto a = LMCAS::detail::make_node<VariableNode>("a");
        auto b = LMCAS::detail::make_node<VariableNode>("b");
        auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));

        auto cond_a = LMCAS::detail::make_node<RelationalNode>(a, zero, RelationalNode::Op::GT);
        auto cond_b = LMCAS::detail::make_node<RelationalNode>(b, zero, RelationalNode::Op::GT);

        // NOT(A OR B)
        auto or_node = LMCAS::detail::make_node<LogicalNode>(cond_a, cond_b, LogicalNode::Op::Or);
        auto not_or = LMCAS::detail::make_node<LogicalNode>(or_node, nullptr, LogicalNode::Op::Not);

        auto result = normalize(not_or);
        // Should become (NOT(A) AND NOT(B))
        auto logical = std::dynamic_pointer_cast<const LogicalNode>(result);
        EXPECT_TRUE(logical != nullptr, "De Morgan OR result is LogicalNode");
        if (logical) {
            EXPECT_TRUE(logical->op() == LogicalNode::Op::And, "De Morgan OR: outer op is And");
            auto left_not = std::dynamic_pointer_cast<const LogicalNode>(logical->left());
            auto right_not = std::dynamic_pointer_cast<const LogicalNode>(logical->right());
            EXPECT_TRUE(left_not != nullptr && left_not->op() == LogicalNode::Op::Not,
                        "De Morgan OR: left is Not");
            EXPECT_TRUE(right_not != nullptr && right_not->op() == LogicalNode::Op::Not,
                        "De Morgan OR: right is Not");
        }
    }

    TEST_CASE("Double negation: NOT(NOT(A)) = A");
    {
        auto a = LMCAS::detail::make_node<VariableNode>("a");
        auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        auto cond_a = LMCAS::detail::make_node<RelationalNode>(a, zero, RelationalNode::Op::GT);

        // NOT(NOT(cond_a))
        auto not_a = LMCAS::detail::make_node<LogicalNode>(cond_a, nullptr, LogicalNode::Op::Not);
        auto not_not_a = LMCAS::detail::make_node<LogicalNode>(not_a, nullptr, LogicalNode::Op::Not);

        auto result = normalize(not_not_a);
        // Should be cond_a (a > 0)
        auto rel = std::dynamic_pointer_cast<const RelationalNode>(result);
        EXPECT_TRUE(rel != nullptr, "Double negation eliminates to RelationalNode");
        if (rel) {
            EXPECT_TRUE(rel->op() == RelationalNode::Op::GT, "Double negation: op is GT");
        }
    }

    TEST_CASE("Implication: A => B = NOT(A) OR B");
    {
        auto a = LMCAS::detail::make_node<VariableNode>("a");
        auto b = LMCAS::detail::make_node<VariableNode>("b");
        auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));

        auto cond_a = LMCAS::detail::make_node<RelationalNode>(a, zero, RelationalNode::Op::GT);
        auto cond_b = LMCAS::detail::make_node<RelationalNode>(b, zero, RelationalNode::Op::GT);

        // A => B
        auto implies = LMCAS::detail::make_node<LogicalNode>(cond_a, cond_b, LogicalNode::Op::Implies);

        auto result = normalize(implies);
        // Should become (NOT(A) OR B)
        auto logical = std::dynamic_pointer_cast<const LogicalNode>(result);
        EXPECT_TRUE(logical != nullptr, "Implication result is LogicalNode");
        if (logical) {
            EXPECT_TRUE(logical->op() == LogicalNode::Op::Or, "Implication: outer op is Or");
            auto left_not = std::dynamic_pointer_cast<const LogicalNode>(logical->left());
            EXPECT_TRUE(left_not != nullptr && left_not->op() == LogicalNode::Op::Not,
                        "Implication: left is Not(A)");
        }
    }

    TEST_CASE("Imaginary unit: i^2 = -1");
    {
        // i = (-1)^(1/2), so i^2 = (-1)^(2/2) = (-1)^1 = -1
        auto neg_one = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
        auto half = LMCAS::detail::make_node<NumberNode>(Rational(1, 2));
        auto i_node = LMCAS::detail::make_node<PowerNode>(neg_one, half);

        // i^2 = ((-1)^(1/2))^2
        auto two = LMCAS::detail::make_node<NumberNode>(BigInt(2));
        auto i_squared = LMCAS::detail::make_node<PowerNode>(i_node, two);

        auto result = normalize(i_squared);
        auto num = std::dynamic_pointer_cast<const NumberNode>(result);
        EXPECT_TRUE(num != nullptr, "i^2 normalizes to NumberNode");
        if (num) {
            bool is_neg_one = false;
            if (std::holds_alternative<BigInt>(num->value()))
                is_neg_one = (std::get<BigInt>(num->value()) == BigInt(-1));
            EXPECT_TRUE(is_neg_one, "i^2 = -1");
        }
    }

    TEST_CASE("Imaginary unit: i^4 = 1");
    {
        // i^4 = ((-1)^(1/2))^4 = (-1)^(4/2) = (-1)^2 = 1
        auto neg_one = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
        auto half = LMCAS::detail::make_node<NumberNode>(Rational(1, 2));
        auto i_node = LMCAS::detail::make_node<PowerNode>(neg_one, half);

        auto four = LMCAS::detail::make_node<NumberNode>(BigInt(4));
        auto i_fourth = LMCAS::detail::make_node<PowerNode>(i_node, four);

        auto result = normalize(i_fourth);
        auto num = std::dynamic_pointer_cast<const NumberNode>(result);
        EXPECT_TRUE(num != nullptr, "i^4 normalizes to NumberNode");
        if (num) {
            bool is_one = false;
            if (std::holds_alternative<BigInt>(num->value()))
                is_one = (std::get<BigInt>(num->value()) == BigInt(1));
            EXPECT_TRUE(is_one, "i^4 = 1");
        }
    }

    TEST_CASE("Imaginary unit: i^3 = -i");
    {
        // i^3 = ((-1)^(1/2))^3 = (-1)^(3/2)
        auto neg_one = LMCAS::detail::make_node<NumberNode>(BigInt(-1));
        auto half = LMCAS::detail::make_node<NumberNode>(Rational(1, 2));
        auto i_node = LMCAS::detail::make_node<PowerNode>(neg_one, half);

        auto three = LMCAS::detail::make_node<NumberNode>(BigInt(3));
        auto i_cubed = LMCAS::detail::make_node<PowerNode>(i_node, three);

        auto result = normalize(i_cubed);
        // Should be -1 * (-1)^(1/2) = -i
        auto mul = std::dynamic_pointer_cast<const MultiplyNode>(result);
        EXPECT_TRUE(mul != nullptr, "i^3 normalizes to MultiplyNode (-i)");
        if (mul) {
            // Check it contains -1 and i
            bool has_neg_one = false;
            bool has_i = false;
            for (const auto& op : mul->operands()) {
                if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                    if (std::holds_alternative<BigInt>(n->value()) &&
                        std::get<BigInt>(n->value()) == BigInt(-1))
                        has_neg_one = true;
                }
                if (auto p = std::dynamic_pointer_cast<const PowerNode>(op)) {
                    if (auto bn = std::dynamic_pointer_cast<const NumberNode>(p->base())) {
                        if (std::holds_alternative<BigInt>(bn->value()) &&
                            std::get<BigInt>(bn->value()) == BigInt(-1)) {
                            if (auto en = std::dynamic_pointer_cast<const NumberNode>(p->exponent())) {
                                if (std::holds_alternative<Rational>(en->value()) &&
                                    std::get<Rational>(en->value()) == Rational(1, 2))
                                    has_i = true;
                            }
                        }
                    }
                }
            }
            EXPECT_TRUE(has_neg_one && has_i, "i^3 = -1 * (-1)^(1/2)");
        }
    }

    TEST_CASE("Quantifier: ForAll x in S: true -> true");
    {
        auto domain = LMCAS::detail::make_node<VariableNode>("S");
        auto pred_true = LMCAS::detail::make_node<NumberNode>(BigInt(1)); // true = 1

        auto forall = LMCAS::detail::make_node<QuantifierNode>(
            QuantifierNode::Type::ForAll, "x", domain, pred_true);

        auto result = normalize(forall);
        auto num = std::dynamic_pointer_cast<const NumberNode>(result);
        EXPECT_TRUE(num != nullptr, "ForAll true simplifies to NumberNode");
        if (num) {
            EXPECT_TRUE(num->is_one(), "ForAll x in S: true = true (1)");
        }
    }

    TEST_CASE("Quantifier: Exists x in S: false -> false");
    {
        auto domain = LMCAS::detail::make_node<VariableNode>("S");
        auto pred_false = LMCAS::detail::make_node<NumberNode>(BigInt(0)); // false = 0

        auto exists = LMCAS::detail::make_node<QuantifierNode>(
            QuantifierNode::Type::Exists, "x", domain, pred_false);

        auto result = normalize(exists);
        auto num = std::dynamic_pointer_cast<const NumberNode>(result);
        EXPECT_TRUE(num != nullptr, "Exists false simplifies to NumberNode");
        if (num) {
            EXPECT_TRUE(num->is_zero(), "Exists x in S: false = false (0)");
        }
    }

    TEST_CASE("SummationNode normalization");
    {
        // Sum(2+3, k=0..n) should normalize body to 5
        auto body = LMCAS::detail::make_node<AddNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
            LMCAS::detail::make_node<NumberNode>(BigInt(2)),
            LMCAS::detail::make_node<NumberNode>(BigInt(3))
        });
        auto lower = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        auto upper = LMCAS::detail::make_node<VariableNode>("n");

        auto sum = LMCAS::detail::make_node<SummationNode>(body, "k", lower, upper);
        auto result = normalize(sum);

        auto sum_result = std::dynamic_pointer_cast<const SummationNode>(result);
        EXPECT_TRUE(sum_result != nullptr, "SummationNode normalization returns SummationNode");
        if (sum_result) {
            auto num = std::dynamic_pointer_cast<const NumberNode>(sum_result->body());
            EXPECT_TRUE(num != nullptr, "SummationNode body normalized");
            if (num) {
                EXPECT_TRUE(std::holds_alternative<BigInt>(num->value()) &&
                            std::get<BigInt>(num->value()) == BigInt(5),
                            "SummationNode body = 5");
            }
        }
    }

    TEST_CASE("ProductNode normalization");
    {
        // Product(1*x, k=1..n) should normalize body to x
        auto x = LMCAS::detail::make_node<VariableNode>("x");
        auto body = LMCAS::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
            LMCAS::detail::make_node<NumberNode>(BigInt(1)), x
        });
        auto lower = LMCAS::detail::make_node<NumberNode>(BigInt(1));
        auto upper = LMCAS::detail::make_node<VariableNode>("n");

        auto prod = LMCAS::detail::make_node<ProductNode>(body, "k", lower, upper);
        auto result = normalize(prod);

        auto prod_result = std::dynamic_pointer_cast<const ProductNode>(result);
        EXPECT_TRUE(prod_result != nullptr, "ProductNode normalization returns ProductNode");
        if (prod_result) {
            auto var = std::dynamic_pointer_cast<const VariableNode>(prod_result->body());
            EXPECT_TRUE(var != nullptr && var->name() == "x", "ProductNode body = x");
        }
    }

    TEST_CASE("TransformNode normalization");
    {
        // Laplace{2+3}(s) should normalize body to 5
        auto body = LMCAS::detail::make_node<AddNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
            LMCAS::detail::make_node<NumberNode>(BigInt(2)),
            LMCAS::detail::make_node<NumberNode>(BigInt(3))
        });

        auto transform = LMCAS::detail::make_node<TransformNode>(
            TransformNode::TransformType::Laplace, body, "t",
            SymbolicFactory::create_variable("s"));
        auto result = normalize(transform);

        auto tr_result = std::dynamic_pointer_cast<const TransformNode>(result);
        EXPECT_TRUE(tr_result != nullptr, "TransformNode normalization returns TransformNode");
        if (tr_result) {
            auto num = std::dynamic_pointer_cast<const NumberNode>(tr_result->body());
            EXPECT_TRUE(num != nullptr, "TransformNode body normalized");
            if (num) {
                EXPECT_TRUE(std::holds_alternative<BigInt>(num->value()) &&
                            std::get<BigInt>(num->value()) == BigInt(5),
                            "TransformNode body = 5");
            }
        }
    }

    TEST_CASE("SetBuilderNode normalization");
    {
        // {x in S | 2+3 > 0} should normalize predicate condition
        auto domain = LMCAS::detail::make_node<VariableNode>("S");
        auto sum_node = LMCAS::detail::make_node<AddNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
            LMCAS::detail::make_node<NumberNode>(BigInt(2)),
            LMCAS::detail::make_node<NumberNode>(BigInt(3))
        });
        auto zero = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        auto pred = LMCAS::detail::make_node<RelationalNode>(sum_node, zero, RelationalNode::Op::GT);

        auto setb = LMCAS::detail::make_node<SetBuilderNode>("x", domain, pred);
        auto result = normalize(setb);

        auto sb_result = std::dynamic_pointer_cast<const SetBuilderNode>(result);
        EXPECT_TRUE(sb_result != nullptr, "SetBuilderNode normalization returns SetBuilderNode");
        if (sb_result) {
            // The predicate should have its left side normalized to 5
            auto rel = std::dynamic_pointer_cast<const RelationalNode>(sb_result->predicate());
            EXPECT_TRUE(rel != nullptr, "SetBuilderNode predicate is RelationalNode");
            if (rel) {
                auto num = std::dynamic_pointer_cast<const NumberNode>(rel->left());
                EXPECT_TRUE(num != nullptr, "SetBuilderNode predicate left normalized");
                if (num) {
                    EXPECT_TRUE(std::holds_alternative<BigInt>(num->value()) &&
                                std::get<BigInt>(num->value()) == BigInt(5),
                                "SetBuilderNode predicate left = 5");
                }
            }
        }
    }

    std::cout << "\n===================================================\n";
    std::cout << "Results: " << g_passes << " passed, " << g_failures << " failed\n";
    return TEST_REPORT();
}
