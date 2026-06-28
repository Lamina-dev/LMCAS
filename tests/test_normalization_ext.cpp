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

static std::string to_str(const std::shared_ptr<SymbolicNode>& node) {
    PrintVisitor pv;
    node->accept(pv);
    return pv.get_result();
}

static std::shared_ptr<SymbolicNode> normalize(const std::shared_ptr<SymbolicNode>& node) {
    NormalizationVisitor v;
    node->accept(v);
    return v.get_result();
}

int main() {
    TEST_CASE("PiecewiseNode normalization");
    {
        // piecewise(2+3 if x>0, 1*y if x<0) should normalize expressions
        auto x = std::make_shared<VariableNode>("x");
        auto y = std::make_shared<VariableNode>("y");
        auto zero = std::make_shared<NumberNode>(BigInt(0));

        auto cond1 = std::make_shared<RelationalNode>(x, zero, RelationalNode::Op::GT);
        auto cond2 = std::make_shared<RelationalNode>(x, zero, RelationalNode::Op::LT);

        // Expression: 2 + 3 (should normalize to 5)
        std::vector<std::shared_ptr<SymbolicNode>> add_ops = {
            std::make_shared<NumberNode>(BigInt(2)),
            std::make_shared<NumberNode>(BigInt(3))
        };
        auto expr1 = std::make_shared<AddNode>(add_ops);

        // Expression: 1 * y (should normalize to y)
        std::vector<std::shared_ptr<SymbolicNode>> mul_ops = {
            std::make_shared<NumberNode>(BigInt(1)),
            y
        };
        auto expr2 = std::make_shared<MultiplyNode>(mul_ops);

        std::vector<PiecewiseNode::Branch> branches;
        PiecewiseNode::Branch b1; b1.expression = expr1; b1.condition = cond1;
        PiecewiseNode::Branch b2; b2.expression = expr2; b2.condition = cond2;
        branches.push_back(b1);
        branches.push_back(b2);
        auto pw = std::make_shared<PiecewiseNode>(branches);
        auto result = normalize(pw);

        auto pw_result = std::dynamic_pointer_cast<PiecewiseNode>(result);
        EXPECT_TRUE(pw_result != nullptr, "PiecewiseNode normalization returns PiecewiseNode");
        if (pw_result) {
            // First branch expression should be 5
            auto num = std::dynamic_pointer_cast<NumberNode>(pw_result->branches[0].expression);
            EXPECT_TRUE(num != nullptr, "First branch expression normalized to number");
            if (num) {
                EXPECT_TRUE(std::holds_alternative<BigInt>(num->value) &&
                            std::get<BigInt>(num->value) == BigInt(5),
                            "First branch expression = 5");
            }
            // Second branch expression should be y
            auto var = std::dynamic_pointer_cast<VariableNode>(pw_result->branches[1].expression);
            EXPECT_TRUE(var != nullptr && var->name == "y", "Second branch expression = y");
        }
    }

    TEST_CASE("De Morgan's law: NOT(A AND B) = NOT(A) OR NOT(B)");
    {
        auto a = std::make_shared<VariableNode>("a");
        auto b = std::make_shared<VariableNode>("b");
        auto zero = std::make_shared<NumberNode>(BigInt(0));

        auto cond_a = std::make_shared<RelationalNode>(a, zero, RelationalNode::Op::GT);
        auto cond_b = std::make_shared<RelationalNode>(b, zero, RelationalNode::Op::GT);

        // NOT(A AND B)
        auto and_node = std::make_shared<LogicalNode>(cond_a, cond_b, LogicalNode::Op::And);
        auto not_and = std::make_shared<LogicalNode>(and_node, nullptr, LogicalNode::Op::Not);

        auto result = normalize(not_and);
        // Should become (NOT(A) OR NOT(B))
        auto logical = std::dynamic_pointer_cast<LogicalNode>(result);
        EXPECT_TRUE(logical != nullptr, "De Morgan AND result is LogicalNode");
        if (logical) {
            EXPECT_TRUE(logical->op == LogicalNode::Op::Or, "De Morgan AND: outer op is Or");
            auto left_not = std::dynamic_pointer_cast<LogicalNode>(logical->left);
            auto right_not = std::dynamic_pointer_cast<LogicalNode>(logical->right);
            EXPECT_TRUE(left_not != nullptr && left_not->op == LogicalNode::Op::Not,
                        "De Morgan AND: left is Not");
            EXPECT_TRUE(right_not != nullptr && right_not->op == LogicalNode::Op::Not,
                        "De Morgan AND: right is Not");
        }
    }

    TEST_CASE("De Morgan's law: NOT(A OR B) = NOT(A) AND NOT(B)");
    {
        auto a = std::make_shared<VariableNode>("a");
        auto b = std::make_shared<VariableNode>("b");
        auto zero = std::make_shared<NumberNode>(BigInt(0));

        auto cond_a = std::make_shared<RelationalNode>(a, zero, RelationalNode::Op::GT);
        auto cond_b = std::make_shared<RelationalNode>(b, zero, RelationalNode::Op::GT);

        // NOT(A OR B)
        auto or_node = std::make_shared<LogicalNode>(cond_a, cond_b, LogicalNode::Op::Or);
        auto not_or = std::make_shared<LogicalNode>(or_node, nullptr, LogicalNode::Op::Not);

        auto result = normalize(not_or);
        // Should become (NOT(A) AND NOT(B))
        auto logical = std::dynamic_pointer_cast<LogicalNode>(result);
        EXPECT_TRUE(logical != nullptr, "De Morgan OR result is LogicalNode");
        if (logical) {
            EXPECT_TRUE(logical->op == LogicalNode::Op::And, "De Morgan OR: outer op is And");
            auto left_not = std::dynamic_pointer_cast<LogicalNode>(logical->left);
            auto right_not = std::dynamic_pointer_cast<LogicalNode>(logical->right);
            EXPECT_TRUE(left_not != nullptr && left_not->op == LogicalNode::Op::Not,
                        "De Morgan OR: left is Not");
            EXPECT_TRUE(right_not != nullptr && right_not->op == LogicalNode::Op::Not,
                        "De Morgan OR: right is Not");
        }
    }

    TEST_CASE("Double negation: NOT(NOT(A)) = A");
    {
        auto a = std::make_shared<VariableNode>("a");
        auto zero = std::make_shared<NumberNode>(BigInt(0));
        auto cond_a = std::make_shared<RelationalNode>(a, zero, RelationalNode::Op::GT);

        // NOT(NOT(cond_a))
        auto not_a = std::make_shared<LogicalNode>(cond_a, nullptr, LogicalNode::Op::Not);
        auto not_not_a = std::make_shared<LogicalNode>(not_a, nullptr, LogicalNode::Op::Not);

        auto result = normalize(not_not_a);
        // Should be cond_a (a > 0)
        auto rel = std::dynamic_pointer_cast<RelationalNode>(result);
        EXPECT_TRUE(rel != nullptr, "Double negation eliminates to RelationalNode");
        if (rel) {
            EXPECT_TRUE(rel->op == RelationalNode::Op::GT, "Double negation: op is GT");
        }
    }

    TEST_CASE("Implication: A => B = NOT(A) OR B");
    {
        auto a = std::make_shared<VariableNode>("a");
        auto b = std::make_shared<VariableNode>("b");
        auto zero = std::make_shared<NumberNode>(BigInt(0));

        auto cond_a = std::make_shared<RelationalNode>(a, zero, RelationalNode::Op::GT);
        auto cond_b = std::make_shared<RelationalNode>(b, zero, RelationalNode::Op::GT);

        // A => B
        auto implies = std::make_shared<LogicalNode>(cond_a, cond_b, LogicalNode::Op::Implies);

        auto result = normalize(implies);
        // Should become (NOT(A) OR B)
        auto logical = std::dynamic_pointer_cast<LogicalNode>(result);
        EXPECT_TRUE(logical != nullptr, "Implication result is LogicalNode");
        if (logical) {
            EXPECT_TRUE(logical->op == LogicalNode::Op::Or, "Implication: outer op is Or");
            auto left_not = std::dynamic_pointer_cast<LogicalNode>(logical->left);
            EXPECT_TRUE(left_not != nullptr && left_not->op == LogicalNode::Op::Not,
                        "Implication: left is Not(A)");
        }
    }

    TEST_CASE("Imaginary unit: i^2 = -1");
    {
        // i = (-1)^(1/2), so i^2 = (-1)^(2/2) = (-1)^1 = -1
        auto neg_one = std::make_shared<NumberNode>(BigInt(-1));
        auto half = std::make_shared<NumberNode>(Rational(1, 2));
        auto i_node = std::make_shared<PowerNode>(neg_one, half);

        // i^2 = ((-1)^(1/2))^2
        auto two = std::make_shared<NumberNode>(BigInt(2));
        auto i_squared = std::make_shared<PowerNode>(i_node, two);

        auto result = normalize(i_squared);
        auto num = std::dynamic_pointer_cast<NumberNode>(result);
        EXPECT_TRUE(num != nullptr, "i^2 normalizes to NumberNode");
        if (num) {
            bool is_neg_one = false;
            if (std::holds_alternative<BigInt>(num->value))
                is_neg_one = (std::get<BigInt>(num->value) == BigInt(-1));
            EXPECT_TRUE(is_neg_one, "i^2 = -1");
        }
    }

    TEST_CASE("Imaginary unit: i^4 = 1");
    {
        // i^4 = ((-1)^(1/2))^4 = (-1)^(4/2) = (-1)^2 = 1
        auto neg_one = std::make_shared<NumberNode>(BigInt(-1));
        auto half = std::make_shared<NumberNode>(Rational(1, 2));
        auto i_node = std::make_shared<PowerNode>(neg_one, half);

        auto four = std::make_shared<NumberNode>(BigInt(4));
        auto i_fourth = std::make_shared<PowerNode>(i_node, four);

        auto result = normalize(i_fourth);
        auto num = std::dynamic_pointer_cast<NumberNode>(result);
        EXPECT_TRUE(num != nullptr, "i^4 normalizes to NumberNode");
        if (num) {
            bool is_one = false;
            if (std::holds_alternative<BigInt>(num->value))
                is_one = (std::get<BigInt>(num->value) == BigInt(1));
            EXPECT_TRUE(is_one, "i^4 = 1");
        }
    }

    TEST_CASE("Imaginary unit: i^3 = -i");
    {
        // i^3 = ((-1)^(1/2))^3 = (-1)^(3/2)
        auto neg_one = std::make_shared<NumberNode>(BigInt(-1));
        auto half = std::make_shared<NumberNode>(Rational(1, 2));
        auto i_node = std::make_shared<PowerNode>(neg_one, half);

        auto three = std::make_shared<NumberNode>(BigInt(3));
        auto i_cubed = std::make_shared<PowerNode>(i_node, three);

        auto result = normalize(i_cubed);
        // Should be -1 * (-1)^(1/2) = -i
        auto mul = std::dynamic_pointer_cast<MultiplyNode>(result);
        EXPECT_TRUE(mul != nullptr, "i^3 normalizes to MultiplyNode (-i)");
        if (mul) {
            // Check it contains -1 and i
            bool has_neg_one = false;
            bool has_i = false;
            for (const auto& op : mul->operands) {
                if (auto n = std::dynamic_pointer_cast<NumberNode>(op)) {
                    if (std::holds_alternative<BigInt>(n->value) &&
                        std::get<BigInt>(n->value) == BigInt(-1))
                        has_neg_one = true;
                }
                if (auto p = std::dynamic_pointer_cast<PowerNode>(op)) {
                    if (auto bn = std::dynamic_pointer_cast<NumberNode>(p->base)) {
                        if (std::holds_alternative<BigInt>(bn->value) &&
                            std::get<BigInt>(bn->value) == BigInt(-1)) {
                            if (auto en = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
                                if (std::holds_alternative<Rational>(en->value) &&
                                    std::get<Rational>(en->value) == Rational(1, 2))
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
        auto domain = std::make_shared<VariableNode>("S");
        auto pred_true = std::make_shared<NumberNode>(BigInt(1)); // true = 1

        auto forall = std::make_shared<QuantifierNode>(
            QuantifierNode::Type::ForAll, "x", domain, pred_true);

        auto result = normalize(forall);
        auto num = std::dynamic_pointer_cast<NumberNode>(result);
        EXPECT_TRUE(num != nullptr, "ForAll true simplifies to NumberNode");
        if (num) {
            EXPECT_TRUE(num->is_one(), "ForAll x in S: true = true (1)");
        }
    }

    TEST_CASE("Quantifier: Exists x in S: false -> false");
    {
        auto domain = std::make_shared<VariableNode>("S");
        auto pred_false = std::make_shared<NumberNode>(BigInt(0)); // false = 0

        auto exists = std::make_shared<QuantifierNode>(
            QuantifierNode::Type::Exists, "x", domain, pred_false);

        auto result = normalize(exists);
        auto num = std::dynamic_pointer_cast<NumberNode>(result);
        EXPECT_TRUE(num != nullptr, "Exists false simplifies to NumberNode");
        if (num) {
            EXPECT_TRUE(num->is_zero(), "Exists x in S: false = false (0)");
        }
    }

    TEST_CASE("SummationNode normalization");
    {
        // Sum(2+3, k=0..n) should normalize body to 5
        auto body = std::make_shared<AddNode>(std::vector<std::shared_ptr<SymbolicNode>>{
            std::make_shared<NumberNode>(BigInt(2)),
            std::make_shared<NumberNode>(BigInt(3))
        });
        auto lower = std::make_shared<NumberNode>(BigInt(0));
        auto upper = std::make_shared<VariableNode>("n");

        auto sum = std::make_shared<SummationNode>(body, "k", lower, upper);
        auto result = normalize(sum);

        auto sum_result = std::dynamic_pointer_cast<SummationNode>(result);
        EXPECT_TRUE(sum_result != nullptr, "SummationNode normalization returns SummationNode");
        if (sum_result) {
            auto num = std::dynamic_pointer_cast<NumberNode>(sum_result->body);
            EXPECT_TRUE(num != nullptr, "SummationNode body normalized");
            if (num) {
                EXPECT_TRUE(std::holds_alternative<BigInt>(num->value) &&
                            std::get<BigInt>(num->value) == BigInt(5),
                            "SummationNode body = 5");
            }
        }
    }

    TEST_CASE("ProductNode_Op normalization");
    {
        // Product(1*x, k=1..n) should normalize body to x
        auto x = std::make_shared<VariableNode>("x");
        auto body = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{
            std::make_shared<NumberNode>(BigInt(1)), x
        });
        auto lower = std::make_shared<NumberNode>(BigInt(1));
        auto upper = std::make_shared<VariableNode>("n");

        auto prod = std::make_shared<ProductNode_Op>(body, "k", lower, upper);
        auto result = normalize(prod);

        auto prod_result = std::dynamic_pointer_cast<ProductNode_Op>(result);
        EXPECT_TRUE(prod_result != nullptr, "ProductNode_Op normalization returns ProductNode_Op");
        if (prod_result) {
            auto var = std::dynamic_pointer_cast<VariableNode>(prod_result->body);
            EXPECT_TRUE(var != nullptr && var->name == "x", "ProductNode_Op body = x");
        }
    }

    TEST_CASE("TransformNode normalization");
    {
        // Laplace{2+3}(s) should normalize body to 5
        auto body = std::make_shared<AddNode>(std::vector<std::shared_ptr<SymbolicNode>>{
            std::make_shared<NumberNode>(BigInt(2)),
            std::make_shared<NumberNode>(BigInt(3))
        });

        auto transform = std::make_shared<TransformNode>(
            TransformNode::TransformType::Laplace, body, "t", "s");
        auto result = normalize(transform);

        auto tr_result = std::dynamic_pointer_cast<TransformNode>(result);
        EXPECT_TRUE(tr_result != nullptr, "TransformNode normalization returns TransformNode");
        if (tr_result) {
            auto num = std::dynamic_pointer_cast<NumberNode>(tr_result->body);
            EXPECT_TRUE(num != nullptr, "TransformNode body normalized");
            if (num) {
                EXPECT_TRUE(std::holds_alternative<BigInt>(num->value) &&
                            std::get<BigInt>(num->value) == BigInt(5),
                            "TransformNode body = 5");
            }
        }
    }

    TEST_CASE("SetBuilderNode normalization");
    {
        // {x in S | 2+3 > 0} should normalize predicate condition
        auto domain = std::make_shared<VariableNode>("S");
        auto sum_node = std::make_shared<AddNode>(std::vector<std::shared_ptr<SymbolicNode>>{
            std::make_shared<NumberNode>(BigInt(2)),
            std::make_shared<NumberNode>(BigInt(3))
        });
        auto zero = std::make_shared<NumberNode>(BigInt(0));
        auto pred = std::make_shared<RelationalNode>(sum_node, zero, RelationalNode::Op::GT);

        auto setb = std::make_shared<SetBuilderNode>("x", domain, pred);
        auto result = normalize(setb);

        auto sb_result = std::dynamic_pointer_cast<SetBuilderNode>(result);
        EXPECT_TRUE(sb_result != nullptr, "SetBuilderNode normalization returns SetBuilderNode");
        if (sb_result) {
            // The predicate should have its left side normalized to 5
            auto rel = std::dynamic_pointer_cast<RelationalNode>(sb_result->predicate);
            EXPECT_TRUE(rel != nullptr, "SetBuilderNode predicate is RelationalNode");
            if (rel) {
                auto num = std::dynamic_pointer_cast<NumberNode>(rel->left);
                EXPECT_TRUE(num != nullptr, "SetBuilderNode predicate left normalized");
                if (num) {
                    EXPECT_TRUE(std::holds_alternative<BigInt>(num->value) &&
                                std::get<BigInt>(num->value) == BigInt(5),
                                "SetBuilderNode predicate left = 5");
                }
            }
        }
    }

    std::cout << "\n===================================================\n";
    std::cout << "Results: " << g_passes << " passed, " << g_failures << " failed\n";
    return g_failures > 0 ? 1 : 0;
}
