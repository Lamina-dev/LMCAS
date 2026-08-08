/**
 * @file test_piecewise_limit.cpp
 * @brief 测试 LimitVisitor 对 PiecewiseNode 的极限处理和方向感知符号求值。
 *
 * 覆盖:
 * - 单侧极限选择正确分支
 * - 双侧极限相等时返回值
 * - 双侧极限不同时返回 nullptr（极限不存在）
 * - 方向感知的 sgn(x) 极限
 * - 方向感知的 |x| 极限
 */
#include "test_common.hpp"
#include "../include/visitors/limit_visitor.hpp"
#include "../include/symbolic_ast.hpp"

int main() {
    auto x = SymbolicExpr::variable("x");
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto neg_one = SymbolicExpr::number(-1);
    auto two = SymbolicExpr::number(2);
    auto three = SymbolicExpr::number(3);


    TEST_CASE("Piecewise: right limit selects x>0 branch");
    {
        // piecewise(x+1 if x>0, x-1 if x<0)
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));

        // Condition: x > 0
        auto cond_gt = lamina::detail::make_node<RelationalNode>(x_node, zero_node, RelationalNode::Op::GT);
        // Condition: x < 0
        auto cond_lt = lamina::detail::make_node<RelationalNode>(x_node, zero_node, RelationalNode::Op::LT);

        // Branch expressions: x+1 and x-1
        std::vector<std::shared_ptr<const SymbolicNode>> add1_ops = {x_node, lamina::detail::make_node<NumberNode>(BigInt(1))};
        auto expr_plus = lamina::detail::make_node<AddNode>(add1_ops);
        std::vector<std::shared_ptr<const SymbolicNode>> add2_ops = {x_node, lamina::detail::make_node<NumberNode>(BigInt(-1))};
        auto expr_minus = lamina::detail::make_node<AddNode>(add2_ops);

        std::vector<PiecewiseNode::Branch> branches;
        branches.push_back({expr_plus, cond_gt});
        branches.push_back({expr_minus, cond_lt});
        auto pw = lamina::detail::make_node<PiecewiseNode>(branches);

        // Right limit at x=0: should select x>0 branch, giving 0+1 = 1
        LimitVisitor rv("x", zero_node, "+");
        pw->accept(rv);
        auto r = rv.get_result();
        EXPECT_TRUE(r != nullptr, "Right limit of piecewise at 0 is not null");
        if (r) {
            NormalizationVisitor norm; r->accept(norm); r = norm.get_result();
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "1", "Right limit = 1");
        }
    }

    TEST_CASE("Piecewise: left limit selects x<0 branch");
    {
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto cond_gt = lamina::detail::make_node<RelationalNode>(x_node, zero_node, RelationalNode::Op::GT);
        auto cond_lt = lamina::detail::make_node<RelationalNode>(x_node, zero_node, RelationalNode::Op::LT);
        std::vector<std::shared_ptr<const SymbolicNode>> add1_ops = {x_node, lamina::detail::make_node<NumberNode>(BigInt(1))};
        auto expr_plus = lamina::detail::make_node<AddNode>(add1_ops);
        std::vector<std::shared_ptr<const SymbolicNode>> add2_ops = {x_node, lamina::detail::make_node<NumberNode>(BigInt(-1))};
        auto expr_minus = lamina::detail::make_node<AddNode>(add2_ops);
        std::vector<PiecewiseNode::Branch> branches;
        branches.push_back({expr_plus, cond_gt});
        branches.push_back({expr_minus, cond_lt});
        auto pw = lamina::detail::make_node<PiecewiseNode>(branches);

        // Left limit at x=0: should select x<0 branch, giving 0-1 = -1
        LimitVisitor lv("x", zero_node, "-");
        pw->accept(lv);
        auto r = lv.get_result();
        EXPECT_TRUE(r != nullptr, "Left limit of piecewise at 0 is not null");
        if (r) {
            NormalizationVisitor norm; r->accept(norm); r = norm.get_result();
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "-1", "Left limit = -1");
        }
    }

    TEST_CASE("Piecewise: two-sided limit DNE when left != right");
    {
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto cond_gt = lamina::detail::make_node<RelationalNode>(x_node, zero_node, RelationalNode::Op::GT);
        auto cond_lt = lamina::detail::make_node<RelationalNode>(x_node, zero_node, RelationalNode::Op::LT);
        std::vector<std::shared_ptr<const SymbolicNode>> add1_ops = {x_node, lamina::detail::make_node<NumberNode>(BigInt(1))};
        auto expr_plus = lamina::detail::make_node<AddNode>(add1_ops);
        std::vector<std::shared_ptr<const SymbolicNode>> add2_ops = {x_node, lamina::detail::make_node<NumberNode>(BigInt(-1))};
        auto expr_minus = lamina::detail::make_node<AddNode>(add2_ops);
        std::vector<PiecewiseNode::Branch> branches;
        branches.push_back({expr_plus, cond_gt});
        branches.push_back({expr_minus, cond_lt});
        auto pw = lamina::detail::make_node<PiecewiseNode>(branches);

        // Two-sided limit at x=0: left=−1, right=1, so DNE (nullptr)
        LimitVisitor tv("x", zero_node, "");
        pw->accept(tv);
        auto r = tv.get_result();
        EXPECT_TRUE(r == nullptr, "Two-sided limit DNE when left != right");
    }

    TEST_CASE("Piecewise: two-sided limit exists when left == right");
    {
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto cond_gt = lamina::detail::make_node<RelationalNode>(x_node, zero_node, RelationalNode::Op::GT);
        auto cond_lt = lamina::detail::make_node<RelationalNode>(x_node, zero_node, RelationalNode::Op::LT);
        // Both branches give x^2, so limit from both sides at 0 is 0
        auto x_sq = lamina::detail::make_node<PowerNode>(x_node, lamina::detail::make_node<NumberNode>(BigInt(2)));
        std::vector<PiecewiseNode::Branch> branches;
        branches.push_back({x_sq, cond_gt});
        branches.push_back({x_sq->clone(), cond_lt});
        auto pw = lamina::detail::make_node<PiecewiseNode>(branches);

        LimitVisitor tv("x", zero_node, "");
        pw->accept(tv);
        auto r = tv.get_result();
        EXPECT_TRUE(r != nullptr, "Two-sided limit exists when branches agree");
        if (r) {
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "0", "Two-sided limit = 0");
        }
    }


    TEST_CASE("sgn(x) right limit at 0 = 1");
    {
        // sgn(x) as x->0+
        auto lim = x->limit("x", zero, "+");
        // We need to construct sgn(x) directly
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        std::vector<std::shared_ptr<const SymbolicNode>> sgn_args = {x_node};
        auto sgn_x = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sgn, sgn_args);

        LimitVisitor rv("x", zero_node, "+");
        sgn_x->accept(rv);
        auto r = rv.get_result();
        EXPECT_TRUE(r != nullptr, "sgn(x) right limit at 0 is not null");
        if (r) {
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "1", "sgn(x) x->0+ = 1");
        }
    }

    TEST_CASE("sgn(x) left limit at 0 = -1");
    {
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        std::vector<std::shared_ptr<const SymbolicNode>> sgn_args = {x_node};
        auto sgn_x = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sgn, sgn_args);

        LimitVisitor lv("x", zero_node, "-");
        sgn_x->accept(lv);
        auto r = lv.get_result();
        EXPECT_TRUE(r != nullptr, "sgn(x) left limit at 0 is not null");
        if (r) {
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "-1", "sgn(x) x->0- = -1");
        }
    }


    TEST_CASE("|x| limit at 0 = 0 (from either side)");
    {
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        std::vector<std::shared_ptr<const SymbolicNode>> abs_args = {x_node};
        auto abs_x = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Abs, abs_args);

        LimitVisitor rv("x", zero_node, "+");
        abs_x->accept(rv);
        auto r = rv.get_result();
        EXPECT_TRUE(r != nullptr, "|x| right limit at 0 is not null");
        if (r) {
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "0", "|x| x->0+ = 0");
        }

        LimitVisitor lv("x", zero_node, "-");
        abs_x->accept(lv);
        r = lv.get_result();
        EXPECT_TRUE(r != nullptr, "|x| left limit at 0 is not null");
        if (r) {
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "0", "|x| x->0- = 0");
        }
    }


    TEST_CASE("1/x right limit at 0 = +inf");
    {
        // x^(-1) as x→0+
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto inv_x = lamina::detail::make_node<PowerNode>(x_node, lamina::detail::make_node<NumberNode>(BigInt(-1)));

        LimitVisitor rv("x", zero_node, "+");
        inv_x->accept(rv);
        auto r = rv.get_result();
        EXPECT_TRUE(r != nullptr, "1/x right limit at 0 is not null");
        if (r) {
            auto func = std::dynamic_pointer_cast<const FunctionNode>(r);
            EXPECT_TRUE(func != nullptr && func->type() == FunctionNode::FuncType::Infinity,
                        "1/x x->0+ = +inf");
        }
    }

    TEST_CASE("1/x left limit at 0 = -inf");
    {
        // x^(-1) as x→0-
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto inv_x = lamina::detail::make_node<PowerNode>(x_node, lamina::detail::make_node<NumberNode>(BigInt(-1)));

        LimitVisitor lv("x", zero_node, "-");
        inv_x->accept(lv);
        auto r = lv.get_result();
        EXPECT_TRUE(r != nullptr, "1/x left limit at 0 is not null");
        if (r) {
            // Should be -1*inf (negative infinity)
            auto mul = std::dynamic_pointer_cast<const MultiplyNode>(r);
            bool is_neg_inf = false;
            if (mul) {
                for (auto& op : mul->operands()) {
                    if (auto f = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                        if (f->type() == FunctionNode::FuncType::Infinity) is_neg_inf = true;
                    }
                }
            }
            EXPECT_TRUE(is_neg_inf, "1/x x->0- = -inf");
        }
    }

    TEST_CASE("1/x^2 right limit at 0 = +inf (even exponent)");
    {
        // x^(-2) as x→0+
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto inv_x2 = lamina::detail::make_node<PowerNode>(x_node, lamina::detail::make_node<NumberNode>(BigInt(-2)));

        LimitVisitor rv("x", zero_node, "+");
        inv_x2->accept(rv);
        auto r = rv.get_result();
        EXPECT_TRUE(r != nullptr, "1/x^2 right limit at 0 is not null");
        if (r) {
            auto func = std::dynamic_pointer_cast<const FunctionNode>(r);
            EXPECT_TRUE(func != nullptr && func->type() == FunctionNode::FuncType::Infinity,
                        "1/x^2 x->0+ = +inf");
        }
    }

    TEST_CASE("1/x^2 left limit at 0 = +inf (even exponent, always positive)");
    {
        // x^(-2) as x→0-
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto zero_node = lamina::detail::make_node<NumberNode>(BigInt(0));
        auto inv_x2 = lamina::detail::make_node<PowerNode>(x_node, lamina::detail::make_node<NumberNode>(BigInt(-2)));

        LimitVisitor lv("x", zero_node, "-");
        inv_x2->accept(lv);
        auto r = lv.get_result();
        EXPECT_TRUE(r != nullptr, "1/x^2 left limit at 0 is not null");
        if (r) {
            auto func = std::dynamic_pointer_cast<const FunctionNode>(r);
            EXPECT_TRUE(func != nullptr && func->type() == FunctionNode::FuncType::Infinity,
                        "1/x^2 x->0- = +inf (even exponent)");
        }
    }


    TEST_CASE("Piecewise with GEQ/LEQ: right limit at boundary");
    {
        // piecewise(x^2 if x>=1, 2*x-1 if x<1)
        // At x=1: right limit uses x>=1 branch → 1^2 = 1
        auto x_node = lamina::detail::make_node<VariableNode>("x");
        auto one_node = lamina::detail::make_node<NumberNode>(BigInt(1));

        auto cond_geq = lamina::detail::make_node<RelationalNode>(x_node, one_node, RelationalNode::Op::GEQ);
        auto cond_lt = lamina::detail::make_node<RelationalNode>(x_node, one_node, RelationalNode::Op::LT);

        auto x_sq = lamina::detail::make_node<PowerNode>(x_node, lamina::detail::make_node<NumberNode>(BigInt(2)));
        std::vector<std::shared_ptr<const SymbolicNode>> lin_ops = {
            lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
                lamina::detail::make_node<NumberNode>(BigInt(2)), x_node}),
            lamina::detail::make_node<NumberNode>(BigInt(-1))};
        auto lin_expr = lamina::detail::make_node<AddNode>(lin_ops);

        std::vector<PiecewiseNode::Branch> branches;
        branches.push_back({x_sq, cond_geq});
        branches.push_back({lin_expr, cond_lt});
        auto pw = lamina::detail::make_node<PiecewiseNode>(branches);

        // Right limit at x=1: x>=1 branch → 1
        LimitVisitor rv("x", one_node, "+");
        pw->accept(rv);
        auto r = rv.get_result();
        EXPECT_TRUE(r != nullptr, "Piecewise GEQ right limit at 1 is not null");
        if (r) {
            NormalizationVisitor norm; r->accept(norm); r = norm.get_result();
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "1",
                          "Right limit at x=1 (x>=1 branch) = 1");
        }

        // Left limit at x=1: x<1 branch → 2(1)-1 = 1
        LimitVisitor lv("x", one_node, "-");
        pw->accept(lv);
        r = lv.get_result();
        EXPECT_TRUE(r != nullptr, "Piecewise LT left limit at 1 is not null");
        if (r) {
            NormalizationVisitor norm; r->accept(norm); r = norm.get_result();
            EXPECT_EQ_STR(lamina::detail::make_expression_ptr(r)->to_string(), "1",
                          "Left limit at x=1 (x<1 branch) = 1");
        }

        // Two-sided limit exists (both = 1)
        LimitVisitor tv("x", one_node, "");
        pw->accept(tv);
        r = tv.get_result();
        EXPECT_TRUE(r != nullptr, "Two-sided limit exists when both sides = 1");
    }

    return TEST_REPORT();
}
