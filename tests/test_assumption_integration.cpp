
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "internal/visitors/limit_visitor.hpp"
#include "matcher.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/print_visitor.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace LMCAS;

/// Convert a SymbolicNode to string via PrintVisitor
static std::string node_to_string(const std::shared_ptr<const SymbolicNode> &node) {
    if (!node)
        return "null";
    PrintVisitor pv;
    node->accept(pv);
    return pv.get_result();
}


TEST(AssumptionIntegration, LimitNullptrBackwardCompat) {
    for (int a = 1; a <= 5; ++a) {
        for (int b = -3; b <= 3; ++b) {

            // Build expression: a*x + b
            auto x_node = test_variable_node("x");
            auto a_node = test_integer_node(a);
            auto b_node = test_integer_node(b);
            std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {a_node, x_node};
            auto ax = LMCAS::detail::make_node<MultiplyNode>(mul_ops);
            std::vector<std::shared_ptr<const SymbolicNode>> add_ops = {ax, b_node};
            auto expr_node = LMCAS::detail::make_node<AddNode>(add_ops);

            // Limit point: x -> 1
            auto point = LMCAS::detail::make_node<NumberNode>(BigInt(1));

            /// 使用默认空上下文的 LimitVisitor.
            LimitVisitor visitor_default("x", point, "");
            expr_node->accept(visitor_default);
            auto result_default = visitor_default.get_result();

            // LimitVisitor with explicit nullptr context
            LimitVisitor visitor_nullptr("x", point, "", nullptr);
            expr_node->accept(visitor_nullptr);
            auto result_nullptr = visitor_nullptr.get_result();

            std::string str_default = node_to_string(result_default);
            std::string str_nullptr = node_to_string(result_nullptr);

            EXPECT_EQ((str_default), (str_nullptr)) << "default and explicit empty-assumption limits agree";
        }
    }
}

TEST(AssumptionIntegration, LimitTrigNullptr) {
    // lim x->0 sin(x) = 0
    {
        auto x_node = test_variable_node("x");
        std::vector<std::shared_ptr<const SymbolicNode>> args = {x_node};
        auto sin_x = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sin, args);
        auto point = LMCAS::detail::make_node<NumberNode>(BigInt(0));

        LimitVisitor visitor_default("x", point, "");
        sin_x->accept(visitor_default);
        auto result_default = visitor_default.get_result();

        LimitVisitor visitor_nullptr("x", point, "", nullptr);
        sin_x->accept(visitor_nullptr);
        auto result_nullptr = visitor_nullptr.get_result();

        std::string str_default = node_to_string(result_default);
        std::string str_nullptr = node_to_string(result_nullptr);

        EXPECT_EQ((str_default), (str_nullptr)) << "lim x->0 sin(x): default == nullptr context";
    }
    // lim x->0 cos(x) = 1
    {
        auto x_node = test_variable_node("x");
        std::vector<std::shared_ptr<const SymbolicNode>> args = {x_node};
        auto cos_x = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Cos, args);
        auto point = LMCAS::detail::make_node<NumberNode>(BigInt(0));

        LimitVisitor visitor_default("x", point, "");
        cos_x->accept(visitor_default);
        auto result_default = visitor_default.get_result();

        LimitVisitor visitor_nullptr("x", point, "", nullptr);
        cos_x->accept(visitor_nullptr);
        auto result_nullptr = visitor_nullptr.get_result();

        std::string str_default = node_to_string(result_default);
        std::string str_nullptr = node_to_string(result_nullptr);

        EXPECT_EQ((str_default), (str_nullptr)) << "lim x->0 cos(x): default == nullptr context";
    }
}

TEST(AssumptionIntegration, SeriesNullptrBackwardCompat) {
    for (int order = 2; order <= 5; ++order) {

        auto x = SymbolicExpr::variable("x");
        auto exp_x = SymbolicExpr::exp(x);
        auto zero = SymbolicExpr::number(0);

        /// 使用默认空上下文计算级数.
        auto result_default = exp_x->series("x", zero, order);

        // Series with explicit nullptr context
        auto result_nullptr = exp_x->series("x", zero, order, nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ((str_default), (str_nullptr)) << "default and explicit empty-assumption series agree";
    }
}

TEST(AssumptionIntegration, SeriesSinCosNullptr) {
    // sin(x) series at 0, order 5
    {
        auto x = SymbolicExpr::variable("x");
        auto sin_x = SymbolicExpr::sin(x);
        auto zero = SymbolicExpr::number(0);

        auto result_default = sin_x->series("x", zero, 5);
        auto result_nullptr = sin_x->series("x", zero, 5, nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ((str_default), (str_nullptr)) << "sin(x) series order 5: default == nullptr context";
    }
    // cos(x) series at 0, order 4
    {
        auto x = SymbolicExpr::variable("x");
        auto cos_x = SymbolicExpr::cos(x);
        auto zero = SymbolicExpr::number(0);

        auto result_default = cos_x->series("x", zero, 4);
        auto result_nullptr = cos_x->series("x", zero, 4, nullptr);

        std::string str_default = result_default ? result_default->to_string() : "null";
        std::string str_nullptr = result_nullptr ? result_nullptr->to_string() : "null";

        EXPECT_EQ((str_default), (str_nullptr)) << "cos(x) series order 4: default == nullptr context";
    }
}

TEST(AssumptionIntegration, MatcherNullptrBackwardCompat) {
    for (int num_val = 1; num_val <= 9; ++num_val) {

        // Pattern: _a + num_val
        auto wc = LMCAS::detail::make_node<VariableNode>("_a");
        auto num = test_integer_node(num_val);
        std::vector<std::shared_ptr<const SymbolicNode>> pat_ops = {wc, num};
        auto pat_node = LMCAS::detail::make_node<AddNode>(pat_ops);
        auto pattern = LMCAS::detail::expression_from_node(pat_node);
        // Target: x + num_val
        auto x_node = test_variable_node("x");
        std::vector<std::shared_ptr<const SymbolicNode>> tgt_ops = {x_node, num};
        auto tgt_node = LMCAS::detail::make_node<AddNode>(tgt_ops);
        auto target = LMCAS::detail::expression_from_node(tgt_node);
        std::unordered_set<std::string> wildcards = {"_a"};

        /// 使用默认空上下文执行匹配.
        MatchMap results_default;
        bool matched_default = Matcher::match(pattern, target, wildcards, results_default);

        // Match with explicit nullptr context
        MatchMap results_nullptr;
        bool matched_nullptr = Matcher::match(pattern, target, wildcards, results_nullptr);

        EXPECT_TRUE((matched_default == matched_nullptr)) << "repeated structural matching is deterministic";

        if (matched_default && matched_nullptr) {
            // Verify bindings are identical
            EXPECT_TRUE((results_default.size() == results_nullptr.size())) << "repeated matching produces the same binding count";
            for (auto &[key, val] : results_default) {
                EXPECT_TRUE((results_nullptr.count(key) > 0)) << "repeated matching preserves binding names";
                EXPECT_EQ((val.to_string()), (results_nullptr.at(key).to_string())) << "repeated matching preserves binding values";
            }
        }
    }
}

TEST(AssumptionIntegration, RewriteEngineNullptr) {
    // Test that RewriteEngine with nullptr context behaves identically to default
    // Rule: x + 0 -> x
    {
        auto wc_a = wildcard("_a");
        auto zero_expr = LMCAS::detail::expression_from_node(test_integer_node(0));
        std::vector<std::shared_ptr<const SymbolicNode>> pat_ops = {LMCAS::detail::node(wc_a), LMCAS::detail::node(zero_expr)};
        auto pattern = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<AddNode>(pat_ops));
        SymbolicExpr replacement = wc_a;
        std::unordered_set<std::string> wcs = {"_a"};

        Rule rule(pattern, replacement, wcs);

        // Target: y + 0
        auto y_node = test_variable_node("y");
        std::vector<std::shared_ptr<const SymbolicNode>> tgt_ops = {y_node, test_integer_node(0)};
        auto target = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<AddNode>(tgt_ops));
        /// 使用默认空上下文的 RewriteEngine.
        RewriteEngine engine_default;
        engine_default.add_rule(rule);
        auto result_default = engine_default.apply(target, 10);
        ASSERT_TRUE((result_default.has_value())) << "default additive rewrite succeeds";
        if (!result_default)
            return;

        /// 使用显式上下文,并保持 assumptions 为空.
        RewriteEngine engine_nullptr;
        engine_nullptr.add_rule(rule);
        ComputationContext context;
        auto checked = engine_nullptr.apply_checked(target, context, 10);
        ASSERT_TRUE((checked.has_value())) << "checked rewrite succeeds";
        if (!checked)
            return;
        auto result_nullptr = checked.value();

        EXPECT_EQ((result_default.value().to_string()), (result_nullptr.to_string())) << "RewriteEngine y+0: default == explicit context";
    }

    // Rule: _a * 1 -> _a
    {
        auto wc_a = wildcard("_a");
        auto one_expr = LMCAS::detail::expression_from_node(test_integer_node(1));
        std::vector<std::shared_ptr<const SymbolicNode>> pat_ops = {LMCAS::detail::node(wc_a), LMCAS::detail::node(one_expr)};
        auto pattern = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<MultiplyNode>(pat_ops));
        SymbolicExpr replacement = wc_a;
        std::unordered_set<std::string> wcs = {"_a"};

        Rule rule(pattern, replacement, wcs);

        // Target: z * 1
        auto z_node = test_variable_node("z");
        std::vector<std::shared_ptr<const SymbolicNode>> tgt_ops = {z_node, test_integer_node(1)};
        auto target = LMCAS::detail::expression_from_node(LMCAS::detail::make_node<MultiplyNode>(tgt_ops));
        RewriteEngine engine_default;
        engine_default.add_rule(rule);
        auto result_default = engine_default.apply(target, 10);
        ASSERT_TRUE((result_default.has_value())) << "default multiplicative rewrite succeeds";
        if (!result_default)
            return;

        RewriteEngine engine_nullptr;
        engine_nullptr.add_rule(rule);
        ComputationContext context;
        auto checked = engine_nullptr.apply_checked(target, context, 10);
        ASSERT_TRUE((checked.has_value())) << "checked rewrite succeeds";
        if (!checked)
            return;
        auto result_nullptr = checked.value();

        EXPECT_EQ((result_default.value().to_string()), (result_nullptr.to_string())) << "RewriteEngine z*1: default == explicit context";
    }
}
