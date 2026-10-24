
#include "test_common.hpp"
#include "expr.hpp"
#include "integration.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/normalization_utils.hpp"
#include "residual_verification.hpp"

#include <memory>
#include <string>
#include <vector>

using namespace LMCAS;

using LMCAS::Integrator;

namespace {

constexpr const char *kVarName = "x";

// Recursively walk the AST and return true if it contains a FunctionNode
// whose type equals `expected`.
bool ast_contains_functype(const std::shared_ptr<const SymbolicNode> &node,
                           FunctionNode::FuncType expected);

static bool children_contain_functype(
    const std::vector<std::shared_ptr<const SymbolicNode>> &children,
    FunctionNode::FuncType expected) {
    for (const auto &child : children) {
        if (ast_contains_functype(child, expected)) {
            return true;
        }
    }
    return false;
}

bool ast_contains_functype(const std::shared_ptr<const SymbolicNode> &node,
                           FunctionNode::FuncType expected) {
    if (!node) {
        return false;
    }
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (fn->type() == expected) {
            return true;
        }
        return children_contain_functype(fn->arguments(), expected);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return children_contain_functype(add->operands(), expected);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return children_contain_functype(mul->operands(), expected);
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (ast_contains_functype(pw->base(), expected)) {
            return true;
        }
        if (ast_contains_functype(pw->exponent(), expected)) {
            return true;
        }
        return false;
    }
    return false;
}

const char *functype_name(FunctionNode::FuncType t) {
    switch (t) {
    case FunctionNode::FuncType::Erf:
        return "Erf";
    case FunctionNode::FuncType::Ei:
        return "Ei";
    case FunctionNode::FuncType::Si:
        return "Si";
    case FunctionNode::FuncType::Ci:
        return "Ci";
    case FunctionNode::FuncType::Li:
        return "Li";
    default:
        return "<other>";
    }
}

// Build the integrand for each named pattern using SymbolicExpr factories.
std::shared_ptr<SymbolicExpr> build_exp_neg_x2() {
    // exp(-x^2)
    auto x = SymbolicExpr::variable(kVarName);
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto neg_x2 = SymbolicExpr::multiply(SymbolicExpr::number(-1), x2);
    return SymbolicExpr::exp(neg_x2);
}

std::shared_ptr<SymbolicExpr> build_exp_neg_2x2() {
    // exp(-2*x^2)
    auto x = SymbolicExpr::variable(kVarName);
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto neg_two_x2 = SymbolicExpr::multiply(SymbolicExpr::number(-2), x2);
    return SymbolicExpr::exp(neg_two_x2);
}

std::shared_ptr<SymbolicExpr> build_exp_over_x() {
    // exp(x)/x  ==  exp(x) * x^(-1)
    auto x = SymbolicExpr::variable(kVarName);
    auto inv_x = SymbolicExpr::power(x, SymbolicExpr::number(-1));
    auto ex = SymbolicExpr::exp(x);
    return SymbolicExpr::multiply(ex, inv_x);
}

std::shared_ptr<SymbolicExpr> build_sin_over_x() {
    // sin(x)/x  ==  sin(x) * x^(-1)
    auto x = SymbolicExpr::variable(kVarName);
    auto inv_x = SymbolicExpr::power(x, SymbolicExpr::number(-1));
    auto sx = SymbolicExpr::sin(x);
    return SymbolicExpr::multiply(sx, inv_x);
}

std::shared_ptr<SymbolicExpr> build_cos_over_x() {
    // cos(x)/x  ==  cos(x) * x^(-1)
    auto x = SymbolicExpr::variable(kVarName);
    auto inv_x = SymbolicExpr::power(x, SymbolicExpr::number(-1));
    auto cx = SymbolicExpr::cos(x);
    return SymbolicExpr::multiply(cx, inv_x);
}

std::shared_ptr<SymbolicExpr> build_inv_ln_x() {
    // 1/ln(x) == ln(x)^(-1)
    auto x = SymbolicExpr::variable(kVarName);
    auto lnx = SymbolicExpr::ln(x);
    return SymbolicExpr::power(lnx, SymbolicExpr::number(-1));
}

struct Case {
    std::string name;
    std::shared_ptr<SymbolicExpr> (*build)();
    FunctionNode::FuncType expected;
};

const std::vector<Case> &cases() {
    static const std::vector<Case> C = {
        {"exp(-x^2) -> Erf", &build_exp_neg_x2, FunctionNode::FuncType::Erf},
        {"exp(-2*x^2) -> Erf", &build_exp_neg_2x2, FunctionNode::FuncType::Erf},
        {"exp(x)/x -> Ei", &build_exp_over_x, FunctionNode::FuncType::Ei},
        {"sin(x)/x -> Si", &build_sin_over_x, FunctionNode::FuncType::Si},
        {"cos(x)/x -> Ci", &build_cos_over_x, FunctionNode::FuncType::Ci},
        {"1/ln(x) -> Li", &build_inv_ln_x, FunctionNode::FuncType::Li},
    };
    return C;
}

void verify_case(const Case &c) {
    auto integrand = c.build();
    if (!integrand) {
        ADD_FAILURE() << c.name + ": failed to build integrand";
        return;
    }

    Integrator integ;
    ComputationContext context;
    auto integrated = integ.integrate_checked(*integrand, kVarName, context);
    if (!integrated) {
        ADD_FAILURE() << c.name + ": integration failed: " + integrated.error().message;
        return;
    }
    auto result = std::make_shared<SymbolicExpr>(std::move(integrated.value()));
    EXPECT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result)));
    auto derivative = differentiate(result, kVarName, context);
    ASSERT_TRUE(derivative) << c.name << ": " << derivative.error().message;
    auto residual = SymbolicExpr::add(derivative.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), integrand));
    auto proof = check_zero_residual(residual, context);
    ASSERT_TRUE(proof) << c.name << ": " << proof.error().message;
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value())) << c.name;

    bool ok = ast_contains_functype(LMCAS::detail::node(result), c.expected);
    if (!ok) {
        std::string msg = c.name + ": expected FunctionNode with FuncType::" + functype_name(c.expected) + " but result was: " + result->to_string();
        ADD_FAILURE() << msg;
        return;
    }

    std::string ok_msg = c.name + " (integrand=" + integrand->to_string() + ", result=" + result->to_string() + ")";
    EXPECT_TRUE(ok) << ok_msg;
}

} // anonymous namespace

TEST(SpecialFunctionAst, IntegrationProducesExpectedFunctionNodes) {
    for (const auto &c : cases()) {
        SCOPED_TRACE(c.name);
        verify_case(c);
    }
}

TEST(SpecialFunctionAst, NearReciprocalExponentDoesNotSnap) {
    auto x = SymbolicExpr::variable("x");
    for (double exponent : {-1.0 + 1e-13, -1.0 - 1e-13}) {
        auto power = SymbolicExpr::number(exponent);
        auto exact = detail::exact_rational_value(detail::node(power));
        ASSERT_TRUE(exact);
        EXPECT_NE(*exact, Rational(-1));
        auto integrand = SymbolicExpr::multiply(SymbolicExpr::exp(x),
            SymbolicExpr::power(x, power));
        Integrator integrator;
        ComputationContext context;
        auto result = integrator.integrate_checked(*integrand, "x", context);
        ASSERT_TRUE(result) << result.error().message;
        if (auto integral = std::dynamic_pointer_cast<const IntegralNode>(detail::node(result.value()))) {
            EXPECT_TRUE(integral->body()->equals(*detail::node(integrand)));
        } else {
            auto derivative = differentiate(detail::make_expression_ptr(result.value()), "x", context);
            ASSERT_TRUE(derivative) << derivative.error().message;
            auto proof = check_equivalent(derivative.value(), integrand, context);
            ASSERT_TRUE(proof) << proof.error().message;
            EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()));
        }
    }
}
