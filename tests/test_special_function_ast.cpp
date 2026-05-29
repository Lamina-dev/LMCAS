// Feature: integration-enhancements, Property 8: Special function AST node correctness
//
// Validates: Requirements 5.7
//
// Property 8: For any recognized special function integrand pattern
// (exp(-x^2), exp(x)/x, sin(x)/x, cos(x)/x, 1/ln(x)), the integration result
// SHALL contain a FunctionNode whose `type` field matches the expected
// FuncType enum value (Erf, Ei, Si, Ci, Li respectively).
//
// Approach
// --------
//   * For each (integrand, expected FuncType) pair listed below, build the
//     integrand expression by AST construction, invoke
//     Integrator::integrate, then walk the resulting AST recursively and
//     check that *some* FunctionNode in the tree has the expected type.
//   * The five required patterns:
//       exp(-x^2)   -> FuncType::Erf
//       exp(x)/x    -> FuncType::Ei
//       sin(x)/x    -> FuncType::Si
//       cos(x)/x    -> FuncType::Ci
//       1/ln(x)     -> FuncType::Li
//   * One scaled bonus pattern is also exercised:
//       exp(-2*x^2) -> FuncType::Erf (with scaled argument sqrt(2)*x)
//   * For each case the test reports both the matched FuncType (success) or
//     the full result AST (failure) so failures are diagnosable.

#include "test_common.hpp"
#include "integration.hpp"
#include "symbolic_ast.hpp"

#include <memory>
#include <string>
#include <vector>

using lamina::Integrator;

namespace {

constexpr const char* kVarName = "x";

// Recursively walk the AST and return true if it contains a FunctionNode
// whose type equals `expected`.
bool ast_contains_functype(const std::shared_ptr<SymbolicNode>& node,
                           FunctionNode::FuncType expected) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (fn->type == expected) return true;
        for (const auto& a : fn->arguments) {
            if (ast_contains_functype(a, expected)) return true;
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (const auto& op : add->operands) {
            if (ast_contains_functype(op, expected)) return true;
        }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (const auto& op : mul->operands) {
            if (ast_contains_functype(op, expected)) return true;
        }
        return false;
    }
    if (auto pw = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (ast_contains_functype(pw->base, expected)) return true;
        if (ast_contains_functype(pw->exponent, expected)) return true;
        return false;
    }
    return false;
}

const char* functype_name(FunctionNode::FuncType t) {
    switch (t) {
        case FunctionNode::FuncType::Erf: return "Erf";
        case FunctionNode::FuncType::Ei:  return "Ei";
        case FunctionNode::FuncType::Si:  return "Si";
        case FunctionNode::FuncType::Ci:  return "Ci";
        case FunctionNode::FuncType::Li:  return "Li";
        default: return "<other>";
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

const std::vector<Case>& cases() {
    static const std::vector<Case> C = {
        {"exp(-x^2) -> Erf",   &build_exp_neg_x2, FunctionNode::FuncType::Erf},
        {"exp(-2*x^2) -> Erf", &build_exp_neg_2x2, FunctionNode::FuncType::Erf},
        {"exp(x)/x -> Ei",     &build_exp_over_x, FunctionNode::FuncType::Ei},
        {"sin(x)/x -> Si",     &build_sin_over_x, FunctionNode::FuncType::Si},
        {"cos(x)/x -> Ci",     &build_cos_over_x, FunctionNode::FuncType::Ci},
        {"1/ln(x) -> Li",      &build_inv_ln_x,   FunctionNode::FuncType::Li},
    };
    return C;
}

void verify_case(const Case& c) {
    auto integrand = c.build();
    if (!integrand) {
        EXPECT_TRUE(false, c.name + ": failed to build integrand");
        return;
    }

    Integrator integ;
    SymbolicExpr result;
    try {
        result = integ.integrate(*integrand, kVarName);
    } catch (const std::exception& e) {
        EXPECT_TRUE(false,
                    c.name + ": exception during integration: " + e.what());
        return;
    }

    bool ok = ast_contains_functype(result.root, c.expected);
    if (!ok) {
        std::string msg = c.name
            + ": expected FunctionNode with FuncType::"
            + functype_name(c.expected)
            + " but result was: " + result.to_string();
        EXPECT_TRUE(false, msg);
        return;
    }

    std::string ok_msg = c.name + " (integrand=" + integrand->to_string()
        + ", result=" + result.to_string() + ")";
    EXPECT_TRUE(true, ok_msg);
}

} // anonymous namespace

int main() {
    TEST_CASE("Property 8: Special function AST node correctness");

    for (const auto& c : cases()) {
        verify_case(c);
    }

    return TEST_REPORT();
}
