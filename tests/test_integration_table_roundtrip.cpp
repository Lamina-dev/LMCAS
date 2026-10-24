#include "expr.hpp"
#include "test_common.hpp"
#include "integration.hpp"
#include "matcher.hpp"
#include "assumption_context.hpp"
#include "residual_verification.hpp"
#include "internal/integration_support.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>
#include <utility>

using namespace LMCAS;

using LMCAS::detail::IntegrationEntry;
using LMCAS::Matcher;
using LMCAS::MatchMap;

namespace {

constexpr const char *kVarName = "x";
constexpr double kTolerance = 1e-10;

std::shared_ptr<SymbolicExpr> concrete_for_wildcard(const std::string &wname) {
    if (wname == "_u") {
        return SymbolicExpr::variable(kVarName);
    }
    if (wname == "_a") {
        return SymbolicExpr::number(2);
    }
    if (wname == "_n") {
        return SymbolicExpr::number(3);
    }
    return SymbolicExpr::number(1);
}

MatchMap make_bindings(const IntegrationEntry &entry) {
    MatchMap m;
    for (const auto &w : entry.wildcards) {
        auto val = concrete_for_wildcard(w);
        m.emplace(w, *val);
    }
    return m;
}

const std::vector<double> &sample_points() {
    static const std::vector<double> pts = {
        0.25, 0.5, 0.75, 1.5, 2.5, 3.0, 3.5, 4.0,
        -0.25, -0.5, -0.75, -1.5, -2.5, -3.0, -3.5, -4.0};
    return pts;
}

struct EntryReport {
    bool checked = false;
    bool failed = false;
    std::string failure_detail;
};

bool evaluate_roundtrip_point(
    const std::shared_ptr<SymbolicExpr> &pattern,
    const std::shared_ptr<SymbolicExpr> &derivative,
    double point, double &pattern_value, double &derivative_value) {
    auto x_value = SymbolicExpr::number(point);
    auto pattern_at = LMCAS::substitute(pattern, kVarName, x_value).value();
    auto derivative_at = LMCAS::substitute(derivative, kVarName, x_value).value();
    if (!pattern_at || !derivative_at) {
        return false;
    }
    pattern_at = LMCAS::simplify(pattern_at).value();
    derivative_at = LMCAS::simplify(derivative_at).value();
    auto pv = test_numeric_eval(pattern_at);
    auto dv = test_numeric_eval(derivative_at);
    if (!pv || !dv || !std::isfinite(*pv) || !std::isfinite(*dv)) {
        return false;
    }
    pattern_value = *pv;
    derivative_value = *dv;
    return true;
}

class AbsoluteArguments final : public detail::RecursiveSymbolicVisitor {
public:
    std::vector<detail::SymbolicNodePtr> arguments;
    void visit(const FunctionNode& node) override {
        if (node.type() == FunctionNode::FuncType::Abs && node.arguments().size() == 1) {
            arguments.push_back(node.arguments()[0]);
        }
        RecursiveSymbolicVisitor::visit(node);
    }
};

class SignedRealValues final : public detail::SymbolicRewriter {
    const std::vector<std::pair<detail::SymbolicNodePtr, bool>>& signs_;

    detail::SymbolicNodePtr signed_value(const detail::SymbolicNodePtr& node) {
        for (const auto& [argument, positive] : signs_) {
            if (argument->equals(*node)) {
                auto value = rewrite(node);
                return positive ? value : SymbolicFactory::create_multiply({
                    detail::node(SymbolicExpr::number(-1)), value});
            }
        }
        return nullptr;
    }

public:
    explicit SignedRealValues(const std::vector<std::pair<detail::SymbolicNodePtr, bool>>& signs)
        : signs_(signs) {}
    void visit(const FunctionNode& node) override {
        if (node.type() == FunctionNode::FuncType::Abs && node.arguments().size() == 1) {
            if (auto value = signed_value(node.arguments()[0])) {
                set_result(std::move(value));
                return;
            }
        }
        SymbolicRewriter::visit(node);
    }
    void visit(const PowerNode& node) override {
        const auto exponent = detail::exact_rational_value(node.exponent());
        const auto* square = dynamic_cast<const PowerNode*>(node.base().get());
        if (exponent && *exponent == Rational(1, 2) && square) {
            const auto degree = detail::exact_rational_value(square->exponent());
            if (degree && *degree == Rational(2)) {
                if (auto value = signed_value(square->base())) {
                    set_result(std::move(value));
                    return;
                }
            }
        }
        SymbolicRewriter::visit(node);
    }
};

ResidualCheckResult verify_on_real_branch(const ExprPtr& primitive, const ExprPtr& pattern) {
    AbsoluteArguments absolute;
    detail::node(primitive)->accept(absolute);
    for (double point : sample_points()) {
        ComputationContext numeric;
        auto value = evaluate_numeric(*primitive, {{kVarName, point}}, numeric);
        if (!value) {
            if (value.error().code == CasErrc::DomainError) { continue; }
            return ResidualCheckResult::failure(value.error());
        }
        if (!value.value().is_finite()) { continue; }
        auto assumptions = std::make_shared<AssumptionContext>();
        auto real = assumptions->assume_domain_checked(kVarName, Domain::Real);
        if (!real) { return ResidualCheckResult::failure(real.error()); }
        std::vector<std::pair<detail::SymbolicNodePtr, bool>> signs;
        bool valid = true;
        for (const auto& argument : absolute.arguments) {
            auto at = evaluate_numeric(detail::expression_from_node(argument),
                {{kVarName, point}}, numeric);
            if (!at) { return ResidualCheckResult::failure(at.error()); }
            if (!at.value().is_finite() || at.value().value == 0.0) {
                valid = false;
                break;
            }
            auto condition = detail::expression_from_node(detail::make_node<RelationalNode>(
                argument, detail::node(SymbolicExpr::number(0)),
                at.value().value > 0.0 ? RelationOp::GT : RelationOp::LT));
            auto assumed = assumptions->assume_checked(condition);
            if (!assumed) { return ResidualCheckResult::failure(assumed.error()); }
            signs.emplace_back(argument, at.value().value > 0.0);
        }
        if (!valid) { continue; }
        ComputationContext branch;
        auto attached = branch.set_assumptions(assumptions);
        if (!attached) { return ResidualCheckResult::failure(attached.error()); }
        SignedRealValues signed_real(signs);
        auto simplified = LMCAS::simplify(detail::make_expression_ptr(
            signed_real.rewrite(detail::node(primitive))), branch);
        if (!simplified) { return ResidualCheckResult::failure(simplified.error()); }
        auto derivative = LMCAS::differentiate(simplified.value(), kVarName, branch);
        if (!derivative) { return ResidualCheckResult::failure(derivative.error()); }
        auto branch_pattern = detail::make_expression_ptr(
            signed_real.rewrite(detail::node(pattern)));
        return check_equivalent(derivative.value(), branch_pattern, branch);
    }
    return ResidualCheckResult::failure(CasErrc::DomainError,
        "no real branch samples for table primitive", "test.integration_table");
}

EntryReport verify_entry(const IntegrationEntry &entry) {
    EntryReport rep;

    MatchMap bindings = make_bindings(entry);

    SymbolicExpr pat_inst = Matcher::replace(entry.pattern, bindings, false);
    auto pattern = LMCAS::detail::make_expression_ptr(pat_inst);
    Integrator integrator;
    ComputationContext context;
    auto assumptions = std::make_shared<AssumptionContext>();
    auto real_variable = assumptions->assume_domain_checked(kVarName, Domain::Real);
    if (!real_variable) {
        rep.failed = true;
        rep.failure_detail = real_variable.error().message;
        return rep;
    }
    auto attached = context.set_assumptions(assumptions);
    if (!attached) {
        rep.failed = true;
        rep.failure_detail = attached.error().message;
        return rep;
    }
    if (entry.condition && !entry.condition(bindings, kVarName, context)) {
        rep.failed = true;
        rep.failure_detail = "concrete rule parameters do not satisfy its conditions";
        return rep;
    }
    auto entry_result = detail::make_expression_ptr(
        Matcher::replace(entry.result, bindings, false));
    auto entry_derivative = LMCAS::differentiate(entry_result, kVarName, context);
    if (!entry_derivative) {
        rep.failed = true;
        rep.failure_detail = entry_derivative.error().message;
        return rep;
    }
    auto entry_proof = verify_on_real_branch(entry_result, pattern);
    if (!entry_proof || !std::holds_alternative<ProvedZeroResidual>(entry_proof.value())) {
        rep.failed = true;
        rep.failure_detail = entry_proof ? "table entry residual was not proved zero" : entry_proof.error().message;
        return rep;
    }
    auto integrated = integrator.integrate_checked(*pattern, kVarName, context);
    if (!integrated) {
        rep.failed = true;
        rep.failure_detail = integrated.error().message;
        return rep;
    }
    auto result_simp = detail::make_expression_ptr(integrated.value());
    if (detail::contains_node_type<IntegralNode>(detail::node(result_simp))) {
        rep.failed = true;
        rep.failure_detail = "integrator left an unevaluated integral: " + result_simp->to_string();
        return rep;
    }
    auto differentiated = LMCAS::differentiate(result_simp, kVarName, context);
    if (!differentiated) {
        rep.failed = true;
        rep.failure_detail = differentiated.error().message;
        return rep;
    }
    auto deriv_simp = differentiated.value();
    auto proof = verify_on_real_branch(result_simp, pattern);
    if (!proof || !std::holds_alternative<ProvedZeroResidual>(proof.value())) {
        rep.failed = true;
        rep.failure_detail = proof ? "derivative residual was not proved zero" : proof.error().message;
        return rep;
    }
    rep.checked = true;
    auto pattern_simp = pattern;

    int checked_points = 0;
    for (double xv : sample_points()) {
        double pattern_value = 0.0;
        double derivative_value = 0.0;
        if (!evaluate_roundtrip_point(pattern_simp, deriv_simp, xv,
                                      pattern_value, derivative_value)) {
            continue;
        }

        ++checked_points;
        double entry_value = 0.0;
        double direct_value = 0.0;
        if (!evaluate_roundtrip_point(pattern_simp, entry_derivative.value(), xv,
                                      entry_value, direct_value) ||
            std::abs(entry_value - direct_value) > kTolerance) {
            rep.failed = true;
            rep.failure_detail = "table primitive derivative disagrees with its pattern";
            break;
        }
        double delta = std::abs(pattern_value - derivative_value);
        if (delta > kTolerance) {
            rep.failed = true;
            std::ostringstream oss;
            oss << "x=" << xv
                << ": d/dx(result)=" << derivative_value
                << " vs pattern=" << pattern_value
                << " |delta|=" << delta
                << " | pattern_expr=" << pattern_simp->to_string()
                << " | result_expr=" << result_simp->to_string()
                << " | derivative_expr=" << deriv_simp->to_string();
            rep.failure_detail = oss.str();
            break;
        }
    }
    if (!rep.failed && checked_points < 3) {
        rep.failed = true;
        rep.failure_detail = "fewer than three real derivative samples were checked";
    }
    return rep;
}


} // anonymous namespace

TEST(IntegrationTableRoundtrip, EntriesDifferentiateToTheirPatterns) {
    const auto entries = detail::make_integration_entries();
    int verified_entries = 0;
    for (const auto &entry : entries) {
        SCOPED_TRACE(entry.name);
        const EntryReport rep = verify_entry(entry);
        EXPECT_FALSE(rep.failed) << rep.failure_detail;
        if (!rep.failed && rep.checked) {
            ++verified_entries;
        }
    }
    EXPECT_GT(verified_entries, 0)
        << "at least one entry verified by exact derivative residual";
}
