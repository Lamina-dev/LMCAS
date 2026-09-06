/**
 * @file assumption_context.cpp
 * @brief Implementation of the AssumptionContext class.
 *
 * Provides scoped push/pop management with read-through query semantics.
 * Each scope has its own PropertyStore and RelationStore. Queries search
 * from the top scope down to root, with child declarations shadowing parent.
 */

#include "assumption_context.hpp"
#include "computation_context.hpp"
#include "inference_engine.hpp"
#include "interval.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <cctype>
#include <cmath>
#include <sstream>

namespace LMCAS {

namespace {


void require_deserialization_update(const Result<void>& result,
                                    int line,
                                    const std::string& keyword) {
    if (result) {
        return;
    }
    if (result.error().code == CasErrc::ResourceLimit) {
        throw std::bad_alloc();
    }
    if (result.error().code == CasErrc::InvalidArgument) {
        throw std::invalid_argument(
            "Line " + std::to_string(line) + ": " + keyword + ": " +
            result.error().message);
    }
    throw std::runtime_error(result.error().message);
}

} // anonymous namespace

// Construction

AssumptionContext::AssumptionContext() {
    // Start with root scope
    scope_stack_.emplace_back();
}

// Scope management

void AssumptionContext::push() {
    scope_stack_.emplace_back();
    ++cache_generation_;
}

AssumptionVoidResult AssumptionContext::pop() {
    if (scope_stack_.size() <= 1) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "cannot pop the root assumption scope", "assumption.pop");
    }
    scope_stack_.pop_back();
    ++cache_generation_;
    return AssumptionVoidResult::success();
}

int AssumptionContext::depth() const {
    return static_cast<int>(scope_stack_.size());
}

// Direct access to current (top) scope stores

PropertyStore& AssumptionContext::current_properties() {
    return scope_stack_.back().properties;
}

const PropertyStore& AssumptionContext::current_properties() const {
    return scope_stack_.back().properties;
}

RelationStore& AssumptionContext::current_relations() {
    return scope_stack_.back().relations;
}

const RelationStore& AssumptionContext::current_relations() const {
    return scope_stack_.back().relations;
}

// Read-through query methods

bool AssumptionContext::has_sign(const std::string& symbol, Sign sign) const {
    // Search from top scope down to root.
    // Return the result from the first scope that has sign info for this symbol.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (!it->properties.get_signs(symbol).empty()) {
            return it->properties.has_sign(symbol, sign);
        }
    }
    // No scope has sign info for this symbol
    return false;
}

bool AssumptionContext::has_domain(const std::string& symbol, Domain domain) const {
    // Search from top scope down to root.
    // Return the result from the first scope that has domain info for this symbol.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_domain(symbol) != Domain::Complex) {
            return it->properties.has_domain(symbol, domain);
        }
    }
    // No scope has domain info — default is Complex.
    // has_domain checks if the symbol has at least the given specificity.
    return domain == Domain::Complex;
}

Domain AssumptionContext::get_domain(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_domain(symbol) != Domain::Complex) {
            return it->properties.get_domain(symbol);
        }
    }
    return Domain::Complex;
}

std::unordered_set<Sign, SignHash> AssumptionContext::get_signs(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (!it->properties.get_signs(symbol).empty()) {
            return it->properties.get_signs(symbol);
        }
    }
    return {};
}

Parity AssumptionContext::get_parity(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_parity(symbol) != Parity::Unknown) {
            return it->properties.get_parity(symbol);
        }
    }
    return Parity::Unknown;
}

Boundedness AssumptionContext::get_boundedness(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_boundedness(symbol) != Boundedness::Unknown) {
            return it->properties.get_boundedness(symbol);
        }
    }
    return Boundedness::Unknown;
}

std::optional<Interval> AssumptionContext::get_bounds(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_bounds(symbol).has_value()) {
            return it->properties.get_bounds(symbol);
        }
    }
    return std::nullopt;
}

// Conditional assumptions

AssumptionVoidResult AssumptionContext::assume_conditional(
    const SymbolicExpr& condition,
    const SymbolicExpr& conclusion) {
    return assume_conditional_checked(condition, conclusion);
}

AssumptionVoidResult AssumptionContext::assume_conditional_checked(
    const SymbolicExpr& condition,
    const SymbolicExpr& conclusion) {
    constexpr const char* operation = "assume_conditional";
    if (!std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(condition))) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "condition expression must be relational", operation);
    }
    if (!std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(conclusion))) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "conclusion expression must be relational", operation);
    }
    try {
        Tribool cond_result = evaluate_condition(condition);
        if (cond_result == Tribool::True &&
            evaluate_condition(conclusion) == Tribool::False) {
            return AssumptionVoidResult::failure(
                CasErrc::InvalidArgument,
                "condition is satisfied but conclusion contradicts the current assumption state",
                operation);
        }
        scope_stack_.back().conditionals.push_back({condition, conclusion});
        ++cache_generation_;
    } catch (const std::bad_alloc&) {
        return AssumptionVoidResult::failure(
            CasErrc::ResourceLimit, "assumption allocation failed", operation);
    } catch (const std::invalid_argument& ex) {
        return AssumptionVoidResult::failure(CasErrc::InvalidArgument, ex.what(), operation);
    } catch (const std::exception& ex) {
        return AssumptionVoidResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
    return AssumptionVoidResult::success();
}

std::vector<AssumptionContext::ConditionalAssumption> AssumptionContext::get_active_conditionals() const {
    std::vector<ConditionalAssumption> result;
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        for (const auto& cond : it->conditionals) {
            result.push_back(cond);
        }
    }
    return result;
}

Tribool AssumptionContext::evaluate_condition(const SymbolicExpr& condition) const {
    if (!LMCAS::detail::node(condition)) {
        return Tribool::Unknown;
    }
    auto rel_node = std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(condition));
    if (!rel_node) {
        return Tribool::Unknown;
    }

    auto lhs = LMCAS::detail::expression_from_node(rel_node->left());
    auto rhs = LMCAS::detail::expression_from_node(rel_node->right());
    RelationalNode::Op op = rel_node->op();
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->relations.has_relation(lhs, rhs, op)) {
            return Tribool::True;
        }
    }
    // For patterns like "variable > 0", check if the variable has the corresponding sign.
    auto lhs_var = std::dynamic_pointer_cast<const VariableNode>(rel_node->left());
    auto rhs_num = std::dynamic_pointer_cast<const NumberNode>(rel_node->right());

    if (lhs_var && rhs_num && rhs_num->is_zero()) {
        const std::string& name = lhs_var->name();
        switch (op) {
            case RelationalNode::Op::GT:
                if (has_sign(name, Sign::Positive)) return Tribool::True;
                if (has_sign(name, Sign::Negative) || has_sign(name, Sign::NonPositive) || has_sign(name, Sign::Zero))
                    return Tribool::False;
                break;
            case RelationalNode::Op::GEQ:
                if (has_sign(name, Sign::NonNegative) || has_sign(name, Sign::Positive) || has_sign(name, Sign::Zero))
                    return Tribool::True;
                if (has_sign(name, Sign::Negative))
                    return Tribool::False;
                break;
            case RelationalNode::Op::LT:
                if (has_sign(name, Sign::Negative)) return Tribool::True;
                if (has_sign(name, Sign::Positive) || has_sign(name, Sign::NonNegative) || has_sign(name, Sign::Zero))
                    return Tribool::False;
                break;
            case RelationalNode::Op::LEQ:
                if (has_sign(name, Sign::NonPositive) || has_sign(name, Sign::Negative) || has_sign(name, Sign::Zero))
                    return Tribool::True;
                if (has_sign(name, Sign::Positive))
                    return Tribool::False;
                break;
            case RelationalNode::Op::NEQ:
                if (has_sign(name, Sign::NonZero) || has_sign(name, Sign::Positive) || has_sign(name, Sign::Negative))
                    return Tribool::True;
                if (has_sign(name, Sign::Zero))
                    return Tribool::False;
                break;
            case RelationalNode::Op::EQ:
                if (has_sign(name, Sign::Zero)) return Tribool::True;
                if (has_sign(name, Sign::Positive) || has_sign(name, Sign::Negative) || has_sign(name, Sign::NonZero))
                    return Tribool::False;
                break;
        }
    }
    auto lhs_num = std::dynamic_pointer_cast<const NumberNode>(rel_node->left());
    auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(rel_node->right());

    if (lhs_num && lhs_num->is_zero() && rhs_var) {
        const std::string& name = rhs_var->name();
        switch (op) {
            case RelationalNode::Op::LT:  // 0 < var → var > 0
                if (has_sign(name, Sign::Positive)) return Tribool::True;
                if (has_sign(name, Sign::Negative) || has_sign(name, Sign::NonPositive) || has_sign(name, Sign::Zero))
                    return Tribool::False;
                break;
            case RelationalNode::Op::LEQ: // 0 <= var → var >= 0
                if (has_sign(name, Sign::NonNegative) || has_sign(name, Sign::Positive) || has_sign(name, Sign::Zero))
                    return Tribool::True;
                if (has_sign(name, Sign::Negative))
                    return Tribool::False;
                break;
            case RelationalNode::Op::GT:  // 0 > var → var < 0
                if (has_sign(name, Sign::Negative)) return Tribool::True;
                if (has_sign(name, Sign::Positive) || has_sign(name, Sign::NonNegative) || has_sign(name, Sign::Zero))
                    return Tribool::False;
                break;
            case RelationalNode::Op::GEQ: // 0 >= var → var <= 0
                if (has_sign(name, Sign::NonPositive) || has_sign(name, Sign::Negative) || has_sign(name, Sign::Zero))
                    return Tribool::True;
                if (has_sign(name, Sign::Positive))
                    return Tribool::False;
                break;
            case RelationalNode::Op::NEQ: // 0 != var → var != 0
                if (has_sign(name, Sign::NonZero) || has_sign(name, Sign::Positive) || has_sign(name, Sign::Negative))
                    return Tribool::True;
                if (has_sign(name, Sign::Zero))
                    return Tribool::False;
                break;
            case RelationalNode::Op::EQ:  // 0 == var → var == 0
                if (has_sign(name, Sign::Zero)) return Tribool::True;
                if (has_sign(name, Sign::Positive) || has_sign(name, Sign::Negative) || has_sign(name, Sign::NonZero))
                    return Tribool::False;
                break;
        }
    }
    return Tribool::Unknown;
}

AssumptionVoidResult AssumptionContext::assume_domain(
    const std::string& variable,
    Domain domain) {
    return assume_domain_checked(variable, domain);
}

AssumptionVoidResult AssumptionContext::assume_domain_checked(
    const std::string& variable,
    Domain domain) {
    constexpr const char* operation = "assume_domain";
    if (variable.empty()) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "variable name must not be empty", operation);
    }
    auto result = scope_stack_.back().properties.declare_domain_checked(variable, domain);
    if (!result) {
        return AssumptionVoidResult::failure(
            result.error().code, result.error().message, operation);
    }
    ++cache_generation_;
    return AssumptionVoidResult::success();
}

AssumptionVoidResult AssumptionContext::assume_sign(
    const std::string& variable,
    Sign sign) {
    return assume_sign_checked(variable, sign);
}

AssumptionVoidResult AssumptionContext::assume_sign_checked(
    const std::string& variable,
    Sign sign) {
    constexpr const char* operation = "assume_sign";
    if (variable.empty()) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "variable name must not be empty", operation);
    }
    auto result = scope_stack_.back().properties.declare_sign_checked(variable, sign);
    if (!result) {
        return AssumptionVoidResult::failure(
            result.error().code, result.error().message, operation);
    }
    ++cache_generation_;
    return AssumptionVoidResult::success();
}

AssumptionVoidResult AssumptionContext::assume(const SymbolicExpr& relation) {
    return assume_checked(relation);
}

AssumptionVoidResult AssumptionContext::assume_checked(const SymbolicExpr& relation) {
    constexpr const char* operation = "assume";
    if (!LMCAS::detail::node(relation)) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "relation expression must not be null", operation);
    }
    if (!std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(relation))) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "relation expression root must be relational", operation);
    }
    try {
        auto rel_node = std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(relation));
        auto lhs = LMCAS::detail::expression_from_node(rel_node->left());
        auto rhs = LMCAS::detail::expression_from_node(rel_node->right());
        auto result = scope_stack_.back().relations.add_relation_checked(
            lhs, rhs, rel_node->op(), scope_stack_.back().properties);
        if (!result.has_value()) {
            return AssumptionVoidResult::failure(result.error());
        }
        ++cache_generation_;
    } catch (const std::bad_alloc&) {
        return AssumptionVoidResult::failure(
            CasErrc::ResourceLimit, "assumption allocation failed", operation);
    } catch (const std::invalid_argument& ex) {
        return AssumptionVoidResult::failure(CasErrc::InvalidArgument, ex.what(), operation);
    } catch (const std::exception& ex) {
        return AssumptionVoidResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
    return AssumptionVoidResult::success();
}

// Expression property queries

namespace {

AssumptionTriboolResult assumption_query_result(InferenceTriboolResult result) {
    if (!result) {
        return AssumptionTriboolResult::failure(result.error());
    }
    return AssumptionTriboolResult::success(result.value());
}


} // anonymous namespace

AssumptionTriboolResult AssumptionContext::is_positive(const SymbolicExpr& expr) const {
    return is_positive_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_positive_checked(const SymbolicExpr& expr) const {
    InferenceEngine inference(*this);
    return assumption_query_result(inference.query_positive_checked(expr));
}

AssumptionTriboolResult AssumptionContext::is_negative(const SymbolicExpr& expr) const {
    return is_negative_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_negative_checked(const SymbolicExpr& expr) const {
    InferenceEngine inference(*this);
    return assumption_query_result(inference.query_negative_checked(expr));
}

AssumptionTriboolResult AssumptionContext::is_nonnegative(const SymbolicExpr& expr) const {
    return is_nonnegative_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_nonnegative_checked(const SymbolicExpr& expr) const {
    InferenceEngine inference(*this);
    return assumption_query_result(inference.query_nonnegative_checked(expr));
}

AssumptionTriboolResult AssumptionContext::is_real(const SymbolicExpr& expr) const {
    return is_real_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_real_checked(const SymbolicExpr& expr) const {
    InferenceEngine inference(*this);
    return assumption_query_result(inference.query_real_checked(expr));
}

AssumptionTriboolResult AssumptionContext::is_integer(const SymbolicExpr& expr) const {
    return is_integer_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_integer_checked(const SymbolicExpr& expr) const {
    InferenceEngine inference(*this);
    return assumption_query_result(inference.query_integer_checked(expr));
}

AssumptionTriboolResult AssumptionContext::is_nonzero(const SymbolicExpr& expr) const {
    return is_nonzero_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_nonzero_checked(const SymbolicExpr& expr) const {
    InferenceEngine inference(*this);
    return assumption_query_result(inference.query_nonzero_checked(expr));
}

// Extended query methods (read-through all scopes)

AssumptionTriboolResult AssumptionContext::is_continuous(
    const std::string& symbol,
    const Interval& interval) const {
    return is_continuous_checked(symbol, interval);
}

AssumptionTriboolResult AssumptionContext::is_continuous_checked(
    const std::string& symbol,
    const Interval& interval) const {
    constexpr const char* operation = "is_continuous";
    if (symbol.empty()) {
        return AssumptionTriboolResult::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto result = it->properties.is_continuous_checked(symbol, interval);
        if (!result) {
            return AssumptionTriboolResult::failure(result.error());
        }
        if (result.value()) {
            return AssumptionTriboolResult::success(Tribool::True);
        }
    }
    return AssumptionTriboolResult::success(Tribool::Unknown);
}

AssumptionTriboolResult AssumptionContext::is_differentiable(
    const std::string& symbol,
    const Interval& interval) const {
    return is_differentiable_checked(symbol, interval);
}

AssumptionTriboolResult AssumptionContext::is_differentiable_checked(
    const std::string& symbol,
    const Interval& interval) const {
    constexpr const char* operation = "is_differentiable";
    if (symbol.empty()) {
        return AssumptionTriboolResult::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto result = it->properties.is_differentiable_checked(symbol, interval);
        if (!result) {
            return AssumptionTriboolResult::failure(result.error());
        }
        if (result.value()) {
            return AssumptionTriboolResult::success(Tribool::True);
        }
    }
    return AssumptionTriboolResult::success(Tribool::Unknown);
}

AssumptionTriboolResult AssumptionContext::is_positive_definite(
    const std::string& symbol) const {
    return is_positive_definite_checked(symbol);
}

AssumptionTriboolResult AssumptionContext::is_positive_definite_checked(
    const std::string& symbol) const {
    constexpr const char* operation = "is_positive_definite";
    if (symbol.empty()) {
        return AssumptionTriboolResult::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        Definiteness d = it->properties.get_definiteness(symbol);
        if (d != Definiteness::Unknown) {
            return AssumptionTriboolResult::success(
                (d == Definiteness::PositiveDefinite) ? Tribool::True : Tribool::False);
        }
    }
    return AssumptionTriboolResult::success(Tribool::Unknown);
}

AssumptionTriboolResult AssumptionContext::is_positive_semidefinite(
    const std::string& symbol) const {
    return is_positive_semidefinite_checked(symbol);
}

AssumptionTriboolResult AssumptionContext::is_positive_semidefinite_checked(
    const std::string& symbol) const {
    constexpr const char* operation = "is_positive_semidefinite";
    if (symbol.empty()) {
        return AssumptionTriboolResult::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        Definiteness d = it->properties.get_definiteness(symbol);
        if (d != Definiteness::Unknown) {
            return AssumptionTriboolResult::success(
                (d == Definiteness::PositiveDefinite ||
                 d == Definiteness::PositiveSemiDefinite) ? Tribool::True : Tribool::False);
        }
    }
    return AssumptionTriboolResult::success(Tribool::Unknown);
}

// Depth limit configuration

void AssumptionContext::set_max_query_depth(int depth) {
    if (depth <= 0) return;
    max_query_depth_ = depth;
}

int AssumptionContext::get_max_query_depth() const {
    return max_query_depth_;
}

// Serialization helpers (file-local)

namespace {

std::string domain_to_string(Domain d) {
    switch (d) {
        case Domain::Complex:     return "Complex";
        case Domain::Real:        return "Real";
        case Domain::Algebraic:   return "Algebraic";
        case Domain::Rational:    return "Rational";
        case Domain::Integer:     return "Integer";
        case Domain::Natural:     return "Natural";
        case Domain::PositiveInt: return "PositiveInt";
    }
    return "Complex";
}

std::string sign_to_string(Sign s) {
    switch (s) {
        case Sign::Positive:    return "Positive";
        case Sign::Negative:    return "Negative";
        case Sign::NonNegative: return "NonNegative";
        case Sign::NonPositive: return "NonPositive";
        case Sign::Zero:        return "Zero";
        case Sign::NonZero:     return "NonZero";
    }
    return "Positive";
}

std::string parity_to_string(Parity p) {
    switch (p) {
        case Parity::Even:    return "Even";
        case Parity::Odd:     return "Odd";
        case Parity::Unknown: return "Unknown";
    }
    return "Unknown";
}

std::string boundedness_to_string(Boundedness b) {
    switch (b) {
        case Boundedness::Bounded:   return "Bounded";
        case Boundedness::Unbounded: return "Unbounded";
        case Boundedness::Unknown:   return "Unknown";
    }
    return "Unknown";
}

std::string finiteness_to_string(Finiteness f) {
    switch (f) {
        case Finiteness::Finite:    return "Finite";
        case Finiteness::Divergent: return "Divergent";
        case Finiteness::Unknown:   return "Unknown";
    }
    return "Unknown";
}

std::string definiteness_to_string(Definiteness d) {
    switch (d) {
        case Definiteness::PositiveDefinite:     return "PositiveDefinite";
        case Definiteness::PositiveSemiDefinite: return "PositiveSemiDefinite";
        case Definiteness::NegativeDefinite:     return "NegativeDefinite";
        case Definiteness::NegativeSemiDefinite: return "NegativeSemiDefinite";
        case Definiteness::Indefinite:           return "Indefinite";
        case Definiteness::Unknown:              return "Unknown";
    }
    return "Unknown";
}

std::string monotonicity_to_string(Monotonicity m) {
    switch (m) {
        case Monotonicity::Increasing:    return "Increasing";
        case Monotonicity::Decreasing:    return "Decreasing";
        case Monotonicity::NonDecreasing: return "NonDecreasing";
        case Monotonicity::NonIncreasing: return "NonIncreasing";
        case Monotonicity::Unknown:       return "Unknown";
    }
    return "Unknown";
}

std::string relop_to_string(RelationalNode::Op op) {
    switch (op) {
        case RelationalNode::Op::GT:  return "GT";
        case RelationalNode::Op::LT:  return "LT";
        case RelationalNode::Op::GEQ: return "GEQ";
        case RelationalNode::Op::LEQ: return "LEQ";
        case RelationalNode::Op::NEQ: return "NEQ";
        case RelationalNode::Op::EQ:  return "EQ";
    }
    return "EQ";
}

} // anonymous namespace

namespace {

Domain parse_domain(const std::string& s, int line_num) {
    if (s == "Complex")     return Domain::Complex;
    if (s == "Real")        return Domain::Real;
    if (s == "Algebraic")   return Domain::Algebraic;
    if (s == "Rational")    return Domain::Rational;
    if (s == "Integer")     return Domain::Integer;
    if (s == "Natural")     return Domain::Natural;
    if (s == "PositiveInt") return Domain::PositiveInt;
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown domain '" + s + "'");
}

Sign parse_sign(const std::string& s, int line_num) {
    if (s == "Positive")    return Sign::Positive;
    if (s == "Negative")    return Sign::Negative;
    if (s == "NonNegative") return Sign::NonNegative;
    if (s == "NonPositive") return Sign::NonPositive;
    if (s == "Zero")        return Sign::Zero;
    if (s == "NonZero")     return Sign::NonZero;
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown sign '" + s + "'");
}

Parity parse_parity(const std::string& s, int line_num) {
    if (s == "Even") return Parity::Even;
    if (s == "Odd")  return Parity::Odd;
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown parity '" + s + "'");
}

Boundedness parse_boundedness(const std::string& s, int line_num) {
    if (s == "Bounded")   return Boundedness::Bounded;
    if (s == "Unbounded") return Boundedness::Unbounded;
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown boundedness '" + s + "'");
}

Finiteness parse_finiteness(const std::string& s, int line_num) {
    if (s == "Finite")    return Finiteness::Finite;
    if (s == "Divergent") return Finiteness::Divergent;
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown finiteness '" + s + "'");
}

Definiteness parse_definiteness(const std::string& s, int line_num) {
    if (s == "PositiveDefinite")     return Definiteness::PositiveDefinite;
    if (s == "PositiveSemiDefinite") return Definiteness::PositiveSemiDefinite;
    if (s == "NegativeDefinite")     return Definiteness::NegativeDefinite;
    if (s == "NegativeSemiDefinite") return Definiteness::NegativeSemiDefinite;
    if (s == "Indefinite")           return Definiteness::Indefinite;
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown definiteness '" + s + "'");
}

Monotonicity parse_monotonicity(const std::string& s, int line_num) {
    if (s == "Increasing")    return Monotonicity::Increasing;
    if (s == "Decreasing")    return Monotonicity::Decreasing;
    if (s == "NonDecreasing") return Monotonicity::NonDecreasing;
    if (s == "NonIncreasing") return Monotonicity::NonIncreasing;
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown monotonicity '" + s + "'");
}

RelationalNode::Op parse_relop(const std::string& s, int line_num) {
    if (s == "GT")  return RelationalNode::Op::GT;
    if (s == "LT")  return RelationalNode::Op::LT;
    if (s == "GEQ") return RelationalNode::Op::GEQ;
    if (s == "LEQ") return RelationalNode::Op::LEQ;
    if (s == "NEQ") return RelationalNode::Op::NEQ;
    if (s == "EQ")  return RelationalNode::Op::EQ;
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown relational operator '" + s + "'");
}

void require_line_end(std::istringstream& input,
                      int line_num,
                      const std::string& keyword) {
    std::string trailing;
    if (input >> trailing) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": " + keyword +
            " has unexpected trailing field '" + trailing + "'");
    }
}

std::shared_ptr<const SymbolicNode> parse_finite_number_node(
    const std::string& token,
    int line_num) {
    const auto invalid = [&]() {
        return std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": invalid finite number '" + token + "'");
    };
    if (token.empty()) {
        throw invalid();
    }

    const auto slash = token.find('/');
    if (slash != std::string::npos) {
        if (slash == 0 || slash + 1 == token.size() ||
            token.find('/', slash + 1) != std::string::npos) {
            throw invalid();
        }
        return LMCAS::detail::make_node<NumberNode>(
            Rational(
                BigInt(token.substr(0, slash)),
                BigInt(token.substr(slash + 1))));
    }
    if (token.find_first_of(".eE") == std::string::npos) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(token));
    }

    std::size_t consumed = 0;
    const double value = std::stod(token, &consumed);
    if (consumed != token.size() || !std::isfinite(value)) {
        throw invalid();
    }
    return LMCAS::detail::make_node<NumberNode>(
        static_cast<lmmc_real_t>(value));
}

std::shared_ptr<const SymbolicNode> parse_serialized_atom_node(
    const std::string& token,
    int line_num) {
    if (token.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": empty expression atom");
    }

    const auto is_identifier_start = [](unsigned char c) {
        return std::isalpha(c) != 0 || c == '_' || c >= 0x80;
    };
    const auto is_identifier_continue = [](unsigned char c) {
        return std::isalnum(c) != 0 || c == '_' || c >= 0x80;
    };
    if (is_identifier_start(static_cast<unsigned char>(token.front()))) {
        for (char c : token) {
            if (!is_identifier_continue(static_cast<unsigned char>(c))) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": unsupported serialized expression '" + token + "'");
            }
        }
        return LMCAS::detail::make_node<VariableNode>(token);
    }
    return parse_finite_number_node(token, line_num);
}

/// Serialize an interval endpoint to a numeric string.
std::string endpoint_to_string(const Endpoint& ep) {
    if (ep.is_neg_infinity) return "-inf";
    if (ep.is_pos_infinity) return "+inf";
    if (ep.value) return ep.value->to_string();
    return "0";
}

/// Serialize an interval to string format: [lo, hi] or (lo, hi) etc.
std::string interval_to_string(const Interval& iv) {
    std::string result;
    result += iv.lower.is_open ? "(" : "[";
    result += endpoint_to_string(iv.lower);
    result += ", ";
    result += endpoint_to_string(iv.upper);
    result += iv.upper.is_open ? ")" : "]";
    return result;
}

/// Parse an interval from the canonical format: [lo, hi], (lo, hi], etc.
Interval parse_interval(const std::string& s, int line_num) {
    const auto malformed = [&]() {
        return std::invalid_argument(
            "Line " + std::to_string(line_num) + ": malformed interval '" + s + "'");
    };
    if (s.size() < 5 ||
        (s.front() != '[' && s.front() != '(') ||
        (s.back() != ']' && s.back() != ')')) {
        throw malformed();
    }

    const bool lower_open = s.front() == '(';
    const bool upper_open = s.back() == ')';
    const std::string inner = s.substr(1, s.size() - 2);
    const auto comma_pos = inner.find(", ");
    if (comma_pos == std::string::npos ||
        inner.find(", ", comma_pos + 2) != std::string::npos) {
        throw malformed();
    }

    const std::string lo_str = inner.substr(0, comma_pos);
    const std::string hi_str = inner.substr(comma_pos + 2);
    if (lo_str.empty() || hi_str.empty()) {
        throw malformed();
    }
    if ((lo_str == "-inf" && !lower_open) ||
        (hi_str == "+inf" && !upper_open)) {
        throw malformed();
    }

    const auto parse_finite_number =
        [&](const std::string& token) -> std::shared_ptr<SymbolicExpr> {
        try {
            return LMCAS::detail::make_expression_ptr(
                parse_finite_number_node(token, line_num));
        } catch (const std::invalid_argument&) {
            throw malformed();
        }
    };

    Interval iv;
    if (lo_str == "-inf") {
        iv.lower = Endpoint::neg_inf();
    } else {
        const auto lo_expr = parse_finite_number(lo_str);
        iv.lower = lower_open ? Endpoint::open(lo_expr) : Endpoint::closed(lo_expr);
    }
    if (hi_str == "+inf") {
        iv.upper = Endpoint::pos_inf();
    } else {
        const auto hi_expr = parse_finite_number(hi_str);
        iv.upper = upper_open ? Endpoint::open(hi_expr) : Endpoint::closed(hi_expr);
    }

    return iv;
}

} // anonymous namespace

// Serialization

std::string AssumptionContext::serialize() const {
    std::ostringstream out;

    for (int scope_idx = 0; scope_idx < static_cast<int>(scope_stack_.size()); ++scope_idx) {
        out << "SCOPE " << scope_idx << "\n";
        const auto& scope = scope_stack_[scope_idx];
        const auto& props = scope.properties;
        auto symbols = props.get_all_symbols();
        for (const auto& sym : symbols) {
            Domain dom = props.get_domain(sym);
            if (dom != Domain::Complex) {
                out << "DOMAIN " << sym << " " << domain_to_string(dom) << "\n";
            }
            auto signs = props.get_signs(sym);
            for (Sign s : signs) {
                out << "SIGN " << sym << " " << sign_to_string(s) << "\n";
            }
            Parity par = props.get_parity(sym);
            if (par != Parity::Unknown) {
                out << "PARITY " << sym << " " << parity_to_string(par) << "\n";
            }
            Boundedness bnd = props.get_boundedness(sym);
            if (bnd != Boundedness::Unknown) {
                out << "BOUNDED " << sym << " " << boundedness_to_string(bnd) << "\n";
            }
            if (props.is_transcendental(sym)) {
                out << "TRANSCENDENTAL " << sym << "\n";
            }
            Finiteness fin = props.get_finiteness(sym);
            if (fin != Finiteness::Unknown) {
                out << "FINITENESS " << sym << " " << finiteness_to_string(fin) << "\n";
            }
            Definiteness def = props.get_definiteness(sym);
            if (def != Definiteness::Unknown) {
                out << "DEFINITENESS " << sym << " " << definiteness_to_string(def) << "\n";
            }
            auto period = props.get_period(sym);
            if (period.has_value() && *period) {
                out << "PERIODIC " << sym << " " << (*period)->to_string() << "\n";
            }
        }
        for (const auto& sym : symbols) {
            auto cont_decls = props.get_continuity_decls(sym);
            for (const auto& decl : cont_decls) {
                if (decl.is_differentiable) {
                    out << "DIFFERENTIABLE " << sym << " "
                        << interval_to_string(decl.interval) << "\n";
                } else {
                    out << "CONTINUOUS " << sym << " "
                        << interval_to_string(decl.interval) << "\n";
                }
            }

            auto mono_decls = props.get_monotonicity_decls(sym);
            for (const auto& decl : mono_decls) {
                out << "MONOTONICITY " << sym << " " << decl.variable << " "
                    << interval_to_string(decl.interval) << " "
                    << monotonicity_to_string(decl.type) << "\n";
            }
        }
        const auto& relations = scope.relations.get_relations();
        for (const auto& rel : relations) {
            std::string lhs_str = rel.lhs.to_string();
            std::string rhs_str = rel.rhs.to_string();
            out << "RELATION " << lhs_str << " " << relop_to_string(rel.op)
                << " " << rhs_str << "\n";
        }
        for (const auto& cond : scope.conditionals) {
            std::string cond_str = cond.condition.to_string();
            std::string concl_str = cond.conclusion.to_string();
            out << "CONDITIONAL (" << cond_str << ") => (" << concl_str << ")\n";
        }
    }

    out << "END\n";
    return out.str();
}

// Deserialization

AssumptionContext AssumptionContext::deserialize_impl(const std::string& data) {
    AssumptionContext ctx;
    std::istringstream input(data);
    std::string line;
    int line_num = 0;
    int current_scope = -1;
    bool ended = false;

    while (std::getline(input, line)) {
        ++line_num;
        if (line.empty()) continue;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) continue;
        if (ended) {
            throw std::invalid_argument(
                "Line " + std::to_string(line_num) + ": data after END");
        }
        if (line == "END") {
            ended = true;
            continue;
        }
        std::istringstream ls(line);
        std::string keyword;
        ls >> keyword;

        if (keyword == "SCOPE") {
            int idx;
            if (!(ls >> idx)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": SCOPE missing index");
            }
            require_line_end(ls, line_num, keyword);
            if (idx != current_scope + 1) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": SCOPE indices must be sequential from zero");
            }
            if (idx > 0) {
                ctx.push();
            }
            current_scope = idx;
        } else if (keyword == "DOMAIN") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": DOMAIN before SCOPE");
            }
            std::string sym, dom_str;
            if (!(ls >> sym >> dom_str)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": DOMAIN requires symbol and domain");
            }
            require_line_end(ls, line_num, keyword);
            Domain dom = parse_domain(dom_str, line_num);
            require_deserialization_update(
                ctx.current_properties().declare_domain_checked(sym, dom), line_num, keyword);
        } else if (keyword == "SIGN") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": SIGN before SCOPE");
            }
            std::string sym, sign_str;
            if (!(ls >> sym >> sign_str)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": SIGN requires symbol and sign");
            }
            require_line_end(ls, line_num, keyword);
            Sign s = parse_sign(sign_str, line_num);
            require_deserialization_update(
                ctx.current_properties().declare_sign_checked(sym, s), line_num, keyword);
        } else if (keyword == "PARITY") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": PARITY before SCOPE");
            }
            std::string sym, par_str;
            if (!(ls >> sym >> par_str)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": PARITY requires symbol and parity");
            }
            require_line_end(ls, line_num, keyword);
            Parity p = parse_parity(par_str, line_num);
            require_deserialization_update(
                ctx.current_properties().declare_parity_checked(sym, p), line_num, keyword);
        } else if (keyword == "BOUNDED") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": BOUNDED before SCOPE");
            }
            std::string sym, bnd_str;
            if (!(ls >> sym >> bnd_str)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": BOUNDED requires symbol and boundedness");
            }
            require_line_end(ls, line_num, keyword);
            Boundedness b = parse_boundedness(bnd_str, line_num);
            require_deserialization_update(
                ctx.current_properties().declare_bounded_checked(sym, b), line_num, keyword);
        } else if (keyword == "TRANSCENDENTAL") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": TRANSCENDENTAL before SCOPE");
            }
            std::string sym;
            if (!(ls >> sym)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": TRANSCENDENTAL requires symbol");
            }
            require_line_end(ls, line_num, keyword);
            require_deserialization_update(
                ctx.current_properties().declare_transcendental_checked(sym), line_num, keyword);
        } else if (keyword == "FINITENESS") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": FINITENESS before SCOPE");
            }
            std::string sym, fin_str;
            if (!(ls >> sym >> fin_str)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": FINITENESS requires symbol and value");
            }
            Finiteness f = parse_finiteness(fin_str, line_num);
            require_line_end(ls, line_num, keyword);
            require_deserialization_update(
                ctx.current_properties().declare_finiteness_checked(sym, f), line_num, keyword);
        } else if (keyword == "DEFINITENESS") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": DEFINITENESS before SCOPE");
            }
            std::string sym, def_str;
            if (!(ls >> sym >> def_str)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": DEFINITENESS requires symbol and value");
            }
            Definiteness d = parse_definiteness(def_str, line_num);
            require_line_end(ls, line_num, keyword);
            require_deserialization_update(
                ctx.current_properties().declare_definiteness_checked(sym, d), line_num, keyword);
        } else if (keyword == "PERIODIC") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": PERIODIC before SCOPE");
            }
            std::string sym;
            if (!(ls >> sym)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": PERIODIC requires symbol");
            }
            std::string period_str;
            std::getline(ls, period_str);
            size_t start = period_str.find_first_not_of(' ');
            if (start == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": PERIODIC requires period value");
            }
            period_str = period_str.substr(start);
            std::shared_ptr<SymbolicExpr> period_expr =
                LMCAS::detail::make_expression_ptr(
                    parse_serialized_atom_node(period_str, line_num));
            require_deserialization_update(
                ctx.current_properties().declare_periodic_checked(sym, period_expr),
                line_num, keyword);
        } else if (keyword == "RELATION") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": RELATION before SCOPE");
            }
            // The rest of the line after "RELATION " is: lhs_str OP rhs_str
            std::string rest;
            std::getline(ls, rest);
            size_t rstart = rest.find_first_not_of(' ');
            if (rstart == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": RELATION requires lhs, op, and rhs");
            }
            rest = rest.substr(rstart);
            std::string op_str;
            size_t op_pos = std::string::npos;
            for (const auto& op_candidate : {"GEQ", "LEQ", "NEQ"}) {
                std::string search = std::string(" ") + op_candidate + " ";
                size_t pos = rest.find(search);
                if (pos != std::string::npos) {
                    op_str = op_candidate;
                    op_pos = pos;
                    break;
                }
            }
            if (op_pos == std::string::npos) {
                for (const auto& op_candidate : {"GT", "LT", "EQ"}) {
                    std::string search = std::string(" ") + op_candidate + " ";
                    size_t pos = rest.find(search);
                    if (pos != std::string::npos) {
                        op_str = op_candidate;
                        op_pos = pos;
                        break;
                    }
                }
            }

            if (op_pos == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": RELATION missing operator");
            }

            std::string lhs_str = rest.substr(0, op_pos);
            std::string rhs_str = rest.substr(op_pos + 1 + op_str.size() + 1);
            RelationalNode::Op op = parse_relop(op_str, line_num);
            auto parse_simple_expr = [&](const std::string& s) -> SymbolicExpr {
                return LMCAS::detail::expression_from_node(
                    parse_serialized_atom_node(s, line_num));
            };

            SymbolicExpr lhs = parse_simple_expr(lhs_str);
            SymbolicExpr rhs = parse_simple_expr(rhs_str);
            auto relation_result = ctx.current_relations().add_relation_checked(
                lhs, rhs, op, ctx.current_properties());
            if (!relation_result.has_value()) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": " +
                    relation_result.error().message);
            }
        } else if (keyword == "CONDITIONAL") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": CONDITIONAL before SCOPE");
            }
            std::string rest;
            std::getline(ls, rest);
            size_t rstart = rest.find_first_not_of(' ');
            if (rstart == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": CONDITIONAL requires condition and conclusion");
            }
            rest = rest.substr(rstart);
            size_t arrow_pos = rest.find(") => (");
            if (arrow_pos == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": CONDITIONAL malformed, expected '(cond) => (concl)'");
            }

            // Extract condition string (between first '(' and the arrow)
            if (rest[0] != '(') {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": CONDITIONAL condition must start with '('");
            }
            std::string cond_str = rest.substr(1, arrow_pos - 1);

            // Extract conclusion string (between '(' after '=>' and final ')')
            std::string concl_part = rest.substr(arrow_pos + 6); // skip ") => ("
            if (concl_part.empty() || concl_part.back() != ')') {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": CONDITIONAL conclusion must end with ')'");
            }
            std::string concl_str = concl_part.substr(0, concl_part.size() - 1);
            // The persistence grammar intentionally supports relational atoms.
            auto parse_relational = [&](const std::string& s) -> SymbolicExpr {
                std::string op_token;
                size_t op_pos_local = std::string::npos;
                for (const auto& candidate : {"GEQ", "LEQ", "NEQ"}) {
                    std::string search = std::string(" ") + candidate + " ";
                    size_t pos = s.find(search);
                    if (pos != std::string::npos) {
                        op_token = candidate;
                        op_pos_local = pos;
                        break;
                    }
                }
                if (op_pos_local == std::string::npos) {
                    for (const auto& candidate : {"GT", "LT", "EQ"}) {
                        std::string search = std::string(" ") + candidate + " ";
                        size_t pos = s.find(search);
                        if (pos != std::string::npos) {
                            op_token = candidate;
                            op_pos_local = pos;
                            break;
                        }
                    }
                }
                if (op_pos_local == std::string::npos) {
                    throw std::invalid_argument(
                        "Line " + std::to_string(line_num) +
                        ": CONDITIONAL cannot parse expression '" + s + "'");
                }
                std::string lhs_s = s.substr(0, op_pos_local);
                std::string rhs_s = s.substr(op_pos_local + 1 + op_token.size() + 1);
                RelationalNode::Op rel_op = parse_relop(op_token, line_num);
                auto parse_token =
                    [&](const std::string& tok) -> std::shared_ptr<const SymbolicNode> {
                    return parse_serialized_atom_node(tok, line_num);
                };

                auto lhs_node = parse_token(lhs_s);
                auto rhs_node = parse_token(rhs_s);
                return LMCAS::detail::expression_from_node(LMCAS::detail::make_node<RelationalNode>(lhs_node, rhs_node, rel_op));
            };

            SymbolicExpr cond_expr = parse_relational(cond_str);
            SymbolicExpr concl_expr = parse_relational(concl_str);
            require_deserialization_update(
                ctx.assume_conditional_checked(cond_expr, concl_expr), line_num, keyword);
        } else if (keyword == "CONTINUOUS") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": CONTINUOUS before SCOPE");
            }
            std::string sym;
            if (!(ls >> sym)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": CONTINUOUS requires symbol");
            }
            std::string iv_str;
            std::getline(ls, iv_str);
            size_t ivstart = iv_str.find_first_not_of(' ');
            if (ivstart == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": CONTINUOUS requires interval");
            }
            iv_str = iv_str.substr(ivstart);
            Interval iv = parse_interval(iv_str, line_num);
            require_deserialization_update(
                ctx.current_properties().declare_continuous_checked(sym, iv),
                line_num, keyword);
        } else if (keyword == "DIFFERENTIABLE") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": DIFFERENTIABLE before SCOPE");
            }
            std::string sym;
            if (!(ls >> sym)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": DIFFERENTIABLE requires symbol");
            }
            std::string iv_str;
            std::getline(ls, iv_str);
            size_t ivstart = iv_str.find_first_not_of(' ');
            if (ivstart == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": DIFFERENTIABLE requires interval");
            }
            iv_str = iv_str.substr(ivstart);
            Interval iv = parse_interval(iv_str, line_num);
            require_deserialization_update(
                ctx.current_properties().declare_differentiable_checked(sym, iv),
                line_num, keyword);
        } else if (keyword == "MONOTONICITY") {
            if (current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": MONOTONICITY before SCOPE");
            }
            std::string sym, var;
            if (!(ls >> sym >> var)) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": MONOTONICITY requires symbol and variable");
            }
            std::string rest_line;
            std::getline(ls, rest_line);
            size_t mstart = rest_line.find_first_not_of(' ');
            if (mstart == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": MONOTONICITY requires interval and type");
            }
            rest_line = rest_line.substr(mstart);
            size_t bracket_end = rest_line.find_first_of("])");
            if (bracket_end == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": MONOTONICITY malformed interval");
            }
            std::string mono_iv_str = rest_line.substr(0, bracket_end + 1);
            Interval iv = parse_interval(mono_iv_str, line_num);

            std::string mono_str = rest_line.substr(bracket_end + 2); // skip bracket + space
            size_t mtstart = mono_str.find_first_not_of(' ');
            if (mtstart == std::string::npos) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": MONOTONICITY requires type after interval");
            }
            mono_str = mono_str.substr(mtstart);
            size_t mend = mono_str.find_last_not_of(" \t\r\n");
            if (mend != std::string::npos) mono_str = mono_str.substr(0, mend + 1);

            Monotonicity m = parse_monotonicity(mono_str, line_num);
            require_deserialization_update(
                ctx.current_properties().declare_monotonicity_checked(sym, var, iv, m),
                line_num, keyword);
        } else if (keyword == "END") {
            throw std::invalid_argument(
                "Line " + std::to_string(line_num) + ": malformed END record");
        } else {
            throw std::invalid_argument(
                "Line " + std::to_string(line_num) + ": unknown keyword '" + keyword + "'");
        }
    }

    if (!ended) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": unexpected end of input (missing END)");
    }

    return ctx;
}

Result<AssumptionContext> AssumptionContext::deserialize(const std::string& data) {
    return deserialize_checked(data);
}

Result<AssumptionContext> AssumptionContext::deserialize_checked(
    const std::string& data) {
    ComputationContext context;
    return deserialize_checked(data, context);
}

Result<AssumptionContext> AssumptionContext::deserialize_checked(
    const std::string& data, ComputationContext& context) {
    constexpr const char* operation = "deserialize";
    auto input_budget = context.require_input_bytes(data.size(), operation);
    if (!input_budget) {
        return Result<AssumptionContext>::failure(input_budget.error());
    }
    try {
        return Result<AssumptionContext>::success(deserialize_impl(data));
    } catch (const std::bad_alloc&) {
        return Result<AssumptionContext>::failure(
            CasErrc::ResourceLimit,
            "assumption deserialization allocation failed", operation);
    } catch (const std::invalid_argument& ex) {
        return Result<AssumptionContext>::failure(
            CasErrc::ParseError, ex.what(), operation);
    } catch (const std::exception& ex) {
        return Result<AssumptionContext>::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

} // namespace LMCAS
