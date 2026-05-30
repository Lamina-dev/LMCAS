/**
 * @file assumption_context.cpp
 * @brief Implementation of the AssumptionContext class.
 *
 * Provides scoped push/pop management with read-through query semantics.
 * Each scope has its own PropertyStore and RelationStore. Queries search
 * from the top scope down to root, with child declarations shadowing parent.
 */

#include "assumption_context.hpp"
#include "query_interface.hpp"
#include "inference_engine.hpp"
#include "interval.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <sstream>

namespace lamina {

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

void AssumptionContext::pop() {
    if (scope_stack_.size() <= 1) {
        throw std::runtime_error("Cannot pop root scope");
    }
    scope_stack_.pop_back();
    ++cache_generation_;
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

void AssumptionContext::assume_conditional(const SymbolicExpr& condition, const SymbolicExpr& conclusion) {
    if (!condition.root) {
        throw std::invalid_argument(
            "assume_conditional: condition expression must not be null/empty");
    }
    if (!conclusion.root) {
        throw std::invalid_argument(
            "assume_conditional: conclusion expression must not be null/empty");
    }

    // Check if the condition is already satisfied. If so, verify the conclusion
    // does not contradict the current state.
    Tribool cond_result = evaluate_condition(condition);
    if (cond_result == Tribool::True) {
        Tribool concl_result = evaluate_condition(conclusion);
        if (concl_result == Tribool::False) {
            throw std::invalid_argument(
                "assume_conditional: contradiction detected — condition '" +
                condition.to_string() + "' is satisfied but conclusion '" +
                conclusion.to_string() + "' contradicts the current assumption state");
        }
    }

    scope_stack_.back().conditionals.push_back({condition, conclusion});
    ++cache_generation_;
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
    if (!condition.root) {
        return Tribool::Unknown;
    }
    auto rel_node = std::dynamic_pointer_cast<RelationalNode>(condition.root);
    if (!rel_node) {
        return Tribool::Unknown;
    }

    SymbolicExpr lhs(rel_node->left);
    SymbolicExpr rhs(rel_node->right);
    RelationalNode::Op op = rel_node->op;
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->relations.has_relation(lhs, rhs, op)) {
            return Tribool::True;
        }
    }
    // For patterns like "variable > 0", check if the variable has the corresponding sign.
    auto lhs_var = std::dynamic_pointer_cast<VariableNode>(rel_node->left);
    auto rhs_num = std::dynamic_pointer_cast<NumberNode>(rel_node->right);

    if (lhs_var && rhs_num && rhs_num->is_zero()) {
        const std::string& name = lhs_var->name;
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
    auto lhs_num = std::dynamic_pointer_cast<NumberNode>(rel_node->left);
    auto rhs_var = std::dynamic_pointer_cast<VariableNode>(rel_node->right);

    if (lhs_num && lhs_num->is_zero() && rhs_var) {
        const std::string& name = rhs_var->name;
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

// Convenience declaration API

void AssumptionContext::assume_domain(const std::string& variable, Domain domain) {
    if (variable.empty()) {
        throw std::invalid_argument("assume_domain: variable name must not be empty");
    }
    scope_stack_.back().properties.declare_domain(variable, domain);
    ++cache_generation_;
}

void AssumptionContext::assume_sign(const std::string& variable, Sign sign) {
    if (variable.empty()) {
        throw std::invalid_argument("assume_sign: variable name must not be empty");
    }
    scope_stack_.back().properties.declare_sign(variable, sign);
    ++cache_generation_;
}

void AssumptionContext::assume(const SymbolicExpr& relation) {
    if (!relation.root) {
        throw std::invalid_argument("assume: expression must not be null/empty");
    }
    auto rel_node = std::dynamic_pointer_cast<RelationalNode>(relation.root);
    if (!rel_node) {
        throw std::invalid_argument("assume: expression root must be a RelationalNode");
    }
    // Extract lhs, rhs, and op from the RelationalNode and store in RelationStore
    SymbolicExpr lhs(rel_node->left);
    SymbolicExpr rhs(rel_node->right);
    scope_stack_.back().relations.add_relation(lhs, rhs, rel_node->op,
                                               scope_stack_.back().properties);
    ++cache_generation_;
}

// Convenience query API (delegates to QueryInterface)

Tribool AssumptionContext::is_positive(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_positive: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_positive(expr);
}

Tribool AssumptionContext::is_negative(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_negative: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_negative(expr);
}

Tribool AssumptionContext::is_nonnegative(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_nonnegative: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_nonnegative(expr);
}

Tribool AssumptionContext::is_real(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_real: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_real(expr);
}

Tribool AssumptionContext::is_integer(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_integer: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_integer(expr);
}

Tribool AssumptionContext::is_nonzero(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_nonzero: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_nonzero(expr);
}

// Extended query methods (read-through all scopes)

Tribool AssumptionContext::is_continuous(const std::string& symbol, const Interval& interval) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.is_continuous(symbol, interval)) {
            return Tribool::True;
        }
    }
    return Tribool::Unknown;
}

Tribool AssumptionContext::is_differentiable(const std::string& symbol, const Interval& interval) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.is_differentiable(symbol, interval)) {
            return Tribool::True;
        }
    }
    return Tribool::Unknown;
}

Tribool AssumptionContext::is_positive_definite(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        Definiteness d = it->properties.get_definiteness(symbol);
        if (d != Definiteness::Unknown) {
            return (d == Definiteness::PositiveDefinite) ? Tribool::True : Tribool::False;
        }
    }
    return Tribool::Unknown;
}

Tribool AssumptionContext::is_positive_semidefinite(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        Definiteness d = it->properties.get_definiteness(symbol);
        if (d != Definiteness::Unknown) {
            return (d == Definiteness::PositiveDefinite ||
                    d == Definiteness::PositiveSemiDefinite) ? Tribool::True : Tribool::False;
        }
    }
    return Tribool::Unknown;
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

/// Parse an interval from string format like [0.000000, 1.000000]
Interval parse_interval(const std::string& s, int line_num) {
    if (s.size() < 5) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": malformed interval '" + s + "'");
    }

    Interval iv;
    bool lower_open = (s[0] == '(');
    bool upper_open = (s.back() == ')');
    std::string inner = s.substr(1, s.size() - 2);
    auto comma_pos = inner.find(", ");
    if (comma_pos == std::string::npos) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": malformed interval '" + s + "'");
    }

    std::string lo_str = inner.substr(0, comma_pos);
    std::string hi_str = inner.substr(comma_pos + 2);
    if (lo_str == "-inf") {
        iv.lower = Endpoint::neg_inf();
    } else {
        double lo_val = std::stod(lo_str);
        auto lo_expr = std::make_shared<SymbolicExpr>(
            std::make_shared<NumberNode>(static_cast<lmmc_real_t>(lo_val)));
        iv.lower = lower_open ? Endpoint::open(lo_expr) : Endpoint::closed(lo_expr);
    }
    if (hi_str == "+inf") {
        iv.upper = Endpoint::pos_inf();
    } else {
        double hi_val = std::stod(hi_str);
        auto hi_expr = std::make_shared<SymbolicExpr>(
            std::make_shared<NumberNode>(static_cast<lmmc_real_t>(hi_val)));
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

AssumptionContext AssumptionContext::deserialize(const std::string& data) {
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
        if (line == "END") {
            ended = true;
            break;
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
            if (idx == 0) {
                current_scope = 0;
            } else {
                ctx.push();
                current_scope = idx;
            }
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
            Domain dom = parse_domain(dom_str, line_num);
            ctx.current_properties().declare_domain(sym, dom);
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
            Sign s = parse_sign(sign_str, line_num);
            ctx.current_properties().declare_sign(sym, s);
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
            Parity p = parse_parity(par_str, line_num);
            ctx.current_properties().declare_parity(sym, p);
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
            Boundedness b = parse_boundedness(bnd_str, line_num);
            ctx.current_properties().declare_bounded(sym, b);
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
            ctx.current_properties().declare_transcendental(sym);
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
            ctx.current_properties().declare_finiteness(sym, f);
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
            ctx.current_properties().declare_definiteness(sym, d);
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
            std::shared_ptr<SymbolicExpr> period_expr;
            try {
                double period_val = std::stod(period_str);
                period_expr = std::make_shared<SymbolicExpr>(
                    std::make_shared<NumberNode>(static_cast<lmmc_real_t>(period_val)));
            } catch (...) {
                period_expr = std::make_shared<SymbolicExpr>(
                    std::make_shared<VariableNode>(period_str));
            }
            ctx.current_properties().declare_periodic(sym, period_expr);
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
                try {
                    double val = std::stod(s);
                    return SymbolicExpr(std::make_shared<NumberNode>(
                        static_cast<lmmc_real_t>(val)));
                } catch (...) {}
                return SymbolicExpr(std::make_shared<VariableNode>(s));
            };

            SymbolicExpr lhs = parse_simple_expr(lhs_str);
            SymbolicExpr rhs = parse_simple_expr(rhs_str);
            ctx.current_relations().add_relation(lhs, rhs, op, ctx.current_properties());
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
            // We parse simple "var OP val" patterns for condition and conclusion.
            // Complex expression deserialization is deferred per task spec.
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
                auto parse_token = [](const std::string& tok) -> std::shared_ptr<SymbolicNode> {
                    try {
                        double val = std::stod(tok);
                        return std::make_shared<NumberNode>(static_cast<lmmc_real_t>(val));
                    } catch (...) {}
                    return std::make_shared<VariableNode>(tok);
                };

                auto lhs_node = parse_token(lhs_s);
                auto rhs_node = parse_token(rhs_s);
                return SymbolicExpr(std::make_shared<RelationalNode>(lhs_node, rhs_node, rel_op));
            };

            SymbolicExpr cond_expr = parse_relational(cond_str);
            SymbolicExpr concl_expr = parse_relational(concl_str);
            ctx.assume_conditional(cond_expr, concl_expr);
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
            ctx.current_properties().declare_continuous(sym, iv);
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
            ctx.current_properties().declare_differentiable(sym, iv);
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
            ctx.current_properties().declare_monotonicity(sym, var, iv, m);
        } else if (keyword == "END") {
            ended = true;
            break;
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

} // namespace lamina
