#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/property_store_support.hpp"
#include <sstream>

namespace LMCAS {

namespace {

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
std::string endpoint_to_string(const Endpoint& ep) {
    if (ep.is_neg_infinity) return "-inf";
    if (ep.is_pos_infinity) return "+inf";
    if (ep.value) return ep.value->to_string();
    return "0";
}
std::string interval_to_string(const Interval& iv) {
    std::string result;
    result += iv.lower.is_open ? "(" : "[";
    result += endpoint_to_string(iv.lower);
    result += ", ";
    result += endpoint_to_string(iv.upper);
    result += iv.upper.is_open ? ")" : "]";
    return result;
}

void serialize_symbol_properties(std::ostream& out, const PropertyStore& props,
                                 const std::string& sym) {
    Domain dom = props.get_domain(sym);
    if (dom != Domain::Complex) {
        out << "DOMAIN " << sym << " " << property_detail::domain_str(dom) << "\n";
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
        out << "BOUNDED " << sym << " " << boundedness_to_string(bnd);
        if (bnd == Boundedness::Bounded) {
            const auto bounds = props.get_bounds(sym);
            if (bounds) {
                out << " " << interval_to_string(*bounds);
            }
        }
        out << "\n";
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
    for (const auto& [variable, period] : props.get_period_decls(sym)) {
        out << "PERIODIC " << sym << " " << variable << " "
            << period.to_string() << "\n";
    }
}

void serialize_symbol_intervals(std::ostream& out, const PropertyStore& props,
                                const std::string& sym) {
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

}

std::string AssumptionContext::serialize() const {
    std::ostringstream out;

    for (int scope_idx = 0; scope_idx < static_cast<int>(scope_stack_.size()); ++scope_idx) {
        out << "SCOPE " << scope_idx << "\n";
        const auto& scope = scope_stack_[scope_idx];
        const auto& props = scope.properties;
        auto symbols = props.get_all_symbols();
        for (const auto& sym : symbols) {
            serialize_symbol_properties(out, props, sym);
        }
        for (const auto& sym : symbols) {
            serialize_symbol_intervals(out, props, sym);
        }
        const auto& relations = scope.relations.get_relations();
        for (const auto i : scope.relations.declaration_order_) {
            const auto& rel = relations[i];
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

}
