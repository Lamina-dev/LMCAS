#include "expr.hpp"

#include <string>

#include "internal/expr_internal.hpp"
#include "internal/expr_common.hpp"
namespace LMCAS {
namespace {

using namespace expr_detail::expr_common;

const char* unit_set_error_name(CasErrc code) noexcept {
    switch (code) {
    case CasErrc::DimensionMismatch: {
        return "DimensionMismatch";
    }
    case CasErrc::UnitInvalid: {
        return "UnitInvalid";
    }
    case CasErrc::UnitStripTypeMismatch: {
        return "UnitStripTypeMismatch";
    }
    case CasErrc::SetElementTypeMismatch: {
        return "SetElementTypeMismatch";
    }
    case CasErrc::SetOperandTypeMismatch: {
        return "SetOperandTypeMismatch";
    }
    case CasErrc::SetElementNotHashable: {
        return "SetElementNotHashable";
    }
    default: {
        return "InternalInvariant";
    }
    }
}

const char* equivalence_error_name(const CasError& error) noexcept {
    if (error.operation == kEquivalentOperation &&
        error.code == CasErrc::ResourceLimit) {
        return "EqvBudgetExceeded";
    }
    if (error.operation == kEquivalentProfileOperation &&
        error.code == CasErrc::UnsupportedExpression) {
        return "EqvRuleDisabled";
    }
    return error_name(error.code);
}

}

const char* error_name(CasErrc code) noexcept {
    switch (code) {
    case CasErrc::InvalidArgument: {
        return "InvalidArgument";
    }
    case CasErrc::ParseError: {
        return "ParseError";
    }
    case CasErrc::UnboundSymbol: {
        return "UnboundSymbol";
    }
    case CasErrc::DomainError: {
        return "DomainError";
    }
    case CasErrc::UnsupportedExpression: {
        return "UnsupportedExpression";
    }
    case CasErrc::Inconclusive: {
        return "Inconclusive";
    }
    case CasErrc::ResourceLimit: {
        return "ResourceLimit";
    }
    case CasErrc::Cancelled: {
        return "Cancelled";
    }
    case CasErrc::NumericFailure: {
        return "NumericFailure";
    }
    case CasErrc::InternalInvariant: {
        return "InternalInvariant";
    }
    default: {
        return unit_set_error_name(code);
    }
    }
}

const char* error_name(const CasError& error) noexcept {
    if (error.operation == kSymOperation &&
        error.code == CasErrc::InvalidArgument &&
        error.message.find("imaginary unit") != std::string::npos) {
        return "ImaginaryUnitReserved";
    }
    if (error.operation == kEvalComplexOperation &&
        error.code == CasErrc::UnboundSymbol) {
        return "ComplexEvalUnboundSymbol";
    }
    if (error.operation == kComplexOperation &&
        error.code == CasErrc::InvalidArgument) {
        return "ComplexTypeMismatch";
    }
    if (error.operation == kExprSetOperation &&
        error.code == CasErrc::InvalidArgument) {
        return "SetElementTypeMismatch";
    }
    if (error.operation == kSolveExprSetOperation &&
        error.code == CasErrc::Inconclusive) {
        return "SetResultInconclusive";
    }
    return equivalence_error_name(error);
}

}
