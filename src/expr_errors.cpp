#include "expr.hpp"

#include <string>

#include "expr_internal.hpp"

namespace LMCAS {
namespace {

constexpr const char* kSymOperation = "LMCAS.sym";
constexpr const char* kComplexOperation = "LMCAS.complex";

} // namespace

const char* error_name(CasErrc code) noexcept {
    switch (code) {
    case CasErrc::InvalidArgument:
        return "InvalidArgument";
    case CasErrc::ParseError:
        return "ParseError";
    case CasErrc::UnboundSymbol:
        return "UnboundSymbol";
    case CasErrc::DomainError:
        return "DomainError";
    case CasErrc::UnsupportedExpression:
        return "UnsupportedExpression";
    case CasErrc::Inconclusive:
        return "Inconclusive";
    case CasErrc::ResourceLimit:
        return "ResourceLimit";
    case CasErrc::Cancelled:
        return "Cancelled";
    case CasErrc::NumericFailure:
        return "NumericFailure";
    case CasErrc::InternalInvariant:
        return "InternalInvariant";
    case CasErrc::DimensionMismatch:
        return "DimensionMismatch";
    case CasErrc::UnitInvalid:
        return "UnitInvalid";
    case CasErrc::UnitStripTypeMismatch:
        return "UnitStripTypeMismatch";
    case CasErrc::SetElementTypeMismatch:
        return "SetElementTypeMismatch";
    case CasErrc::SetOperandTypeMismatch:
        return "SetOperandTypeMismatch";
    case CasErrc::SetElementNotHashable:
        return "SetElementNotHashable";
    }
    return "InternalInvariant";
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

} // namespace LMCAS
