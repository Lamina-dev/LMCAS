#include "lsr_expr.hpp"

#include <string>

#include "lsr_expr_internal.hpp"

namespace lamina::lsr {
namespace {

constexpr const char* kSymOperation = "lsr.sym";
constexpr const char* kIntegerOperation = "lsr.integer";
constexpr const char* kRationalOperation = "lsr.rational";
constexpr const char* kApproxOperation = "lsr.approx_real";
constexpr const char* kConstantOperation = "lsr.constant";
constexpr const char* kComplexOperation = "lsr.complex";
constexpr const char* kExprOperation = "lsr.expr_op";
constexpr const char* kMathOperation = "lsr.math";
constexpr const char* kRealOperation = "lsr.real";
constexpr const char* kImagOperation = "lsr.imag";
constexpr const char* kConjOperation = "lsr.conj";
constexpr const char* kAbsOperation = "lsr.abs";
constexpr const char* kSimplifyOperation = "lsr.simplify";
constexpr const char* kExpandOperation = "lsr.expand";
constexpr const char* kDifferentiateOperation = "lsr.differentiate";
constexpr const char* kSubstituteOperation = "lsr.substitute";
constexpr const char* kExprMatchOperation = "lsr.expr_match";

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

} // namespace lamina::lsr
