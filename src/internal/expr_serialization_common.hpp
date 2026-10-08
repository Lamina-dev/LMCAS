#pragma once

#include "computation_context.hpp"
#include "conditional_result.hpp"
#include "internal/symbolic_ast.hpp"

#include <array>
#include <string_view>

namespace LMCAS::serialization_detail {

inline constexpr std::string_view kHeader = "LMCAS_EXPR/1\n";
inline constexpr const char* kEncode = "LMCAS.serialize_expr";
inline constexpr const char* kDecode = "LMCAS.parse_serialized_expr";
using NodePtr = detail::SymbolicNodePtr;

// These names, not the C++ enum ordinals, are the wire values.
inline constexpr std::array<std::string_view, 36> kFunctions = {
    "sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "atan",
    "sinh", "cosh", "tanh", "ln", "log", "abs", "sqrt", "exp",
    "lambert_w", "atan2", "infinity", "erf", "ei", "si", "ci", "li",
    "max", "min", "sgn", "floor", "ceil", "round", "real", "imag",
    "conjugate", "complex_abs", "complex_arg"
};
inline constexpr std::array<std::string_view, 6> kRelations = {
    "eq", "neq", "lt", "gt", "leq", "geq"
};
inline constexpr std::array<std::string_view, 4> kLogic = {
    "and", "or", "not", "implies"
};
inline constexpr std::array<std::string_view, 5> kTransforms = {
    "laplace", "inverse_laplace", "fourier", "inverse_fourier", "z"
};
inline constexpr std::array<std::string_view, 2> kQuantifiers = {"forall", "exists"};
inline constexpr std::array<std::string_view, 3> kDirections = {
    "both", "below", "above"
};

template <typename Enum, std::size_t N>
std::string_view enum_name(Enum value, const std::array<std::string_view, N>& names,
                           const char* operation) {
    auto index = static_cast<std::size_t>(value);
    if (index >= names.size()) {
        throw CasError{CasErrc::UnsupportedExpression, "unsupported expression operator", operation};
    }
    return names[index];
}

inline void check(Result<void> result) {
    if (!result) throw result.error();
}

class Recursion {
public:
    Recursion(ComputationContext& context, const char* operation)
        : context_(context) {
        check(context_.enter_recursion(operation));
    }
    ~Recursion() { context_.leave_recursion(); }
    Recursion(const Recursion&) = delete;
    Recursion& operator=(const Recursion&) = delete;
private:
    ComputationContext& context_;
};

} // namespace LMCAS::serialization_detail
