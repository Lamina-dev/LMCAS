#include "internal/numeric_evaluation_support.hpp"
#include "lmmc/numeric.h"
#include <cmath>
#include <optional>

namespace LMCAS::detail {
namespace {

std::optional<Result<ApproxReal>> trigonometric(FunctionNode::FuncType type, double x) {
    switch (type) {
        case FunctionNode::FuncType::Sin: {
            return numeric_approximation(std::sin(x));
        }
        case FunctionNode::FuncType::Cos: {
            return numeric_approximation(std::cos(x));
        }
        case FunctionNode::FuncType::Tan: {
            return numeric_approximation(std::tan(x));
        }
        case FunctionNode::FuncType::Cot: {
            if (std::sin(x) == 0.0) {
                return numeric_failure(CasErrc::DomainError, "cot is undefined");
            }
            return numeric_approximation(1.0 / std::tan(x));
        }
        case FunctionNode::FuncType::Sec: {
            if (std::cos(x) == 0.0) {
                return numeric_failure(CasErrc::DomainError, "sec is undefined");
            }
            return numeric_approximation(1.0 / std::cos(x));
        }
        case FunctionNode::FuncType::Csc: {
            if (std::sin(x) == 0.0) {
                return numeric_failure(CasErrc::DomainError, "csc is undefined");
            }
            return numeric_approximation(1.0 / std::sin(x));
        }
        default: {
            return std::nullopt;
        }
    }
}

std::optional<Result<ApproxReal>> inverse_trigonometric(FunctionNode::FuncType type, double x) {
    switch (type) {
        case FunctionNode::FuncType::ArcSin: {
            if (x < -1.0 || x > 1.0) {
                return numeric_failure(CasErrc::DomainError, "asin real domain is [-1, 1]");
            }
            return numeric_approximation(std::asin(x));
        }
        case FunctionNode::FuncType::ArcCos: {
            if (x < -1.0 || x > 1.0) {
                return numeric_failure(CasErrc::DomainError, "acos real domain is [-1, 1]");
            }
            return numeric_approximation(std::acos(x));
        }
        case FunctionNode::FuncType::ArcTan: {
            return numeric_approximation(std::atan(x));
        }
        default: {
            return std::nullopt;
        }
    }
}

std::optional<Result<ApproxReal>> hyperbolic(FunctionNode::FuncType type, double x) {
    switch (type) {
        case FunctionNode::FuncType::Sinh: {
            return numeric_approximation(std::sinh(x));
        }
        case FunctionNode::FuncType::Cosh: {
            return numeric_approximation(std::cosh(x));
        }
        case FunctionNode::FuncType::Tanh: {
            return numeric_approximation(std::tanh(x));
        }
        default: {
            return std::nullopt;
        }
    }
}

std::optional<Result<ApproxReal>> elementary(FunctionNode::FuncType type, double x) {
    switch (type) {
        case FunctionNode::FuncType::Ln:
        case FunctionNode::FuncType::Log: {
            if (x <= 0.0) {
                return numeric_failure(CasErrc::DomainError, "logarithm requires a positive real argument");
            }
            return numeric_approximation(std::log(x));
        }
        case FunctionNode::FuncType::Abs: {
            return numeric_approximation(std::abs(x));
        }
        case FunctionNode::FuncType::Sqrt: {
            if (x < 0.0) {
                return numeric_failure(CasErrc::DomainError, "real square root requires a non-negative argument");
            }
            return numeric_approximation(std::sqrt(x));
        }
        case FunctionNode::FuncType::Exp: {
            return numeric_approximation(std::exp(x));
        }
        case FunctionNode::FuncType::Sgn: {
            return numeric_approximation((x > 0.0) - (x < 0.0));
        }
        case FunctionNode::FuncType::Floor: {
            return numeric_approximation(std::floor(x));
        }
        case FunctionNode::FuncType::Ceil: {
            return numeric_approximation(std::ceil(x));
        }
        case FunctionNode::FuncType::Round: {
            return numeric_approximation(std::round(x));
        }
        default: {
            return std::nullopt;
        }
    }
}

std::optional<Result<ApproxReal>> special(FunctionNode::FuncType type, double x) {
    switch (type) {
        case FunctionNode::FuncType::Erf: {
            return numeric_approximation(std::erf(x));
        }
        case FunctionNode::FuncType::LambertW: {
            lmmc_real_t value = 0.0;
            const auto status = lmmc_lambertw(static_cast<lmmc_real_t>(x), &value);
            if (status == LMMC_STATUS_INVALID_ARGUMENT) {
                return numeric_failure(CasErrc::DomainError, "LambertW argument is outside the real branch");
            }
            if (status == LMMC_STATUS_CONVERGENCE_FAILED) {
                return numeric_failure(CasErrc::NumericFailure, "LambertW failed to certify convergence");
            }
            if (status == LMMC_STATUS_NUMERICAL_FAILURE) {
                return numeric_failure(CasErrc::NumericFailure, "LambertW encountered a nonfinite numerical value");
            }
            if (status != LMMC_STATUS_OK) {
                return numeric_failure(CasErrc::InternalInvariant, "Unexpected LambertW status");
            }
            return numeric_approximation(static_cast<double>(value));
        }
        case FunctionNode::FuncType::RealPart:
        case FunctionNode::FuncType::Conjugate: {
            return numeric_approximation(x);
        }
        case FunctionNode::FuncType::ComplexArg: {
            return numeric_approximation(std::atan2(0.0, x));
        }
        case FunctionNode::FuncType::ImagPart: {
            return numeric_approximation(0.0);
        }
        case FunctionNode::FuncType::ComplexAbs: {
            return numeric_approximation(std::abs(x));
        }
        default: {
            return std::nullopt;
        }
    }
}

}

Result<ApproxReal> evaluate_unary_numeric(FunctionNode::FuncType type, double x) {
    auto evaluated = trigonometric(type, x);
    if (!evaluated) {
        evaluated = inverse_trigonometric(type, x);
    }
    if (!evaluated) {
        evaluated = hyperbolic(type, x);
    }
    if (!evaluated) {
        evaluated = elementary(type, x);
    }
    if (!evaluated) {
        evaluated = special(type, x);
    }
    if (!evaluated) {
        return numeric_failure(CasErrc::UnsupportedExpression,
                               "function is not supported by real numeric evaluation");
    }
    return std::move(*evaluated);
}

}
