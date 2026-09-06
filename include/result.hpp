#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <type_traits>
#include <variant>

namespace LMCAS {

enum class CasErrc {
    InvalidArgument,
    ParseError,
    UnboundSymbol,
    DomainError,
    DimensionMismatch,
    UnitInvalid,
    UnitStripTypeMismatch,
    SetElementTypeMismatch,
    SetOperandTypeMismatch,
    SetElementNotHashable,
    UnsupportedExpression,
    Inconclusive,
    ResourceLimit,
    Cancelled,
    NumericFailure,
    InternalInvariant
};

struct CasError {
    CasErrc code = CasErrc::InternalInvariant;
    std::string message;
    std::string operation;
};


template <typename T>
class [[nodiscard]] Result {
public:
    /** Allow direct return of a success value from a Result-returning function. */
    Result(T value) : storage_(std::move(value)) {}

    template <typename U = T,
              std::enable_if_t<
                  std::is_constructible_v<U, std::nullptr_t>, int> = 0>
    Result(std::nullptr_t) : storage_(T(nullptr)) {}
    static Result success(T value) {
        return Result(std::move(value));
    }

    static Result failure(CasErrc code, std::string message,
                          std::string operation = {}) {
        return Result(CasError{code, std::move(message), std::move(operation)});
    }

    static Result failure(CasError error) {
        return Result(std::move(error));
    }

    bool has_value() const noexcept {
        return std::holds_alternative<T>(storage_);
    }

    explicit operator bool() const noexcept { return has_value(); }

    T& value() {
        if (!has_value()) throw std::logic_error("Result does not contain a value");
        return std::get<T>(storage_);
    }

    const T& value() const {
        if (!has_value()) throw std::logic_error("Result does not contain a value");
        return std::get<T>(storage_);
    }

    CasError& error() {
        if (has_value()) throw std::logic_error("Result does not contain an error");
        return std::get<CasError>(storage_);
    }

    const CasError& error() const {
        if (has_value()) throw std::logic_error("Result does not contain an error");
        return std::get<CasError>(storage_);
    }

private:
    explicit Result(CasError error) : storage_(std::move(error)) {}

    std::variant<T, CasError> storage_;
};


template <>
class [[nodiscard]] Result<void> {
public:
    static Result success() { return Result(); }

    static Result failure(CasErrc code, std::string message,
                          std::string operation = {}) {
        return Result(CasError{code, std::move(message), std::move(operation)});
    }

    static Result failure(CasError error) { return Result(std::move(error)); }

    bool has_value() const noexcept { return !has_error_; }
    explicit operator bool() const noexcept { return has_value(); }

    const CasError& error() const {
        if (!has_error_) throw std::logic_error("Result does not contain an error");
        return error_;
    }

private:
    Result() = default;
    explicit Result(CasError error) : has_error_(true), error_(std::move(error)) {}

    bool has_error_ = false;
    CasError error_{};
};

enum class Completeness { Complete, Inconclusive };

template <typename T>
struct MathResult {
    T value;
    Completeness completeness = Completeness::Complete;
    std::string reason;
};

} // namespace LMCAS
