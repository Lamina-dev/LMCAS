#include "internal/symbolic_ast.hpp"

namespace LMCAS {
namespace {

template <typename Number>
int compare_values(const Number& left, const Number& right) {
    if (left < right) {
        return -1;
    }
    if (left == right) {
        return 0;
    }
    return 1;
}

int compare_exact_to_rational(const NumberNode& left, const Rational& right) {
    if (const auto* rational = std::get_if<Rational>(&left.value())) {
        return compare_values(*rational, right);
    }
    return compare_values(Rational(std::get<BigInt>(left.value())), right);
}

int compare_exact_nodes(const NumberNode& left, const NumberNode& right) {
    if (const auto* rational = std::get_if<Rational>(&right.value())) {
        return compare_exact_to_rational(left, *rational);
    }
    if (const auto* rational = std::get_if<Rational>(&left.value())) {
        return -compare_exact_to_rational(right, *rational);
    }
    return compare_values(std::get<BigInt>(left.value()), std::get<BigInt>(right.value()));
}

}

int NumberNode::compare_same_type(const SymbolicNode& other) const {
    const auto& right = static_cast<const NumberNode&>(other);
    return compare_exact_nodes(*this, right);
}

bool NumberNode::is_negative_one() const {
    if (const auto* integer = std::get_if<BigInt>(&value_)) {
        const auto narrowed = integer->try_to_int64();
        return narrowed && *narrowed == -1;
    }
    return std::get<Rational>(value_) == Rational(-1);
}

}
