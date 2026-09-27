#include "internal/property_store_support.hpp"

namespace LMCAS::property_detail {

std::string domain_str(Domain domain) {
    static const std::pair<Domain, const char*> names[] = {
        {Domain::Complex, "Complex"},
        {Domain::Real, "Real"},
        {Domain::Algebraic, "Algebraic"},
        {Domain::Rational, "Rational"},
        {Domain::Integer, "Integer"},
        {Domain::Natural, "Natural"},
        {Domain::PositiveInt, "PositiveInt"}
    };
    for (const auto& name : names) {
        if (domain == name.first) return name.second;
    }
    return "Complex";
}
PropertyStoreResult invalid_empty_symbol(const std::string& operation) {
    return PropertyStoreResult::failure(
        CasErrc::InvalidArgument, "symbol name must not be empty", operation);
}

}
