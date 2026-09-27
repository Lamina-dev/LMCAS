#include "internal/visitors/expand_visitor.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {

void ExpandVisitor::visit(const PowerNode& node) {
    std::size_t size = 1;
    auto expanded_base = expand_child(node.base(), size);
    auto expanded_exponent = expand_child(node.exponent(), size);
    const auto number = std::dynamic_pointer_cast<const NumberNode>(expanded_exponent);
    BigInt exponent;
    if (try_get_integer_value(number, exponent)) {
        if (exponent.is_zero() &&
            normalization_can_discard(expanded_base, detail::no_facts(), Domain::Real, context_) &&
            normalization_proved(detail::query_nonzero_value(expanded_base, detail::no_facts(), Domain::Complex, context_))) {
            set_result(detail::make_node<NumberNode>(BigInt(1)), 1);
            return;
        }
        const auto count = exponent.try_to_uint64();
        if (count && *count > 0 && *count < 20) {
            NormalizationVisitor normalizer(context_, detail::no_facts(), Domain::Real, rewrite_budget());
            set_result(expanded_base);
            for (std::uint64_t i = 1; i < *count; ++i) {
                set_result(normalizer.expand_product(result, expanded_base));
            }
            return;
        }
    }
    set_result(detail::make_node<PowerNode>(expanded_base, expanded_exponent), size);
}

}
