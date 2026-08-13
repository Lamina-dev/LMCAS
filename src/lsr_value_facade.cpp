#include "lsr_expr.hpp"

#include <utility>

namespace lamina::lsr {

ExprResult with_unit(const ExprPtr& value, const std::string& unit,
                     ComputationContext& context) {
    return attach_unit(value, unit, context);
}

ExprResult convert_to_unit(const ExprPtr& quantity, const std::string& unit,
                           ComputationContext& context) {
    return convert_unit(quantity, unit, context);
}

ExprResult strip_to_base_value(const ExprPtr& quantity,
                               ComputationContext& context) {
    return strip_unit(quantity, UnitStripMode::BaseValue, context);
}

ExprResult strip_to_display_value(const ExprPtr& quantity,
                                  ComputationContext& context) {
    return strip_unit(quantity, UnitStripMode::DisplayValue, context);
}

ExprResult finite_set(std::vector<ExprPtr> elements,
                      ComputationContext& context) {
    return make_finite_set(std::move(elements), context);
}

ExprResult interval(const ExprPtr& lower, const ExprPtr& upper,
                    bool lower_closed, bool upper_closed,
                    ComputationContext& context) {
    return make_interval(lower, upper, lower_closed, upper_closed, context);
}

ExprResult member(const ExprPtr& element, const ExprPtr& set,
                  ComputationContext& context) {
    return make_membership(element, set, context);
}

} // namespace lamina::lsr
