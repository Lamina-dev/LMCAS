#pragma once

#include "symbolic.hpp"
#include <string>
#include <vector>

namespace LMCAS::detail {

ExpressionResult simplify_expression(const ExprPtr& expression,
                                     ComputationContext& context);

ExpressionResult expand_expression(const ExprPtr& expression,
                                   ComputationContext& context);
ExpressionResult differentiate_expression(const ExprPtr& expression,
                                          const std::string& variable,
                                          ComputationContext& context);
ExpressionResult substitute_expression(const ExprPtr& expression,
                                       const std::string& variable,
                                       const ExprPtr& replacement,
                                       ComputationContext& context);
ExpressionResult substitute_expression(
    const ExprPtr& expression,
    const std::vector<Binding>& replacements,
    ComputationContext& context);

}
