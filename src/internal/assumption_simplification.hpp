#pragma once

#include "symbolic.hpp"

namespace LMCAS::detail {

ExpressionResult simplify_expression(const ExprPtr& expression,
                                     ComputationContext& context);

}
