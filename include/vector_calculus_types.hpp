#pragma once

#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <utility>

namespace LMCAS {

/** @brief 向量场由标量表达式组成。 */
using VectorField = std::vector<std::shared_ptr<SymbolicExpr>>;
using VectorCalculusExprResult = Result<std::shared_ptr<SymbolicExpr>>;
using VectorCalculusFieldResult = Result<VectorField>;

}
