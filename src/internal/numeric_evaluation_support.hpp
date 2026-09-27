#pragma once

#include "numeric_evaluation.hpp"
#include "internal/symbolic_ast.hpp"
#include <array>

namespace LMCAS::detail {

Result<ApproxReal> numeric_failure(CasErrc code, std::string message);
Result<ApproxReal> numeric_approximation(double value);
Result<ApproxReal> evaluate_unary_numeric(FunctionNode::FuncType type, double value);

class NumericEvaluator {
public:
    NumericEvaluator(const NumericBindings& bindings, ComputationContext& context)
        : bindings_(bindings), context_(context) {}

    Result<ApproxReal> evaluate(const std::shared_ptr<const SymbolicNode>& node);

private:
    const NumericBindings& bindings_;
    ComputationContext& context_;

    Result<ApproxReal> number(const NumberNode& node);
    Result<ApproxReal> variable(const VariableNode& node);
    Result<ApproxReal> sum(const AddNode& node);
    Result<ApproxReal> product(const MultiplyNode& node);
    Result<ApproxReal> norm(const std::array<const PowerNode*, 2>& squares);
    Result<ApproxReal> power(const PowerNode& node);
    Result<ApproxReal> function(const FunctionNode& node);
    Result<ApproxReal> atan2(const FunctionNode& node);
    Result<ApproxReal> extremum(const FunctionNode& node);
    Result<ApproxReal> root(const std::shared_ptr<const SymbolicNode>& node);
};

}
