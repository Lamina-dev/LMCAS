#pragma once

#include "internal/solver_groebner_internal.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS::groebner_detail {

    class PolyBuilder : public LMCAS::detail::SymbolicVisitor {
        std::vector<std::string> vars;

        std::vector<std::string>& ext_vars;

        std::unordered_map<std::string, size_t>& transcendental_map;

        std::unordered_map<size_t, std::shared_ptr<const SymbolicNode>>* aux_to_node;
        Poly result;
        bool strict_mode;

        size_t get_or_create_aux_var(const std::shared_ptr<const SymbolicNode>& node);

        void represent_as_aux_or_fail(const std::shared_ptr<const SymbolicNode>& node);

        bool use_numeric_function(const FunctionNode& node);
        bool function_depends_on_variables(const FunctionNode& node) const;

    public:

        PolyBuilder(const std::vector<std::string>& v, PolyContext& ctx, bool strict = false);


        PolyBuilder(const std::vector<std::string>& v,
                    std::vector<std::string>& ext_v,
                    std::unordered_map<std::string, size_t>& trans_map,
                    std::unordered_map<size_t, std::shared_ptr<const SymbolicNode>>* aux_map,
                    bool strict = false)
            : vars(v), ext_vars(ext_v), transcendental_map(trans_map),
              aux_to_node(aux_map), result(ext_v.size()), strict_mode(strict) {}

        Poly get_result() const { return result; }
        bool failed = false;

        void visit(const NumberNode& node) override;

        void visit(const VariableNode& node) override;

        void visit(const AddNode& node) override;

        void visit(const MultiplyNode& node) override;

        void visit(const PowerNode& node) override;

        void visit(const FunctionNode& node) override;
        void visit(const UninterpretedFunctionNode& node) override;

        void visit(const MatrixNode& node) override;
        void visit(const RelationalNode& node) override;
        void visit(const LogicalNode& node) override;
        void visit(const PiecewiseNode& node) override;
        void visit(const SummationNode& node) override;
        void visit(const ProductNode& node) override;
        void visit(const TransformNode& node) override;
        void visit(const QuantifierNode& node) override;
        void visit(const SetBuilderNode& node) override;
        void visit(const ComplexNode& node) override;
        void visit(const FiniteSetNode& node) override;
        void visit(const IntervalNode& node) override;
        void visit(const MembershipNode& node) override;
        void visit(const QuantityNode& node) override;
        void visit(const IntegralNode& node) override;
        void visit(const LimitNode& node) override;
        void visit(const RootOfNode& node) override;
    };

}
