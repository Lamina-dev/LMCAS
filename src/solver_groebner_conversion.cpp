#include "internal/solver_groebner_builder.hpp"
#include "internal/expression_analysis.hpp"

namespace LMCAS::groebner_detail {

        size_t PolyBuilder::get_or_create_aux_var(const std::shared_ptr<const SymbolicNode>& node) {

            auto tmp = LMCAS::detail::expression_from_node(node);
            std::string key = tmp.to_string();

            auto it = transcendental_map.find(key);
            if (it != transcendental_map.end()) {
                return it->second;
            }

            size_t idx = ext_vars.size();
            std::string aux_name = "__aux_" + std::to_string(idx) + "_";
            ext_vars.push_back(aux_name);
            transcendental_map[key] = idx;

            if (aux_to_node) {
                (*aux_to_node)[idx] = node;
            }
            return idx;
        }

        void PolyBuilder::represent_as_aux_or_fail(const std::shared_ptr<const SymbolicNode>& node) {
            if (strict_mode) {
                failed = true;
                return;
            }
            size_t idx = get_or_create_aux_var(node);
            result = Poly(ext_vars.size());
            Monomial m(ext_vars.size(), 0);
            m[idx] = 1;
            result.add_term(m, Rational(1));
        }

bool PolyBuilder::function_depends_on_variables(const FunctionNode& node) const {
    for (const auto& argument : node.arguments()) {
        for (const auto& variable : vars) {
            if (expression_depends_on_variable(argument, variable)) return true;
        }
    }
    return false;
}

bool PolyBuilder::use_numeric_function(const FunctionNode& node) {
    auto expression = LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<FunctionNode>(node.type(), node.arguments()));
    auto simplified = expression.simplify();
    if (!simplified || !simplified->is_number()) return false;
    auto number = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(simplified));
    if (!number) return false;
    if (std::holds_alternative<Rational>(number->value())) {
        result = Poly(ext_vars.size());
        result.add_term(Monomial(ext_vars.size(), 0), std::get<Rational>(number->value()));
        return true;
    }
    if (std::holds_alternative<BigInt>(number->value())) {
        result = Poly(ext_vars.size());
        result.add_term(Monomial(ext_vars.size(), 0), Rational(std::get<BigInt>(number->value())));
        return true;
    }
    return false;
}

void PolyBuilder::visit(const FunctionNode& node) {
    if (strict_mode) {
        failed = true;
        return;
    }
    if (!function_depends_on_variables(node) && use_numeric_function(node)) return;
    represent_as_aux_or_fail(
        LMCAS::detail::make_node<FunctionNode>(node.type(), node.arguments()));
}

        void PolyBuilder::visit(const UninterpretedFunctionNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const MatrixNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const RelationalNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const LogicalNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const PiecewiseNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const SummationNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const ProductNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const TransformNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const QuantifierNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const SetBuilderNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const ComplexNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const FiniteSetNode& node) { represent_as_aux_or_fail(node.clone()); }

        void PolyBuilder::visit(const IntervalNode& node) { represent_as_aux_or_fail(node.clone()); }

        void PolyBuilder::visit(const MembershipNode& node) { represent_as_aux_or_fail(node.clone()); }

        void PolyBuilder::visit(const QuantityNode& node) { represent_as_aux_or_fail(node.clone()); }

        void PolyBuilder::visit(const IntegralNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const LimitNode& node) {
            represent_as_aux_or_fail(node.clone());
        }

        void PolyBuilder::visit(const RootOfNode& node) {
            represent_as_aux_or_fail(node.clone());
        }


    Poly to_poly(const SymbolicExpr& expr, PolyContext& ctx) {
        PolyBuilder b(ctx.ext_vars, ctx);
        LMCAS::detail::node(expr)->accept(b);
        return b.get_result();
    }

}
