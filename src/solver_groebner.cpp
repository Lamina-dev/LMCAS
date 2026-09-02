#include "solver.hpp"
#include "solve_strategies.hpp"
#include "symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "assumption_context.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <queue>
#include <unordered_map>
#include <optional>
#include <limits>
#include "internal/solver_support.hpp"

namespace lamina {
using namespace solver_detail;
namespace {

    using Monomial = std::vector<int>;

    struct MonomialLess {
        bool operator()(const Monomial& a, const Monomial& b) const {

            for (size_t i = 0; i < a.size(); ++i) {
                if (a[i] != b[i]) {
                    return a[i] > b[i];
                }
            }
            return false;
        }
    };

    struct Term {
        Rational coeff;
        Monomial mono;
    };

    using PolyTerms = std::map<Monomial, Rational, MonomialLess>;

    struct Poly {
        PolyTerms terms;
        size_t num_vars;
        int sugar = 0;

        Poly() : num_vars(0) {}
        Poly(size_t n) : num_vars(n) {}

        bool is_zero() const { return terms.empty(); }

        Monomial LM() const {
            if (terms.empty()) return std::vector<int>(num_vars, 0);
            return terms.begin()->first;
        }

        Rational LC() const {
            if (terms.empty()) return Rational(0);
            return terms.begin()->second;
        }

        Term LT() const {
             if (terms.empty()) return {Rational(0), std::vector<int>(num_vars, 0)};
             return {terms.begin()->second, terms.begin()->first};
        }

        void add_term(const Monomial& m, const Rational& c) {
            if (c == Rational(0)) return;
            auto it = terms.find(m);
            if (it != terms.end()) {

                Rational new_c = it->second + c;
                if (new_c == Rational(0)) {
                    terms.erase(it);
                } else {
                    it->second = new_c;
                }
            } else {
                terms[m] = c;
            }
        }
    };

    Monomial mul_mono(const Monomial& a, const Monomial& b) {
        Monomial res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = a[i] + b[i];
        return res;
    }

    Monomial lcm_mono(const Monomial& a, const Monomial& b) {
        Monomial res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = std::max(a[i], b[i]);
        return res;
    }

    bool divides_mono(const Monomial& a, const Monomial& b) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] > b[i]) return false;
        }
        return true;
    }

    Monomial div_mono(const Monomial& num, const Monomial& den) {
        Monomial res(num.size());
        for (size_t i = 0; i < num.size(); ++i) res[i] = num[i] - den[i];
        return res;
    }

    Poly add_poly(const Poly& a, const Poly& b) {
        Poly res = a;
        for (auto const& [m, c] : b.terms) {
            res.add_term(m, c);
        }
        return res;
    }

    Poly sub_poly(const Poly& a, const Poly& b) {
        Poly res = a;
        for (auto const& [m, c] : b.terms) {

            Rational c_neg = c * Rational(-1);
            res.add_term(m, c_neg);
        }
        return res;
    }

    Poly mul_poly_term(const Poly& p, const Term& t) {
        Poly res(p.num_vars);
        if (t.coeff == Rational(0)) return res;
        for (auto const& [m, c] : p.terms) {
            res.add_term(mul_mono(m, t.mono), c * t.coeff);
        }
        return res;
    }

    Poly mul_poly(const Poly& a, const Poly& b) {
        Poly res(a.num_vars);
        for (auto const& [ma, ca] : a.terms) {
            for (auto const& [mb, cb] : b.terms) {
                res.add_term(mul_mono(ma, mb), ca * cb);
            }
        }
        return res;
    }

    struct PolyContext;

    class PolyBuilder : public lamina::detail::SymbolicVisitor {
        std::vector<std::string> vars;

        std::vector<std::string>& ext_vars;

        std::unordered_map<std::string, size_t>& transcendental_map;

        std::unordered_map<size_t, std::shared_ptr<const SymbolicNode>>* aux_to_node;
        Poly result;
        bool strict_mode;

        size_t get_or_create_aux_var(const std::shared_ptr<const SymbolicNode>& node) {

            auto tmp = lamina::detail::expression_from_node(node);
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

        void represent_as_aux_or_fail(const std::shared_ptr<const SymbolicNode>& node) {
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

    public:

        PolyBuilder(const std::vector<std::string>& v, PolyContext& ctx, bool strict = false);

        PolyBuilder(const std::vector<std::string>& v,
                    std::vector<std::string>& ext_v,
                    std::unordered_map<std::string, size_t>& trans_map,
                    bool strict = false)
            : vars(v), ext_vars(ext_v), transcendental_map(trans_map),
              aux_to_node(nullptr), result(ext_v.size()), strict_mode(strict) {}

        PolyBuilder(const std::vector<std::string>& v,
                    std::vector<std::string>& ext_v,
                    std::unordered_map<std::string, size_t>& trans_map,
                    std::unordered_map<size_t, std::shared_ptr<const SymbolicNode>>* aux_map,
                    bool strict = false)
            : vars(v), ext_vars(ext_v), transcendental_map(trans_map),
              aux_to_node(aux_map), result(ext_v.size()), strict_mode(strict) {}

        Poly get_result() const { return result; }
        bool failed = false;

        void visit(const NumberNode& node) override {
            result = Poly(ext_vars.size());

            if (std::holds_alternative<Rational>(node.value())) {
                result.add_term(std::vector<int>(ext_vars.size(), 0), std::get<Rational>(node.value()));
            } else if (std::holds_alternative<BigInt>(node.value())) {
                result.add_term(std::vector<int>(ext_vars.size(), 0), Rational(std::get<BigInt>(node.value())));
            } else if (std::holds_alternative<lmmc_real_t>(node.value())) {
                result.add_term(std::vector<int>(ext_vars.size(), 0), Rational((long long)std::get<lmmc_real_t>(node.value())));
            }
        }

        void visit(const VariableNode& node) override {
            result = Poly(ext_vars.size());

            auto it = std::find(ext_vars.begin(), ext_vars.end(), node.name());
            if (it != ext_vars.end()) {
                Monomial m(ext_vars.size(), 0);
                m[std::distance(ext_vars.begin(), it)] = 1;
                result.add_term(m, Rational(1));
            } else {

                if (strict_mode) {
                    failed = true;
                    return;
                }
                size_t idx = get_or_create_aux_var(lamina::detail::make_node<VariableNode>(node.name()));

                result = Poly(ext_vars.size());
                Monomial m(ext_vars.size(), 0);
                m[idx] = 1;
                result.add_term(m, Rational(1));
            }
        }

        void visit(const AddNode& node) override {
            Poly sum(ext_vars.size());
            for (auto& op : node.operands()) {
                PolyBuilder b(vars, ext_vars, transcendental_map, aux_to_node, strict_mode);
                op->accept(b);
                if (b.failed) { failed = true; return; }

                Poly b_res = b.get_result();
                if (b_res.num_vars > sum.num_vars) {

                    Poly new_sum(b_res.num_vars);
                    for (auto& [m, c] : sum.terms) {
                        Monomial extended = m;
                        extended.resize(b_res.num_vars, 0);
                        new_sum.add_term(extended, c);
                    }
                    sum = new_sum;
                } else if (sum.num_vars > b_res.num_vars) {
                    Poly new_b(sum.num_vars);
                    for (auto& [m, c] : b_res.terms) {
                        Monomial extended = m;
                        extended.resize(sum.num_vars, 0);
                        new_b.add_term(extended, c);
                    }
                    b_res = new_b;
                }
                sum = add_poly(sum, b_res);
            }
            result = sum;
        }

        void visit(const MultiplyNode& node) override {
            Poly prod(ext_vars.size());
            prod.add_term(std::vector<int>(ext_vars.size(), 0), Rational(1));

            for (auto& op : node.operands()) {
                PolyBuilder b(vars, ext_vars, transcendental_map, aux_to_node, strict_mode);
                op->accept(b);
                if (b.failed) { failed = true; return; }
                Poly b_res = b.get_result();

                if (b_res.num_vars > prod.num_vars) {
                    Poly new_prod(b_res.num_vars);
                    for (auto& [m, c] : prod.terms) {
                        Monomial extended = m;
                        extended.resize(b_res.num_vars, 0);
                        new_prod.add_term(extended, c);
                    }
                    prod = new_prod;
                } else if (prod.num_vars > b_res.num_vars) {
                    Poly new_b(prod.num_vars);
                    for (auto& [m, c] : b_res.terms) {
                        Monomial extended = m;
                        extended.resize(prod.num_vars, 0);
                        new_b.add_term(extended, c);
                    }
                    b_res = new_b;
                }
                prod = mul_poly(prod, b_res);
            }
            result = prod;
        }

        void visit(const PowerNode& node) override {
            PolyBuilder b_base(vars, ext_vars, transcendental_map, aux_to_node, strict_mode);
            node.base()->accept(b_base);
            if (b_base.failed) { failed = true; return; }
            Poly base = b_base.get_result();

            long long exp = 0;
            bool is_integer_exp = false;
            if (node.exponent()->is_number()) {
                auto num_node = std::dynamic_pointer_cast<const NumberNode>(node.exponent());
                if (num_node) {
                    const auto& val = num_node->value();
                    if (std::holds_alternative<BigInt>(val)) {
                        exp = std::get<BigInt>(val).to_int();
                        is_integer_exp = true;
                    } else if (std::holds_alternative<Rational>(val)) {
                        const auto& r = std::get<Rational>(val);
                        if (r.is_integer()) {
                            exp = r.to_BigInt().to_int();
                            is_integer_exp = true;
                        }
                    } else if (std::holds_alternative<lmmc_real_t>(val)) {
                        lmmc_real_t v = std::get<lmmc_real_t>(val);
                        lmmc_real_t rounded = std::round(v);
                        int eq;
                        lmmc_double_nearly_equal_tol(v, rounded, 1e-12, 1e-12, &eq);
                        if (eq) {
                            exp = static_cast<long long>(rounded);
                            is_integer_exp = true;
                        }
                    }
                }
            }

            if (!is_integer_exp) {

                if (strict_mode) { failed = true; return; }
                size_t idx = get_or_create_aux_var(
                    lamina::detail::make_node<PowerNode>(node.base(), node.exponent()));
                result = Poly(ext_vars.size());
                Monomial m(ext_vars.size(), 0);
                m[idx] = 1;
                result.add_term(m, Rational(1));
                return;
            }

            if (exp == 0) {
                result = Poly(ext_vars.size());
                result.add_term(std::vector<int>(ext_vars.size(), 0), Rational(1));
            } else if (exp > 0) {
                Poly res(ext_vars.size());
                res.add_term(std::vector<int>(ext_vars.size(), 0), Rational(1));
                for (long long i = 0; i < exp; ++i) {

                    if (base.num_vars > res.num_vars) {
                        Poly new_res(base.num_vars);
                        for (auto& [m, c] : res.terms) {
                            Monomial extended = m;
                            extended.resize(base.num_vars, 0);
                            new_res.add_term(extended, c);
                        }
                        res = new_res;
                    }
                    res = mul_poly(res, base);
                }
                result = res;
            } else {

                if (strict_mode) { failed = true; return; }
                size_t idx = get_or_create_aux_var(
                    lamina::detail::make_node<PowerNode>(node.base(), node.exponent()));
                result = Poly(ext_vars.size());
                Monomial m(ext_vars.size(), 0);
                m[idx] = 1;
                result.add_term(m, Rational(1));
            }
        }

        void visit(const FunctionNode& node) override {

            if (strict_mode) { failed = true; return; }

            bool expression_depends_on_variables = false;
            for (const auto& arg : node.arguments()) {
                for (const auto& v : vars) {
                    if (expression_depends_on_variable(arg, v)) {
                        expression_depends_on_variables = true;
                        break;
                    }
                }
                if (expression_depends_on_variables) break;
            }

            if (!expression_depends_on_variables) {

                auto func_expr = lamina::detail::expression_from_node(lamina::detail::make_node<FunctionNode>(node.type(), node.arguments()));
                auto simplified = func_expr.simplify();
                if (simplified && simplified->is_number()) {
                    auto nn = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified));
                    if (nn) {

                        if (std::holds_alternative<Rational>(nn->value())) {
                            result = Poly(ext_vars.size());
                            result.add_term(std::vector<int>(ext_vars.size(), 0),
                                            std::get<Rational>(nn->value()));
                            return;
                        }
                        if (std::holds_alternative<BigInt>(nn->value())) {
                            result = Poly(ext_vars.size());
                            result.add_term(std::vector<int>(ext_vars.size(), 0),
                                            Rational(std::get<BigInt>(nn->value())));
                            return;
                        }

                    }
                }
            }

            size_t idx = get_or_create_aux_var(
                lamina::detail::make_node<FunctionNode>(node.type(), node.arguments()));
            result = Poly(ext_vars.size());
            Monomial m(ext_vars.size(), 0);
            m[idx] = 1;
            result.add_term(m, Rational(1));
        }
        void visit(const UninterpretedFunctionNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }

        void visit(const MatrixNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const RelationalNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const LogicalNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const PiecewiseNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const SummationNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const ProductNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const TransformNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const QuantifierNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const SetBuilderNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const ComplexNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const FiniteSetNode& node) override { represent_as_aux_or_fail(node.clone()); }
        void visit(const IntervalNode& node) override { represent_as_aux_or_fail(node.clone()); }
        void visit(const MembershipNode& node) override { represent_as_aux_or_fail(node.clone()); }
        void visit(const QuantityNode& node) override { represent_as_aux_or_fail(node.clone()); }
        void visit(const IntegralNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const LimitNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
        void visit(const RootOfNode& node) override {
            represent_as_aux_or_fail(node.clone());
        }
    };

    struct PolyContext {
        std::vector<std::string> ext_vars;
        std::unordered_map<std::string, size_t> transcendental_map;
        std::unordered_map<size_t, std::shared_ptr<const SymbolicNode>> aux_to_node;
        size_t num_original_vars;

        PolyContext(const std::vector<std::string>& vars)
            : ext_vars(vars), num_original_vars(vars.size()) {}
    };

    inline PolyBuilder::PolyBuilder(const std::vector<std::string>& v, PolyContext& ctx, bool strict)
        : vars(v), ext_vars(ctx.ext_vars), transcendental_map(ctx.transcendental_map),
          aux_to_node(&ctx.aux_to_node), result(ctx.ext_vars.size()), strict_mode(strict) {}

    Poly to_poly(const SymbolicExpr& expr, PolyContext& ctx) {
        PolyBuilder b(ctx.ext_vars, ctx);
        lamina::detail::node(expr)->accept(b);
        return b.get_result();
    }

    [[maybe_unused]] Poly to_poly(const SymbolicExpr& expr, const std::vector<std::string>& vars) {

        std::vector<std::string> ext_vars = vars;
        std::unordered_map<std::string, size_t> trans_map;
        PolyBuilder b(vars, ext_vars, trans_map);
        lamina::detail::node(expr)->accept(b);
        return b.get_result();
    }

    [[maybe_unused]] SymbolicExpr from_poly(const Poly& p, const std::vector<std::string>& vars) {
        if (p.terms.empty()) return lamina::detail::expression_from_node(SymbolicFactory::create_number(BigInt(0)));

        std::vector<std::shared_ptr<const SymbolicNode>> add_ops;

        for (auto const& [m, c] : p.terms) {

            std::vector<std::shared_ptr<const SymbolicNode>> mul_ops;

            if (c.get_denominator() == BigInt(1)) {
                mul_ops.push_back(SymbolicFactory::create_number(c.get_numerator()));
            } else {
                mul_ops.push_back(SymbolicFactory::create_number(c));
            }

            for (size_t i = 0; i < m.size() && i < vars.size(); ++i) {
                if (m[i] > 0) {
                    auto var = SymbolicFactory::create_variable(vars[i]);
                    if (m[i] == 1) {
                         mul_ops.push_back(var);
                    } else {
                         auto pow = SymbolicFactory::create_power(var, SymbolicFactory::create_number(BigInt(m[i])));
                         mul_ops.push_back(pow);
                    }
                }
            }

            if (mul_ops.size() == 0) {

            } else if (mul_ops.size() == 1) {
                add_ops.push_back(mul_ops[0]);
            } else {
                add_ops.push_back(SymbolicFactory::create_multiply(mul_ops));
            }
        }

        if (add_ops.empty()) return lamina::detail::expression_from_node(SymbolicFactory::create_number(BigInt(0)));
        if (add_ops.size() == 1) return lamina::detail::expression_from_node(add_ops[0]);
        return lamina::detail::expression_from_node(SymbolicFactory::create_add(add_ops));
    }

    SymbolicExpr from_poly_ext(const Poly& p, const PolyContext& ctx,
                               const std::vector<std::string>&) {
        if (p.terms.empty()) return lamina::detail::expression_from_node(SymbolicFactory::create_number(BigInt(0)));

        std::vector<std::shared_ptr<const SymbolicNode>> add_ops;

        for (auto const& [m, c] : p.terms) {
            std::vector<std::shared_ptr<const SymbolicNode>> mul_ops;

            if (c.get_denominator() == BigInt(1)) {
                mul_ops.push_back(SymbolicFactory::create_number(c.get_numerator()));
            } else {
                mul_ops.push_back(SymbolicFactory::create_number(c));
            }

            for (size_t i = 0; i < m.size() && i < ctx.ext_vars.size(); ++i) {
                if (m[i] > 0) {

                    std::shared_ptr<const SymbolicNode> var_node;
                    if (i < ctx.num_original_vars) {
                        var_node = SymbolicFactory::create_variable(ctx.ext_vars[i]);
                    } else {
                        auto it = ctx.aux_to_node.find(i);
                        if (it != ctx.aux_to_node.end()) {

                            var_node = it->second;
                        } else {

                            var_node = SymbolicFactory::create_variable(ctx.ext_vars[i]);
                        }
                    }

                    if (m[i] == 1) {
                        mul_ops.push_back(var_node);
                    } else {
                        auto pow = SymbolicFactory::create_power(
                            var_node,
                            SymbolicFactory::create_number(BigInt(m[i])));
                        mul_ops.push_back(pow);
                    }
                }
            }

            if (mul_ops.empty()) {
            } else if (mul_ops.size() == 1) {
                add_ops.push_back(mul_ops[0]);
            } else {
                add_ops.push_back(SymbolicFactory::create_multiply(mul_ops));
            }
        }

        if (add_ops.empty()) return lamina::detail::expression_from_node(SymbolicFactory::create_number(BigInt(0)));
        if (add_ops.size() == 1) return lamina::detail::expression_from_node(add_ops[0]);
        return lamina::detail::expression_from_node(SymbolicFactory::create_add(add_ops));
    }

    Poly reduce(Poly p, const std::vector<Poly>& G) {
        Poly r(p.num_vars);

        while (!p.is_zero()) {
            bool reduced = false;
            Monomial p_lm = p.LM();

            for (const auto& g : G) {
                if (divides_mono(g.LM(), p_lm)) {

                    Term factor;
                    factor.coeff = p.LC() / g.LC();
                    factor.mono = div_mono(p_lm, g.LM());

                    Poly term_g = mul_poly_term(g, factor);
                    p = sub_poly(p, term_g);
                    reduced = true;
                    break;
                }
            }

            if (!reduced) {
                Term lt = p.LT();
                r.add_term(lt.mono, lt.coeff);

                p.terms.erase(p.terms.begin());
            }
        }
        return r;
    }

    Poly s_poly(const Poly& f, const Poly& g) {
        Monomial L = lcm_mono(f.LM(), g.LM());

        Term t1;
        t1.mono = div_mono(L, f.LM());
        t1.coeff = Rational(1) / f.LC();

        Term t2;
        t2.mono = div_mono(L, g.LM());
        t2.coeff = Rational(1) / g.LC();

        return sub_poly(mul_poly_term(f, t1), mul_poly_term(g, t2));
    }

    bool coprime_leading_monomials(const Poly& f, const Poly& g) {
        if (f.is_zero() || g.is_zero()) return false;
        const Monomial& lm_f = f.LM();
        const Monomial& lm_g = g.LM();
        for (size_t i = 0; i < lm_f.size() && i < lm_g.size(); ++i) {
            if (lm_f[i] > 0 && lm_g[i] > 0) return false;
        }
        return true;
    }

    bool chain_criterion(const std::vector<Poly>& G, size_t i, size_t j,
                         const std::set<std::pair<size_t,size_t>>& processed_pairs) {
        Monomial L = lcm_mono(G[i].LM(), G[j].LM());
        for (size_t k = 0; k < G.size(); ++k) {
            if (k == i || k == j) continue;
            if (G[k].is_zero()) continue;
            if (divides_mono(G[k].LM(), L)) {
                auto pair_ik = (i < k) ? std::make_pair(i, k) : std::make_pair(k, i);
                auto pair_kj = (k < j) ? std::make_pair(k, j) : std::make_pair(j, k);
                if (processed_pairs.count(pair_ik) && processed_pairs.count(pair_kj)) {
                    return true;
                }
            }
        }
        return false;
    }

    struct SugarPair {
        size_t i, j;
        int sugar_degree;

        bool operator>(const SugarPair& other) const {
            return sugar_degree > other.sugar_degree;
        }
    };

    int compute_spoly_sugar(const Poly& f, const Poly& g) {
        Monomial L = lcm_mono(f.LM(), g.LM());

        int deg_L = 0;
        for (int e : L) deg_L += e;

        int deg_lm_f = 0;
        for (int e : f.LM()) deg_lm_f += e;

        int deg_lm_g = 0;
        for (int e : g.LM()) deg_lm_g += e;

        int sugar_f = f.sugar + deg_L - deg_lm_f;
        int sugar_g = g.sugar + deg_L - deg_lm_g;

        return std::max(sugar_f, sugar_g);
    }

}

std::vector<SymbolicExpr> Solver::groebner_basis(
    const std::vector<SymbolicExpr>& polynomials,
    const std::vector<std::string>& variables)
{

    PolyContext ctx(variables);

    std::vector<Poly> G;
    for (const auto& expr : polynomials) {
        Poly p = to_poly(expr, ctx);
        if (!p.is_zero()) G.push_back(p);
    }

    size_t max_vars = ctx.ext_vars.size();
    for (auto& p : G) {
        if (p.num_vars < max_vars) {
            Poly resized(max_vars);
            for (auto& [m, c] : p.terms) {
                Monomial extended = m;
                extended.resize(max_vars, 0);
                resized.add_term(extended, c);
            }
            p = resized;
        }
    }

    for (auto& p : G) {
        if (!p.is_zero()) {
            int max_deg = 0;
            for (const auto& [m, c] : p.terms) {
                int d = 0;
                for (int e : m) d += e;
                max_deg = std::max(max_deg, d);
            }
            p.sugar = max_deg;
        }
    }

    std::priority_queue<SugarPair, std::vector<SugarPair>, std::greater<SugarPair>> pairs;
    for (size_t i = 0; i < G.size(); ++i) {
        for (size_t j = i + 1; j < G.size(); ++j) {
            int sugar = compute_spoly_sugar(G[i], G[j]);
            pairs.push({i, j, sugar});
        }
    }

    std::set<std::pair<size_t, size_t>> processed_pairs;

    while (!pairs.empty()) {
        auto [i, j, sugar_deg] = pairs.top();
        pairs.pop();

        auto canonical_pair = (i < j) ? std::make_pair(i, j) : std::make_pair(j, i);
        processed_pairs.insert(canonical_pair);

        if (coprime_leading_monomials(G[i], G[j])) {
            continue;
        }

        if (chain_criterion(G, i, j, processed_pairs)) {
            continue;
        }

        Poly S = s_poly(G[i], G[j]);

        if (S.num_vars < max_vars) {
            Poly resized(max_vars);
            for (auto& [m, c] : S.terms) {
                Monomial extended = m;
                extended.resize(max_vars, 0);
                resized.add_term(extended, c);
            }
            S = resized;
        }

        Poly r = reduce(S, G);

        if (!r.is_zero()) {

            r.sugar = sugar_deg;
            size_t new_idx = G.size();
            G.push_back(r);
            for (size_t k = 0; k < new_idx; ++k) {
                int new_sugar = compute_spoly_sugar(G[k], G[new_idx]);
                pairs.push({k, new_idx, new_sugar});
            }
        }
    }

    std::vector<SymbolicExpr> result;
    for (const auto& p : G) {
        result.push_back(from_poly_ext(p, ctx, variables));
    }
    return result;
}
static std::vector<std::map<std::string, SymbolicExpr>>
solve_polynomial_system_impl(
    const std::vector<SymbolicExpr>& equations,
    const std::vector<std::string>& variables,
    ComputationContext& context)
{
    std::vector<SymbolicExpr> cleared_equations;
    std::vector<std::shared_ptr<SymbolicExpr>> denom_constraints;
    cleared_equations.reserve(equations.size());
    for (const auto& eq : equations) {
        if (!lamina::detail::node(eq)) return {};
        std::vector<std::shared_ptr<const SymbolicNode>> den_factors;
        std::vector<std::shared_ptr<SymbolicExpr>> den_local_constraints;
        if (!collect_denominator_factors(lamina::detail::node(eq), den_factors, den_local_constraints)) return {};
        auto denom_expr = multiply_factors(den_factors);
        auto cleared = to_ptr(eq);
        if (!den_factors.empty()) {
            if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(cleared))) {
                std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
                new_ops.reserve(add->operands().size());
                for (const auto& op : add->operands()) {
                    auto prod = multiply_no_expand(op, den_factors);
                    new_ops.push_back(lamina::detail::node(prod));
                }
                cleared = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(new_ops));
            } else {
                cleared = multiply_no_expand(lamina::detail::node(cleared), den_factors);
            }
        } else {
            cleared = cleared->simplify();
        }

        if (!cleared || !lamina::detail::node(cleared) || !is_polynomial_node(lamina::detail::node(cleared))) return {};
        cleared_equations.push_back(*cleared);

        for (const auto& c : den_local_constraints) {
            denom_constraints.push_back(c);
        }
    }
    if (cleared_equations.size() == 1 && variables.size() == 1) {
        auto roots = detail::propagate_result(solve_finite_checked(
            detail::make_expression_ptr(cleared_equations[0]),
            variables[0], context, SolveOptions{}));
        std::vector<std::map<std::string, SymbolicExpr>> single_solutions;
        for (const auto& r : roots) {
            single_solutions.push_back({{variables[0], *r}});
        }
        if (denom_constraints.empty()) return single_solutions;

        std::vector<std::map<std::string, SymbolicExpr>> filtered;
        filtered.reserve(single_solutions.size());
        for (const auto& sol : single_solutions) {
            bool ok = true;
            for (const auto& den : denom_constraints) {
                auto sub = den;
                for (const auto& [name, val] : sol) {
                    sub = sub->substitute(name, lamina::detail::make_expression_ptr(val));
                    if (!sub) break;
                }
                if (!sub) continue;
                sub = sub->simplify();
                if (sub && sub->is_zero()) {
                    ok = false;
                    break;
                }
            }
            if (ok) filtered.push_back(sol);
        }
        return filtered;
    }

    auto G_basis = Solver::groebner_basis(cleared_equations, variables);
    std::vector<std::shared_ptr<SymbolicExpr>> basis;
    basis.reserve(G_basis.size());
    for (const auto& g : G_basis) {
        auto g_ptr = lamina::detail::make_expression_ptr(g);
        auto simp = g_ptr->simplify();
        if (simp && !simp->is_zero()) {
            basis.push_back(simp);
        }
    }
    if (variables.empty()) {
        for (const auto& p : basis) {
            if (p->is_number() && !p->is_zero()) return {};
        }
        return {{} };
    }
    auto substitute_all = [&](const std::shared_ptr<SymbolicExpr>& expr,
                              const std::map<std::string, SymbolicExpr>& subs) {
        auto res = expr;
        for (const auto& [name, val] : subs) {
            res = res->substitute(name, lamina::detail::make_expression_ptr(val));
            if (!res) return std::shared_ptr<SymbolicExpr>(nullptr);
        }
        return res->simplify();
    };
    auto solve_rec = [&](auto&& self, int var_pos,
                         const std::map<std::string, SymbolicExpr>& partial)
        -> std::vector<std::map<std::string, SymbolicExpr>> {
        std::vector<std::shared_ptr<SymbolicExpr>> reduced;
        reduced.reserve(basis.size());

        for (const auto& p : basis) {
            auto r = substitute_all(p, partial);
            if (!r) continue;
            if (r->is_zero()) continue;

            bool depends = false;
            for (int i = 0; i <= var_pos && i < (int)variables.size(); ++i) {
                if (contains(*r, variables[i])) {
                    depends = true;
                    break;
                }
            }

            if (!depends) {
                if (r->is_number() && !r->is_zero()) return std::vector<std::map<std::string, SymbolicExpr>>{};
                continue;
            }

            reduced.push_back(r);
        }
        if (var_pos < 0) {
            return {partial};
        }
        const auto& curr_var = variables[var_pos];
        bool curr_var_appears = false;
        std::shared_ptr<SymbolicExpr> target = nullptr;
        int best_deg = std::numeric_limits<int>::max();

        for (const auto& r : reduced) {
            if (!contains(*r, curr_var)) continue;
            curr_var_appears = true;

            bool has_other = false;
            for (int i = 0; i < var_pos; ++i) {
                if (contains(*r, variables[i])) {
                    has_other = true;
                    break;
                }
            }
            if (has_other) continue;

            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(r, curr_var);
            int deg = poly.degree();
            if (deg >= 1 && deg < best_deg) {
                best_deg = deg;
                target = r;
            }
        }
        if (!curr_var_appears) {
            auto next_partial = partial;
            next_partial.insert_or_assign(curr_var, *SymbolicExpr::variable(curr_var));
            return self(self, var_pos - 1, next_partial);
        }
        if (!target) return {};

        auto roots = detail::propagate_result(solve_finite_checked(
            target, curr_var, context, SolveOptions{}));
        if (roots.empty()) return {};
        std::vector<std::map<std::string, SymbolicExpr>> results;
        for (const auto& r : roots) {
            auto next_partial = partial;
            next_partial.insert_or_assign(curr_var, *r);
            auto sub_res = self(self, var_pos - 1, next_partial);
            results.insert(results.end(), sub_res.begin(), sub_res.end());
        }
        return results;
    };
    std::map<std::string, SymbolicExpr> empty;
    auto candidates = solve_rec(solve_rec, static_cast<int>(variables.size()) - 1, empty);
    if (denom_constraints.empty()) return candidates;
    std::vector<std::map<std::string, SymbolicExpr>> filtered;
    filtered.reserve(candidates.size());
    for (const auto& sol : candidates) {
        bool ok = true;
        for (const auto& den : denom_constraints) {
            auto sub = den;
            for (const auto& [name, val] : sol) {
                sub = sub->substitute(name, lamina::detail::make_expression_ptr(val));
                if (!sub) break;
            }
            if (!sub) continue;
            sub = sub->simplify();
            if (sub && sub->is_zero()) {
                ok = false;
                break;
            }
        }
        if (ok) filtered.push_back(sol);
    }
    return filtered;
}
PolynomialSystemResult Solver::solve_polynomial_system_checked(
    const std::vector<SymbolicExpr>& equations,
    const std::vector<std::string>& variables,
    ComputationContext& context)
{
    constexpr const char* operation = "solve_polynomial_system";
    if (equations.empty() || variables.empty()) {
        return PolynomialSystemResult::failure(
            CasErrc::InvalidArgument,
            "polynomial system requires equations and variables", operation);
    }
    auto budget = context.consume_steps(
        equations.size() * variables.size() + 1, operation);
    if (!budget) return PolynomialSystemResult::failure(budget.error());
    try {
        return PolynomialSystemResult::success(
            solve_polynomial_system_impl(equations, variables, context));
    } catch (const detail::ResultPropagation& propagation) {
        return PolynomialSystemResult::failure(propagation.error());
    } catch (const std::bad_alloc&) {
        return PolynomialSystemResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving polynomial system", operation);
    } catch (const std::exception& ex) {
        return PolynomialSystemResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}
PolynomialSystemResult Solver::solve_polynomial_system_checked(
    const std::vector<SymbolicExpr>& equations,
    const std::vector<std::string>& variables)
{
    ComputationContext context;
    return solve_polynomial_system_checked(equations, variables, context);
}
std::vector<SymbolicExpr> Solver::reduced_groebner_basis(
    const std::vector<SymbolicExpr>& polynomials,
    const std::vector<std::string>& variables)
{

    auto gb_exprs = groebner_basis(polynomials, variables);
    if (gb_exprs.empty()) return {};

    PolyContext ctx(variables);
    std::vector<Poly> G;
    for (const auto& expr : gb_exprs) {
        Poly p = to_poly(expr, ctx);
        if (!p.is_zero()) G.push_back(p);
    }

    size_t max_vars = ctx.ext_vars.size();
    for (auto& p : G) {
        if (p.num_vars < max_vars) {
            Poly resized(max_vars);
            for (auto& [m, c] : p.terms) {
                Monomial extended = m;
                extended.resize(max_vars, 0);
                resized.add_term(extended, c);
            }
            p = resized;
        }
    }

    std::vector<bool> marked(G.size(), false);
    for (size_t i = 0; i < G.size(); ++i) {
        if (marked[i]) continue;
        for (size_t j = 0; j < G.size(); ++j) {
            if (i == j || marked[j]) continue;
            if (divides_mono(G[j].LM(), G[i].LM())) {
                marked[i] = true;
                break;
            }
        }
    }
    std::vector<Poly> minimal;
    for (size_t i = 0; i < G.size(); ++i) {
        if (!marked[i]) minimal.push_back(G[i]);
    }

    for (size_t i = 0; i < minimal.size(); ++i) {
        std::vector<Poly> others;
        for (size_t j = 0; j < minimal.size(); ++j) {
            if (j != i) others.push_back(minimal[j]);
        }
        minimal[i] = reduce(minimal[i], others);
    }

    std::vector<Poly> reduced;
    for (auto& p : minimal) {
        if (!p.is_zero()) reduced.push_back(p);
    }

    for (auto& p : reduced) {
        Rational lc = p.LC();
        if (lc != Rational(0) && lc != Rational(1)) {
            Rational inv = Rational(1) / lc;
            Poly monic(p.num_vars);
            for (auto& [m, c] : p.terms) {
                monic.add_term(m, c * inv);
            }
            p = monic;
        }
    }

    std::vector<SymbolicExpr> result;
    for (const auto& p : reduced) {
        result.push_back(from_poly_ext(p, ctx, variables));
    }
    return result;
}

bool Solver::ideal_membership(
    const SymbolicExpr& polynomial,
    const std::vector<SymbolicExpr>& basis,
    const std::vector<std::string>& variables)
{

    PolyContext ctx(variables);

    std::vector<Poly> basis_polys;
    for (const auto& b : basis) {
        Poly p = to_poly(b, ctx);
        if (!p.is_zero()) basis_polys.push_back(p);
    }

    Poly poly = to_poly(polynomial, ctx);

    size_t max_vars = ctx.ext_vars.size();
    for (auto& p : basis_polys) {
        if (p.num_vars < max_vars) {
            Poly resized(max_vars);
            for (auto& [m, c] : p.terms) {
                Monomial extended = m;
                extended.resize(max_vars, 0);
                resized.add_term(extended, c);
            }
            p = resized;
        }
    }
    if (poly.num_vars < max_vars) {
        Poly resized(max_vars);
        for (auto& [m, c] : poly.terms) {
            Monomial extended = m;
            extended.resize(max_vars, 0);
            resized.add_term(extended, c);
        }
        poly = resized;
    }

    Poly remainder = reduce(poly, basis_polys);
    return remainder.is_zero();
}

std::vector<SymbolicExpr> Solver::elimination_ideal(
    const std::vector<SymbolicExpr>& basis,
    const std::vector<std::string>& variables,
    int elim_count)
{
    if (elim_count <= 0) {
        return basis;
    }
    if (elim_count >= static_cast<int>(variables.size())) {
        return {};
    }

    PolyContext ctx(variables);
    std::vector<SymbolicExpr> result;

    for (const auto& expr : basis) {
        Poly p = to_poly(expr, ctx);
        if (p.is_zero()) continue;

        bool involves_eliminated = false;
        for (const auto& [m, c] : p.terms) {
            for (int i = 0; i < elim_count && i < static_cast<int>(m.size()); ++i) {
                if (m[i] != 0) {
                    involves_eliminated = true;
                    break;
                }
            }
            if (involves_eliminated) break;
        }

        if (!involves_eliminated) {
            result.push_back(expr);
        }
    }

    return result;
}

} // namespace lamina
