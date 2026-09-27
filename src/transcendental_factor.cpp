#include "transcendental_factor.hpp"
#include "internal/symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/transcendental_support.hpp"
#include "internal/multivariate_factor_support.hpp"

#include <exception>
#include <new>
#include <string>
#include <vector>

namespace LMCAS {

static std::vector<std::shared_ptr<SymbolicExpr>> tf_restore_factors(
    const std::vector<Polynomial<Rational>>& true_factors,
    const TransSubstitutionResult& sub_result,
    const std::shared_ptr<SymbolicExpr>& pyth_simplified) {
    /// --- Phase 6: 逆换元 ---
    std::vector<std::shared_ptr<SymbolicExpr>> symbolic_factors;
    symbolic_factors.reserve(true_factors.size());

    for (const auto& poly_factor : true_factors) {
        /// 将多项式因子转回符号表达式
        auto sym_factor = poly_to_symbolic(poly_factor);
        if (!sym_factor) continue;

        /// 逆换元:将不定元替换回原始超越表达式
        auto back_sub = tf_back_substitute(sym_factor, sub_result.mappings);
        if (back_sub) {
            symbolic_factors.push_back(back_sub);
        }
    }

    /// 若逆换元后因子数 <= 1,分解无效
    if (symbolic_factors.size() <= 1) {
        return {pyth_simplified};
    }

    /// --- Phase 6.2: 化简与常数提取 ---
    auto final_factors = tf_simplify_factors(symbolic_factors);

    /// 若化简后仅剩一个非常数因子,分解无效
    size_t non_const_count = 0;
    for (const auto& f : final_factors) {
        if (f && !f->is_number()) non_const_count++;
    }
    if (non_const_count <= 1 && final_factors.size() <= 1) {
        return {pyth_simplified};
    }

    return final_factors;
}

static Result<std::vector<std::shared_ptr<SymbolicExpr>>> tf_factor_substitution(
    const TransSubstitutionResult& sub_result, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& pyth_simplified,
    ComputationContext& context) {
    using FactorsResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;
    std::vector<std::string> indeterminates;
    indeterminates.reserve(sub_result.mappings.size());
    for (const auto& m : sub_result.mappings) {
        indeterminates.push_back(m.indeterminate);
    }
    auto built = tf_build_polynomial(
        sub_result.poly_expr, indeterminates, var);
    if (!built) { return FactorsResult::failure(built.error()); }
    const auto& poly_result = built.value();
    if (!poly_result.success) {
        return FactorsResult::success({pyth_simplified});
    }
    if (poly_result.poly.degree() <= 1) {
        return FactorsResult::success({pyth_simplified});
    }
    auto factored =
        factor_univariate_bridge_checked(poly_result.poly, context);
    if (!factored) {
        return FactorsResult::failure(factored.error());
    }
    const auto& factorization = factored.value();
    if (factorization.completeness != Completeness::Complete ||
        factorization.value.size() <= 1) {
        return FactorsResult::success({pyth_simplified});
    }
    return tf_restore_factors(
        factorization.value, sub_result, pyth_simplified);
}

Result<std::vector<std::shared_ptr<SymbolicExpr>>> factor_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context) try {

    using FactorsResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;
    if (!expr || !LMCAS::detail::node(expr) || var.empty()) {
        return FactorsResult::failure(CasErrc::InvalidArgument,
            "expression and variable must be nonempty", "factor_transcendental");
    }
    auto step = context.consume_steps(1, "factor_transcendental");
    if (!step) {
        return FactorsResult::failure(step.error());
    }
    auto pyth_simplified = tf_simplify_pythagorean(expr, var);
    auto mult_factors = tf_detect_multiplicative_structure(pyth_simplified);
    if (!mult_factors.empty()) {
        return mult_factors;
    }
    auto exp_factors = tf_detect_exponential_separation(pyth_simplified, var);
    if (!exp_factors.empty()) {
        return exp_factors;
    }

    if (!tf_contains_transcendental(LMCAS::detail::node(pyth_simplified), var)) {
        return FactorsResult::success({pyth_simplified});
    }

    TransSubstitutionResult sub_result = detect_trans_substitutions(pyth_simplified, var);
    if (sub_result.mappings.empty()) {
        return FactorsResult::success({pyth_simplified});
    }
    if (tf_is_linear_irreducible(sub_result, var)) {
        return FactorsResult::success({pyth_simplified});
    }

    return tf_factor_substitution(
        sub_result, var, pyth_simplified, context);
}
catch (const std::bad_alloc&) {
    return Result<std::vector<std::shared_ptr<SymbolicExpr>>>::failure(
        CasErrc::ResourceLimit, "allocation failed while factoring expression",
        "factor_transcendental");
} catch (const std::exception& ex) {
    return Result<std::vector<std::shared_ptr<SymbolicExpr>>>::failure(
        CasErrc::InternalInvariant, ex.what(), "factor_transcendental");
}

Result<std::vector<std::shared_ptr<SymbolicExpr>>> factor_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {
    ComputationContext context;
    return factor_transcendental(expr, var, context);
}

}
