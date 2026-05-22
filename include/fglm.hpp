#pragma once

// FGLM Algorithm for Gröbner Basis Conversion
// Converts a Gröbner basis from one monomial ordering to another.
// Only works for zero-dimensional ideals (finitely many solutions).
//
// Reference: Faugère, Gianni, Lazard, Mora (1993)
// "Efficient Computation of Zero-Dimensional Gröbner Bases by Change of Ordering"

#include "monomial_order.hpp"
#include "rational.hpp"
#include <vector>
#include <map>
#include <algorithm>
#include <queue>
#include <set>
#include <stdexcept>
#include <cassert>
#include <functional>

namespace lamina {

// ============================================================================
// FGLMPoly: A multivariate polynomial with configurable monomial ordering.
// Terms are stored as (monomial, coefficient) pairs sorted by the given order.
// ============================================================================
struct FGLMPoly {
    std::vector<std::pair<Monomial, Rational>> terms; // sorted largest-first
    size_t num_vars;

    FGLMPoly() : num_vars(0) {}
    explicit FGLMPoly(size_t n) : num_vars(n) {}

    bool is_zero() const { return terms.empty(); }

    /// Leading monomial (largest term under the current sort order)
    Monomial LM() const {
        if (terms.empty()) return Monomial();
        return terms.front().first;
    }

    /// Leading coefficient
    Rational LC() const {
        if (terms.empty()) return Rational(0);
        return terms.front().second;
    }

    /// Add a term. Does NOT automatically sort or combine duplicates.
    void add_term(const Monomial& m, const Rational& c) {
        if (!c.is_zero()) {
            terms.emplace_back(m, c);
        }
    }

    /// Sort terms by the given monomial order (largest first).
    void sort_terms(const MonomialOrder& order) {
        std::sort(terms.begin(), terms.end(),
            [&order](const std::pair<Monomial, Rational>& a,
                     const std::pair<Monomial, Rational>& b) {
                return order(a.first, b.first);
            });
    }

    /// Remove zero-coefficient terms and combine duplicates.
    void normalize() {
        // Combine duplicates (assumes sorted)
        std::vector<std::pair<Monomial, Rational>> cleaned;
        for (auto& [m, c] : terms) {
            if (c.is_zero()) continue;
            if (!cleaned.empty() && cleaned.back().first == m) {
                cleaned.back().second = cleaned.back().second + c;
                if (cleaned.back().second.is_zero()) {
                    cleaned.pop_back();
                }
            } else {
                cleaned.emplace_back(m, c);
            }
        }
        terms = std::move(cleaned);
    }

    /// Polynomial reduction: compute normal form of this poly w.r.t. a basis.
    FGLMPoly reduce(const std::vector<FGLMPoly>& basis,
                    const MonomialOrder& order) const {
        FGLMPoly r(num_vars);
        FGLMPoly f = *this;
        f.sort_terms(order);
        f.normalize();

        int max_steps = 10000; // safety limit
        while (!f.is_zero() && max_steps-- > 0) {
            bool reduced = false;
            for (const auto& g : basis) {
                if (g.is_zero()) continue;
                Monomial lm_f = f.LM();
                Monomial lm_g = g.LM();
                if (divides_monomial(lm_g, lm_f)) {
                    // Compute quotient monomial and coefficient
                    Monomial quot_mon(num_vars, 0);
                    for (size_t i = 0; i < num_vars; ++i) {
                        int fi = (i < lm_f.size()) ? lm_f[i] : 0;
                        int gi = (i < lm_g.size()) ? lm_g[i] : 0;
                        quot_mon[i] = fi - gi;
                    }
                    Rational quot_coeff = f.LC() / g.LC();

                    // f = f - quot_coeff * x^quot_mon * g
                    // Build the subtracted polynomial directly
                    FGLMPoly subtracted(num_vars);
                    for (const auto& [gm, gc] : g.terms) {
                        Monomial product_mon(num_vars, 0);
                        for (size_t i = 0; i < num_vars; ++i) {
                            int qi = (i < quot_mon.size()) ? quot_mon[i] : 0;
                            int gmi = (i < gm.size()) ? gm[i] : 0;
                            product_mon[i] = qi + gmi;
                        }
                        subtracted.add_term(product_mon, quot_coeff * gc);
                    }
                    
                    // Subtract: f = f - subtracted
                    for (const auto& [sm, sc] : subtracted.terms) {
                        f.add_term(sm, -sc);
                    }
                    f.sort_terms(order);
                    f.normalize();
                    reduced = true;
                    break;
                }
            }
            if (!reduced) {
                // Leading term is irreducible; move to remainder
                r.add_term(f.LM(), f.LC());
                f.terms.erase(f.terms.begin());
                f.normalize();
            }
        }
        r.sort_terms(order);
        r.normalize();
        return r;
    }

    /// Create a polynomial consisting of a single monomial with coefficient 1.
    static FGLMPoly from_monomial(const Monomial& m, size_t num_vars) {
        FGLMPoly p(num_vars);
        p.add_term(m, Rational(1));
        return p;
    }

    /// Create the zero polynomial.
    static FGLMPoly zero(size_t num_vars) {
        return FGLMPoly(num_vars);
    }
};

// ============================================================================
// Helper: Compute normal form of a polynomial w.r.t. a Gröbner basis.
// This is multivariate polynomial division returning the remainder.
// ============================================================================
inline FGLMPoly normal_form(const FGLMPoly& f,
                            const std::vector<FGLMPoly>& basis,
                            const MonomialOrder& order) {
    return f.reduce(basis, order);
}

// ============================================================================
// Helper: Check if an ideal is zero-dimensional.
// An ideal I in Q[x1,...,xn] is zero-dimensional iff for each variable xi,
// there exists a basis element whose leading monomial is a pure power of xi.
// ============================================================================
inline bool is_zero_dimensional(const std::vector<FGLMPoly>& basis, size_t num_vars) {
    if (basis.empty()) return false;

    // For each variable, check if some LM is a pure power of that variable
    for (size_t var = 0; var < num_vars; ++var) {
        bool found = false;
        for (const auto& g : basis) {
            if (g.is_zero()) continue;
            Monomial lm = g.LM();
            // Check if lm is a pure power of variable var
            bool is_pure = true;
            for (size_t i = 0; i < num_vars; ++i) {
                int exp_i = (i < lm.size()) ? lm[i] : 0;
                if (i == var) {
                    if (exp_i == 0) { is_pure = false; break; }
                } else {
                    if (exp_i != 0) { is_pure = false; break; }
                }
            }
            if (is_pure) { found = true; break; }
        }
        if (!found) return false;
    }
    return true;
}

// ============================================================================
// Helper: Compute the dimension of Q[x1,...,xn]/I.
// Counts standard monomials (those not divisible by any LM in the basis)
// up to a maximum degree bound.
// ============================================================================
inline int quotient_dimension(const std::vector<FGLMPoly>& basis,
                              size_t num_vars,
                              int max_degree = 50) {
    if (basis.empty()) return -1; // infinite or unknown

    // Collect leading monomials
    std::vector<Monomial> lead_mons;
    for (const auto& g : basis) {
        if (!g.is_zero()) {
            lead_mons.push_back(g.LM());
        }
    }

    // Count monomials of degree 0..max_degree not divisible by any LM
    int count = 0;
    const int limit = 1000; // safety cap

    // Generate monomials degree by degree
    // For degree d, enumerate all monomials of total degree d in num_vars variables
    for (int d = 0; d <= max_degree && count <= limit; ++d) {
        // Generate all monomials of degree d using stars-and-bars
        // We use a recursive approach via a stack
        std::vector<Monomial> degree_mons;
        Monomial current(num_vars, 0);

        // Helper lambda to generate monomials of exact degree d
        std::function<void(size_t, int)> generate =
            [&](size_t var_idx, int remaining_deg) {
                if (var_idx == num_vars - 1) {
                    current[var_idx] = remaining_deg;
                    degree_mons.push_back(current);
                    current[var_idx] = 0;
                    return;
                }
                for (int e = 0; e <= remaining_deg; ++e) {
                    current[var_idx] = e;
                    generate(var_idx + 1, remaining_deg - e);
                }
                current[var_idx] = 0;
            };

        generate(0, d);

        for (const auto& m : degree_mons) {
            // Check if m is divisible by any leading monomial
            bool divisible = false;
            for (const auto& lm : lead_mons) {
                if (divides_monomial(lm, m)) {
                    divisible = true;
                    break;
                }
            }
            if (!divisible) {
                ++count;
                if (count > limit) return count; // too large
            }
        }

        // If no standard monomials found at this degree and d > 0,
        // all higher degrees will also be zero (for zero-dim ideals)
        if (degree_mons.size() > 0) {
            bool any_standard = false;
            for (const auto& m : degree_mons) {
                bool divisible = false;
                for (const auto& lm : lead_mons) {
                    if (divides_monomial(lm, m)) { divisible = true; break; }
                }
                if (!divisible) { any_standard = true; break; }
            }
            if (!any_standard && d > 0) break; // no more standard monomials
        }
    }
    return count;
}

// ============================================================================
// Internal helpers for the FGLM algorithm
// ============================================================================
namespace detail {

/// Monomial comparator for use in ordered containers (less-than under order).
/// Returns true if a < b (i.e., b > a under the given ordering).
struct MonomialLessUnder {
    const MonomialOrder* order;
    MonomialLessUnder(const MonomialOrder* o) : order(o) {}
    bool operator()(const Monomial& a, const Monomial& b) const {
        return order->operator()(b, a); // a < b means b > a
    }
};

/// Enumerate monomials in ascending order under the given ordering.
/// Uses a priority queue (min-heap) to generate monomials smallest-first.
/// Has a maximum degree bound to prevent unbounded memory growth.
class MonomialEnumerator {
public:
    MonomialEnumerator(size_t num_vars, const MonomialOrder& order, int max_deg = 100)
        : num_vars_(num_vars), order_(order), max_degree_(max_deg) {
        // Start with the constant monomial (all zeros)
        Monomial one(num_vars, 0);
        push(one);
    }

    /// Get the next smallest monomial. Returns false if exhausted.
    bool next(Monomial& out) {
        if (heap_.empty()) return false;
        out = heap_.front();
        std::pop_heap(heap_.begin(), heap_.end(), cmp_);
        heap_.pop_back();

        // Generate successors: multiply by each variable (if within degree bound)
        int out_deg = total_degree(out);
        if (out_deg < max_degree_) {
            for (size_t i = 0; i < num_vars_; ++i) {
                Monomial succ = out;
                succ[i] += 1;
                push(succ);
            }
        }
        return true;
    }

private:
    size_t num_vars_;
    MonomialOrder order_;
    int max_degree_;
    std::vector<Monomial> heap_;
    std::set<Monomial> visited_;

    // Comparator for max-heap that gives us smallest first
    // We want a min-heap: the "greatest" element in the heap comparator
    // should be the smallest monomial under our ordering.
    struct HeapCmp {
        const MonomialOrder* order;
        HeapCmp() : order(nullptr) {}
        HeapCmp(const MonomialOrder* o) : order(o) {}
        bool operator()(const Monomial& a, const Monomial& b) const {
            // For std::pop_heap to give us the smallest, we need
            // the comparator to return true when a > b
            return order->operator()(a, b);
        }
    };
    HeapCmp cmp_{&order_};

    void push(const Monomial& m) {
        if (visited_.count(m)) return;
        visited_.insert(m);
        heap_.push_back(m);
        std::push_heap(heap_.begin(), heap_.end(), cmp_);
    }
};

/// Incremental Gaussian elimination over Rational vectors.
/// Used to detect linear dependence as we add normal form vectors.
class GaussianEliminator {
public:
    GaussianEliminator() {}

    /// Try to add a vector. Returns true if it was linearly independent
    /// (and has been added to the basis). Returns false if dependent.
    /// If dependent, `combination` is filled with the coefficients expressing
    /// the input as a linear combination of previously added basis vectors.
    bool add_vector(const std::vector<Rational>& v,
                    std::vector<Rational>& combination) {
        size_t n = v.size();
        // Extend all existing rows if needed
        for (auto& row : rows_) {
            row.resize(n, Rational(0));
        }

        // Copy v and track the linear combination
        std::vector<Rational> working = v;
        // combination[i] = coefficient of the i-th basis vector
        combination.assign(basis_count_, Rational(0));

        // Reduce using existing pivots
        for (size_t i = 0; i < pivots_.size(); ++i) {
            size_t col = pivots_[i];
            if (col >= working.size()) continue;
            if (working[col].is_zero()) continue;

            Rational factor = working[col] / rows_[i][col];
            for (size_t j = 0; j < n; ++j) {
                if (j < working.size() && j < rows_[i].size()) {
                    working[j] = working[j] - factor * rows_[i][j];
                }
            }
            // Track combination
            for (size_t j = 0; j < combination.size(); ++j) {
                combination[j] = combination[j] - factor * combinations_[i][j];
            }
        }

        // Find pivot in the reduced vector
        size_t pivot_col = n; // sentinel
        for (size_t j = 0; j < n; ++j) {
            if (!working[j].is_zero()) {
                pivot_col = j;
                break;
            }
        }

        if (pivot_col == n) {
            // Linearly dependent — combination holds the coefficients
            // Negate to get: v = sum(combination[i] * basis_i)
            for (auto& c : combination) {
                c = -c;
            }
            return false;
        }

        // Linearly independent — normalize and store
        Rational pivot_val = working[pivot_col];
        for (size_t j = 0; j < n; ++j) {
            working[j] = working[j] / pivot_val;
        }

        // Normalize the combination tracking
        std::vector<Rational> comb_row(basis_count_ + 1, Rational(0));
        for (size_t j = 0; j < combination.size(); ++j) {
            comb_row[j] = combination[j] / pivot_val;
        }
        // The new basis vector's own coefficient is 1/pivot_val
        comb_row[basis_count_] = Rational(1) / pivot_val;

        rows_.push_back(working);
        pivots_.push_back(pivot_col);
        combinations_.push_back(comb_row);
        ++basis_count_;

        return true;
    }

    size_t rank() const { return basis_count_; }

private:
    std::vector<std::vector<Rational>> rows_;
    std::vector<size_t> pivots_;
    std::vector<std::vector<Rational>> combinations_;
    size_t basis_count_ = 0;
};

/// Convert a polynomial (normal form) to a coefficient vector
/// w.r.t. a given list of monomials. Extends monomial_index as needed.
inline std::vector<Rational> poly_to_vector(
    const FGLMPoly& p,
    std::vector<Monomial>& all_monomials,
    std::map<Monomial, size_t>& monomial_index) {

    // First, ensure all monomials in p are in the index
    for (const auto& [m, c] : p.terms) {
        if (monomial_index.find(m) == monomial_index.end()) {
            size_t idx = all_monomials.size();
            all_monomials.push_back(m);
            monomial_index[m] = idx;
        }
    }

    // Build the vector
    std::vector<Rational> v(all_monomials.size(), Rational(0));
    for (const auto& [m, c] : p.terms) {
        v[monomial_index[m]] = c;
    }
    return v;
}

} // namespace detail

// ============================================================================
// FGLM Conversion Algorithm
//
// Given a Gröbner basis G under source_order, compute the Gröbner basis
// under target_order. The ideal must be zero-dimensional.
//
// Algorithm:
// 1. Enumerate monomials in ascending target_order.
// 2. For each monomial t:
//    a. Compute NF(t, G_source) under source_order.
//    b. Express NF(t) as a coefficient vector.
//    c. If linearly dependent on previous NFs → new GB element found.
//    d. Otherwise, add to the linear algebra basis.
// 3. Terminate when we have found enough generators or exhausted the
//    quotient space dimension.
// ============================================================================
inline std::vector<FGLMPoly> fglm_convert(
    const std::vector<FGLMPoly>& source_basis,
    const MonomialOrder& source_order,
    const MonomialOrder& target_order,
    size_t num_vars) {

    if (source_basis.empty()) {
        return {};
    }

    // Verify zero-dimensionality
    if (!is_zero_dimensional(source_basis, num_vars)) {
        throw std::runtime_error(
            "FGLM: ideal is not zero-dimensional. "
            "FGLM only works for zero-dimensional ideals.");
    }

    // Compute expected quotient dimension
    int dim = quotient_dimension(source_basis, num_vars);
    if (dim <= 0 || dim > 1000) {
        throw std::runtime_error(
            "FGLM: quotient dimension is too large or could not be determined. "
            "dim = " + std::to_string(dim));
    }

    std::vector<FGLMPoly> target_basis; // result GB under target_order

    // Monomials that form the basis of the quotient ring (standard monomials)
    std::vector<Monomial> basis_monomials;

    // Pre-compute the set of standard monomials (those not divisible by any LM)
    // These form the basis of the quotient ring and bound the vector space dimension.
    std::vector<Monomial> standard_monomials;
    std::map<Monomial, size_t> monomial_index;
    {
        std::vector<Monomial> source_lms;
        for (const auto& g : source_basis) {
            if (!g.is_zero()) source_lms.push_back(g.LM());
        }
        
        // Generate standard monomials degree by degree
        for (int d = 0; d <= dim + 2 && (int)standard_monomials.size() < dim; ++d) {
            Monomial current(num_vars, 0);
            std::function<void(size_t, int)> gen = [&](size_t var_idx, int remaining) {
                if (var_idx == num_vars - 1) {
                    current[var_idx] = remaining;
                    // Check if standard
                    bool is_std = true;
                    for (const auto& lm : source_lms) {
                        if (divides_monomial(lm, current)) { is_std = false; break; }
                    }
                    if (is_std) standard_monomials.push_back(current);
                    current[var_idx] = 0;
                    return;
                }
                for (int e = 0; e <= remaining; ++e) {
                    current[var_idx] = e;
                    gen(var_idx + 1, remaining - e);
                }
                current[var_idx] = 0;
            };
            gen(0, d);
        }
        
        // Build index
        for (size_t i = 0; i < standard_monomials.size(); ++i) {
            monomial_index[standard_monomials[i]] = i;
        }
    }
    
    // Helper: convert a normal form polynomial to a fixed-size vector
    // indexed by standard monomials
    auto nf_to_vector = [&](const FGLMPoly& nf) -> std::vector<Rational> {
        std::vector<Rational> v(standard_monomials.size(), Rational(0));
        for (const auto& [m, c] : nf.terms) {
            auto it = monomial_index.find(m);
            if (it != monomial_index.end()) {
                v[it->second] = c;
            }
        }
        return v;
    };

    // Gaussian eliminator for linear dependence detection
    detail::GaussianEliminator gauss;

    // Monomial enumerator in ascending target_order (with degree bound)
    int max_mono_degree = dim + 5; // generous bound: max degree in GB shouldn't exceed dim much
    detail::MonomialEnumerator enumerator(num_vars, target_order, max_mono_degree);

    // Leading monomials of the target basis (to skip reducible monomials)
    std::vector<Monomial> target_lms;

    // Safety limit on iterations
    const int max_iterations = dim * 10 + 50;
    int iterations = 0;

    Monomial t;
    while (enumerator.next(t) && iterations < max_iterations) {
        ++iterations;

        // Skip monomials divisible by a leading monomial of the target basis.
        // These cannot contribute new basis elements.
        bool skip = false;
        for (const auto& lm : target_lms) {
            if (divides_monomial(lm, t)) {
                skip = true;
                break;
            }
        }
        if (skip) continue;

        // Compute NF(t, G_source) under source_order
        FGLMPoly t_poly = FGLMPoly::from_monomial(t, num_vars);
        FGLMPoly nf = normal_form(t_poly, source_basis, source_order);

        // Convert normal form to a coefficient vector (fixed dimension)
        std::vector<Rational> v = nf_to_vector(nf);

        // Check linear independence
        std::vector<Rational> combination;
        bool independent = gauss.add_vector(v, combination);

        if (independent) {
            // t is a new standard monomial (basis of quotient ring)
            basis_monomials.push_back(t);

            // Termination: found all standard monomials
            if (static_cast<int>(basis_monomials.size()) >= dim) {
                // We've found the full quotient basis; any further
                // monomials will be dependent. But we continue to find
                // all generators of the target GB.
                // Actually, once we have dim standard monomials,
                // every subsequent monomial must be dependent.
                // We can stop if we also have enough generators.
            }
        } else {
            // t is linearly dependent on previous standard monomials.
            // Construct the new GB element: t - linear_combination
            FGLMPoly new_elem(num_vars);
            new_elem.add_term(t, Rational(1)); // leading term is t

            // Subtract the linear combination of basis monomials
            for (size_t i = 0; i < combination.size() && i < basis_monomials.size(); ++i) {
                if (!combination[i].is_zero()) {
                    // combination[i] * NF(basis_monomials[i])
                    // But we need: t = sum(combination[i] * basis_monomials[i])
                    // in the quotient ring. The GB element is:
                    // t - sum(combination[i] * basis_monomials[i])
                    new_elem.add_term(basis_monomials[i], -combination[i]);
                }
            }

            new_elem.sort_terms(target_order);
            new_elem.normalize();

            // Make monic (leading coefficient = 1)
            if (!new_elem.is_zero()) {
                Rational lc = new_elem.LC();
                if (!(lc == Rational(1))) {
                    for (auto& [m, c] : new_elem.terms) {
                        c = c / lc;
                    }
                }
            }

            if (!new_elem.is_zero()) {
                target_lms.push_back(new_elem.LM());
                target_basis.push_back(std::move(new_elem));
            }
        }

        // Heuristic termination: if we have found n generators and
        // all standard monomials, we can stop.
        if (static_cast<int>(basis_monomials.size()) >= dim) {
            // All remaining monomials will be reducible or dependent.
            // Check if we have generators for all variables.
            bool complete = true;
            for (size_t var = 0; var < num_vars; ++var) {
                bool has_gen = false;
                for (const auto& lm : target_lms) {
                    // Check if some LM involves only variable var
                    // (or at least has a pure power of var)
                    bool is_pure = true;
                    for (size_t i = 0; i < num_vars; ++i) {
                        int exp_i = (i < lm.size()) ? lm[i] : 0;
                        if (i == var) {
                            if (exp_i == 0) { is_pure = false; break; }
                        } else {
                            if (exp_i != 0) { is_pure = false; break; }
                        }
                    }
                    if (is_pure) { has_gen = true; break; }
                }
                if (!has_gen) { complete = false; break; }
            }
            if (complete) break;
        }
    }

    return target_basis;
}

// ============================================================================
// Convenience: Convert grevlex GB to lex GB (most common use case)
// ============================================================================
inline std::vector<FGLMPoly> grevlex_to_lex(
    const std::vector<FGLMPoly>& grevlex_basis,
    size_t num_vars) {
    return fglm_convert(
        grevlex_basis,
        MonomialOrder::grevlex(),
        MonomialOrder::lex(),
        num_vars);
}

} // namespace lamina
