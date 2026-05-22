#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

namespace lamina {

using Monomial = std::vector<int>;

/// Supported monomial ordering types for multivariate polynomial rings.
enum class MonomialOrderType {
    Lex,        // Pure lexicographic
    GrevLex,    // Graded reverse lexicographic (same as DegRevLex)
    DegLex,     // Graded lexicographic (total degree first, then lex)
    DegRevLex   // Degree reverse lexicographic (standard name for GrevLex)
};

/// Returns the total degree of a monomial (sum of all exponents).
inline int total_degree(const Monomial& m) {
    return std::accumulate(m.begin(), m.end(), 0);
}

/// Returns the component-wise maximum (LCM of monomials viewed as divisors).
inline Monomial lcm_monomial(const Monomial& a, const Monomial& b) {
    Monomial result(std::max(a.size(), b.size()), 0);
    for (size_t i = 0; i < a.size(); ++i)
        result[i] = std::max(result[i], a[i]);
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = std::max(result[i], b[i]);
    return result;
}

/// Returns true if divisor divides target component-wise (each exponent <=).
inline bool divides_monomial(const Monomial& divisor, const Monomial& target) {
    if (divisor.size() > target.size()) {
        for (size_t i = target.size(); i < divisor.size(); ++i)
            if (divisor[i] != 0) return false;
    }
    for (size_t i = 0; i < std::min(divisor.size(), target.size()); ++i) {
        if (divisor[i] > target[i]) return false;
    }
    return true;
}

/// Configurable monomial ordering comparator.
/// All orderings return true if a > b (i.e., a is "larger" and comes first).
class MonomialOrder {
public:
    explicit MonomialOrder(MonomialOrderType type) : type_(type) {}

    /// Compare two monomials. Returns true if a > b under this ordering.
    inline bool operator()(const Monomial& a, const Monomial& b) const {
        switch (type_) {
            case MonomialOrderType::Lex:
                return compare_lex(a, b);
            case MonomialOrderType::GrevLex:
            case MonomialOrderType::DegRevLex:
                return compare_grevlex(a, b);
            case MonomialOrderType::DegLex:
                return compare_deglex(a, b);
        }
        return false;
    }

    MonomialOrderType type() const { return type_; }

    /// Factory: pure lexicographic ordering.
    /// Compares left-to-right; first nonzero difference decides.
    static MonomialOrder lex() { return MonomialOrder(MonomialOrderType::Lex); }

    /// Factory: graded reverse lexicographic ordering.
    /// Compares total degree first; ties broken right-to-left in reverse.
    static MonomialOrder grevlex() { return MonomialOrder(MonomialOrderType::GrevLex); }

    /// Factory: graded lexicographic ordering.
    /// Compares total degree first; ties broken by lex.
    static MonomialOrder deglex() { return MonomialOrder(MonomialOrderType::DegLex); }

    /// Factory: degree reverse lexicographic (standard name for grevlex).
    static MonomialOrder degrevlex() { return MonomialOrder(MonomialOrderType::DegRevLex); }

private:
    MonomialOrderType type_;

    /// Lex: a > b iff the leftmost nonzero entry of (a - b) is positive.
    static inline bool compare_lex(const Monomial& a, const Monomial& b) {
        size_t n = std::max(a.size(), b.size());
        for (size_t i = 0; i < n; ++i) {
            int ai = (i < a.size()) ? a[i] : 0;
            int bi = (i < b.size()) ? b[i] : 0;
            if (ai != bi) return ai > bi;
        }
        return false; // equal
    }

    /// GrevLex / DegRevLex: compare total degree first.
    /// If equal, compare from right to left; a > b iff the rightmost
    /// nonzero entry of (a - b) is negative.
    static inline bool compare_grevlex(const Monomial& a, const Monomial& b) {
        int deg_a = total_degree(a);
        int deg_b = total_degree(b);
        if (deg_a != deg_b) return deg_a > deg_b;

        // Same total degree: scan from right to left
        size_t n = std::max(a.size(), b.size());
        for (size_t i = n; i > 0; --i) {
            int ai = (i - 1 < a.size()) ? a[i - 1] : 0;
            int bi = (i - 1 < b.size()) ? b[i - 1] : 0;
            if (ai != bi) return ai < bi; // reverse: a > b when a[i] < b[i]
        }
        return false; // equal
    }

    /// DegLex: compare total degree first; ties broken by lex.
    static inline bool compare_deglex(const Monomial& a, const Monomial& b) {
        int deg_a = total_degree(a);
        int deg_b = total_degree(b);
        if (deg_a != deg_b) return deg_a > deg_b;
        return compare_lex(a, b);
    }
};

} // namespace lamina
