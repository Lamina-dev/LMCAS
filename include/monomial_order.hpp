#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

namespace lamina {

using Monomial = std::vector<int>;

enum class MonomialOrderType {
    Lex,
    GrevLex,
    DegLex,
    DegRevLex
};

inline int total_degree(const Monomial& m) {
    return std::accumulate(m.begin(), m.end(), 0);
}

inline Monomial lcm_monomial(const Monomial& a, const Monomial& b) {
    Monomial result(std::max(a.size(), b.size()), 0);
    for (size_t i = 0; i < a.size(); ++i)
        result[i] = std::max(result[i], a[i]);
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = std::max(result[i], b[i]);
    return result;
}

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

class MonomialOrder {
public:
    explicit MonomialOrder(MonomialOrderType type) : type_(type) {}

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

    static MonomialOrder lex() { return MonomialOrder(MonomialOrderType::Lex); }

    static MonomialOrder grevlex() { return MonomialOrder(MonomialOrderType::GrevLex); }

    static MonomialOrder deglex() { return MonomialOrder(MonomialOrderType::DegLex); }

    static MonomialOrder degrevlex() { return MonomialOrder(MonomialOrderType::DegRevLex); }

private:
    MonomialOrderType type_;

    static inline bool compare_lex(const Monomial& a, const Monomial& b) {
        size_t n = std::max(a.size(), b.size());
        for (size_t i = 0; i < n; ++i) {
            int ai = (i < a.size()) ? a[i] : 0;
            int bi = (i < b.size()) ? b[i] : 0;
            if (ai != bi) return ai > bi;
        }
        return false;
    }

    static inline bool compare_grevlex(const Monomial& a, const Monomial& b) {
        int deg_a = total_degree(a);
        int deg_b = total_degree(b);
        if (deg_a != deg_b) return deg_a > deg_b;

        size_t n = std::max(a.size(), b.size());
        for (size_t i = n; i > 0; --i) {
            int ai = (i - 1 < a.size()) ? a[i - 1] : 0;
            int bi = (i - 1 < b.size()) ? b[i - 1] : 0;
            if (ai != bi) return ai < bi;
        }
        return false;
    }

    static inline bool compare_deglex(const Monomial& a, const Monomial& b) {
        int deg_a = total_degree(a);
        int deg_b = total_degree(b);
        if (deg_a != deg_b) return deg_a > deg_b;
        return compare_lex(a, b);
    }
};

}
