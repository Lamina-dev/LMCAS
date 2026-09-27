#include "internal/fglm_internal.hpp"
#include <algorithm>

namespace LMCAS::fglm_detail {

    MonomialEnumerator::MonomialEnumerator(size_t num_vars, const MonomialOrder& order, int max_deg)
        : num_vars_(num_vars), order_(order), max_degree_(max_deg) {

        Monomial one(num_vars, 0);
        push(one);
    }

    bool MonomialEnumerator::next(Monomial& out) {
        if (heap_.empty()) return false;
        out = heap_.front();
        std::pop_heap(heap_.begin(), heap_.end(), cmp_);
        heap_.pop_back();

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

    void MonomialEnumerator::push(const Monomial& m) {
        if (visited_.count(m)) return;
        visited_.insert(m);
        heap_.push_back(m);
        std::push_heap(heap_.begin(), heap_.end(), cmp_);
    }

    bool GaussianEliminator::add_vector(const std::vector<Rational>& v,
                    std::vector<Rational>& combination) {
        size_t n = v.size();

        for (auto& row : rows_) {
            row.resize(n, Rational(0));
        }

        std::vector<Rational> working = v;

        combination.assign(basis_count_, Rational(0));

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

            for (size_t j = 0; j < combinations_[i].size(); ++j) {
                combination[j] =
                    combination[j] - factor * combinations_[i][j];
            }
        }

        size_t pivot_col = n;
        for (size_t j = 0; j < n; ++j) {
            if (!working[j].is_zero()) {
                pivot_col = j;
                break;
            }
        }

        if (pivot_col == n) {

            for (auto& c : combination) {
                c = -c;
            }
            return false;
        }

        Rational pivot_val = working[pivot_col];
        for (size_t j = 0; j < n; ++j) {
            working[j] = working[j] / pivot_val;
        }

        std::vector<Rational> comb_row(basis_count_ + 1, Rational(0));
        for (size_t j = 0; j < combination.size(); ++j) {
            comb_row[j] = combination[j] / pivot_val;
        }

        comb_row[basis_count_] = Rational(1) / pivot_val;

        rows_.push_back(working);
        pivots_.push_back(pivot_col);
        combinations_.push_back(comb_row);
        ++basis_count_;

        return true;
    }

}
