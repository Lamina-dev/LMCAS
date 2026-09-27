#include "internal/mixed_transcendental_support.hpp"

#include <algorithm>
#include <cmath>

namespace LMCAS::detail {
namespace {

constexpr lmmc_real_t minimum_width = 1e-6;
constexpr int initial_divisions = 64;

struct SampledInterval {
    lmmc_real_t lo;
    lmmc_real_t hi;
    lmmc_real_t f_lo;
    lmmc_real_t f_hi;
};

class RootIsolation {
public:
    RootIsolation(const std::shared_ptr<SymbolicExpr>& expr,
                  const std::shared_ptr<SymbolicExpr>& derivative,
                  const std::string& var)
        : expr_(expr), derivative_(derivative), var_(var) {}

    bool initialize(const SearchInterval& interval, ComputationContext* context,
                    std::optional<CasError>* failure, bool* complete) {
        const auto step = (interval.hi - interval.lo) / initial_divisions;
        std::vector<lmmc_real_t> samples(initial_divisions + 1);
        for (int i = 0; i <= initial_divisions; ++i) {
            if (!consume_mixed_step(context, failure, complete)) {
                return false;
            }
            samples[i] = sample(interval.lo + i * step, step * 0.5);
            if (samples[i] == 0.0) {
                const auto x = interval.lo + i * step;
                if (mixed_evaluate_at(expr_, var_, x) == 0.0)
                    result_.push_back({x, x, true});
            }
        }
        for (int i = 0; i < initial_divisions; ++i) {
            const auto fa = samples[i];
            const auto fb = samples[i + 1];
            if (std::isnan(fa) && std::isnan(fb)) {
                continue;
            }
            if (std::isnan(fa) || std::isnan(fb) || fa * fb < 0.0) {
                work_.push_back({interval.lo + i * step,
                                 interval.lo + (i + 1) * step, fa, fb});
            }
        }
        return true;
    }

    void run(const SolveOptions& opts, ComputationContext* context,
             std::optional<CasError>* failure, bool* complete) {
        while (!work_.empty()) {
            if (!consume_mixed_step(context, failure, complete)) {
                break;
            }
            if (opts.max_roots > 0 &&
                result_.size() >= static_cast<std::size_t>(opts.max_roots)) {
                if (complete) {
                    *complete = false;
                }
                break;
            }
            auto current = work_.back();
            work_.pop_back();
            process(current);
        }
    }

    std::vector<IsolatedInterval> finish(int max_roots, bool* complete) {
        std::sort(result_.begin(), result_.end(),
                  [](const IsolatedInterval& a, const IsolatedInterval& b) {
                      return a.lo < b.lo;
                  });
        if (max_roots > 0 && result_.size() > static_cast<std::size_t>(max_roots)) {
            if (complete) {
                *complete = false;
            }
            result_.resize(static_cast<std::size_t>(max_roots));
        }
        return std::move(result_);
    }

private:
    lmmc_real_t sample(lmmc_real_t x, lmmc_real_t offset) const {
        return mixed_evaluate_with_retry(expr_, var_, x, offset);
    }

    void accept(const SampledInterval& interval, bool confirmed) {
        result_.push_back({interval.lo, interval.hi, confirmed});
    }

    bool resample_endpoints(SampledInterval& interval) {
        const auto width = interval.hi - interval.lo;
        if (std::isnan(interval.f_lo)) {
            interval.f_lo = sample(interval.lo, width * 0.25);
            if (std::isnan(interval.f_lo)) {
                if (width > minimum_width * 2.0) {
                    const auto mid = (interval.lo + interval.hi) * 0.5;
                    const auto fm = sample(mid, width * 0.125);
                    work_.push_back({mid, interval.hi, fm, interval.f_hi});
                }
                return false;
            }
        }
        if (std::isnan(interval.f_hi)) {
            interval.f_hi = sample(interval.hi, width * 0.25);
            if (std::isnan(interval.f_hi)) {
                if (width > minimum_width * 2.0) {
                    const auto mid = (interval.lo + interval.hi) * 0.5;
                    const auto fm = sample(mid, width * 0.125);
                    work_.push_back({interval.lo, mid, interval.f_lo, fm});
                }
                return false;
            }
        }
        return true;
    }

    void subdivide(const SampledInterval& interval) {
        const auto mid = (interval.lo + interval.hi) * 0.5;
        const auto fm = sample(mid, (interval.hi - interval.lo) * 0.125);
        if (std::isnan(fm)) {
            accept(interval, false);
            return;
        }
        if (fm == 0.0 && mixed_evaluate_at(expr_, var_, mid) == 0.0) {
            result_.push_back({mid, mid, true});
        }
        if (interval.f_lo * fm < 0.0) {
            work_.push_back({interval.lo, mid, interval.f_lo, fm});
        }
        if (fm * interval.f_hi < 0.0) {
            work_.push_back({mid, interval.hi, fm, interval.f_hi});
        }
        if (interval.f_lo * fm >= 0.0 && fm * interval.f_hi >= 0.0) {
            accept(interval, false);
        }
    }

    void process(const SampledInterval& current) {
        auto interval = current;
        if (!resample_endpoints(interval) || interval.f_lo * interval.f_hi >= 0.0) {
            return;
        }
        const auto width = interval.hi - interval.lo;
        if (!derivative_) {
            if (width > minimum_width * 4.0) {
                subdivide(interval);
            } else {
                accept(interval, false);
            }
            return;
        }
        const auto da = mixed_evaluate_at(derivative_, var_, interval.lo);
        const auto db = mixed_evaluate_at(derivative_, var_, interval.hi);
        const bool finite = !std::isnan(da) && !std::isnan(db);
        const bool monotone = finite && da * db > 0.0;
        if (width <= minimum_width || monotone) {
            accept(interval, monotone);
            return;
        }
        if (finite && da * db < 0.0 && width > minimum_width * 2.0) {
            subdivide(interval);
        } else {
            accept(interval, false);
        }
    }

    const std::shared_ptr<SymbolicExpr>& expr_;
    const std::shared_ptr<SymbolicExpr>& derivative_;
    const std::string& var_;
    std::vector<SampledInterval> work_;
    std::vector<IsolatedInterval> result_;
};

}

std::vector<IsolatedInterval> isolate_roots_with_context(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var, const SearchInterval& interval,
    const SolveOptions& opts, ComputationContext* context,
    std::optional<CasError>* failure, bool* complete) {
    if (interval.lo >= interval.hi) {
        return {};
    }
    RootIsolation isolation(expr, derivative, var);
    if (!isolation.initialize(interval, context, failure, complete)) {
        return {};
    }
    isolation.run(opts, context, failure, complete);
    return isolation.finish(opts.max_roots, complete);
}

}

namespace LMCAS {

std::vector<IsolatedInterval> isolate_roots(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var, const SearchInterval& interval,
    const SolveOptions& opts) {
    return detail::isolate_roots_with_context(
        expr, derivative, var, interval, opts, nullptr, nullptr, nullptr);
}

}
