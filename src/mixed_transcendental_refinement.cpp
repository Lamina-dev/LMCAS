#include "internal/mixed_transcendental_support.hpp"
#include "internal/solver_support.hpp"

#include <cmath>
#include <limits>

namespace LMCAS::detail {
namespace {

class RootRefinement {
public:
    RootRefinement(const std::shared_ptr<SymbolicExpr>& expr,
                   const std::shared_ptr<SymbolicExpr>& derivative,
                   const std::string& var, const IsolatedInterval& interval,
                   const SolveOptions& opts)
        : expr_(expr), derivative_(derivative), var_(var), opts_(opts),
          lo_(interval.lo), hi_(interval.hi), x_(midpoint()),
          f_lo_(evaluate(lo_)), f_hi_(evaluate(hi_)),
          best_{x_, std::numeric_limits<lmmc_real_t>::max(), 0} {}

    std::optional<NumericRoot> endpoint_root() const {
        auto lower = root_at(lo_, f_lo_, 0);
        if (lower) {
            return lower;
        }
        return root_at(hi_, f_hi_, 0);
    }

    std::optional<NumericRoot> bisect(
        ComputationContext* context, std::optional<CasError>* failure,
        bool* complete) {
        if (std::isnan(f_lo_) || std::isnan(f_hi_) || f_lo_ * f_hi_ > 0.0) {
            return root_at(x_, evaluate(x_), 1);
        }
        for (int i = 1; i <= opts_.max_newton_iterations; ++i) {
            if (!consume_mixed_step(context, failure, complete)) {
                return std::nullopt;
            }
            const auto mid = midpoint();
            const auto f_mid = evaluate(mid);
            if (std::isnan(f_mid)) {
                hi_ = mid;
                continue;
            }
            remember(mid, f_mid, i);
            auto root = root_at(mid, f_mid, i);
            if (root) {
                return root;
            }
            if (f_lo_ * f_mid < 0.0) {
                hi_ = mid;
            } else {
                lo_ = mid;
                f_lo_ = f_mid;
            }
            if (hi_ - lo_ < opts_.tolerance) {
                const auto final_x = midpoint();
                remember(final_x, evaluate(final_x), i);
                break;
            }
        }
        return best_root();
    }

    std::optional<NumericRoot> newton(
        ComputationContext* context, std::optional<CasError>* failure,
        bool* complete) {
        if (!resample_endpoints()) {
            return std::nullopt;
        }
        for (int i = 1; i <= opts_.max_newton_iterations; ++i) {
            if (!consume_mixed_step(context, failure, complete)) {
                return std::nullopt;
            }
            const auto fx = evaluate(x_);
            if (std::isnan(fx)) {
                x_ = midpoint();
                continue;
            }
            remember(x_, fx, i);
            auto root = root_at(x_, fx, i);
            if (root) {
                return root;
            }
            const auto dfx = mixed_evaluate_at(derivative_, var_, x_);
            if (std::isnan(dfx) || std::fabs(dfx) < 1e-15) {
                bisection_step();
                continue;
            }
            if (accept_newton_step(fx, dfx)) {
                continue;
            }
            bisection_step();
            if (hi_ - lo_ < opts_.tolerance) {
                const auto final_x = midpoint();
                root = root_at(final_x, evaluate(final_x), i);
                if (root) {
                    return root;
                }
                break;
            }
        }
        return best_root();
    }

private:
    lmmc_real_t midpoint() const {
        return solver_detail::finite_midpoint(lo_, hi_);
    }

    lmmc_real_t evaluate(lmmc_real_t x) const {
        return mixed_evaluate_at(expr_, var_, x);
    }

    std::optional<NumericRoot> root_at(
        lmmc_real_t x, lmmc_real_t fx, int iterations) const {
        if (!std::isnan(fx) && std::fabs(fx) < opts_.tolerance) {
            return NumericRoot{x, std::fabs(fx), iterations};
        }
        return std::nullopt;
    }

    void remember(lmmc_real_t x, lmmc_real_t fx, int iterations) {
        const auto residual = std::fabs(fx);
        if (!std::isnan(residual) && residual < best_.residual) {
            best_ = {x, residual, iterations};
        }
    }

    std::optional<NumericRoot> best_root() const {
        if (best_.residual < opts_.tolerance) {
            return best_;
        }
        return std::nullopt;
    }

    bool resample_endpoints() {
        if (std::isnan(f_lo_)) {
            f_lo_ = evaluate(lo_ + (hi_ - lo_) * 0.01);
            lo_ = lo_ + (hi_ - lo_) * 0.01;
        }
        if (std::isnan(f_hi_)) {
            f_hi_ = evaluate(hi_ - (hi_ - lo_) * 0.01);
            hi_ = hi_ - (hi_ - lo_) * 0.01;
        }
        return !std::isnan(f_lo_) && !std::isnan(f_hi_);
    }

    void bisection_step() {
        const auto mid = midpoint();
        const auto f_mid = evaluate(mid);
        if (!std::isnan(f_mid)) {
            if (!std::isnan(f_lo_) && f_lo_ * f_mid < 0.0) {
                hi_ = mid;
                f_hi_ = f_mid;
            } else {
                lo_ = mid;
                f_lo_ = f_mid;
            }
        }
        x_ = midpoint();
    }

    bool accept_newton_step(lmmc_real_t fx, lmmc_real_t dfx) {
        const auto next = x_ - fx / dfx;
        if (next < lo_ || next > hi_) {
            return false;
        }
        const auto f_next = evaluate(next);
        if (std::isnan(f_next) || std::fabs(f_next) > std::fabs(fx)) {
            return false;
        }
        if (!std::isnan(f_lo_) && fx * f_lo_ < 0.0) {
            hi_ = x_;
            f_hi_ = fx;
        } else if (!std::isnan(f_hi_) && fx * f_hi_ < 0.0) {
            lo_ = x_;
            f_lo_ = fx;
        }
        x_ = next;
        return true;
    }

    const std::shared_ptr<SymbolicExpr>& expr_;
    const std::shared_ptr<SymbolicExpr>& derivative_;
    const std::string& var_;
    const SolveOptions& opts_;
    lmmc_real_t lo_;
    lmmc_real_t hi_;
    lmmc_real_t x_;
    lmmc_real_t f_lo_;
    lmmc_real_t f_hi_;
    NumericRoot best_;
};

}

std::optional<NumericRoot> refine_root_with_context(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var, const IsolatedInterval& interval,
    const SolveOptions& opts, ComputationContext* context,
    std::optional<CasError>* failure, bool* complete) {
    RootRefinement refinement(expr, derivative, var, interval, opts);
    auto endpoint = refinement.endpoint_root();
    if (endpoint) {
        return endpoint;
    }
    if (!derivative) {
        return refinement.bisect(context, failure, complete);
    }
    return refinement.newton(context, failure, complete);
}

}

namespace LMCAS {

std::optional<NumericRoot> refine_root(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& var, const IsolatedInterval& interval,
    const SolveOptions& opts) {
    return detail::refine_root_with_context(
        expr, derivative, var, interval, opts, nullptr, nullptr, nullptr);
}

}
