#include "solve_mixed_transcendental.hpp"
#include <algorithm>
#include <cmath>

namespace LMCAS {

std::vector<lmmc_real_t> deduplicate_roots(
    std::vector<NumericRoot>& roots,
    lmmc_real_t tolerance,
    int max_roots)
{
    if (roots.empty()) {
        return {};
    }

    std::sort(roots.begin(), roots.end(),
        [](const NumericRoot& a, const NumericRoot& b) {
            return a.value < b.value;
        });
    lmmc_real_t threshold = 10.0 * tolerance;
    std::vector<NumericRoot> survivors;
    survivors.reserve(roots.size());
    survivors.push_back(roots[0]);

    for (size_t i = 1; i < roots.size(); ++i) {
        NumericRoot& last = survivors.back();
        if (std::fabs(roots[i].value - last.value) < threshold) {
            if (roots[i].residual < last.residual) {
                last = roots[i];
            }
        } else {
            survivors.push_back(roots[i]);
        }
    }
    std::vector<lmmc_real_t> result;
    result.reserve(survivors.size());
    for (const auto& r : survivors) {
        result.push_back(r.value);
    }

    if (max_roots > 0 && result.size() > static_cast<std::size_t>(max_roots)) {
        result.resize(static_cast<size_t>(max_roots));
    }

    return result;
}

}
