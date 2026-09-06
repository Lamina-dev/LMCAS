#include "internal/lmmc_lifecycle.hpp"

#include "lmmc/init.h"

#include <exception>

namespace LMCAS::detail {
namespace {

class LmmcLifecycle final {
public:
    LmmcLifecycle() {
        if (lmmc_init() != LMMC_STATUS_OK) {
            std::terminate();
        }
    }

    ~LmmcLifecycle() {
        if (lmmc_deinit() != LMMC_STATUS_OK) {
            std::terminate();
        }
    }

    LmmcLifecycle(const LmmcLifecycle&) = delete;
    LmmcLifecycle& operator=(const LmmcLifecycle&) = delete;
};

} // namespace

void ensure_lmmc_lifecycle() noexcept {
    static thread_local LmmcLifecycle lifecycle;
    (void)lifecycle;
}

} // namespace LMCAS::detail
