#include "lmmc/init.h"

namespace {

struct LmmcLifecycle {
    LmmcLifecycle() {
        lmmc_init();
    }

    ~LmmcLifecycle() {
        lmmc_deinit();
    }
};

LmmcLifecycle lmmc_lifecycle;

} // namespace
