#include "internal/lmmc_lifecycle.hpp"

#include "lmmc/init.h"

#include <cstdint>
#include <exception>

#ifdef _WIN32
#include <windows.h>
#endif

namespace LMCAS::detail {
namespace {

#ifdef _WIN32
/**
 * @brief 使用原生 FLS，在 LMMP/LMMC 线程状态有效时清理，兼容内部测试链接。
 * @note 可执行文件的 C++ TLS 析构可能晚于依赖的 MinGW DLL 仿真 TLS 销毁。
 */
// FlsFree can invoke callbacks for other threads. Only thread exit or the
// unloading thread may release its lease; live foreign threads are reclaimed
// by process exit. Explicit unload requires all other CAS threads to be joined.
void CALLBACK release_lmmc_lifecycle(void* owner) {
    if (reinterpret_cast<std::uintptr_t>(owner) == GetCurrentThreadId() &&
        lmmc_deinit() != LMMC_STATUS_OK) {
        std::terminate();
    }
}

class LmmcLifecycleSlot final {
public:
    LmmcLifecycleSlot() : index_(FlsAlloc(release_lmmc_lifecycle)) {
        if (index_ == FLS_OUT_OF_INDEXES) {
            std::terminate();
        }
    }

    ~LmmcLifecycleSlot() {
        if (void* owner = FlsGetValue(index_)) {
            if (!FlsSetValue(index_, nullptr)) {
                std::terminate();
            }
            release_lmmc_lifecycle(owner);
        }
        if (!FlsFree(index_)) {
            std::terminate();
        }
    }

    DWORD index() const noexcept { return index_; }

    LmmcLifecycleSlot(const LmmcLifecycleSlot&) = delete;
    LmmcLifecycleSlot& operator=(const LmmcLifecycleSlot&) = delete;

private:
    DWORD index_;
};
#else
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
#endif

}

void ensure_lmmc_lifecycle() noexcept {
#ifdef _WIN32
    static LmmcLifecycleSlot slot;
    if (FlsGetValue(slot.index())) {
        return;
    }
    if (lmmc_init() != LMMC_STATUS_OK) {
        std::terminate();
    }
    void* owner = reinterpret_cast<void*>(static_cast<std::uintptr_t>(GetCurrentThreadId()));
    if (!FlsSetValue(slot.index(), owner)) {
        if (lmmc_deinit() != LMMC_STATUS_OK) {
            std::terminate();
        }
        std::terminate();
    }
#else
    static thread_local LmmcLifecycle lifecycle;
    (void)lifecycle;
#endif
}

}
