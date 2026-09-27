#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include "bigint.hpp"
#include "lmmc/init.h"
#include <lmmp.h>

#include <atomic>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cwchar>
#include <exception>
#include <limits>
#include <new>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "%s (Win32 error %lu)\n", message,
                     static_cast<unsigned long>(GetLastError()));
        ExitProcess(1);
    }
}

struct Accounting {
    DWORD owner = 0;
    std::atomic<std::size_t> allocations{0};
    std::atomic<std::size_t> outstanding{0};
    std::atomic<bool> cross_thread{false};
};

struct alignas(std::max_align_t) Allocation {
    Accounting* account;
};

thread_local Accounting* current_account = nullptr;

void* counted_alloc(std::size_t bytes) {
    if (bytes > std::numeric_limits<std::size_t>::max() - sizeof(Allocation)) {
        return nullptr;
    }
    void* memory = std::malloc(sizeof(Allocation) + bytes);
    if (!memory) {
        return nullptr;
    }
    require(current_account != nullptr, "Allocation without a thread account");
    auto* allocation = new (memory) Allocation{current_account};
    ++current_account->allocations;
    ++current_account->outstanding;
    return allocation + 1;
}

void counted_free(void* pointer) {
    if (!pointer) {
        return;
    }
    auto* allocation = static_cast<Allocation*>(pointer) - 1;
    Accounting& account = *allocation->account;
    if (account.owner != GetCurrentThreadId()) {
        account.cross_thread = true;
    }
    --account.outstanding;
    std::free(allocation);
}

void* counted_realloc(void* pointer, std::size_t bytes) {
    if (!pointer) {
        return counted_alloc(bytes);
    }
    if (bytes > std::numeric_limits<std::size_t>::max() - sizeof(Allocation)) {
        return nullptr;
    }
    auto* allocation = static_cast<Allocation*>(pointer) - 1;
    Accounting* account = allocation->account;
    if (account->owner != GetCurrentThreadId()) {
        account->cross_thread = true;
    }
    void* memory = std::realloc(allocation, sizeof(Allocation) + bytes);
    if (!memory) {
        return nullptr;
    }
    auto* resized = new (memory) Allocation{account};
    ++account->allocations;
    return resized + 1;
}

void install_allocator(Accounting& account) {
    account.owner = GetCurrentThreadId();
    current_account = &account;
    const lmmp_heap_allocator_t allocator{counted_alloc, counted_free, counted_realloc};
    lmmp_set_heap_allocator(&allocator);
}

void check_released(const Accounting& account) {
    require(account.allocations != 0, "Allocator did not observe this CAS thread");
    require(account.outstanding == 0, "Thread lease leaked LMMP allocations");
    require(!account.cross_thread, "Allocation released by a different thread");
}

void restore_allocator() {
    const lmmp_heap_allocator_t allocator{std::malloc, std::free, std::realloc};
    lmmp_set_heap_allocator(&allocator);
    lmmp_global_deinit();
    current_account = nullptr;
}

// Registered before the function-local lifecycle slot, so its destructor observes
// the main thread's actual lease release rather than releasing it for the test.
Accounting main_account;
struct MainThreadObservation {
    bool active = false;
    ~MainThreadObservation() {
        if (active) {
            check_released(main_account);
            restore_allocator();
        }
    }
} main_observation;

int internal_probe() {
    const LMCAS::BigInt value("12345678901234567890");
    const LMCAS::BigInt product = value * LMCAS::BigInt(9);
    return product.to_string() == "111111110111111111010" ? 0 : 1;
}

using Probe = int (*)();

struct Worker {
    Probe probe;
    Accounting* account;
    HANDLE ready;
    HANDLE release;
    HANDLE thread;
    int result = 0;

    static DWORD WINAPI run(void* parameter) {
        auto& worker = *static_cast<Worker*>(parameter);
        if (worker.account) {
            install_allocator(*worker.account);
        }
        if (worker.probe) {
            worker.result = worker.probe();
        }
        if (worker.account) {
            require(worker.account->allocations != 0, "Worker CAS allocation was not observed");
        }
        require(SetEvent(worker.ready) != FALSE, "Cannot signal worker readiness");
        require(WaitForSingleObject(worker.release, INFINITE) == WAIT_OBJECT_0,
                "Worker wait failed");
        return static_cast<DWORD>(worker.result);
    }

    explicit Worker(Probe function, Accounting* observation = nullptr)
        : probe(function), account(observation),
          ready(CreateEventW(nullptr, TRUE, FALSE, nullptr)),
          release(CreateEventW(nullptr, TRUE, FALSE, nullptr)), thread(nullptr) {
        require(ready && release, "Cannot create worker events");
        thread = CreateThread(nullptr, 0, run, this, 0, nullptr);
        require(thread != nullptr, "Cannot create worker thread");
        require(WaitForSingleObject(ready, 10000) == WAIT_OBJECT_0,
                "Worker did not complete its CAS computation");
        require(result == 0, "Worker arithmetic probe failed");
    }

    void join() {
        require(SetEvent(release) != FALSE, "Cannot release worker");
        require(WaitForSingleObject(thread, 10000) == WAIT_OBJECT_0, "Worker did not exit");
        DWORD exit_code = 1;
        require(GetExitCodeThread(thread, &exit_code) != FALSE && exit_code == 0,
                "Worker exited with an error");
        require(CloseHandle(thread) != FALSE, "Cannot close worker thread");
        require(CloseHandle(release) != FALSE, "Cannot close worker release event");
        require(CloseHandle(ready) != FALSE, "Cannot close worker ready event");
        if (account) {
            check_released(*account);
        }
    }
};

struct Module {
    HMODULE handle;
    Probe probe;

    explicit Module(const wchar_t* path) : handle(LoadLibraryW(path)), probe(nullptr) {
        require(handle != nullptr, "Cannot load lifecycle module");
        const FARPROC address = GetProcAddress(handle, "lmcas_lifecycle_probe");
        require(address != nullptr, "Lifecycle module has no probe");
        static_assert(sizeof(address) == sizeof(probe));
        std::memcpy(&probe, &address, sizeof(probe));
    }

    void unload(const wchar_t* lmcas_path) {
        require(FreeLibrary(handle) != FALSE, "Cannot unload lifecycle module");
        require(GetModuleHandleW(lmcas_path) == nullptr, "LMCAS remains loaded after FreeLibrary");
    }
};

void unload_round(const wchar_t* module_path, const wchar_t* lmcas_path,
                  bool use_main, bool non_cas_worker) {
    Module module(module_path);
    Accounting main_thread;
    Accounting worker_thread;
    if (use_main) {
        install_allocator(main_thread);
        require(module.probe() == 0, "Main module arithmetic probe failed");
    }
    Worker worker(non_cas_worker ? nullptr : module.probe,
                  non_cas_worker ? nullptr : &worker_thread);
    if (!non_cas_worker) {
        worker.join();
    }
    // Every CAS participant except this unloading thread has exited. A worker
    // that has never touched CAS is allowed to remain alive during FreeLibrary.
    module.unload(lmcas_path);
    if (use_main) {
        check_released(main_thread);
        restore_allocator();
    }
    if (non_cas_worker) {
        worker.join();
    }
}

int child(const wchar_t* scenario, const wchar_t* module_path, const wchar_t* lmcas_path) {
    std::set_terminate([] { ExitProcess(77); });
    SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
    if (std::wcscmp(scenario, L"internal-single") == 0) {
        main_observation.active = true;
        install_allocator(main_account);
        for (int iteration = 0; iteration < 2; ++iteration) {
            require(internal_probe() == 0, "Internal arithmetic probe failed");
            require(lmmc_stack_reset(327680) == LMMC_STATUS_OK,
                    "Internal reentry acquired more than one lease");
        }
        return 0;
    }
    if (std::wcscmp(scenario, L"internal-joined") == 0) {
        Accounting account;
        Worker worker(internal_probe, &account);
        worker.join();
        main_observation.active = true;
        install_allocator(main_account);
        require(internal_probe() == 0, "Main arithmetic after worker exit failed");
        require(lmmc_stack_reset(327680) == LMMC_STATUS_OK, "Main lease is not unique");
        return 0;
    }
    if (std::wcscmp(scenario, L"internal-main-live-worker") == 0 ||
        std::wcscmp(scenario, L"internal-worker-only") == 0) {
        if (std::wcscmp(scenario, L"internal-main-live-worker") == 0) {
            require(internal_probe() == 0, "Main internal arithmetic probe failed");
        }
        Worker worker(internal_probe);
        // Preserve the worker's wait state until normal CRT process teardown.
        // Its handles and thread-local resources belong to process reclamation.
        std::exit(0);
    }
    if (std::wcscmp(scenario, L"module-joined-unload") == 0) {
        unload_round(module_path, lmcas_path, false, false);
        return 0;
    }
    if (std::wcscmp(scenario, L"module-live-worker") == 0) {
        Module module(module_path);
        Worker worker(module.probe);
        std::exit(0);
    }
    if (std::wcscmp(scenario, L"module-reload") == 0) {
        for (int iteration = 0; iteration < 256; ++iteration) {
            unload_round(module_path, lmcas_path, true, false);
        }
        return 0;
    }
    if (std::wcscmp(scenario, L"module-non-cas-worker") == 0) {
        unload_round(module_path, lmcas_path, true, true);
        return 0;
    }
    return 1;
}

std::wstring quote_argument(const std::wstring& argument) {
    std::wstring quoted = L"\"";
    std::size_t backslashes = 0;
    for (wchar_t character : argument) {
        if (character == L'\\') {
            ++backslashes;
        } else {
            quoted.append(backslashes * (character == L'"' ? 2 : 1), L'\\');
            backslashes = 0;
            if (character == L'"') {
                quoted += L'\\';
            }
            quoted += character;
        }
    }
    quoted.append(backslashes * 2, L'\\');
    quoted += L'"';
    return quoted;
}

bool run_child(const std::wstring& executable, const wchar_t* scenario,
               const wchar_t* module_path, const wchar_t* lmcas_path) {
    std::wstring command = quote_argument(executable) + L" --child " +
        quote_argument(scenario) + L" " + quote_argument(module_path) + L" " +
        quote_argument(lmcas_path);
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process{};
    if (!CreateProcessW(executable.c_str(), command.data(), nullptr, nullptr, FALSE,
                        0, nullptr, nullptr, &startup, &process)) {
        std::fwprintf(stderr, L"Cannot start scenario %ls: %lu\n", scenario,
                      static_cast<unsigned long>(GetLastError()));
        return false;
    }
    CloseHandle(process.hThread);
    const DWORD wait = WaitForSingleObject(process.hProcess, 30000);
    DWORD exit_code = 1;
    const bool completed = wait == WAIT_OBJECT_0 &&
        GetExitCodeProcess(process.hProcess, &exit_code) != FALSE;
    if (!completed) {
        require(TerminateProcess(process.hProcess, 1) != FALSE, "Cannot terminate stalled child");
        require(WaitForSingleObject(process.hProcess, 5000) == WAIT_OBJECT_0,
                "Terminated child did not exit");
    }
    CloseHandle(process.hProcess);
    std::fwprintf(completed && exit_code == 0 ? stdout : stderr,
                  L"%ls: %ls (exit %lu)\n", scenario,
                  completed ? (exit_code == 0 ? L"passed" : L"failed") : L"timed out or wait failed",
                  static_cast<unsigned long>(exit_code));
    // 77 means terminate in a Windows child, never a CTest platform skip.
    return completed && exit_code == 0;
}

} // namespace

int wmain(int argc, wchar_t** argv) {
    if (argc == 5 && std::wcscmp(argv[1], L"--child") == 0) {
        return child(argv[2], argv[3], argv[4]);
    }
    if (argc != 3) {
        return 1;
    }
    std::vector<wchar_t> path(32768);
    const DWORD length = GetModuleFileNameW(nullptr, path.data(), static_cast<DWORD>(path.size()));
    require(length != 0 && length < path.size(), "Cannot find lifecycle driver executable");
    const std::wstring executable(path.data(), length);
    const wchar_t* scenarios[] = {
        L"internal-single", L"internal-joined", L"internal-main-live-worker",
        L"internal-worker-only", L"module-joined-unload", L"module-live-worker",
        L"module-reload", L"module-non-cas-worker"
    };
    for (const wchar_t* scenario : scenarios) {
        if (!run_child(executable, scenario, argv[1], argv[2])) {
            return 1;
        }
    }
    return 0;
}
#else
int main() {
    return 77;
}
#endif
