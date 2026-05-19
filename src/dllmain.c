/**
 * @file dllmain.c
 * @brief DLL entry point for the lmcas shared library.
 *
 * On process attach/detach we initialise / tear down the LAMMP
 * global resources via the LMMC lifecycle functions.  This keeps
 * the LAMMP initialisation concern inside the library boundary
 * (LMMC) rather than leaking into application or CAS code.
 */

#ifdef _WIN32
#include <windows.h>
#include "lmmc/init.h"

BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
    (void)hinstDLL;
    (void)lpvReserved;

    switch (fdwReason)
    {
    case DLL_PROCESS_ATTACH:
        lmmc_init();
        break;

    case DLL_PROCESS_DETACH:
        lmmc_deinit();
        break;

    case DLL_THREAD_ATTACH:
        /* new LAMMP may need per-thread lmmp_global_init_() here */
        break;

    case DLL_THREAD_DETACH:
        /* new LAMMP may need per-thread lmmp_global_deinit() here */
        break;
    }

    return TRUE;
}

#else
/* On non-Windows platforms we rely on explicit lmmc_init() calls
   at the application level or __attribute__((constructor)) in the
   future.  Add support as needed. */
#endif
