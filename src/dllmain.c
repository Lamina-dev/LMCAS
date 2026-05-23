#ifdef _WIN32
#include <windows.h>
#include "lmmc/init.h"

BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
    (void)hinstDLL;

    switch (fdwReason)
    {
    case DLL_PROCESS_ATTACH:

        lmmc_init();
        break;

    case DLL_PROCESS_DETACH:

        if (lpvReserved == NULL) {
            lmmc_deinit();
        }
        break;

    case DLL_THREAD_ATTACH:

        break;

    case DLL_THREAD_DETACH:

        break;
    }

    return TRUE;
}

#else

#endif
