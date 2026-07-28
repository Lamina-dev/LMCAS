#pragma once

#ifndef LAMINA_API
#if defined(_WIN32) || defined(__CYGWIN__)
#if defined(LAMINA_CORE_EXPORTS)
#define LAMINA_API __declspec(dllexport)
#else
#define LAMINA_API __declspec(dllimport)
#endif
#else
#define LAMINA_API
#endif
#endif
