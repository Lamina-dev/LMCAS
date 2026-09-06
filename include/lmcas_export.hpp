#pragma once

#ifndef LMCAS_API
#if defined(_WIN32) || defined(__CYGWIN__)
#if defined(LMCAS_EXPORTS)
#define LMCAS_API __declspec(dllexport)
#else
#define LMCAS_API __declspec(dllimport)
#endif
#else
#define LMCAS_API __attribute__((visibility("default")))
#endif
#endif
