#include "bigint.hpp"

#ifdef _WIN32
#define LIFECYCLE_EXPORT __declspec(dllexport)
#else
#define LIFECYCLE_EXPORT __attribute__((visibility("default")))
#endif

extern "C" LIFECYCLE_EXPORT int lmcas_lifecycle_probe(void) {
    const LMCAS::BigInt value("12345678901234567890");
    const LMCAS::BigInt product = value * LMCAS::BigInt(9);
    return product.to_string() == "111111110111111111010" ? 0 : 1;
}
