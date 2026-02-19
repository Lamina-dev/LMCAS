#include "lammp/lmmp.h"
#include <string.h>

// Mock implementation of LAMMP functions for testing PrintVisitor

void lmmp_mul_(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn) {
    // Check for null pointers if needed, but for now do nothing
}

void lmmp_div_(mp_ptr q, mp_ptr r, mp_srcptr n, mp_size_t nn, mp_srcptr d, mp_size_t dn) {
    // Do nothing
}

// mp_byte_t is uint8_t
mp_size_t lmmp_from_str_(mp_ptr res, const uint8_t* str, mp_size_t len, int base) {
    if (res) res[0] = 1;
    return 1;
}

// Mock implementation that simply returns last limb value % 10 as single digit for verification
// This allows differentiating numbers like 1, 2, ...
mp_size_t lmmp_to_str_(uint8_t* str, mp_srcptr n, mp_size_t nn, int base) {
    if (str && n && nn > 0) {
        str[0] = (uint8_t)(n[0] % 10); // Simple hack: last digit of first limb
        return 1;
    }
    if (str) str[0] = 0;
    return 1;
}

mp_limb_t lmmp_add_n_(mp_ptr res, mp_srcptr s1, mp_srcptr s2, mp_size_t n) {
    return 0;
}

mp_limb_t lmmp_sub_n_(mp_ptr res, mp_srcptr s1, mp_srcptr s2, mp_size_t n) {
    return 0;
}

mp_limb_t lmmp_add_nc_(mp_ptr res, mp_srcptr s1, mp_srcptr s2, mp_size_t n, mp_limb_t carry) {
    return 0;
}

mp_size_t lmmp_pow_size_(mp_srcptr b, mp_size_t bn, unsigned long exp) {
    return 1;
}

mp_size_t lmmp_pow_(mp_ptr dst, mp_size_t rn, mp_srcptr base, mp_size_t n, unsigned long exp) {
    if (dst && rn > 0) dst[0] = 1; 
    return 1;
}

mp_size_t lmmp_gcd_lehmer_(mp_ptr dst, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn) {
    if (dst) dst[0] = 1;
    return 1;
}
