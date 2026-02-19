#include <iostream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstdint>

// Minimal mock of LAMMP/GMP for testing BigInt logic
// Only supports single-limb operations correctly for now.

extern "C" {

// Typedefs usually in lmmp.h
typedef uint64_t mp_limb_t;
typedef int64_t mp_size_t;
typedef mp_limb_t* mp_ptr;
typedef const mp_limb_t* mp_srcptr;
typedef enum { LMMP_ERROR_NONE = 0 } lmmp_error_t;

// r = a * b
void lmmp_mul_(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn) {
    if (!r) return;
    // Zero initialize result area
    for(mp_size_t i=0; i<an+bn; ++i) r[i] = 0;
    
    // Simple O(N*M) schoolbook multiplication
    // But since we only support single-limb inputs properly in this mock:
    unsigned long long v_a = (an > 0) ? a[0] : 0;
    unsigned long long v_b = (bn > 0) ? b[0] : 0;
    
    // 128-bit result of 64*64 mul
    unsigned long long res_lo = v_a * v_b; // simple 64-bit mul for now
    if (an > 0 && bn > 0) {
        r[0] = (mp_limb_t)res_lo;
    }
}

// q = n / d, r = n % d
void lmmp_div_(mp_ptr q, mp_ptr r, mp_srcptr n, mp_size_t nn, mp_srcptr d, mp_size_t dn) {
    unsigned long long v_n = (nn > 0) ? n[0] : 0;
    unsigned long long v_d = (dn > 0) ? d[0] : 0; 
    
    if (v_d == 0) return; 

    if (q) {
        for(mp_size_t i=0; i<nn-dn+1; ++i) q[i] = 0;
        q[0] = (mp_limb_t)(v_n / v_d);
    }
    if (r) {
        for(mp_size_t i=0; i<dn; ++i) r[i] = 0;
        r[0] = (mp_limb_t)(v_n % v_d);
    }
}

// Convert number to string (base 10)
mp_size_t lmmp_to_str_(uint8_t* str, mp_srcptr n, mp_size_t nn, int base) {
    if (!str || !n || nn == 0) return 0;
    
    unsigned long long val = n[0]; 
    if (val == 0) {
        str[0] = 0;
        return 1;
    }
    
    int len = 0;
    while(val > 0) {
        str[len++] = (uint8_t)(val % 10);
        val /= 10;
    }
    return len;
}

// Convert string to number
mp_size_t lmmp_from_str_(mp_ptr res, const uint8_t* str, mp_size_t len, int base) {
    if (res && str) {
       unsigned long long val = 0;
       unsigned long long mult = 1;
       for(int i=0; i<len; ++i) {
           val += str[i] * mult;
           mult *= 10;
       }
       res[0] = (mp_limb_t)val;
       return 1;
    }
    return 0;
}

// res = s1 + s2 (n limbs)
// Returns carry
mp_limb_t lmmp_add_n_(mp_ptr res, mp_srcptr s1, mp_srcptr s2, mp_size_t n) {
    mp_limb_t carry = 0;
    for(mp_size_t i=0; i<n; ++i) {
        unsigned long long a = s1[i];
        unsigned long long b = s2[i];
        unsigned long long sum = a + b + carry;
        res[i] = (mp_limb_t)sum;
        
        // Carry if wraparound occurred
        carry = (sum < a) || (carry && sum == a) ? 1 : 0; 
    }
    return carry;
}

// res = s1 - s2 (n limbs)
// Returns borrow
mp_limb_t lmmp_sub_n_(mp_ptr res, mp_srcptr s1, mp_srcptr s2, mp_size_t n) {
    mp_limb_t borrow = 0;
    for(mp_size_t i=0; i<n; ++i) {
        unsigned long long a = s1[i];
        unsigned long long b = s2[i];
        
        // a - b - borrow
        unsigned long long diff = a - b - borrow;
        res[i] = (mp_limb_t)diff;
        
        // Borrow logic: 
        // if a < b: borrow=1
        // if a == b: borrow unchanged? No.
        // Actually: if borrow_in=1:
        //    if a > b: no borrow out (e.g. 10 - 5 - 1 = 4)
        //    if a <= b: borrow out (e.g. 5 - 5 - 1 = -1 -> borrow)
        // if borrow_in=0:
        //    if a < b: borrow out
        bool new_borrow;
        if (borrow) new_borrow = (a <= b);
        else new_borrow = (a < b);
        
        borrow = new_borrow ? 1 : 0;
    }
    return borrow;
}

mp_limb_t lmmp_add_nc_(mp_ptr res, mp_srcptr s1, mp_srcptr s2, mp_size_t n, mp_limb_t carry) {
     for(mp_size_t i=0; i<n; ++i) {
        unsigned long long a = s1[i];
        unsigned long long b = s2[i];
        unsigned long long sum = a + b + carry;
        res[i] = (mp_limb_t)sum;
        carry = (sum < a) || (carry && sum == a) ? 1 : 0; 
    }
    return carry;
}

mp_size_t lmmp_pow_size_(mp_srcptr b, mp_size_t bn, unsigned long exp) { return 1; }

mp_size_t lmmp_pow_(mp_ptr dst, mp_size_t rn, mp_srcptr base, mp_size_t n, unsigned long exp) {
     unsigned long long b_val = (n > 0) ? base[0] : 0;
     unsigned long long res = 1;
     for(unsigned long i=0; i<exp; ++i) res *= b_val;
     if(dst) dst[0] = (mp_limb_t)res;
     return 1;
}

mp_size_t lmmp_gcd_lehmer_(mp_ptr dst, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn) {
    unsigned long long u = (un > 0) ? up[0] : 0;
    unsigned long long v = (vn > 0) ? vp[0] : 0;
    while(v) {
        unsigned long long t = v;
        v = u % v;
        u = t;
    }
    if(dst) dst[0] = (mp_limb_t)u;
    return 1;
}

void lmmp_abort(lmmp_error_t type, const char* msg, const char* file, int line) {
    // printf("Mock Abort: %s\n", msg);
    // exit(1);
}

}
