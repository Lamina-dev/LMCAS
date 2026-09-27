#include "bigint.hpp"
#include "internal/lmmc_lifecycle.hpp"

namespace LMCAS {
namespace {
void bigint_mul(mp_ptr destination, mp_srcptr lhs, mp_size_t lhs_size,
                 mp_srcptr rhs, mp_size_t rhs_size) {
    LMCAS::detail::ensure_lmmc_lifecycle();
    ::lmmp_mul_(destination, lhs, lhs_size, rhs, rhs_size);
}
void bigint_div(mp_ptr quotient, mp_ptr remainder, mp_srcptr numerator,
                 mp_size_t numerator_size, mp_srcptr denominator,
                 mp_size_t denominator_size) {
    LMCAS::detail::ensure_lmmc_lifecycle();
    ::lmmp_div_(quotient, remainder, numerator, numerator_size, denominator,
                denominator_size);
}
}

void BigInt::add_abs(BigInt& dst, const BigInt& a, const BigInt& b) {
        mp_size_t n = std::max(a.size_, b.size_);
        dst.realloc_to(n + 1);

        mp_limb_t cy = 0;

        mp_srcptr ap = a.data_;
        mp_srcptr bp = b.data_;
        mp_size_t na = a.size_;
        mp_size_t nb = b.size_;

        if (na < nb) { std::swap(ap, bp); std::swap(na, nb); }

        if (nb > 0) {
             cy = lmmp_add_n_(dst.data_, ap, bp, nb);
        }

        if (na > nb) {
             if (dst.data_ != ap) std::memcpy(dst.data_ + nb, ap + nb, (na - nb) * sizeof(mp_limb_t));

             mp_size_t k = nb;
             while (cy && k < na) {
                 dst.data_[k]++;
                 cy = (dst.data_[k] == 0);
                 k++;
             }
        }
        dst.data_[na] = cy;
        dst.size_ = na + cy;
        dst.normalize();
    }

void BigInt::sub_abs(BigInt& dst, const BigInt& a, const BigInt& b) {
        mp_size_t na = a.size_;
        mp_size_t nb = b.size_;
        dst.realloc_to(na);

        mp_limb_t bw = 0;
        if (nb > 0)
            bw = lmmp_sub_n_(dst.data_, a.data_, b.data_, nb);

        if (na > nb) {
             if (dst.data_ != a.data_) std::memcpy(dst.data_ + nb, a.data_ + nb, (na - nb) * sizeof(mp_limb_t));
             mp_size_t k = nb;
             while (bw && k < na) {
                 dst.data_[k]--;
                 bw = (dst.data_[k] == (mp_limb_t)-1);
                 k++;
             }
        }
        dst.size_ = na;
        dst.normalize();
    }

BigInt BigInt::operator+(const BigInt& other) const {
        if (sign_ == ZERO) return other;
        if (other.sign_ == ZERO) return *this;

        BigInt res;
        if (sign_ == other.sign_) {
            add_abs(res, *this, other);
            res.sign_ = sign_;
        } else {

            int cmp = cmp_abs(*this, other);
            if (cmp == 0) {
                return BigInt(0);
            }
            if (cmp > 0) {
                sub_abs(res, *this, other);
                res.sign_ = sign_;
            } else {
                sub_abs(res, other, *this);
                res.sign_ = other.sign_;
            }
        }
        return res;
    }

BigInt BigInt::operator-(const BigInt& other) const {
        return *this + (-other);
    }

BigInt BigInt::operator*(const BigInt& other) const {
        if (size_ == 0 || other.size_ == 0) return BigInt(0);

        BigInt res;
        mp_size_t na = size_;
        mp_size_t nb = other.size_;
        res.realloc_to(na + nb);

        if (na >= nb) {
            bigint_mul(res.data_, data_, na, other.data_, nb);
        } else {
            bigint_mul(res.data_, other.data_, nb, data_, na);
        }

        res.size_ = na + nb;
        res.sign_ = (sign_ == other.sign_) ? POSITIVE : NEGATIVE;
        res.normalize();
        return res;
    }

BigInt BigInt::operator/(const BigInt& other) const {
        if (other.size_ == 0) throw std::domain_error("Division by zero");
        if (size_ < other.size_) return BigInt(0);

        BigInt q, r;
        mp_size_t na = size_;
        mp_size_t nb = other.size_;

        q.realloc_to(na - nb + 1);
        r.realloc_to(nb);

        bigint_div(q.data_, r.data_, data_, na, other.data_, nb);

        q.size_ = na - nb + 1;
        q.sign_ = (sign_ == other.sign_) ? POSITIVE : NEGATIVE;
        q.normalize();
        return q;
    }

BigInt BigInt::operator%(const BigInt& other) const {
        if (other.size_ == 0) throw std::domain_error("Division by zero");
        if (size_ < other.size_) return *this;

        BigInt q;
        mp_size_t na = size_;
        mp_size_t nb = other.size_;
        q.realloc_to(na - nb + 1);

        BigInt r;
        r.realloc_to(nb);

        bigint_div(q.data_, r.data_, data_, na, other.data_, nb);

        r.size_ = nb;
        r.sign_ = sign_;
        r.normalize();
        return r;
    }

BigInt& BigInt::operator+=(const BigInt& other) { *this = *this + other; return *this; }

BigInt& BigInt::operator-=(const BigInt& other) { *this = *this - other; return *this; }

BigInt& BigInt::operator*=(const BigInt& other) { *this = *this * other; return *this; }

BigInt& BigInt::operator/=(const BigInt& other) { *this = *this / other; return *this; }

BigInt& BigInt::operator%=(const BigInt& other) { *this = *this % other; return *this; }

BigInt BigInt::power(unsigned long exp) const {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (exp == 0) return BigInt(1);
        if (size_ == 0) return BigInt(0);

        BigInt res;
        mp_size_t needed = lmmp_pow_size_(data_, size_, exp);
        res.realloc_to(needed);

        res.size_ = lmmp_pow_(res.data_, needed, data_, size_, exp);

        if (sign_ == NEGATIVE && (exp & 1)) {
            res.sign_ = NEGATIVE;
        } else {
            res.sign_ = POSITIVE;
        }
        res.normalize();
        return res;
    }

BigInt BigInt::power(BigInt exp) const {
        if (exp.sign_ == NEGATIVE) throw std::domain_error("Negative exponent in integer power");
        if (exp.size_ == 0) return BigInt(1);

        if (exp.size_ <= 1) {
            mp_limb_t lo = (exp.size_ == 0) ? 0 : exp.data_[0];
            if (lo <= static_cast<mp_limb_t>(std::numeric_limits<unsigned long>::max())) {
                return power(static_cast<unsigned long>(lo));
            }
        }

        BigInt base = *this;
        BigInt res(1);

        while (!exp.is_zero()) {
            if (exp.data_[0] & 1) {
                res = res * base;
            }
            base = base * base;

             mp_limb_t carry = 0;
             for (mp_size_t i = exp.size_; i > 0; --i) {
                 mp_limb_t cur = exp.data_[i-1];
                 mp_limb_t next_carry = cur & 1;
                 exp.data_[i-1] = (cur >> 1) | (carry << (LIMB_BITS - 1));
                 carry = next_carry;
             }
             exp.normalize();
        }
        return res;
    }

BigInt BigInt::sqrt() const {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (sign_ == NEGATIVE) throw std::domain_error("Sqrt of negative number");
        if (size_ == 0) return BigInt(0);
        if (*this == BigInt(1)) return BigInt(1);

        BigInt res;

        const mp_size_t root_size = (size_ + 1) / 2;
        res.realloc_to(root_size + 1);

        lmmp_sqrt_(res.data_, nullptr, data_, size_, 0);

        res.size_ = root_size;
        res.sign_ = POSITIVE;
        res.normalize();

        return res;
    }

bool BigInt::is_odd() const {
        if (size_ == 0) return false;
        return (data_[0] & 1);
    }

bool BigInt::is_even() const {
        return !is_odd();
    }

mp_size_t BigInt::trailing_zeros() const {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (is_zero()) return 0;
        mp_size_t count = 0;
        for (mp_size_t i = 0; i < size_; ++i) {
            if (data_[i] == 0) {
                count += LIMB_BITS;
            } else {
                count += lmmp_tailing_zeros_(data_[i]);
                break;
            }
        }
        return count;
    }

BigInt& BigInt::operator>>=(mp_size_t shift) {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (shift == 0) return *this;
        if (is_zero()) return *this;

        mp_size_t limb_shift = shift / LIMB_BITS;
        mp_size_t bit_shift = shift % LIMB_BITS;

        if (limb_shift >= size_) {
            zero();
            return *this;
        }

        if (limb_shift > 0) {
            std::memmove(data_, data_ + limb_shift, (size_ - limb_shift) * sizeof(mp_limb_t));
            size_ -= limb_shift;
        }

        if (bit_shift > 0) {
             lmmp_shr_(data_, data_, size_, bit_shift);
             if (size_ > 0 && data_[size_-1] == 0) size_--;
        }
        normalize();
        return *this;
    }

BigInt& BigInt::operator<<=(mp_size_t shift) {
        if (shift == 0) return *this;
        if (is_zero()) return *this;

        mp_size_t limb_shift = shift / LIMB_BITS;
        mp_size_t bit_shift = shift % LIMB_BITS;

        mp_size_t old_size = size_;
        mp_size_t needed = old_size + limb_shift + (bit_shift > 0 ? 1 : 0);
        realloc_to(needed);

        if (bit_shift > 0) {
            mp_limb_t carry = lmmp_shl_(data_ + limb_shift, data_, old_size, bit_shift);
            if (carry) {
                data_[old_size + limb_shift] = carry;
                size_ = old_size + limb_shift + 1;
            } else {
                 size_ = old_size + limb_shift;
            }
        } else {
             std::memmove(data_ + limb_shift, data_, old_size * sizeof(mp_limb_t));
             size_ += limb_shift;
        }

        if (limb_shift > 0) {
            std::memset(data_, 0, limb_shift * sizeof(mp_limb_t));
        }

        normalize();
        return *this;
    }

BigInt BigInt::operator>>(mp_size_t shift) const {
        BigInt res = *this;
        res >>= shift;
        return res;
    }

BigInt BigInt::operator<<(mp_size_t shift) const {
        BigInt res = *this;
        res <<= shift;
        return res;
    }

}
