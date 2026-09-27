#include "bigint.hpp"
#include "internal/lmmc_lifecycle.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace LMCAS {


void BigInt::realloc_to(mp_size_t new_alloc) {
        if (new_alloc <= capacity_) return;
        const std::size_t requested = static_cast<std::size_t>(new_alloc);
        if (requested > std::numeric_limits<std::size_t>::max() - 3) {
            throw std::length_error("BigInt allocation size overflow");
        }
        const std::size_t rounded = (requested + 3) & ~std::size_t(3);
        if (rounded > std::numeric_limits<std::size_t>::max() / sizeof(mp_limb_t) ||
            rounded > static_cast<std::size_t>(std::numeric_limits<mp_size_t>::max())) {
            throw std::length_error("BigInt allocation byte size overflow");
        }

        mp_ptr new_data =
            static_cast<mp_ptr>(::operator new[](rounded * sizeof(mp_limb_t)));

        if (size_ > 0 && data_) {
             std::memcpy(new_data, data_, size_ * sizeof(mp_limb_t));
        }
        if (data_) ::operator delete[](data_);
        data_ = new_data;
        capacity_ = static_cast<mp_size_t>(rounded);
    }

void BigInt::normalize() {
        if (size_ > 0 && !data_) {
            throw std::logic_error("BigInt invariant violated: nonzero size without storage");
        }
        while (size_ > 0 && data_[size_ - 1] == 0) {
            size_--;
        }
        if (size_ == 0) {
            sign_ = ZERO;
        } else if (sign_ == ZERO) {
            sign_ = POSITIVE;
        }
    }

void BigInt::zero() {
        size_ = 0;
        sign_ = ZERO;
    }

BigInt::BigInt() = default;

BigInt::~BigInt() {
        if (data_) ::operator delete[](data_);
    }

BigInt::BigInt(const BigInt& other) {
        if (other.size_ > 0) {
            realloc_to(other.size_);
            std::memcpy(data_, other.data_, other.size_ * sizeof(mp_limb_t));
            size_ = other.size_;
            sign_ = other.sign_;
        } else {
            zero();
        }
    }

BigInt::BigInt(BigInt&& other) noexcept
        : data_(other.data_), size_(other.size_), capacity_(other.capacity_), sign_(other.sign_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
        other.sign_ = ZERO;
    }

BigInt& BigInt::operator=(const BigInt& other) {
        if (this != &other) {
            if (other.size_ > capacity_) {
                realloc_to(other.size_);
            }
            if (other.size_ > 0)
                std::memcpy(data_, other.data_, other.size_ * sizeof(mp_limb_t));
            size_ = other.size_;
            sign_ = other.sign_;
        }
        return *this;
    }

BigInt& BigInt::operator=(BigInt&& other) noexcept {
        if (this != &other) {
            if (data_) ::operator delete[](data_);
            data_ = other.data_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            sign_ = other.sign_;

            other.data_ = nullptr;
            other.capacity_ = 0;
            other.zero();
        }
        return *this;
    }

void BigInt::import_unsigned(std::uint64_t magnitude) {
    zero();
    if (magnitude == 0) {
        return;
    }
    constexpr auto limb_bits = std::numeric_limits<mp_limb_t>::digits;
    constexpr auto value_bits = std::numeric_limits<std::uint64_t>::digits;
    realloc_to((value_bits + limb_bits - 1) / limb_bits);
    do {
        data_[size_++] = static_cast<mp_limb_t>(magnitude);
        if constexpr (limb_bits < value_bits) {
            magnitude >>= limb_bits % value_bits;
        } else {
            magnitude = 0;
        }
    } while (magnitude != 0);
    sign_ = POSITIVE;
}

BigInt::BigInt(std::int64_t value) {
    const auto bits = static_cast<std::uint64_t>(value);
    import_unsigned(value < 0 ? std::uint64_t{0} - bits : bits);
    if (value < 0) {
        sign_ = NEGATIVE;
    }
}

BigInt::BigInt(std::uint64_t value) {
    import_unsigned(value);
}

BigInt::BigInt(int value) : BigInt(static_cast<std::int64_t>(value)) {}

BigInt::BigInt(unsigned int value) : BigInt(static_cast<std::uint64_t>(value)) {}

BigInt::BigInt(const std::string& str) {
        if (str.empty()) {
            throw std::invalid_argument("BigInt: decimal string must contain digits");
        }
        size_t start = 0;
        int sign = POSITIVE;
        if (str[0] == '-') {
            sign = NEGATIVE;
            start = 1;
        } else if (str[0] == '+') {
            start = 1;
        }

        if (start == str.length()) {
            throw std::invalid_argument("BigInt: sign must be followed by digits");
        }

        size_t len = str.length() - start;

        std::vector<mp_byte_t> digit_buf(len);
        for (size_t i = 0; i < len; ++i) {
            char c = str[start + i];
            if (c < '0' || c > '9') {
                throw std::invalid_argument(
                    "BigInt: invalid character '" + std::string(1, c) +
                    "' at position " + std::to_string(start + i) +
                    " in input \"" + str + "\"");
            }
            uint8_t d = static_cast<uint8_t>(c - '0');

            digit_buf[len - 1 - i] = d;
        }

        mp_size_t needed = len / 19 + 2;
        realloc_to(needed);
        LMCAS::detail::ensure_lmmc_lifecycle();

        size_ = lmmp_from_str_(data_, digit_buf.data(), len, 10);

        if (size_ == 0) {
            zero();
        } else {
            sign_ = sign;
        }
        normalize();
    }

std::string BigInt::to_string() const {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (size_ == 0) return "0";

        size_t len_needed = size_ * 20 + 5;
        std::vector<mp_byte_t> buf(len_needed);

        mp_size_t str_len = lmmp_to_str_((mp_byte_t*)buf.data(), data_, size_, 10);

        if (str_len == 0) return "0";

        std::string res;
        res.reserve(str_len + 2);
        if (sign_ == NEGATIVE) res += '-';

        for(mp_size_t i = str_len; i > 0; --i) {
            res += (char)(buf[i - 1] + '0');
        }
        return res;
    }


std::size_t BigInt::bit_length() const noexcept {
    if (size_ == 0) {
        return 0;
    }
    std::size_t top_bits = 0;
    for (mp_limb_t top = data_[size_ - 1]; top != 0; top >>= 1) {
        ++top_bits;
    }
    constexpr auto limb_bits = std::numeric_limits<mp_limb_t>::digits;
    return (static_cast<std::size_t>(size_) - 1) * limb_bits + top_bits;
}

std::optional<std::uint64_t> BigInt::magnitude_uint64() const noexcept {
    constexpr auto value_bits = std::numeric_limits<std::uint64_t>::digits;
    constexpr auto limb_bits = std::numeric_limits<mp_limb_t>::digits;
    if (bit_length() > value_bits) {
        return std::nullopt;
    }
    std::uint64_t magnitude = 0;
    for (mp_size_t i = size_; i > 0; --i) {
        if constexpr (limb_bits < value_bits) {
            magnitude <<= limb_bits % value_bits;
        }
        magnitude |= static_cast<std::uint64_t>(data_[i - 1]);
    }
    return magnitude;
}

std::optional<std::uint64_t> BigInt::try_to_uint64() const noexcept {
    if (sign_ == NEGATIVE) {
        return std::nullopt;
    }
    return magnitude_uint64();
}

std::optional<std::int64_t> BigInt::try_to_int64() const noexcept {
    const auto magnitude = magnitude_uint64();
    if (!magnitude) {
        return std::nullopt;
    }
    constexpr auto maximum =
        static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
    if (sign_ == NEGATIVE) {
        if (*magnitude > maximum + std::uint64_t{1}) {
            return std::nullopt;
        }
        if (*magnitude == maximum + std::uint64_t{1}) {
            return std::numeric_limits<std::int64_t>::min();
        }
        return -static_cast<std::int64_t>(*magnitude);
    }
    if (*magnitude > maximum) {
        return std::nullopt;
    }
    return static_cast<std::int64_t>(*magnitude);
}

int BigInt::to_int() const {
    const auto value = try_to_int64();
    if (!value) {
        return sign_ == NEGATIVE ? std::numeric_limits<int>::min()
                                 : std::numeric_limits<int>::max();
    }
    if (*value < std::numeric_limits<int>::min()) {
        return std::numeric_limits<int>::min();
    }
    if (*value > std::numeric_limits<int>::max()) {
        return std::numeric_limits<int>::max();
    }
    return static_cast<int>(*value);
}

lmmc_real_t BigInt::to_double() const {
        if (size_ == 0) return 0.0;
        lmmc_real_t res = 0.0;
        const lmmc_real_t base_mul =
            std::ldexp(lmmc_real_t{1}, std::numeric_limits<mp_limb_t>::digits);

        for (mp_size_t i = size_; i > 0; --i) {
            LMMC_REAL_MUL(&res, &res, &base_mul);
            lmmc_real_t val = (lmmc_real_t)data_[i-1];
            LMMC_REAL_ADD(&res, &res, &val);
        }
        if (sign_ == NEGATIVE) {
            lmmc_real_t zero = 0.0;
            LMMC_REAL_SUB(&res, &zero, &res);
        }
        return res;
    }

int BigInt::cmp_abs(const BigInt& a, const BigInt& b) {
        if (a.size_ != b.size_) {
            return a.size_ > b.size_ ? 1 : -1;
        }
        if (a.size_ == 0) {
            return 0;
        }
        for (mp_size_t i = a.size_; i > 0; --i) {
            if (a.data_[i - 1] != b.data_[i - 1]) {
                return a.data_[i - 1] > b.data_[i - 1] ? 1 : -1;
            }
        }
        return 0;
    }

bool BigInt::operator==(const BigInt& other) const {
        return sign_ == other.sign_ && cmp_abs(*this, other) == 0;
    }

bool BigInt::operator!=(const BigInt& other) const { return !(*this == other); }

bool BigInt::operator<(const BigInt& other) const {
        if (sign_ != other.sign_) return sign_ < other.sign_;
        if (sign_ == ZERO) return false;
        int cmp = cmp_abs(*this, other);
        return sign_ == POSITIVE ? cmp < 0 : cmp > 0;
    }

bool BigInt::operator>(const BigInt& other) const { return other < *this; }

bool BigInt::operator<=(const BigInt& other) const { return !(*this > other); }

bool BigInt::operator>=(const BigInt& other) const { return !(*this < other); }

bool BigInt::operator!() const { return size_ == 0; }

BigInt::operator bool() const { return size_ != 0; }

bool BigInt::is_zero() const { return size_ == 0; }

std::vector<uint64_t> BigInt::get_digits() const {
        std::vector<uint64_t> d;
        d.reserve(size_);
        for(mp_size_t i = 0; i < size_; ++i) {
            d.push_back(data_[i]);
        }
        return d;
    }

std::size_t BigInt::hash() const {
        std::size_t seed = 0;
        for (mp_size_t i = 0; i < size_; ++i) {
            seed ^= std::hash<uint64_t>{}(data_[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        if (sign_ == NEGATIVE) { seed ^= std::hash<int>{}(-1) + 0x9e3779b9 + (seed << 6) + (seed >> 2); }
        return seed;
    }

BigInt BigInt::abs() const {
        BigInt ret = *this;
        if (ret.size_ > 0) {
            ret.sign_ = POSITIVE;
        }
        return ret;
    }

bool BigInt::is_negative() const { return sign_ == NEGATIVE; }

BigInt BigInt::negate() const {
        BigInt ret = *this;
        if (ret.size_ > 0) {
            ret.sign_ = (ret.sign_ == POSITIVE) ? NEGATIVE : POSITIVE;
        }
        return ret;
    }

BigInt BigInt::operator-() const { return negate(); }

} // namespace LMCAS
