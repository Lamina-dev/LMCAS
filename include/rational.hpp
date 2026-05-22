#pragma once
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "bigint.hpp"
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <map>
#include <stdint.h>

class Rational {
private:
    BigInt numerator;
    BigInt denominator;

    void simplify() {
        if (!denominator) {
            throw std::runtime_error("Denominator cannot be zero");
        }

        if (denominator.IsNegative()) {
            numerator = -numerator;
            denominator = denominator.Abs();
        }

        BigInt g = BigInt::gcd(numerator, denominator);
        if (g != 0 && g != 1) {
            numerator = numerator / g;
            denominator = denominator / g;
        }
    }

public:

    Rational() : numerator(0), denominator(1) {}

    Rational(const BigInt& num) : numerator(num), denominator(1) {}

    Rational(const BigInt& num, const BigInt& den) : numerator(num), denominator(den) {
        simplify();
    }

    Rational(int num) : numerator(num), denominator(1) {}

    Rational(int num, int den) : numerator(num), denominator(den) {
        simplify();
    }

    Rational(const std::string& num):Rational() {
        if (num.empty() || num == "0") return;
        BigInt &up = numerator, &down = denominator;
        BigInt BigInt0to10[11];
        for (int i = 0; i < 11; i++) BigInt0to10[i] = BigInt(i);
        uint64_t i = 0;
        if (num[i] == '-') i++;
        while (i < num.size() && num[i] != '.' && num[i] != 'e' && num[i] != 'E') {
            up *= BigInt0to10[10];
            up += BigInt0to10[num[i] - '0'];
            i++;
        }
        if (i == num.size() || num[i] == 'e' || num[i] == 'E') {
            if (i != num.size() && (num[i] == 'e' || num[i] == 'E')) {
                i++;
                if (num[i] == '-') {
                    i++;
                    down *= BigInt0to10[10].power(BigInt(std::string(num.begin() + i, num.end())));
                } else {
                    up *= BigInt0to10[10].power(BigInt(std::string(num.begin() + i, num.end())));
                }
            }
            return ;
        }
        i++;
        std::vector<short> n;
        while (i < num.size() && num[i] != '.' && num[i] != 'e' && num[i] != 'E') {
            n.emplace_back(num[i] - '0');
            i++;
        }
        if (i == num.size() || ((num[i] == 'e' || num[i] == 'E'))) {

            for (short& i: n) {
                up *= BigInt0to10[10];
                down *= BigInt0to10[10];
                up += BigInt0to10[i];
            }
        } else {

            std::pair<uint64_t, uint64_t> xun = jiance(n);

            BigInt xup = BigInt(0), xdown = BigInt(0);
            for (uint64_t i = xun.first; i <= xun.second; i++) {
                xup *= BigInt0to10[10];
                xup += BigInt0to10[n[i]];
                xdown *= BigInt0to10[10];
                xdown += BigInt0to10[9];
            }

            for (uint64_t i = 0; i < xun.first; i++) {
                up *= BigInt0to10[10];
                up += BigInt0to10[n[i]];
                down *= BigInt0to10[10];
            }

            up = up * xdown;
            down = down * xdown;

            up = up + xup;
        }

        if (i != num.size() && (num[i] == 'e' || num[i] == 'E')) {
            i++;
            if (num[i] == '-') {
                i++;
                down *= BigInt0to10[10].power(BigInt(std::string(num.begin() + i, num.end())));
            } else {
                up *= BigInt0to10[10].power(BigInt(std::string(num.begin() + i, num.end())));
            }
        }

        if (num[0] == '-') up = BigInt(0) - up;

        simplify();
    }

private:
    std::pair<uint64_t, uint64_t> jiance(std::vector<short>& n) {

        uint64_t h1[10] = {};
        uint64_t h2[10] = {};

        for (short& i: n) h1[i]++;
        bool flag = true;

        for (uint64_t i = 0; i < n.size(); i++) {
            flag = true;
            for (int j = 0; j < 10; j++) {
                if ((h1[j] - h2[j]) & 1) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                uint64_t i1 = i + 1, i2 = n.size();
                uint64_t ti1, ti2;
                bool ok = true;
                while (!((i2 - i1) & 1) && ok) {

                    ti1 = i1;
                    ti2 = i1 + ((i2 - i1) >> 1);
                    while (ti2 < i2) {
                        if (n[ti1] != n[ti2]) {
                            ok = false;
                            break;
                        }
                        ti1++;
                        ti2++;
                    }
                    if (ok) {
                        i2 = i1 + ((i2 - i1) >> 1);
                    } else if (i2 == n.size()) {
                        break;
                    } else {
                        return {i1, i2 - 1};
                    }
                }
                if (ok) {
                    return {i1 - 1, i2 - 1};
                }
            }
            h2[n[i]]++;
        }
        return {n.size() - 1, n.size() - 1};
    }

public:

    static Rational from_double(double value) {
        if (value == 0.0) {
            return Rational();
        }
        if (std::floor(value) == value) {
            return Rational(BigInt(std::to_string(value)));
        }
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(15) << value;
        std::string str = oss.str();

        size_t e_pos = str.find('e');
        std::string base = str.substr(0, e_pos);
        bool negative = false;
        if (base[0] == '-') {
            negative = true;
            base.erase(0, 1);
        }
        base.erase(1, 1);
        while (base.back() == '0') {
            base.pop_back();
        }
        int exponent = std::stoi(str.substr(e_pos + 1));

        size_t decimal_places = std::max(0, static_cast<int>(base.length()) - exponent - 1);

        BigInt num(base);
        if (negative) num = -num;
        BigInt den("1" + std::string(decimal_places, '0'));

        Rational r(num, den);
        r.simplify();
        return r;
    }

    BigInt get_numerator() const { return numerator; }
    BigInt get_denominator() const { return denominator; }

    bool is_integer() const {
        return denominator == BigInt(1);
    }

    bool is_zero() const {
        return !numerator;
    }

    std::size_t hash() const {
        std::size_t seed = numerator.hash();
        std::size_t d_hash = denominator.hash();
        seed ^= d_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }

    std::string to_string() const {
        if (is_integer()) {
            return numerator.ToString();
        }
        return numerator.ToString() + "/" + denominator.ToString();
    }

    inline std::string to_float_string() {
        if (numerator % denominator == BigInt(0)) return (numerator / denominator).ToString();
        std::string re;
        const BigInt ten(10);
        const BigInt zero(0);
        if ((numerator < zero) ^ (denominator < zero)) re += '-';
        BigInt num = numerator.Abs(), den = denominator.Abs();
        re += (num / den).ToString();
        num %= den;
        re += '.';

        std::map<BigInt, size_t> num_index;

        while (num != zero && num_index.find(num) == num_index.end()) {
            num_index[num] = re.size();
            re += '0';
            num *= ten;
            while (num >= den) {
                num -= den;
                re.back()++;
            }
        }
        if (num_index.find(num) != num_index.end()) {
            re.insert(re.begin() + num_index[num], '(');
            re += ")...";
        }
        return re;

    }

    inline std::string to_float_string(int64_t n) const{
        if (is_zero()) return "0";
        if (n == 0) return (numerator / denominator).ToString();
        if (n > 0) {
            BigInt pow10n(10);
            pow10n = pow10n.power(BigInt(n));
            BigInt scaled = numerator * pow10n / denominator;
            bool negative = scaled < BigInt(0);
            if (negative) scaled = scaled * BigInt(-1);
            std::string re = scaled.ToString();
            if ((int64_t)re.size() <= n) re.insert(re.begin(), n - re.size() + 1, '0');
            re.insert(re.end() - n,'.');
            while (re.back() == '0') re.pop_back();
            if (re.back() == '.') re.pop_back();
            if (negative) re.insert(re.begin(), '-');
            return re;
        } else {
            BigInt pow10n(10);
            pow10n = pow10n.power(BigInt(n).Abs());
            std::string re = (numerator / (denominator * pow10n)).ToString();
            if (re.size() == 1 && re[0] == '0') return re;
            for (; n < 0; n++) re.push_back('0');
            return re;
        }
    }

    inline void floor(const int64_t& n) {
        floor_without_sim(n);
        simplify();
    }

private:

    inline void floor_without_sim(const int64_t& n) {
        if (n == 0) return;
        if (n > 0) {
            BigInt pow10n(10);
            pow10n = pow10n.power(BigInt(n));
            numerator *= pow10n;
            numerator /= denominator;
            denominator = pow10n;
        } else {
            BigInt pow10n(10);
            pow10n = pow10n.power(BigInt(n).Abs());
            denominator *= pow10n;
            numerator /= denominator;
            numerator *= pow10n;
            denominator = BigInt(1);
        }
    }

public:

    inline void sqrt_self(int64_t n) {
        if (numerator < BigInt(0)) throw std::runtime_error("Sqrt negative number");
        Rational t(numerator * denominator);

        const int64_t nadd1 = n + 1;
        const BigInt ten(10);
        Rational ans( (t * t + Rational(6) * t + Rational(1))  /  (Rational(4) * (t + Rational(1))) );
        Rational temp;
        while (temp.numerator / ten != ans.numerator / ten) {
            temp = (ans + t / ans) / 2;
            temp.floor_without_sim(nadd1);
            ans = (temp + t / temp) / 2;
            ans.floor_without_sim(nadd1);
        }
        numerator = ans.numerator;
        denominator *= ans.denominator;
        floor(n);
    }

    inline Rational sqrt(int64_t n) const{
        Rational re = *this;
        re.sqrt_self(n);
        return re;
    }

    inline void radicand_self(const BigInt& radical, int64_t n) {
        if (numerator < BigInt(0) && (radical % BigInt(2) == BigInt(0))) throw std::runtime_error("Radicand negative number");
        Rational t(numerator * denominator.power(radical - BigInt(1)));

        const int64_t nadd1 = n + 1;
        const Rational rat_radical(radical);
        const BigInt ten(10);
        Rational ans(Rational(t.numerator,radical) * (Rational(radical - BigInt(1)) + Rational(t.numerator, t.numerator.power(radical))));

        Rational temp;
        while (temp.numerator / ten != ans.numerator / ten) {
            temp = ans / rat_radical * ((rat_radical - Rational(1)) + (t / ans.power(radical)));
            temp.floor_without_sim(nadd1);
            ans = temp / rat_radical * ((rat_radical - Rational(1)) + (t / temp.power(radical)));
            ans.floor_without_sim(nadd1);
        }

        numerator = ans.numerator;
        denominator *= ans.denominator;
        floor(n);
    }

    inline Rational radicand(const BigInt& radical, int64_t n) {
        Rational re = *this;
        re.radicand_self(radical, n);
        return re;
    }

    BigInt to_BigInt() const {
        if (!is_integer()) {
            throw std::runtime_error("Cannot convert non-integer fraction to BigInt");
        }
        return numerator;
    }

    lmmc_real_t to_double() const {
        if (is_zero()) return 0.0;
        lmmc_real_t res = 0.0;
#ifdef LMMC_SCN_REAL
        sscanf(to_float_string(15).c_str(), "%" LMMC_SCN_REAL, &res);
#else
        res = std::stod(to_float_string(15));
#endif
        return res;
    }

    Rational operator+(const Rational& other) const {
        BigInt new_num = numerator * other.denominator + other.numerator * denominator;
        BigInt new_den = denominator * other.denominator;
        return Rational(new_num, new_den);
    }

    Rational operator-(const Rational& other) const {
        BigInt new_num = numerator * other.denominator - other.numerator * denominator;
        BigInt new_den = denominator * other.denominator;
        return Rational(new_num, new_den);
    }

    Rational operator*(const Rational& other) const {
        BigInt new_num = numerator * other.numerator;
        BigInt new_den = denominator * other.denominator;
        return Rational(new_num, new_den);
    }

    Rational operator/(const Rational& other) const {
        if (other.is_zero()) {
            throw std::runtime_error("Division by zero");
        }
        BigInt new_num = numerator * other.denominator;
        BigInt new_den = denominator * other.numerator;
        return Rational(new_num, new_den);
    }

    Rational power(const BigInt& exponent) const {
        if (exponent < BigInt(0)) {

            if (is_zero()) {
                throw std::runtime_error("Cannot raise zero to negative power");
            }
            BigInt pos_exp = -exponent;
            return Rational(denominator.power(pos_exp), numerator.power(pos_exp));
        }

        if (!exponent) {
            if (is_zero()) {
                throw std::runtime_error("0^0 is undefined");
            }
            return Rational(1);
        }

        return Rational(numerator.power(exponent), denominator.power(exponent));
    }

    inline void power_self(const Rational& exponent,size_t n) {
        *this = power(exponent.numerator);
        radicand_self(exponent.denominator,n);
    }

    inline Rational power(const Rational& exponent, size_t n) {
        Rational re = *this;
        re.power_self(exponent, n);
        return re;
    }

    bool operator==(const Rational& other) const {
        return numerator * other.denominator == other.numerator * denominator;
    }

    bool operator!=(const Rational& other) const {
        return !(*this == other);
    }

    bool operator<(const Rational& other) const {
        BigInt left = numerator * other.denominator;
        BigInt right = other.numerator * denominator;
        return left < right;
    }

    bool operator<=(const Rational& other) const {
        return *this < other || *this == other;
    }

    bool operator>(const Rational& other) const {
        return !(*this <= other);
    }

    bool operator>=(const Rational& other) const {
        return !(*this < other);
    }

    Rational reciprocal() const {
        if (is_zero()) {
            throw std::runtime_error("Cannot take reciprocal of zero");
        }
        return Rational(denominator, numerator);
    }

    Rational abs() const {
        Rational result = *this;
        result.numerator = result.numerator.Abs();
        return result;
    }

    Rational operator-() const {
        Rational result = *this;
        result.numerator = -result.numerator;
        return result;
    }
};
