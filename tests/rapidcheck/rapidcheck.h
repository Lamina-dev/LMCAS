/**
 * @file rapidcheck.h
 * @brief Minimal header-only property-based testing framework.
 *
 * Provides a lightweight rc::check() API for property-based testing,
 * compatible with the project's custom test harness. Generates random
 * inputs and verifies that properties hold across many iterations.
 *
 * Usage:
 *   rc::check("description", [](int a, int b) { RC_ASSERT(a + b == b + a); });
 */
#pragma once

#include <functional>
#include <random>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <type_traits>
#include <cstdint>

namespace rc {

/// Exception thrown when a property assertion fails.
struct PropertyFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

/// Global random engine for property generation.
inline std::mt19937& global_rng() {
    static std::mt19937 rng(42); // Fixed seed for reproducibility
    return rng;
}

/// Reset the RNG seed (call between test runs if needed).
inline void seed(uint32_t s) { global_rng().seed(s); }

// ============================================================
// Generators: produce random values of various types
// ============================================================

namespace gen {

/// Generate a random int in [lo, hi].
inline int inRange(int lo, int hi) {
    std::uniform_int_distribution<int> dist(lo, hi);
    return dist(global_rng());
}

/// Generate a random size_t in [lo, hi].
inline size_t inRangeSizeT(size_t lo, size_t hi) {
    std::uniform_int_distribution<size_t> dist(lo, hi);
    return dist(global_rng());
}

/// Generate a random double in [lo, hi].
inline double inRangeDouble(double lo, double hi) {
    std::uniform_real_distribution<double> dist(lo, hi);
    return dist(global_rng());
}

/// Generate a random bool.
inline bool boolean() {
    return inRange(0, 1) == 1;
}

/// Pick a random element from a vector.
template<typename T>
T elementOf(const std::vector<T>& vec) {
    if (vec.empty()) throw std::runtime_error("elementOf: empty vector");
    return vec[inRangeSizeT(0, vec.size() - 1)];
}

} // namespace gen

// ============================================================
// Arbitrary: type-based random generation
// ============================================================

template<typename T>
struct Arbitrary;

template<>
struct Arbitrary<int> {
    static int generate() { return gen::inRange(-100, 100); }
};

template<>
struct Arbitrary<unsigned int> {
    static unsigned int generate() { return static_cast<unsigned int>(gen::inRange(0, 100)); }
};

template<>
struct Arbitrary<double> {
    static double generate() { return gen::inRangeDouble(-100.0, 100.0); }
};

template<>
struct Arbitrary<bool> {
    static bool generate() { return gen::boolean(); }
};

template<>
struct Arbitrary<std::string> {
    static std::string generate() {
        int len = gen::inRange(1, 8);
        std::string s;
        for (int i = 0; i < len; ++i) {
            s += static_cast<char>('a' + gen::inRange(0, 25));
        }
        return s;
    }
};

// ============================================================
// RC_ASSERT macro
// ============================================================

#define RC_ASSERT(expr) \
    do { \
        if (!(expr)) { \
            std::ostringstream oss_; \
            oss_ << "Property assertion failed: " #expr \
                 << " at " << __FILE__ << ":" << __LINE__; \
            throw ::rc::PropertyFailure(oss_.str()); \
        } \
    } while(0)

#define RC_CLASSIFY(cond, label) (void)(cond)

// ============================================================
// check: run a property test with N iterations
// ============================================================

namespace detail {

/// Helper to invoke a callable with generated arguments.
template<typename F, typename... Args>
void invoke_with_generated(F&& f, std::tuple<Args...>*) {
    f(Arbitrary<std::decay_t<Args>>::generate()...);
}

/// Extract function argument types from a lambda/function.
template<typename Ret, typename Class, typename... Args>
std::tuple<Args...>* arg_types(Ret(Class::*)(Args...) const) { return nullptr; }

template<typename Ret, typename Class, typename... Args>
std::tuple<Args...>* arg_types(Ret(Class::*)(Args...)) { return nullptr; }

template<typename Ret, typename... Args>
std::tuple<Args...>* arg_types(Ret(*)(Args...)) { return nullptr; }

} // namespace detail

/// Number of iterations per property check (configurable).
inline int& num_iterations() {
    static int n = 100;
    return n;
}

/**
 * @brief Run a property-based test.
 *
 * Generates random inputs and verifies the property holds for all of them.
 * Reports pass/fail using the project's test harness globals.
 *
 * @param description Human-readable property description
 * @param property A callable (lambda) whose arguments are randomly generated
 * @return true if all iterations passed, false otherwise
 */
template<typename F>
bool check(const std::string& description, F&& property) {
    for (int i = 0; i < num_iterations(); ++i) {
        try {
            using FType = std::decay_t<F>;
            detail::invoke_with_generated(
                std::forward<F>(property),
                detail::arg_types(&FType::operator()));
        } catch (const PropertyFailure& e) {
            std::cerr << "[FAIL] Property: " << description
                      << "\n  Iteration " << i << ": " << e.what() << std::endl;
            ::g_failures++;
            return false;
        } catch (const std::exception& e) {
            std::cerr << "[FAIL] Property: " << description
                      << "\n  Iteration " << i << " threw: " << e.what() << std::endl;
            ::g_failures++;
            return false;
        }
    }
    std::cout << "[PASS] Property (" << num_iterations() << " iterations): "
              << description << std::endl;
    ::g_passes++;
    return true;
}

/**
 * @brief Run a property-based test with a no-argument callable.
 *
 * The callable is responsible for generating its own random inputs
 * using rc::gen:: functions.
 */
inline bool check(const std::string& description, std::function<void()> property) {
    for (int i = 0; i < num_iterations(); ++i) {
        try {
            property();
        } catch (const PropertyFailure& e) {
            std::cerr << "[FAIL] Property: " << description
                      << "\n  Iteration " << i << ": " << e.what() << std::endl;
            ::g_failures++;
            return false;
        } catch (const std::exception& e) {
            std::cerr << "[FAIL] Property: " << description
                      << "\n  Iteration " << i << " threw: " << e.what() << std::endl;
            ::g_failures++;
            return false;
        }
    }
    std::cout << "[PASS] Property (" << num_iterations() << " iterations): "
              << description << std::endl;
    ::g_passes++;
    return true;
}

} // namespace rc
