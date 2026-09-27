#include "lmmc/stdlib.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

static bool num_set_operations_failed(
    lmmc_std_num_set_t& set_a,
    lmmc_std_num_set_t& set_b,
    lmmc_std_num_set_t& set_union,
    lmmc_std_num_set_t& set_difference) {
    int subset = 0;
    if (lmmc_std_num_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
        set_union.size != 4 ||
        lmmc_std_num_set_difference(&set_union,
                                        &set_b,
                                        &set_difference) !=
                LMMC_STATUS_OK ||
        set_difference.size != 1) {
        return true;
    }
    return lmmc_std_num_set_subset(&set_difference, &set_union, &subset) !=
                LMMC_STATUS_OK ||
        !subset ||
        lmmc_std_num_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT;
}

static int check_num_sets() {
    lmmc_real_t set_a_values[] = {1.0, 2.0, 2.0, -0.0};
    lmmc_real_t set_b_values[] = {2.0, 3.0, 0.0};
    lmmc_std_num_set_t set_a = {};
    lmmc_std_num_set_t set_b = {};
    lmmc_std_num_set_t set_union = {};
    lmmc_std_num_set_t set_difference = {};
    int contains = 0;
    auto cleanup_sets = [&]() {
        lmmc_std_num_set_destroy(&set_a);
        lmmc_std_num_set_destroy(&set_b);
        lmmc_std_num_set_destroy(&set_union);
        lmmc_std_num_set_destroy(&set_difference);
    };
    if (lmmc_std_num_set_make(set_a_values, 4, &set_a) !=
            LMMC_STATUS_OK ||
        lmmc_std_num_set_make(set_b_values, 3, &set_b) !=
            LMMC_STATUS_OK ||
        set_a.size != 3 ||
        lmmc_std_num_set_contains(&set_a, 0.0, &contains) !=
            LMMC_STATUS_OK ||
        !contains ||
        num_set_operations_failed(set_a, set_b, set_union, set_difference)) {
        cleanup_sets();
        std::cerr << "failed to call installed LMMC standard library num set adapters\n";
        return 15;
    }
    cleanup_sets();
    return 0;
}

static bool text_set_operations_failed(
    lmmc_std_text_set_t& set_a,
    lmmc_std_text_set_t& set_b,
    lmmc_std_text_set_t& set_union,
    lmmc_std_text_set_t& set_difference) {
    int subset = 0;
    if (lmmc_std_text_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
        set_union.size != 3 ||
        lmmc_std_text_set_difference(&set_union,
                                         &set_b,
                                         &set_difference) !=
                LMMC_STATUS_OK ||
        set_difference.size != 1) {
        return true;
    }
    return lmmc_std_text_set_subset(&set_difference,
                                     &set_union,
                                     &subset) != LMMC_STATUS_OK ||
        !subset ||
        lmmc_std_text_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT;
}

static int check_text_sets() {
    const char* text_a_values[] = {"alpha", "beta", "alpha"};
    const char* text_b_values[] = {"beta", "gamma"};
    lmmc_std_text_set_t set_a = {};
    lmmc_std_text_set_t set_b = {};
    lmmc_std_text_set_t set_union = {};
    lmmc_std_text_set_t set_difference = {};
    int contains = 0;
    auto cleanup_sets = [&]() {
        lmmc_std_text_set_destroy(&set_a);
        lmmc_std_text_set_destroy(&set_b);
        lmmc_std_text_set_destroy(&set_union);
        lmmc_std_text_set_destroy(&set_difference);
    };
    if (lmmc_std_text_set_make(text_a_values, 3, &set_a) !=
            LMMC_STATUS_OK ||
        lmmc_std_text_set_make(text_b_values, 2, &set_b) !=
            LMMC_STATUS_OK ||
        set_a.size != 2 ||
        lmmc_std_text_set_contains(&set_a, "alpha", &contains) !=
            LMMC_STATUS_OK ||
        !contains ||
        text_set_operations_failed(set_a, set_b, set_union, set_difference)) {
        cleanup_sets();
        std::cerr << "failed to call installed LMMC standard library text set adapters\n";
        return 15;
    }
    cleanup_sets();
    return 0;
}

static int check_bool_sets() {
    int bool_a_values[] = {1, 1, 0};
    int bool_b_values[] = {0};
    lmmc_std_bool_set_t set_a = {};
    lmmc_std_bool_set_t set_b = {};
    lmmc_std_bool_set_t set_difference = {};
    int contains = 0;
    int subset = 0;
    auto cleanup_sets = [&]() {
        lmmc_std_bool_set_destroy(&set_a);
        lmmc_std_bool_set_destroy(&set_b);
        lmmc_std_bool_set_destroy(&set_difference);
    };
    if (lmmc_std_bool_set_make(bool_a_values, 3, &set_a) !=
            LMMC_STATUS_OK ||
        lmmc_std_bool_set_make(bool_b_values, 1, &set_b) !=
            LMMC_STATUS_OK ||
        set_a.size != 2 ||
        lmmc_std_bool_set_contains(&set_a, 1, &contains) !=
            LMMC_STATUS_OK ||
        !contains) {
        cleanup_sets();
        std::cerr << "failed to call installed LMMC standard library bool set adapters\n";
        return 15;
    }
    if (lmmc_std_bool_set_difference(&set_a,
                                     &set_b,
                                     &set_difference) !=
            LMMC_STATUS_OK ||
        set_difference.size != 1 ||
        lmmc_std_bool_set_subset(&set_difference,
                                 &set_a,
                                 &subset) != LMMC_STATUS_OK ||
        !subset ||
        lmmc_std_bool_set_make(nullptr, 1, &set_difference) !=
            LMMC_STATUS_INVALID_ARGUMENT) {
        cleanup_sets();
        std::cerr << "failed to call installed LMMC standard library bool set adapters\n";
        return 15;
    }
    cleanup_sets();
    return 0;
}

static bool complex_set_operations_failed(
    lmmc_std_complex_set_t& set_a,
    lmmc_std_complex_set_t& set_b,
    lmmc_std_complex_set_t& set_union,
    lmmc_std_complex_set_t& set_difference) {
    int subset = 0;
    if (lmmc_std_complex_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
        set_union.size != 3 ||
        lmmc_std_complex_set_difference(&set_union,
                                            &set_b,
                                            &set_difference) !=
                LMMC_STATUS_OK ||
        set_difference.size != 1) {
        return true;
    }
    return lmmc_std_complex_set_subset(&set_difference,
                                        &set_union,
                                        &subset) != LMMC_STATUS_OK ||
        !subset ||
        lmmc_std_complex_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT;
}

static int check_complex_sets() {
    lmmc_complex_t complex_a_values[] = {
        {1.0, 2.0}, {1.0, 2.0}, {-0.0, 0.0}};
    lmmc_complex_t complex_b_values[] = {{0.0, -0.0}, {3.0, 4.0}};
    lmmc_std_complex_set_t set_a = {};
    lmmc_std_complex_set_t set_b = {};
    lmmc_std_complex_set_t set_union = {};
    lmmc_std_complex_set_t set_difference = {};
    lmmc_complex_t zero_complex = {0.0, 0.0};
    int contains = 0;
    auto cleanup_sets = [&]() {
        lmmc_std_complex_set_destroy(&set_a);
        lmmc_std_complex_set_destroy(&set_b);
        lmmc_std_complex_set_destroy(&set_union);
        lmmc_std_complex_set_destroy(&set_difference);
    };
    if (lmmc_std_complex_set_make(complex_a_values, 3, &set_a) !=
            LMMC_STATUS_OK ||
        lmmc_std_complex_set_make(complex_b_values, 2, &set_b) !=
            LMMC_STATUS_OK ||
        set_a.size != 2 ||
        lmmc_std_complex_set_contains(&set_a,
                                      &zero_complex,
                                      &contains) != LMMC_STATUS_OK ||
        !contains ||
        complex_set_operations_failed(set_a, set_b, set_union, set_difference)) {
        cleanup_sets();
        std::cerr
            << "failed to call installed LMMC standard library complex set adapters\n";
        return 15;
    }
    cleanup_sets();
    return 0;
}

static int check_num_keys() {
    uint64_t lmmc_hash_z = 0;
    uint64_t lmmc_hash_w = 0;
    int lmmc_num_equal = 0;
    if (lmmc_std_num_equal(2.0, 2.0, &lmmc_num_equal) !=
            LMMC_STATUS_OK ||
        !lmmc_num_equal ||
        lmmc_std_num_hash(0.0, &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_std_num_hash(-0.0, &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC standard library num key adapters\n";
        return 15;
    }
    return 0;
}

static int check_bool_text_keys() {
    uint64_t lmmc_hash_z = 0;
    uint64_t lmmc_hash_w = 0;
    int lmmc_bool_equal = 0;
    int lmmc_text_equal = 0;
    if (lmmc_std_bool_equal(1, 1, &lmmc_bool_equal) != LMMC_STATUS_OK ||
        !lmmc_bool_equal) {
        std::cerr << "failed to call installed LMMC standard library bool/text key adapters\n";
        return 15;
    }
    if (lmmc_std_bool_hash(1, &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_std_bool_hash(0, &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z == lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC standard library bool/text key adapters\n";
        return 15;
    }
    if (lmmc_std_text_equal("alpha", "alpha", &lmmc_text_equal) !=
            LMMC_STATUS_OK ||
        !lmmc_text_equal) {
        std::cerr << "failed to call installed LMMC standard library bool/text key adapters\n";
        return 15;
    }
    if (lmmc_std_text_hash("alpha", &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_std_text_hash("alpha", &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC standard library bool/text key adapters\n";
        return 15;
    }
    return 0;
}

static int check_complex_keys() {
    lmmc_complex_t lmmc_z = {};
    lmmc_complex_t lmmc_w = {};
    uint64_t lmmc_hash_z = 0;
    uint64_t lmmc_hash_w = 0;
    int lmmc_complex_equal = 0;
    if (lmmc_std_math_complex(3.0, 4.0, &lmmc_z) != LMMC_STATUS_OK ||
        lmmc_std_math_complex(3.0, 4.0, &lmmc_w) != LMMC_STATUS_OK ||
        lmmc_std_math_complex_equal(&lmmc_z, &lmmc_w,
                                    &lmmc_complex_equal) != LMMC_STATUS_OK ||
        !lmmc_complex_equal) {
        std::cerr << "failed to call installed LMMC standard library complex key adapters\n";
        return 15;
    }
    if (lmmc_std_math_complex_hash(&lmmc_z, &lmmc_hash_z) !=
            LMMC_STATUS_OK ||
        lmmc_std_math_complex_hash(&lmmc_w, &lmmc_hash_w) !=
            LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC standard library complex key adapters\n";
        return 15;
    }
    return 0;
}

int run_lmmc_set_consumer_checks() {
    const auto checks = {
        check_num_sets,
        check_text_sets,
        check_bool_sets,
        check_complex_sets,
        check_num_keys,
        check_bool_text_keys,
        check_complex_keys,
    };
    for (const auto check : checks) {
        if (const int status = check(); status != 0) {
            return status;
        }
    }
    return 0;
}
