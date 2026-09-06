#include "lmmc/stdlib.h"

#include <cmath>
#include <cstdint>
#include <iostream>

int run_lmmc_stdlib_consumer_checks() {
    lmmc_real_t lmmc_out = 0.0;
    lmmc_complex_t lmmc_z = {};
    lmmc_complex_t lmmc_w = {};
    uint64_t lmmc_hash_z = 0;
    uint64_t lmmc_hash_w = 0;
    int lmmc_complex_equal = 0;
    int lmmc_num_equal = 0;
    int lmmc_bool_equal = 0;
    int lmmc_text_equal = 0;
    const char* lmmc_constant_name = nullptr;
    const char* lmmc_constant_unit = nullptr;
    lmmc_real_t lmmc_constant_value = 0.0;
    if (lmmc_std_num_equal(2.0, 2.0, &lmmc_num_equal) !=
            LMMC_STATUS_OK ||
        !lmmc_num_equal ||
        lmmc_std_num_hash(0.0, &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_std_num_hash(-0.0, &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC standard library num key adapters\n";
        return 15;
    }
    {
        lmmc_real_t set_a_values[] = {1.0, 2.0, 2.0, -0.0};
        lmmc_real_t set_b_values[] = {2.0, 3.0, 0.0};
        lmmc_std_num_set_t set_a = {};
        lmmc_std_num_set_t set_b = {};
        lmmc_std_num_set_t set_union = {};
        lmmc_std_num_set_t set_difference = {};
        int contains = 0;
        int subset = 0;
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
            lmmc_std_num_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
            set_union.size != 4 ||
            lmmc_std_num_set_difference(&set_union,
                                        &set_b,
                                        &set_difference) !=
                LMMC_STATUS_OK ||
            set_difference.size != 1 ||
            lmmc_std_num_set_subset(&set_difference, &set_union, &subset) !=
                LMMC_STATUS_OK ||
            !subset ||
            lmmc_std_num_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_sets();
            std::cerr << "failed to call installed LMMC standard library num set adapters\n";
            return 15;
        }
        cleanup_sets();
    }
    if (lmmc_std_bool_equal(1, 1, &lmmc_bool_equal) != LMMC_STATUS_OK ||
        !lmmc_bool_equal ||
        lmmc_std_bool_hash(1, &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_std_bool_hash(0, &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z == lmmc_hash_w ||
        lmmc_std_text_equal("alpha", "alpha", &lmmc_text_equal) !=
            LMMC_STATUS_OK ||
        !lmmc_text_equal ||
        lmmc_std_text_hash("alpha", &lmmc_hash_z) != LMMC_STATUS_OK ||
        lmmc_std_text_hash("alpha", &lmmc_hash_w) != LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC standard library bool/text key adapters\n";
        return 15;
    }
    {
        const char* text_a_values[] = {"alpha", "beta", "alpha"};
        const char* text_b_values[] = {"beta", "gamma"};
        lmmc_std_text_set_t set_a = {};
        lmmc_std_text_set_t set_b = {};
        lmmc_std_text_set_t set_union = {};
        lmmc_std_text_set_t set_difference = {};
        int contains = 0;
        int subset = 0;
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
            lmmc_std_text_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
            set_union.size != 3 ||
            lmmc_std_text_set_difference(&set_union,
                                         &set_b,
                                         &set_difference) !=
                LMMC_STATUS_OK ||
            set_difference.size != 1 ||
            lmmc_std_text_set_subset(&set_difference,
                                     &set_union,
                                     &subset) != LMMC_STATUS_OK ||
            !subset ||
            lmmc_std_text_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_sets();
            std::cerr << "failed to call installed LMMC standard library text set adapters\n";
            return 15;
        }
        cleanup_sets();
    }
    {
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
            !contains ||
            lmmc_std_bool_set_difference(&set_a,
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
    }
    if (lmmc_std_math_complex(3.0, 4.0, &lmmc_z) != LMMC_STATUS_OK ||
        lmmc_std_math_complex(3.0, 4.0, &lmmc_w) != LMMC_STATUS_OK ||
        lmmc_std_math_complex_equal(&lmmc_z, &lmmc_w,
                                    &lmmc_complex_equal) != LMMC_STATUS_OK ||
        !lmmc_complex_equal ||
        lmmc_std_math_complex_hash(&lmmc_z, &lmmc_hash_z) !=
            LMMC_STATUS_OK ||
        lmmc_std_math_complex_hash(&lmmc_w, &lmmc_hash_w) !=
            LMMC_STATUS_OK ||
        lmmc_hash_z != lmmc_hash_w) {
        std::cerr << "failed to call installed LMMC standard library complex key adapters\n";
        return 15;
    }
    {
        lmmc_complex_t complex_a_values[] = {
            {1.0, 2.0}, {1.0, 2.0}, {-0.0, 0.0}};
        lmmc_complex_t complex_b_values[] = {{0.0, -0.0}, {3.0, 4.0}};
        lmmc_std_complex_set_t set_a = {};
        lmmc_std_complex_set_t set_b = {};
        lmmc_std_complex_set_t set_union = {};
        lmmc_std_complex_set_t set_difference = {};
        lmmc_complex_t zero_complex = {0.0, 0.0};
        int contains = 0;
        int subset = 0;
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
            lmmc_std_complex_set_union(&set_a, &set_b, &set_union) !=
                LMMC_STATUS_OK ||
            set_union.size != 3 ||
            lmmc_std_complex_set_difference(&set_union,
                                            &set_b,
                                            &set_difference) !=
                LMMC_STATUS_OK ||
            set_difference.size != 1 ||
            lmmc_std_complex_set_subset(&set_difference,
                                        &set_union,
                                        &subset) != LMMC_STATUS_OK ||
            !subset ||
            lmmc_std_complex_set_make(nullptr, 1, &set_difference) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_sets();
            std::cerr
                << "failed to call installed LMMC standard library complex set adapters\n";
            return 15;
        }
        cleanup_sets();
    }
    {
        lmmc_real_t diagonal_values[] = {1.0, 2.0, 3.0};
        lmmc_vec_t diagonal = {3, diagonal_values, 0};
        lmmc_mat_t eye = {};
        lmmc_mat_t diag = {};
        auto cleanup_mats = [&]() {
            lmmc_mat_destroy(&eye);
            lmmc_mat_destroy(&diag);
        };
        if (lmmc_std_linalg_eye(3, &eye) != LMMC_STATUS_OK ||
            eye.rows != 3 || eye.cols != 3 ||
            std::abs(eye.data[0 * eye.stride + 0] - 1.0) > 1e-12 ||
            std::abs(eye.data[1 * eye.stride + 1] - 1.0) > 1e-12 ||
            std::abs(eye.data[0 * eye.stride + 2]) > 1e-12 ||
            lmmc_std_linalg_diag(&diagonal, &diag) != LMMC_STATUS_OK ||
            diag.rows != 3 || diag.cols != 3 ||
            std::abs(diag.data[0 * diag.stride + 0] - 1.0) > 1e-12 ||
            std::abs(diag.data[1 * diag.stride + 1] - 2.0) > 1e-12 ||
            std::abs(diag.data[2 * diag.stride + 2] - 3.0) > 1e-12 ||
            std::abs(diag.data[0 * diag.stride + 1]) > 1e-12 ||
            lmmc_std_linalg_eye(0, &eye) != LMMC_STATUS_INVALID_ARGUMENT ||
            lmmc_std_linalg_diag(nullptr, &diag) !=
                LMMC_STATUS_INVALID_ARGUMENT) {
            cleanup_mats();
            std::cerr
                << "failed to call installed LMMC standard library matrix constructors\n";
            return 15;
        }
        cleanup_mats();
    }
    if (lmmc_std_constants_count() == 0 ||
        lmmc_std_constants_get("C", &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.99792458e8) > 1e-3 ||
        std::string(lmmc_std_constants_unit("C")) != "m*s^-1" ||
        lmmc_std_constants_entry(7,
                                 &lmmc_constant_name,
                                 &lmmc_constant_value,
                                 &lmmc_constant_unit) != LMMC_STATUS_OK ||
        std::string(lmmc_constant_name ? lmmc_constant_name : "") != "C" ||
        std::abs(lmmc_constant_value - 2.99792458e8) > 1e-3 ||
        std::string(lmmc_constant_unit ? lmmc_constant_unit : "") !=
            "m*s^-1" ||
        lmmc_std_constants_get("NO_SUCH_CONSTANT", &lmmc_out) !=
            LMMC_STATUS_INVALID_ARGUMENT ||
        lmmc_std_constants_unit("NO_SUCH_CONSTANT") != nullptr ||
        lmmc_std_constants_entry(lmmc_std_constants_count(),
                                 &lmmc_constant_name,
                                 &lmmc_constant_value,
                                 &lmmc_constant_unit) !=
            LMMC_STATUS_INVALID_ARGUMENT) {
        std::cerr << "failed to call installed LMMC standard library constants adapters\n";
        return 15;
    }

    lmmc_real_t stats_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t stats_scaled[] = {2.0, 4.0, 6.0, 8.0};
    if (lmmc_std_stats_mean(stats_values, 4, &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12 ||
        lmmc_std_stats_median(stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12 ||
        lmmc_std_stats_quantile(stats_values, 4, 0.5, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12 ||
        lmmc_std_stats_cov(stats_values, stats_scaled, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 3.333333333333333) > 1e-12 ||
        lmmc_std_stats_corr(stats_values, stats_scaled, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 1.0) > 1e-12 ||
        lmmc_std_stats_normal_pdf(0.0, 0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 0.3989422804014327) > 1e-12 ||
        lmmc_std_stats_binomial_pmf(2, 4, 0.5, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 0.375) > 1e-12 ||
        lmmc_std_stats_mean(stats_values, 0, &lmmc_out) !=
            LMMC_STATUS_EMPTY_INPUT) {
        std::cerr << "failed to call installed LMMC standard library stats adapters\n";
        return 15;
    }

    lmmc_rng_t* lmmc_rng = nullptr;
    lmmc_real_t lmmc_rand_a = 0.0;
    lmmc_real_t lmmc_rand_b = 0.0;
    int64_t lmmc_rand_int = 0;
    if (lmmc_rng_create(&lmmc_rng) != LMMC_STATUS_OK ||
        lmmc_std_random_seed(lmmc_rng, 1234) != LMMC_STATUS_OK ||
        lmmc_std_random_rand(lmmc_rng, &lmmc_rand_a) != LMMC_STATUS_OK ||
        lmmc_rand_a < 0.0 || lmmc_rand_a >= 1.0 ||
        lmmc_std_random_seed(lmmc_rng, 1234) != LMMC_STATUS_OK ||
        lmmc_std_random_rand(lmmc_rng, &lmmc_rand_b) != LMMC_STATUS_OK ||
        lmmc_rand_a != lmmc_rand_b ||
        lmmc_std_random_randint(lmmc_rng, 1, 3, &lmmc_rand_int) !=
            LMMC_STATUS_OK ||
        lmmc_rand_int < 1 || lmmc_rand_int > 3 ||
        lmmc_std_random_normal(lmmc_rng, 0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        !std::isfinite(lmmc_out) ||
        lmmc_std_random_choice(lmmc_rng, stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        (lmmc_out != 1.0 && lmmc_out != 2.0 && lmmc_out != 3.0 &&
         lmmc_out != 4.0)) {
        if (lmmc_rng) lmmc_rng_destroy(lmmc_rng);
        std::cerr << "failed to call installed LMMC standard library random adapters\n";
        return 15;
    }
    lmmc_rng_destroy(lmmc_rng);
    if (lmmc_std_random_default_seed(4321) != LMMC_STATUS_OK ||
        lmmc_std_random_default_rand(&lmmc_rand_a) != LMMC_STATUS_OK ||
        lmmc_std_random_default_seed(4321) != LMMC_STATUS_OK ||
        lmmc_std_random_default_rand(&lmmc_rand_b) != LMMC_STATUS_OK ||
        lmmc_rand_a != lmmc_rand_b ||
        lmmc_std_random_default_randint(2, 4, &lmmc_rand_int) !=
            LMMC_STATUS_OK ||
        lmmc_rand_int < 2 || lmmc_rand_int > 4 ||
        lmmc_std_random_default_normal(0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        !std::isfinite(lmmc_out) ||
        lmmc_std_random_default_choice(stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        (lmmc_out != 1.0 && lmmc_out != 2.0 && lmmc_out != 3.0 &&
         lmmc_out != 4.0)) {
        lmmc_std_random_default_deinit();
        std::cerr << "failed to call installed LMMC standard library default random adapters\n";
        return 15;
    }
    lmmc_std_random_default_deinit();

    int lmmc_dimensionless = 0;
    if (lmmc_std_units_strip_num(10.0, "km", &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10000.0) > 1e-12 ||
        lmmc_std_units_convert_from_si(10000.0, "km", &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10.0) > 1e-12 ||
        lmmc_std_units_convert_num(10000.0, "km", &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10.0) > 1e-12 ||
        lmmc_std_units_strip_scalar(10.0, &lmmc_out) != LMMC_STATUS_OK ||
        lmmc_out != 10.0 ||
        lmmc_std_units_is_dimensionless_num(10.0, &lmmc_dimensionless) !=
            LMMC_STATUS_OK ||
        !lmmc_dimensionless) {
        std::cerr << "failed to call installed LMMC standard library unit adapters\n";
        return 15;
    }

    return 0;
}
