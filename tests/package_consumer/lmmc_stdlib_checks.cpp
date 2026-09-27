#include "lmmc/stdlib.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

static int check_constant_lookup() {
    lmmc_real_t lmmc_out = 0.0;
    if (lmmc_std_constants_get("C", &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.99792458e8) > 1e-3 ||
        std::string(lmmc_std_constants_unit("C")) != "m*s^-1" ||
        lmmc_std_constants_get("NO_SUCH_CONSTANT", &lmmc_out) !=
            LMMC_STATUS_INVALID_ARGUMENT ||
        lmmc_std_constants_unit("NO_SUCH_CONSTANT") != nullptr) {
        std::cerr << "failed to call installed LMMC standard library constants adapters\n";
        return 15;
    }
    return 0;
}

static int check_constant_entries() {
    const char* lmmc_constant_name = nullptr;
    const char* lmmc_constant_unit = nullptr;
    lmmc_real_t lmmc_constant_value = 0.0;
    if (lmmc_std_constants_entry(7,
                                 &lmmc_constant_name,
                                 &lmmc_constant_value,
                                 &lmmc_constant_unit) != LMMC_STATUS_OK) {
        std::cerr << "failed to call installed LMMC standard library constants adapters\n";
        return 15;
    }
    if (!lmmc_constant_name ||
        std::strcmp(lmmc_constant_name, "C") != 0 ||
        std::abs(lmmc_constant_value - 2.99792458e8) > 1e-3) {
        std::cerr << "failed to call installed LMMC standard library constants adapters\n";
        return 15;
    }
    if (!lmmc_constant_unit ||
        std::strcmp(lmmc_constant_unit, "m*s^-1") != 0 ||
        lmmc_std_constants_entry(lmmc_std_constants_count(),
                                 &lmmc_constant_name,
                                 &lmmc_constant_value,
                                 &lmmc_constant_unit) !=
            LMMC_STATUS_INVALID_ARGUMENT) {
        std::cerr << "failed to call installed LMMC standard library constants adapters\n";
        return 15;
    }
    return 0;
}

static int check_central_tendency() {
    lmmc_real_t lmmc_out = 0.0;
    lmmc_real_t stats_values[] = {1.0, 2.0, 3.0, 4.0};
    if (lmmc_std_stats_mean(stats_values, 4, &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12 ||
        lmmc_std_stats_median(stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12) {
        std::cerr << "failed to call installed LMMC standard library stats adapters\n";
        return 15;
    }
    if (lmmc_std_stats_quantile(stats_values, 4, 0.5, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 2.5) > 1e-12 ||
        lmmc_std_stats_mean(stats_values, 0, &lmmc_out) !=
            LMMC_STATUS_EMPTY_INPUT) {
        std::cerr << "failed to call installed LMMC standard library stats adapters\n";
        return 15;
    }
    return 0;
}

static int check_covariance_correlation() {
    lmmc_real_t lmmc_out = 0.0;
    lmmc_real_t stats_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t stats_scaled[] = {2.0, 4.0, 6.0, 8.0};
    if (lmmc_std_stats_cov(stats_values, stats_scaled, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 3.333333333333333) > 1e-12 ||
        lmmc_std_stats_corr(stats_values, stats_scaled, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 1.0) > 1e-12) {
        std::cerr << "failed to call installed LMMC standard library stats adapters\n";
        return 15;
    }
    return 0;
}

static int check_probability_distributions() {
    lmmc_real_t lmmc_out = 0.0;
    if (lmmc_std_stats_normal_pdf(0.0, 0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 0.3989422804014327) > 1e-12 ||
        lmmc_std_stats_binomial_pmf(2, 4, 0.5, &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 0.375) > 1e-12) {
        std::cerr << "failed to call installed LMMC standard library stats adapters\n";
        return 15;
    }
    return 0;
}

static int check_units() {
    lmmc_real_t lmmc_out = 0.0;
    int lmmc_dimensionless = 0;
    if (lmmc_std_units_strip_num(10.0, "km", &lmmc_out) != LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10000.0) > 1e-12) {
        std::cerr << "failed to call installed LMMC standard library unit adapters\n";
        return 15;
    }
    if (lmmc_std_units_convert_from_si(10000.0, "km", &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10.0) > 1e-12) {
        std::cerr << "failed to call installed LMMC standard library unit adapters\n";
        return 15;
    }
    if (lmmc_std_units_convert_num(10000.0, "km", &lmmc_out) !=
            LMMC_STATUS_OK ||
        std::abs(lmmc_out - 10.0) > 1e-12) {
        std::cerr << "failed to call installed LMMC standard library unit adapters\n";
        return 15;
    }
    if (lmmc_std_units_strip_scalar(10.0, &lmmc_out) != LMMC_STATUS_OK ||
        lmmc_out != 10.0) {
        std::cerr << "failed to call installed LMMC standard library unit adapters\n";
        return 15;
    }
    if (lmmc_std_units_is_dimensionless_num(10.0, &lmmc_dimensionless) !=
            LMMC_STATUS_OK ||
        !lmmc_dimensionless) {
        std::cerr << "failed to call installed LMMC standard library unit adapters\n";
        return 15;
    }
    return 0;
}

static bool seeded_repeatability_failed(lmmc_rng_t* lmmc_rng) {
    lmmc_real_t lmmc_rand_a = 0.0;
    lmmc_real_t lmmc_rand_b = 0.0;
    if (lmmc_std_random_seed(lmmc_rng, 1234) != LMMC_STATUS_OK ||
        lmmc_std_random_rand(lmmc_rng, &lmmc_rand_a) != LMMC_STATUS_OK ||
        lmmc_rand_a < 0.0 ||
        lmmc_rand_a >= 1.0) {
        return true;
    }
    return lmmc_std_random_seed(lmmc_rng, 1234) != LMMC_STATUS_OK ||
        lmmc_std_random_rand(lmmc_rng, &lmmc_rand_b) != LMMC_STATUS_OK ||
        lmmc_rand_a != lmmc_rand_b;
}

static bool seeded_distributions_failed(lmmc_rng_t* lmmc_rng) {
    lmmc_real_t lmmc_out = 0.0;
    lmmc_real_t stats_values[] = {1.0, 2.0, 3.0, 4.0};
    int64_t lmmc_rand_int = 0;
    if (lmmc_std_random_randint(lmmc_rng, 1, 3, &lmmc_rand_int) !=
            LMMC_STATUS_OK ||
        lmmc_rand_int < 1 ||
        lmmc_rand_int > 3) {
        return true;
    }
    if (lmmc_std_random_normal(lmmc_rng, 0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        !std::isfinite(lmmc_out)) {
        return true;
    }
    return lmmc_std_random_choice(lmmc_rng, stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        (lmmc_out != 1.0 && lmmc_out != 2.0 && lmmc_out != 3.0 &&
         lmmc_out != 4.0);
}

static bool default_repeatability_failed() {
    lmmc_real_t lmmc_rand_a = 0.0;
    lmmc_real_t lmmc_rand_b = 0.0;
    return lmmc_std_random_default_seed(4321) != LMMC_STATUS_OK ||
        lmmc_std_random_default_rand(&lmmc_rand_a) != LMMC_STATUS_OK ||
        lmmc_std_random_default_seed(4321) != LMMC_STATUS_OK ||
        lmmc_std_random_default_rand(&lmmc_rand_b) != LMMC_STATUS_OK ||
        lmmc_rand_a != lmmc_rand_b;
}

static bool default_distributions_failed() {
    lmmc_real_t lmmc_out = 0.0;
    lmmc_real_t stats_values[] = {1.0, 2.0, 3.0, 4.0};
    int64_t lmmc_rand_int = 0;
    if (lmmc_std_random_default_randint(2, 4, &lmmc_rand_int) !=
            LMMC_STATUS_OK ||
        lmmc_rand_int < 2 ||
        lmmc_rand_int > 4) {
        return true;
    }
    if (lmmc_std_random_default_normal(0.0, 1.0, &lmmc_out) !=
            LMMC_STATUS_OK ||
        !std::isfinite(lmmc_out)) {
        return true;
    }
    return lmmc_std_random_default_choice(stats_values, 4, &lmmc_out) !=
            LMMC_STATUS_OK ||
        (lmmc_out != 1.0 && lmmc_out != 2.0 && lmmc_out != 3.0 &&
         lmmc_out != 4.0);
}

static int check_seeded_random() {
    lmmc_rng_t* lmmc_rng = nullptr;
    const bool failed = lmmc_rng_create(&lmmc_rng) != LMMC_STATUS_OK ||
        seeded_repeatability_failed(lmmc_rng) ||
        seeded_distributions_failed(lmmc_rng);
    if (lmmc_rng) {
        lmmc_rng_destroy(lmmc_rng);
    }
    if (failed) {
        std::cerr << "failed installed LMMC seeded random checks\n";
        return 15;
    }
    return 0;
}

static int check_default_random() {
    const bool failed = default_repeatability_failed() ||
        default_distributions_failed();
    lmmc_std_random_default_deinit();
    if (failed) {
        std::cerr << "failed installed LMMC default random checks\n";
        return 15;
    }
    return 0;
}

static int check_identity_matrix() {
    lmmc_mat_t eye = {};
    bool failed = lmmc_std_linalg_eye(3, &eye) != LMMC_STATUS_OK ||
        eye.rows != 3 ||
        eye.cols != 3;
    if (!failed) {
        failed = std::abs(eye.data[0 * eye.stride + 0] - 1.0) > 1e-12 ||
        std::abs(eye.data[1 * eye.stride + 1] - 1.0) > 1e-12 ||
        std::abs(eye.data[0 * eye.stride + 2]) > 1e-12 ||
        lmmc_std_linalg_eye(0, &eye) != LMMC_STATUS_INVALID_ARGUMENT;
    }
    lmmc_mat_destroy(&eye);
    if (failed) {
        std::cerr << "failed installed LMMC identity_matrix checks\n";
        return 15;
    }
    return 0;
}

static int check_diagonal_matrix() {
    lmmc_real_t diagonal_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t diagonal = {3, diagonal_values, 0};
    lmmc_mat_t diag = {};
    bool failed = lmmc_std_linalg_diag(&diagonal, &diag) != LMMC_STATUS_OK;
    if (!failed) {
        failed = diag.rows != 3 ||
        diag.cols != 3;
    }
    if (!failed) {
        failed = std::abs(diag.data[0 * diag.stride + 0] - 1.0) > 1e-12 ||
        std::abs(diag.data[1 * diag.stride + 1] - 2.0) > 1e-12 ||
        std::abs(diag.data[2 * diag.stride + 2] - 3.0) > 1e-12 ||
        std::abs(diag.data[0 * diag.stride + 1]) > 1e-12 ||
        lmmc_std_linalg_diag(nullptr, &diag) !=
                LMMC_STATUS_INVALID_ARGUMENT;
    }
    lmmc_mat_destroy(&diag);
    if (failed) {
        std::cerr << "failed installed LMMC diagonal_matrix checks\n";
        return 15;
    }
    return 0;
}

int run_lmmc_scalar_consumer_checks() {
    const auto checks = {
        check_constant_lookup,
        check_constant_entries,
        check_central_tendency,
        check_covariance_correlation,
        check_probability_distributions,
        check_units,
        check_seeded_random,
        check_default_random,
        check_identity_matrix,
        check_diagonal_matrix,
    };
    for (const auto check : checks) {
        if (const int status = check(); status != 0) {
            return status;
        }
    }
    return 0;
}

int run_lmmc_set_consumer_checks();

int run_lmmc_stdlib_consumer_checks() {
    if (const int status = run_lmmc_set_consumer_checks(); status != 0) {
        return status;
    }
    return run_lmmc_scalar_consumer_checks();
}
