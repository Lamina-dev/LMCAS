#include "lmmc/stdlib.h"
#include "lmmc/itersolve.h"

#include <cmath>
#include <iostream>
#include <string>

static int check_dot() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_real_t dot = 0.0;
    const bool failed = lmmc_std_linalg_dot(&a_vec, &b_vec, &dot) != LMMC_STATUS_OK ||
        dot != 32.0;
    if (failed) {
        std::cerr << "failed installed LMMC dot check\n";
        return 15;
    }
    return 0;
}

static int check_norm() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_real_t norm = 0.0;
    const bool failed = lmmc_std_linalg_norm(&a_vec, &norm) != LMMC_STATUS_OK ||
        std::abs(norm - std::sqrt(14.0)) > 1e-12;
    if (failed) {
        std::cerr << "failed installed LMMC norm check\n";
        return 15;
    }
    return 0;
}

static int check_cross() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_vec_t cross = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_cross(&a_vec, &b_vec, &cross) != LMMC_STATUS_OK ||
        cross.size != 3 ||
        !cross.data ||
        cross.data[0] != -3.0 ||
        cross.data[1] != 6.0 ||
        cross.data[2] != -3.0;
    lmmc_vec_destroy(&cross);
    if (failed) {
        std::cerr << "failed installed LMMC cross check\n";
        return 15;
    }
    return 0;
}

static int check_vec_add() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_vec_t vec_add = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_add(&a_vec, &b_vec, &vec_add) !=
            LMMC_STATUS_OK ||
        vec_add.size != 3 ||
        !vec_add.data ||
        vec_add.data[0] != 5.0 ||
        vec_add.data[1] != 7.0 ||
        vec_add.data[2] != 9.0;
    lmmc_vec_destroy(&vec_add);
    if (failed) {
        std::cerr << "failed installed LMMC vec_add check\n";
        return 15;
    }
    return 0;
}

static int check_vec_add_scalar() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t vec_add_scalar = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_add_scalar(&a_vec, 10.0, &vec_add_scalar) !=
            LMMC_STATUS_OK ||
        vec_add_scalar.size != 3 ||
        !vec_add_scalar.data ||
        vec_add_scalar.data[0] != 11.0 ||
        vec_add_scalar.data[1] != 12.0 ||
        vec_add_scalar.data[2] != 13.0;
    lmmc_vec_destroy(&vec_add_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC vec_add_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_vec_sub() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_vec_t vec_sub = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_sub(&b_vec, &a_vec, &vec_sub) !=
            LMMC_STATUS_OK ||
        vec_sub.size != 3 ||
        !vec_sub.data ||
        vec_sub.data[0] != 3.0 ||
        vec_sub.data[1] != 3.0 ||
        vec_sub.data[2] != 3.0;
    lmmc_vec_destroy(&vec_sub);
    if (failed) {
        std::cerr << "failed installed LMMC vec_sub check\n";
        return 15;
    }
    return 0;
}

static int check_vec_sub_scalar() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t vec_sub_scalar = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_sub_scalar(&a_vec, 1.0, &vec_sub_scalar) !=
            LMMC_STATUS_OK ||
        vec_sub_scalar.size != 3 ||
        !vec_sub_scalar.data ||
        vec_sub_scalar.data[0] != 0.0 ||
        vec_sub_scalar.data[1] != 1.0 ||
        vec_sub_scalar.data[2] != 2.0;
    lmmc_vec_destroy(&vec_sub_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC vec_sub_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_scalar_sub_vec() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t scalar_sub_vec = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_scalar_sub_vec(10.0, &a_vec, &scalar_sub_vec) !=
            LMMC_STATUS_OK ||
        scalar_sub_vec.size != 3 ||
        !scalar_sub_vec.data ||
        scalar_sub_vec.data[0] != 9.0 ||
        scalar_sub_vec.data[1] != 8.0 ||
        scalar_sub_vec.data[2] != 7.0;
    lmmc_vec_destroy(&scalar_sub_vec);
    if (failed) {
        std::cerr << "failed installed LMMC scalar_sub_vec check\n";
        return 15;
    }
    return 0;
}

static int check_vec_mul() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_vec_t vec_mul = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_mul(&a_vec, &b_vec, &vec_mul) !=
            LMMC_STATUS_OK ||
        vec_mul.size != 3 ||
        !vec_mul.data ||
        vec_mul.data[0] != 4.0 ||
        vec_mul.data[1] != 10.0 ||
        vec_mul.data[2] != 18.0;
    lmmc_vec_destroy(&vec_mul);
    if (failed) {
        std::cerr << "failed installed LMMC vec_mul check\n";
        return 15;
    }
    return 0;
}

static int check_vec_mul_scalar() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t vec_mul_scalar = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_mul_scalar(&a_vec, 2.0, &vec_mul_scalar) !=
            LMMC_STATUS_OK ||
        vec_mul_scalar.size != 3 ||
        !vec_mul_scalar.data ||
        vec_mul_scalar.data[0] != 2.0 ||
        vec_mul_scalar.data[1] != 4.0 ||
        vec_mul_scalar.data[2] != 6.0;
    lmmc_vec_destroy(&vec_mul_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC vec_mul_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_vec_div() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_vec_t vec_div = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_div(&b_vec, &a_vec, &vec_div) !=
            LMMC_STATUS_OK ||
        vec_div.size != 3 ||
        !vec_div.data ||
        vec_div.data[0] != 4.0 ||
        vec_div.data[1] != 2.5 ||
        vec_div.data[2] != 2.0;
    lmmc_vec_destroy(&vec_div);
    if (failed) {
        std::cerr << "failed installed LMMC vec_div check\n";
        return 15;
    }
    return 0;
}

static int check_vec_div_scalar() {
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_vec_t vec_div_scalar = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_div_scalar(&b_vec, 2.0, &vec_div_scalar) !=
            LMMC_STATUS_OK ||
        vec_div_scalar.size != 3 ||
        !vec_div_scalar.data ||
        vec_div_scalar.data[0] != 2.0 ||
        vec_div_scalar.data[1] != 2.5 ||
        vec_div_scalar.data[2] != 3.0;
    lmmc_vec_destroy(&vec_div_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC vec_div_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_scalar_div_vec() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t scalar_div_vec = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_scalar_div_vec(12.0, &a_vec, &scalar_div_vec) !=
            LMMC_STATUS_OK ||
        scalar_div_vec.size != 3 ||
        !scalar_div_vec.data ||
        scalar_div_vec.data[0] != 12.0 ||
        scalar_div_vec.data[1] != 6.0 ||
        scalar_div_vec.data[2] != 4.0;
    lmmc_vec_destroy(&scalar_div_vec);
    if (failed) {
        std::cerr << "failed installed LMMC scalar_div_vec check\n";
        return 15;
    }
    return 0;
}

static int check_vec_pow() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t vec_pow = {0, nullptr, 0};
    lmmc_real_t short_values[] = {1.0, 2.0};
    lmmc_vec_t short_vec = {2, short_values, 0};
    const bool failed = lmmc_std_linalg_vec_pow(&a_vec, &short_vec, &vec_pow) !=
            LMMC_STATUS_DIMENSION_MISMATCH;
    lmmc_vec_destroy(&vec_pow);
    if (failed) {
        std::cerr << "failed installed LMMC vec_pow check\n";
        return 15;
    }
    return 0;
}

static int check_vec_pow_scalar() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t vec_pow_scalar = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_pow_scalar(&a_vec, 2.0, &vec_pow_scalar) !=
            LMMC_STATUS_OK ||
        vec_pow_scalar.size != 3 ||
        !vec_pow_scalar.data ||
        vec_pow_scalar.data[0] != 1.0 ||
        vec_pow_scalar.data[1] != 4.0 ||
        vec_pow_scalar.data[2] != 9.0;
    lmmc_vec_destroy(&vec_pow_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC vec_pow_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_vec_compare_scalar() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_std_bool_vec_t vec_cmp = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_compare_scalar(&a_vec,
                                           LMMC_STD_COMPARE_GT,
                                           1.0,
                                           &vec_cmp) != LMMC_STATUS_OK ||
        vec_cmp.size != 3 ||
        !vec_cmp.data ||
        vec_cmp.data[0] != 0 ||
        vec_cmp.data[1] != 1 ||
        vec_cmp.data[2] != 1;
    lmmc_std_bool_vec_destroy(&vec_cmp);
    if (failed) {
        std::cerr << "failed installed LMMC vec_compare_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_vec_scale() {
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t vec_scale = {0, nullptr, 0};
    const bool failed = lmmc_std_linalg_vec_scale(&a_vec, 2.0, &vec_scale) !=
            LMMC_STATUS_OK ||
        vec_scale.size != 3 ||
        !vec_scale.data ||
        vec_scale.data[0] != 2.0 ||
        vec_scale.data[1] != 4.0 ||
        vec_scale.data[2] != 6.0;
    lmmc_vec_destroy(&vec_scale);
    if (failed) {
        std::cerr << "failed installed LMMC vec_scale check\n";
        return 15;
    }
    return 0;
}

static int check_matvec() {
    lmmc_vec_t matvec = {0, nullptr, 0};
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t short_values[] = {1.0, 2.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_vec_t short_vec = {2, short_values, 0};
    const bool failed = lmmc_std_linalg_matvec(&matrix, &short_vec, &matvec) !=
            LMMC_STATUS_OK ||
        matvec.size != 2 ||
        !matvec.data ||
        matvec.data[0] != 5.0 ||
        matvec.data[1] != 11.0;
    lmmc_vec_destroy(&matvec);
    if (failed) {
        std::cerr << "failed installed LMMC matvec check\n";
        return 15;
    }
    return 0;
}

static int check_shape_vec() {
    lmmc_vec_t shape_vec = {0, nullptr, 0};
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    const bool failed = lmmc_std_linalg_shape_vec(&matrix, &shape_vec) !=
            LMMC_STATUS_OK ||
        shape_vec.size != 2 ||
        !shape_vec.data ||
        shape_vec.data[0] != 2.0 ||
        shape_vec.data[1] != 2.0;
    lmmc_vec_destroy(&shape_vec);
    if (failed) {
        std::cerr << "failed installed LMMC shape_vec check\n";
        return 15;
    }
    return 0;
}

static int check_dimension_mismatch_name() {
    const bool failed = std::string(lmmc_std_error_name(LMMC_STATUS_DIMENSION_MISMATCH)) !=
            "DimensionMismatch";
    if (failed) {
        std::cerr << "failed installed LMMC dimension_mismatch_name check\n";
        return 15;
    }
    return 0;
}

static int check_unit_strip_error_name() {
    const bool failed = std::string(lmmc_std_error_name(
            LMMC_STATUS_UNIT_STRIP_INVALID)) != "UnitStripInvalid";
    if (failed) {
        std::cerr << "failed installed LMMC unit_strip_error_name check\n";
        return 15;
    }
    return 0;
}

static int check_unknown_strip_unit() {
    lmmc_real_t lmmc_out = 0.0;
    const bool failed = lmmc_std_units_strip_num(1.0, "unknown", &lmmc_out) !=
            LMMC_STATUS_UNIT_STRIP_INVALID;
    if (failed) {
        std::cerr << "failed installed LMMC unknown_strip_unit check\n";
        return 15;
    }
    return 0;
}

static int check_legacy_strip_syntax() {
    lmmc_real_t lmmc_out = 0.0;
    const bool failed = lmmc_std_units_strip_num(1.0, "num<m>", &lmmc_out) !=
            LMMC_STATUS_UNIT_STRIP_LEGACY_SYNTAX;
    if (failed) {
        std::cerr << "failed installed LMMC legacy_strip_syntax check\n";
        return 15;
    }
    return 0;
}

static int check_unknown_conversion_unit() {
    lmmc_real_t lmmc_out = 0.0;
    const bool failed = lmmc_std_units_convert_from_si(1.0, "unknown", &lmmc_out) !=
            LMMC_STATUS_INVALID_ARGUMENT;
    if (failed) {
        std::cerr << "failed installed LMMC unknown_conversion_unit check\n";
        return 15;
    }
    return 0;
}

static lmmc_status_t lsqr_forward(
    const lmmc_vec_t* input, lmmc_vec_t* output, void*) {
    if (input->size != 2 || output->size != 3) {
        return LMMC_STATUS_DIMENSION_MISMATCH;
    }
    output->data[0] = input->data[0];
    output->data[1] = input->data[1];
    output->data[2] = input->data[0] + input->data[1];
    return LMMC_STATUS_OK;
}

static lmmc_status_t lsqr_transpose(
    const lmmc_vec_t* input, lmmc_vec_t* output, void*) {
    if (input->size != 3 || output->size != 2) {
        return LMMC_STATUS_DIMENSION_MISMATCH;
    }
    output->data[0] = input->data[0] + input->data[2];
    output->data[1] = input->data[1] + input->data[2];
    return LMMC_STATUS_OK;
}

static int check_matrix_free_lsqr() {
    lmmc_real_t rhs_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t solution_values[] = {0.0, 0.0};
    lmmc_vec_t rhs = {3, rhs_values, 0};
    lmmc_vec_t solution = {2, solution_values, 0};
    lmmc_itersolve_config_t config{};
    lmmc_itersolve_result_t result{};
    if (lmmc_itersolve_default_config(2, &config) != LMMC_STATUS_OK) {
        std::cerr << "failed installed LMMC LSQR config check\n";
        return 15;
    }
    config.abs_tol = 1e-12;
    config.rel_tol = 1e-12;
    config.apply_op = lsqr_forward;
    config.apply_transpose_op = lsqr_transpose;
    const bool failed =
        lmmc_lsqr_solve(nullptr, &rhs, &config, &solution, &result) !=
            LMMC_STATUS_OK ||
        !result.converged ||
        std::abs(solution.data[0] - 1.0) > 1e-10 ||
        std::abs(solution.data[1] - 2.0) > 1e-10;
    if (failed) {
        std::cerr << "failed installed LMMC matrix-free LSQR check\n";
        return 15;
    }
    return 0;
}

int run_lmmc_vector_consumer_checks() {
    const auto checks = {
        check_dot,
        check_norm,
        check_cross,
        check_vec_add,
        check_vec_add_scalar,
        check_vec_sub,
        check_vec_sub_scalar,
        check_scalar_sub_vec,
        check_vec_mul,
        check_vec_mul_scalar,
        check_vec_div,
        check_vec_div_scalar,
        check_scalar_div_vec,
        check_vec_pow,
        check_vec_pow_scalar,
        check_vec_compare_scalar,
        check_vec_scale,
        check_matvec,
        check_shape_vec,
        check_dimension_mismatch_name,
        check_unit_strip_error_name,
        check_unknown_strip_unit,
        check_legacy_strip_syntax,
        check_unknown_conversion_unit,
        check_matrix_free_lsqr,
    };
    for (const auto check : checks) {
        if (const int status = check(); status != 0) {
            return status;
        }
    }
    return 0;
}

int run_lmmc_matrix_consumer_checks();

int run_lmmc_linalg_consumer_checks() {
    if (const int status = run_lmmc_vector_consumer_checks(); status != 0) {
        return status;
    }
    return run_lmmc_matrix_consumer_checks();
}
