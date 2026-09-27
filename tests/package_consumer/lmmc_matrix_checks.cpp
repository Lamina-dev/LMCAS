#include "lmmc/stdlib.h"

#include <cmath>
#include <iostream>
#include <string>

static int check_matmul() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_mat_t matmul = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_matmul(&matrix, &rhs_matrix, &matmul) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = matmul.rows != 2 ||
        matmul.cols != 2 ||
        !matmul.data;
    }
    if (!failed) {
        failed = matmul.data[0] != 19.0 ||
        matmul.data[1] != 22.0 ||
        matmul.data[matmul.stride] != 43.0 ||
        matmul.data[matmul.stride + 1] != 50.0;
    }
    lmmc_mat_destroy(&matmul);
    if (failed) {
        std::cerr << "failed installed LMMC matmul check\n";
        return 15;
    }
    return 0;
}

static int check_mat_add() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_mat_t mat_add = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_add(&matrix, &rhs_matrix, &mat_add) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_add.rows != 2 ||
        mat_add.cols != 2 ||
        !mat_add.data;
    }
    if (!failed) {
        failed = mat_add.data[0] != 6.0 ||
        mat_add.data[1] != 8.0 ||
        mat_add.data[mat_add.stride] != 10.0 ||
        mat_add.data[mat_add.stride + 1] != 12.0;
    }
    lmmc_mat_destroy(&mat_add);
    if (failed) {
        std::cerr << "failed installed LMMC mat_add check\n";
        return 15;
    }
    return 0;
}

static int check_mat_add_scalar() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t mat_add_scalar = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_add_scalar(&matrix, 10.0, &mat_add_scalar) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_add_scalar.rows != 2 ||
        mat_add_scalar.cols != 2 ||
        !mat_add_scalar.data;
    }
    if (!failed) {
        failed = mat_add_scalar.data[0] != 11.0 ||
        mat_add_scalar.data[1] != 12.0 ||
        mat_add_scalar.data[mat_add_scalar.stride] != 13.0 ||
        mat_add_scalar.data[mat_add_scalar.stride + 1] != 14.0;
    }
    lmmc_mat_destroy(&mat_add_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC mat_add_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_mat_sub() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_mat_t mat_sub = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_sub(&rhs_matrix, &matrix, &mat_sub) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_sub.rows != 2 ||
        mat_sub.cols != 2 ||
        !mat_sub.data;
    }
    if (!failed) {
        failed = mat_sub.data[0] != 4.0 ||
        mat_sub.data[1] != 4.0 ||
        mat_sub.data[mat_sub.stride] != 4.0 ||
        mat_sub.data[mat_sub.stride + 1] != 4.0;
    }
    lmmc_mat_destroy(&mat_sub);
    if (failed) {
        std::cerr << "failed installed LMMC mat_sub check\n";
        return 15;
    }
    return 0;
}

static int check_mat_sub_scalar() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t mat_sub_scalar = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_sub_scalar(&matrix, 1.0, &mat_sub_scalar) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_sub_scalar.rows != 2 ||
        mat_sub_scalar.cols != 2 ||
        !mat_sub_scalar.data;
    }
    if (!failed) {
        failed = mat_sub_scalar.data[0] != 0.0 ||
        mat_sub_scalar.data[1] != 1.0 ||
        mat_sub_scalar.data[mat_sub_scalar.stride] != 2.0 ||
        mat_sub_scalar.data[mat_sub_scalar.stride + 1] != 3.0;
    }
    lmmc_mat_destroy(&mat_sub_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC mat_sub_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_scalar_sub_mat() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t scalar_sub_mat = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_scalar_sub_mat(10.0, &matrix, &scalar_sub_mat) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = scalar_sub_mat.rows != 2 ||
        scalar_sub_mat.cols != 2 ||
        !scalar_sub_mat.data;
    }
    if (!failed) {
        failed = scalar_sub_mat.data[0] != 9.0 ||
        scalar_sub_mat.data[1] != 8.0 ||
        scalar_sub_mat.data[scalar_sub_mat.stride] != 7.0 ||
        scalar_sub_mat.data[scalar_sub_mat.stride + 1] != 6.0;
    }
    lmmc_mat_destroy(&scalar_sub_mat);
    if (failed) {
        std::cerr << "failed installed LMMC scalar_sub_mat check\n";
        return 15;
    }
    return 0;
}

static int check_mat_mul_elem() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_mat_t mat_mul_elem = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_mul_elem(&matrix, &rhs_matrix, &mat_mul_elem) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_mul_elem.rows != 2 ||
        mat_mul_elem.cols != 2 ||
        !mat_mul_elem.data;
    }
    if (!failed) {
        failed = mat_mul_elem.data[0] != 5.0 ||
        mat_mul_elem.data[1] != 12.0 ||
        mat_mul_elem.data[mat_mul_elem.stride] != 21.0 ||
        mat_mul_elem.data[mat_mul_elem.stride + 1] != 32.0;
    }
    lmmc_mat_destroy(&mat_mul_elem);
    if (failed) {
        std::cerr << "failed installed LMMC mat_mul_elem check\n";
        return 15;
    }
    return 0;
}

static int check_mat_mul_scalar() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t mat_mul_scalar = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_mul_scalar(&matrix, 2.0, &mat_mul_scalar) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_mul_scalar.rows != 2 ||
        mat_mul_scalar.cols != 2 ||
        !mat_mul_scalar.data;
    }
    if (!failed) {
        failed = mat_mul_scalar.data[0] != 2.0 ||
        mat_mul_scalar.data[1] != 4.0 ||
        mat_mul_scalar.data[mat_mul_scalar.stride] != 6.0 ||
        mat_mul_scalar.data[mat_mul_scalar.stride + 1] != 8.0;
    }
    lmmc_mat_destroy(&mat_mul_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC mat_mul_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_mat_div() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_mat_t mat_div = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_div(&rhs_matrix, &matrix, &mat_div) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_div.rows != 2 ||
        mat_div.cols != 2 ||
        !mat_div.data;
    }
    if (!failed) {
        failed = mat_div.data[0] != 5.0 ||
        mat_div.data[1] != 3.0 ||
        std::abs(mat_div.data[mat_div.stride] - 7.0 / 3.0) > 1e-12 ||
        mat_div.data[mat_div.stride + 1] != 2.0;
    }
    lmmc_mat_destroy(&mat_div);
    if (failed) {
        std::cerr << "failed installed LMMC mat_div check\n";
        return 15;
    }
    return 0;
}

static int check_mat_div_scalar() {
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_mat_t mat_div_scalar = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_div_scalar(&rhs_matrix, 2.0, &mat_div_scalar) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_div_scalar.rows != 2 ||
        mat_div_scalar.cols != 2 ||
        !mat_div_scalar.data;
    }
    if (!failed) {
        failed = mat_div_scalar.data[0] != 2.5 ||
        mat_div_scalar.data[1] != 3.0 ||
        mat_div_scalar.data[mat_div_scalar.stride] != 3.5 ||
        mat_div_scalar.data[mat_div_scalar.stride + 1] != 4.0;
    }
    lmmc_mat_destroy(&mat_div_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC mat_div_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_scalar_div_mat() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t scalar_div_mat = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_scalar_div_mat(12.0, &matrix, &scalar_div_mat) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = scalar_div_mat.rows != 2 ||
        scalar_div_mat.cols != 2 ||
        !scalar_div_mat.data;
    }
    if (!failed) {
        failed = scalar_div_mat.data[0] != 12.0 ||
        scalar_div_mat.data[1] != 6.0 ||
        scalar_div_mat.data[scalar_div_mat.stride] != 4.0 ||
        scalar_div_mat.data[scalar_div_mat.stride + 1] != 3.0;
    }
    lmmc_mat_destroy(&scalar_div_mat);
    if (failed) {
        std::cerr << "failed installed LMMC scalar_div_mat check\n";
        return 15;
    }
    return 0;
}

static int check_mat_pow_elem() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_mat_t mat_pow_elem = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_pow_elem(&matrix, &rhs_matrix, &mat_pow_elem) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_pow_elem.rows != 2 ||
        mat_pow_elem.cols != 2 ||
        !mat_pow_elem.data;
    }
    if (!failed) {
        failed = mat_pow_elem.data[0] != 1.0 ||
        mat_pow_elem.data[1] != 64.0 ||
        mat_pow_elem.data[mat_pow_elem.stride] != 2187.0 ||
        mat_pow_elem.data[mat_pow_elem.stride + 1] != 65536.0;
    }
    lmmc_mat_destroy(&mat_pow_elem);
    if (failed) {
        std::cerr << "failed installed LMMC mat_pow_elem check\n";
        return 15;
    }
    return 0;
}

static int check_mat_pow_scalar() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t mat_pow_scalar = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_pow_scalar(&matrix, 2.0, &mat_pow_scalar) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_pow_scalar.rows != 2 ||
        mat_pow_scalar.cols != 2 ||
        !mat_pow_scalar.data;
    }
    if (!failed) {
        failed = mat_pow_scalar.data[0] != 1.0 ||
        mat_pow_scalar.data[1] != 4.0 ||
        mat_pow_scalar.data[mat_pow_scalar.stride] != 9.0 ||
        mat_pow_scalar.data[mat_pow_scalar.stride + 1] != 16.0;
    }
    lmmc_mat_destroy(&mat_pow_scalar);
    if (failed) {
        std::cerr << "failed installed LMMC mat_pow_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_mat_pow_int() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t mat_pow_int = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_pow_int(&matrix, 2, &mat_pow_int) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_pow_int.rows != 2 ||
        mat_pow_int.cols != 2 ||
        !mat_pow_int.data;
    }
    if (!failed) {
        failed = mat_pow_int.data[0] != 7.0 ||
        mat_pow_int.data[1] != 10.0 ||
        mat_pow_int.data[mat_pow_int.stride] != 15.0 ||
        mat_pow_int.data[mat_pow_int.stride + 1] != 22.0;
    }
    lmmc_mat_destroy(&mat_pow_int);
    if (failed) {
        std::cerr << "failed installed LMMC mat_pow_int check\n";
        return 15;
    }
    return 0;
}

static int check_mat_compare_scalar() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_std_bool_mat_t mat_cmp = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_compare_scalar(&matrix,
                                           LMMC_STD_COMPARE_GE,
                                           3.0,
                                           &mat_cmp) != LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_cmp.rows != 2 ||
        mat_cmp.cols != 2 ||
        !mat_cmp.data;
    }
    if (!failed) {
        failed = mat_cmp.data[0] != 0 ||
        mat_cmp.data[1] != 0 ||
        mat_cmp.data[mat_cmp.stride] != 1 ||
        mat_cmp.data[mat_cmp.stride + 1] != 1;
    }
    lmmc_std_bool_mat_destroy(&mat_cmp);
    if (failed) {
        std::cerr << "failed installed LMMC mat_compare_scalar check\n";
        return 15;
    }
    return 0;
}

static int check_mat_scale() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t mat_scale = {0, 0, 0, nullptr, 0};
    bool failed = lmmc_std_linalg_mat_scale(&matrix, 2.0, &mat_scale) !=
            LMMC_STATUS_OK;
    if (!failed) {
        failed = mat_scale.rows != 2 ||
        mat_scale.cols != 2 ||
        !mat_scale.data;
    }
    if (!failed) {
        failed = mat_scale.data[0] != 2.0 ||
        mat_scale.data[1] != 4.0 ||
        mat_scale.data[mat_scale.stride] != 6.0 ||
        mat_scale.data[mat_scale.stride + 1] != 8.0;
    }
    lmmc_mat_destroy(&mat_scale);
    if (failed) {
        std::cerr << "failed installed LMMC mat_scale check\n";
        return 15;
    }
    return 0;
}

static int check_eig_table() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_std_eig_table_t eig_table = {};
    bool failed = lmmc_std_linalg_eig_table(&matrix, &eig_table) !=
            LMMC_STATUS_OK ||
        lmmc_std_eig_table_count(&eig_table) != 4;
    if (!failed) {
        failed = std::string(lmmc_std_eig_table_key(&eig_table, 0)
                        ? lmmc_std_eig_table_key(&eig_table, 0)
                        : "") != "values_real";
    }
    if (!failed) {
        failed = std::string(lmmc_std_eig_table_key(&eig_table, 3)
                        ? lmmc_std_eig_table_key(&eig_table, 3)
                        : "") != "vectors_imag" ||
        lmmc_std_eig_table_key(&eig_table, 4) != nullptr ||
        lmmc_std_eig_table_get(&eig_table, "values_real") == nullptr;
    }
    lmmc_std_eig_table_destroy(&eig_table);
    if (failed) {
        std::cerr << "failed installed LMMC eig_table check\n";
        return 15;
    }
    return 0;
}

static int check_svd_table() {
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_std_svd_table_t svd_table = {};
    bool failed = lmmc_std_linalg_svd_table(&matrix, &svd_table) !=
            LMMC_STATUS_OK ||
        lmmc_std_svd_table_count(&svd_table) != 3;
    if (!failed) {
        failed = std::string(lmmc_std_svd_table_key(&svd_table, 0)
                        ? lmmc_std_svd_table_key(&svd_table, 0)
                        : "") != "U";
    }
    if (!failed) {
        failed = std::string(lmmc_std_svd_table_key(&svd_table, 2)
                        ? lmmc_std_svd_table_key(&svd_table, 2)
                        : "") != "Vt" ||
        lmmc_std_svd_table_key(&svd_table, 3) != nullptr ||
        lmmc_std_svd_table_get(&svd_table, "S") == nullptr;
    }
    lmmc_std_svd_table_destroy(&svd_table);
    if (failed) {
        std::cerr << "failed installed LMMC svd_table check\n";
        return 15;
    }
    return 0;
}

int run_lmmc_matrix_consumer_checks() {
    const auto checks = {
        check_matmul,
        check_mat_add,
        check_mat_add_scalar,
        check_mat_sub,
        check_mat_sub_scalar,
        check_scalar_sub_mat,
        check_mat_mul_elem,
        check_mat_mul_scalar,
        check_mat_div,
        check_mat_div_scalar,
        check_scalar_div_mat,
        check_mat_pow_elem,
        check_mat_pow_scalar,
        check_mat_pow_int,
        check_mat_compare_scalar,
        check_mat_scale,
        check_eig_table,
        check_svd_table,
    };
    for (const auto check : checks) {
        if (const int status = check(); status != 0) {
            return status;
        }
    }
    return 0;
}
