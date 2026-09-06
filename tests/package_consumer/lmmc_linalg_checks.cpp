#include "lmmc/stdlib.h"

#include <cmath>
#include <iostream>
#include <string>

int run_lmmc_linalg_consumer_checks()
{
    lmmc_real_t lmmc_out = 0.0;
    lmmc_real_t a_values[] = {1.0, 2.0, 3.0};
    lmmc_real_t b_values[] = {4.0, 5.0, 6.0};
    lmmc_vec_t a_vec = {3, a_values, 0};
    lmmc_vec_t b_vec = {3, b_values, 0};
    lmmc_vec_t cross = {0, nullptr, 0};
    lmmc_vec_t matvec = {0, nullptr, 0};
    lmmc_vec_t vec_add = {0, nullptr, 0};
    lmmc_vec_t vec_add_scalar = {0, nullptr, 0};
    lmmc_vec_t vec_sub = {0, nullptr, 0};
    lmmc_vec_t vec_sub_scalar = {0, nullptr, 0};
    lmmc_vec_t scalar_sub_vec = {0, nullptr, 0};
    lmmc_vec_t vec_mul = {0, nullptr, 0};
    lmmc_vec_t vec_mul_scalar = {0, nullptr, 0};
    lmmc_vec_t vec_div = {0, nullptr, 0};
    lmmc_vec_t vec_div_scalar = {0, nullptr, 0};
    lmmc_vec_t scalar_div_vec = {0, nullptr, 0};
    lmmc_vec_t vec_pow = {0, nullptr, 0};
    lmmc_vec_t vec_pow_scalar = {0, nullptr, 0};
    lmmc_vec_t vec_scale = {0, nullptr, 0};
    lmmc_vec_t shape_vec = {0, nullptr, 0};
    lmmc_std_bool_vec_t vec_cmp = {0, nullptr, 0};
    lmmc_real_t matrix_values[] = {1.0, 2.0, 3.0, 4.0};
    lmmc_real_t rhs_values[] = {5.0, 6.0, 7.0, 8.0};
    lmmc_real_t short_values[] = {1.0, 2.0};
    lmmc_mat_t matrix = {2, 2, 2, matrix_values, 0};
    lmmc_mat_t rhs_matrix = {2, 2, 2, rhs_values, 0};
    lmmc_vec_t short_vec = {2, short_values, 0};
    lmmc_mat_t matmul = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_add = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_add_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_sub = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_sub_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t scalar_sub_mat = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_mul_elem = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_mul_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_div = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_div_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t scalar_div_mat = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_pow_elem = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_pow_scalar = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_pow_int = {0, 0, 0, nullptr, 0};
    lmmc_mat_t mat_scale = {0, 0, 0, nullptr, 0};
    lmmc_std_eig_table_t eig_table = {};
    lmmc_std_svd_table_t svd_table = {};
    lmmc_std_bool_mat_t mat_cmp = {0, 0, 0, nullptr, 0};
    lmmc_real_t dot = 0.0;
    lmmc_real_t norm = 0.0;
    if (lmmc_std_linalg_dot(&a_vec, &b_vec, &dot) != LMMC_STATUS_OK ||
        dot != 32.0 ||
        lmmc_std_linalg_norm(&a_vec, &norm) != LMMC_STATUS_OK ||
        std::abs(norm - std::sqrt(14.0)) > 1e-12 ||
        lmmc_std_linalg_cross(&a_vec, &b_vec, &cross) != LMMC_STATUS_OK ||
        cross.size != 3 || !cross.data ||
        cross.data[0] != -3.0 || cross.data[1] != 6.0 ||
        cross.data[2] != -3.0 ||
        lmmc_std_linalg_vec_add(&a_vec, &b_vec, &vec_add) !=
            LMMC_STATUS_OK ||
        vec_add.size != 3 || !vec_add.data ||
        vec_add.data[0] != 5.0 || vec_add.data[1] != 7.0 ||
        vec_add.data[2] != 9.0 ||
        lmmc_std_linalg_vec_add_scalar(&a_vec, 10.0, &vec_add_scalar) !=
            LMMC_STATUS_OK ||
        vec_add_scalar.size != 3 || !vec_add_scalar.data ||
        vec_add_scalar.data[0] != 11.0 ||
        vec_add_scalar.data[1] != 12.0 ||
        vec_add_scalar.data[2] != 13.0 ||
        lmmc_std_linalg_vec_sub(&b_vec, &a_vec, &vec_sub) !=
            LMMC_STATUS_OK ||
        vec_sub.size != 3 || !vec_sub.data ||
        vec_sub.data[0] != 3.0 || vec_sub.data[1] != 3.0 ||
        vec_sub.data[2] != 3.0 ||
        lmmc_std_linalg_vec_sub_scalar(&a_vec, 1.0, &vec_sub_scalar) !=
            LMMC_STATUS_OK ||
        vec_sub_scalar.size != 3 || !vec_sub_scalar.data ||
        vec_sub_scalar.data[0] != 0.0 ||
        vec_sub_scalar.data[1] != 1.0 ||
        vec_sub_scalar.data[2] != 2.0 ||
        lmmc_std_linalg_scalar_sub_vec(10.0, &a_vec, &scalar_sub_vec) !=
            LMMC_STATUS_OK ||
        scalar_sub_vec.size != 3 || !scalar_sub_vec.data ||
        scalar_sub_vec.data[0] != 9.0 ||
        scalar_sub_vec.data[1] != 8.0 ||
        scalar_sub_vec.data[2] != 7.0 ||
        lmmc_std_linalg_vec_mul(&a_vec, &b_vec, &vec_mul) !=
            LMMC_STATUS_OK ||
        vec_mul.size != 3 || !vec_mul.data ||
        vec_mul.data[0] != 4.0 || vec_mul.data[1] != 10.0 ||
        vec_mul.data[2] != 18.0 ||
        lmmc_std_linalg_vec_mul_scalar(&a_vec, 2.0, &vec_mul_scalar) !=
            LMMC_STATUS_OK ||
        vec_mul_scalar.size != 3 || !vec_mul_scalar.data ||
        vec_mul_scalar.data[0] != 2.0 ||
        vec_mul_scalar.data[1] != 4.0 ||
        vec_mul_scalar.data[2] != 6.0 ||
        lmmc_std_linalg_vec_div(&b_vec, &a_vec, &vec_div) !=
            LMMC_STATUS_OK ||
        vec_div.size != 3 || !vec_div.data ||
        vec_div.data[0] != 4.0 || vec_div.data[1] != 2.5 ||
        vec_div.data[2] != 2.0 ||
        lmmc_std_linalg_vec_div_scalar(&b_vec, 2.0, &vec_div_scalar) !=
            LMMC_STATUS_OK ||
        vec_div_scalar.size != 3 || !vec_div_scalar.data ||
        vec_div_scalar.data[0] != 2.0 ||
        vec_div_scalar.data[1] != 2.5 ||
        vec_div_scalar.data[2] != 3.0 ||
        lmmc_std_linalg_scalar_div_vec(12.0, &a_vec, &scalar_div_vec) !=
            LMMC_STATUS_OK ||
        scalar_div_vec.size != 3 || !scalar_div_vec.data ||
        scalar_div_vec.data[0] != 12.0 ||
        scalar_div_vec.data[1] != 6.0 ||
        scalar_div_vec.data[2] != 4.0 ||
        lmmc_std_linalg_vec_pow(&a_vec, &short_vec, &vec_pow) !=
            LMMC_STATUS_DIMENSION_MISMATCH ||
        lmmc_std_linalg_vec_pow_scalar(&a_vec, 2.0, &vec_pow_scalar) !=
            LMMC_STATUS_OK ||
        vec_pow_scalar.size != 3 || !vec_pow_scalar.data ||
        vec_pow_scalar.data[0] != 1.0 ||
        vec_pow_scalar.data[1] != 4.0 ||
        vec_pow_scalar.data[2] != 9.0 ||
        lmmc_std_linalg_vec_compare_scalar(&a_vec,
                                           LMMC_STD_COMPARE_GT,
                                           1.0,
                                           &vec_cmp) != LMMC_STATUS_OK ||
        vec_cmp.size != 3 || !vec_cmp.data ||
        vec_cmp.data[0] != 0 || vec_cmp.data[1] != 1 ||
        vec_cmp.data[2] != 1 ||
        lmmc_std_linalg_vec_scale(&a_vec, 2.0, &vec_scale) !=
            LMMC_STATUS_OK ||
        vec_scale.size != 3 || !vec_scale.data ||
        vec_scale.data[0] != 2.0 || vec_scale.data[1] != 4.0 ||
        vec_scale.data[2] != 6.0 ||
        lmmc_std_linalg_matvec(&matrix, &short_vec, &matvec) !=
            LMMC_STATUS_OK ||
        matvec.size != 2 || !matvec.data ||
        matvec.data[0] != 5.0 || matvec.data[1] != 11.0 ||
        lmmc_std_linalg_shape_vec(&matrix, &shape_vec) !=
            LMMC_STATUS_OK ||
        shape_vec.size != 2 || !shape_vec.data ||
        shape_vec.data[0] != 2.0 || shape_vec.data[1] != 2.0 ||
        lmmc_std_linalg_matmul(&matrix, &rhs_matrix, &matmul) !=
            LMMC_STATUS_OK ||
        matmul.rows != 2 || matmul.cols != 2 || !matmul.data ||
        matmul.data[0] != 19.0 || matmul.data[1] != 22.0 ||
        matmul.data[matmul.stride] != 43.0 ||
        matmul.data[matmul.stride + 1] != 50.0 ||
        lmmc_std_linalg_mat_add(&matrix, &rhs_matrix, &mat_add) !=
            LMMC_STATUS_OK ||
        mat_add.rows != 2 || mat_add.cols != 2 || !mat_add.data ||
        mat_add.data[0] != 6.0 || mat_add.data[1] != 8.0 ||
        mat_add.data[mat_add.stride] != 10.0 ||
        mat_add.data[mat_add.stride + 1] != 12.0 ||
        lmmc_std_linalg_mat_add_scalar(&matrix, 10.0, &mat_add_scalar) !=
            LMMC_STATUS_OK ||
        mat_add_scalar.rows != 2 || mat_add_scalar.cols != 2 ||
        !mat_add_scalar.data ||
        mat_add_scalar.data[0] != 11.0 ||
        mat_add_scalar.data[1] != 12.0 ||
        mat_add_scalar.data[mat_add_scalar.stride] != 13.0 ||
        mat_add_scalar.data[mat_add_scalar.stride + 1] != 14.0 ||
        lmmc_std_linalg_mat_sub(&rhs_matrix, &matrix, &mat_sub) !=
            LMMC_STATUS_OK ||
        mat_sub.rows != 2 || mat_sub.cols != 2 || !mat_sub.data ||
        mat_sub.data[0] != 4.0 || mat_sub.data[1] != 4.0 ||
        mat_sub.data[mat_sub.stride] != 4.0 ||
        mat_sub.data[mat_sub.stride + 1] != 4.0 ||
        lmmc_std_linalg_mat_sub_scalar(&matrix, 1.0, &mat_sub_scalar) !=
            LMMC_STATUS_OK ||
        mat_sub_scalar.rows != 2 || mat_sub_scalar.cols != 2 ||
        !mat_sub_scalar.data ||
        mat_sub_scalar.data[0] != 0.0 ||
        mat_sub_scalar.data[1] != 1.0 ||
        mat_sub_scalar.data[mat_sub_scalar.stride] != 2.0 ||
        mat_sub_scalar.data[mat_sub_scalar.stride + 1] != 3.0 ||
        lmmc_std_linalg_scalar_sub_mat(10.0, &matrix, &scalar_sub_mat) !=
            LMMC_STATUS_OK ||
        scalar_sub_mat.rows != 2 || scalar_sub_mat.cols != 2 ||
        !scalar_sub_mat.data ||
        scalar_sub_mat.data[0] != 9.0 ||
        scalar_sub_mat.data[1] != 8.0 ||
        scalar_sub_mat.data[scalar_sub_mat.stride] != 7.0 ||
        scalar_sub_mat.data[scalar_sub_mat.stride + 1] != 6.0 ||
        lmmc_std_linalg_mat_mul_elem(&matrix, &rhs_matrix, &mat_mul_elem) !=
            LMMC_STATUS_OK ||
        mat_mul_elem.rows != 2 || mat_mul_elem.cols != 2 ||
        !mat_mul_elem.data ||
        mat_mul_elem.data[0] != 5.0 || mat_mul_elem.data[1] != 12.0 ||
        mat_mul_elem.data[mat_mul_elem.stride] != 21.0 ||
        mat_mul_elem.data[mat_mul_elem.stride + 1] != 32.0 ||
        lmmc_std_linalg_mat_mul_scalar(&matrix, 2.0, &mat_mul_scalar) !=
            LMMC_STATUS_OK ||
        mat_mul_scalar.rows != 2 || mat_mul_scalar.cols != 2 ||
        !mat_mul_scalar.data ||
        mat_mul_scalar.data[0] != 2.0 ||
        mat_mul_scalar.data[1] != 4.0 ||
        mat_mul_scalar.data[mat_mul_scalar.stride] != 6.0 ||
        mat_mul_scalar.data[mat_mul_scalar.stride + 1] != 8.0 ||
        lmmc_std_linalg_mat_div(&rhs_matrix, &matrix, &mat_div) !=
            LMMC_STATUS_OK ||
        mat_div.rows != 2 || mat_div.cols != 2 || !mat_div.data ||
        mat_div.data[0] != 5.0 || mat_div.data[1] != 3.0 ||
        std::abs(mat_div.data[mat_div.stride] - 7.0 / 3.0) > 1e-12 ||
        mat_div.data[mat_div.stride + 1] != 2.0 ||
        lmmc_std_linalg_mat_div_scalar(&rhs_matrix, 2.0, &mat_div_scalar) !=
            LMMC_STATUS_OK ||
        mat_div_scalar.rows != 2 || mat_div_scalar.cols != 2 ||
        !mat_div_scalar.data ||
        mat_div_scalar.data[0] != 2.5 ||
        mat_div_scalar.data[1] != 3.0 ||
        mat_div_scalar.data[mat_div_scalar.stride] != 3.5 ||
        mat_div_scalar.data[mat_div_scalar.stride + 1] != 4.0 ||
        lmmc_std_linalg_scalar_div_mat(12.0, &matrix, &scalar_div_mat) !=
            LMMC_STATUS_OK ||
        scalar_div_mat.rows != 2 || scalar_div_mat.cols != 2 ||
        !scalar_div_mat.data ||
        scalar_div_mat.data[0] != 12.0 ||
        scalar_div_mat.data[1] != 6.0 ||
        scalar_div_mat.data[scalar_div_mat.stride] != 4.0 ||
        scalar_div_mat.data[scalar_div_mat.stride + 1] != 3.0 ||
        lmmc_std_linalg_mat_pow_elem(&matrix, &rhs_matrix, &mat_pow_elem) !=
            LMMC_STATUS_OK ||
        mat_pow_elem.rows != 2 || mat_pow_elem.cols != 2 ||
        !mat_pow_elem.data ||
        mat_pow_elem.data[0] != 1.0 ||
        mat_pow_elem.data[1] != 64.0 ||
        mat_pow_elem.data[mat_pow_elem.stride] != 2187.0 ||
        mat_pow_elem.data[mat_pow_elem.stride + 1] != 65536.0 ||
        lmmc_std_linalg_mat_pow_scalar(&matrix, 2.0, &mat_pow_scalar) !=
            LMMC_STATUS_OK ||
        mat_pow_scalar.rows != 2 || mat_pow_scalar.cols != 2 ||
        !mat_pow_scalar.data ||
        mat_pow_scalar.data[0] != 1.0 ||
        mat_pow_scalar.data[1] != 4.0 ||
        mat_pow_scalar.data[mat_pow_scalar.stride] != 9.0 ||
        mat_pow_scalar.data[mat_pow_scalar.stride + 1] != 16.0 ||
        lmmc_std_linalg_mat_pow_int(&matrix, 2, &mat_pow_int) !=
            LMMC_STATUS_OK ||
        mat_pow_int.rows != 2 || mat_pow_int.cols != 2 ||
        !mat_pow_int.data ||
        mat_pow_int.data[0] != 7.0 ||
        mat_pow_int.data[1] != 10.0 ||
        mat_pow_int.data[mat_pow_int.stride] != 15.0 ||
        mat_pow_int.data[mat_pow_int.stride + 1] != 22.0 ||
        lmmc_std_linalg_mat_compare_scalar(&matrix,
                                           LMMC_STD_COMPARE_GE,
                                           3.0,
                                           &mat_cmp) != LMMC_STATUS_OK ||
        mat_cmp.rows != 2 || mat_cmp.cols != 2 || !mat_cmp.data ||
        mat_cmp.data[0] != 0 || mat_cmp.data[1] != 0 ||
        mat_cmp.data[mat_cmp.stride] != 1 ||
        mat_cmp.data[mat_cmp.stride + 1] != 1 ||
        lmmc_std_linalg_mat_scale(&matrix, 2.0, &mat_scale) !=
            LMMC_STATUS_OK ||
        mat_scale.rows != 2 || mat_scale.cols != 2 || !mat_scale.data ||
        mat_scale.data[0] != 2.0 || mat_scale.data[1] != 4.0 ||
        mat_scale.data[mat_scale.stride] != 6.0 ||
        mat_scale.data[mat_scale.stride + 1] != 8.0 ||
        lmmc_std_linalg_eig_table(&matrix, &eig_table) !=
            LMMC_STATUS_OK ||
        lmmc_std_eig_table_count(&eig_table) != 4 ||
        std::string(lmmc_std_eig_table_key(&eig_table, 0)
                        ? lmmc_std_eig_table_key(&eig_table, 0)
                        : "") != "values_real" ||
        std::string(lmmc_std_eig_table_key(&eig_table, 3)
                        ? lmmc_std_eig_table_key(&eig_table, 3)
                        : "") != "vectors_imag" ||
        lmmc_std_eig_table_key(&eig_table, 4) != nullptr ||
        lmmc_std_eig_table_get(&eig_table, "values_real") == nullptr ||
        lmmc_std_linalg_svd_table(&matrix, &svd_table) !=
            LMMC_STATUS_OK ||
        lmmc_std_svd_table_count(&svd_table) != 3 ||
        std::string(lmmc_std_svd_table_key(&svd_table, 0)
                        ? lmmc_std_svd_table_key(&svd_table, 0)
                        : "") != "U" ||
        std::string(lmmc_std_svd_table_key(&svd_table, 2)
                        ? lmmc_std_svd_table_key(&svd_table, 2)
                        : "") != "Vt" ||
        lmmc_std_svd_table_key(&svd_table, 3) != nullptr ||
        lmmc_std_svd_table_get(&svd_table, "S") == nullptr ||
        std::string(lmmc_std_error_name(LMMC_STATUS_DIMENSION_MISMATCH)) !=
            "DimensionMismatch" ||
        std::string(lmmc_std_error_name(
            LMMC_STATUS_UNIT_STRIP_INVALID)) != "UnitStripInvalid" ||
        lmmc_std_units_strip_num(1.0, "unknown", &lmmc_out) !=
            LMMC_STATUS_UNIT_STRIP_INVALID ||
        lmmc_std_units_strip_num(1.0, "num<m>", &lmmc_out) !=
            LMMC_STATUS_UNIT_STRIP_LEGACY_SYNTAX ||
        lmmc_std_units_convert_from_si(1.0, "unknown", &lmmc_out) !=
            LMMC_STATUS_INVALID_ARGUMENT) {
        lmmc_vec_destroy(&cross);
        lmmc_vec_destroy(&matvec);
        lmmc_vec_destroy(&vec_add);
        lmmc_vec_destroy(&vec_add_scalar);
        lmmc_vec_destroy(&vec_sub);
        lmmc_vec_destroy(&vec_sub_scalar);
        lmmc_vec_destroy(&scalar_sub_vec);
        lmmc_vec_destroy(&vec_mul);
        lmmc_vec_destroy(&vec_mul_scalar);
        lmmc_vec_destroy(&vec_div);
        lmmc_vec_destroy(&vec_div_scalar);
        lmmc_vec_destroy(&scalar_div_vec);
        lmmc_vec_destroy(&vec_pow);
        lmmc_vec_destroy(&vec_pow_scalar);
        lmmc_vec_destroy(&vec_scale);
        lmmc_vec_destroy(&shape_vec);
        lmmc_std_bool_vec_destroy(&vec_cmp);
        lmmc_mat_destroy(&matmul);
        lmmc_mat_destroy(&mat_add);
        lmmc_mat_destroy(&mat_add_scalar);
        lmmc_mat_destroy(&mat_sub);
        lmmc_mat_destroy(&mat_sub_scalar);
        lmmc_mat_destroy(&scalar_sub_mat);
        lmmc_mat_destroy(&mat_mul_elem);
        lmmc_mat_destroy(&mat_mul_scalar);
        lmmc_mat_destroy(&mat_div);
        lmmc_mat_destroy(&mat_div_scalar);
        lmmc_mat_destroy(&scalar_div_mat);
        lmmc_mat_destroy(&mat_pow_elem);
        lmmc_mat_destroy(&mat_pow_scalar);
        lmmc_mat_destroy(&mat_pow_int);
        lmmc_mat_destroy(&mat_scale);
        lmmc_std_eig_table_destroy(&eig_table);
        lmmc_std_svd_table_destroy(&svd_table);
        lmmc_std_bool_mat_destroy(&mat_cmp);
        std::cerr << "failed to call installed LMMC standard library linalg adapters\n";
        return 15;
    }
    lmmc_vec_destroy(&cross);
    lmmc_vec_destroy(&matvec);
    lmmc_vec_destroy(&vec_add);
    lmmc_vec_destroy(&vec_add_scalar);
    lmmc_vec_destroy(&vec_sub);
    lmmc_vec_destroy(&vec_sub_scalar);
    lmmc_vec_destroy(&scalar_sub_vec);
    lmmc_vec_destroy(&vec_mul);
    lmmc_vec_destroy(&vec_mul_scalar);
    lmmc_vec_destroy(&vec_div);
    lmmc_vec_destroy(&vec_div_scalar);
    lmmc_vec_destroy(&scalar_div_vec);
    lmmc_vec_destroy(&vec_pow);
    lmmc_vec_destroy(&vec_pow_scalar);
    lmmc_vec_destroy(&vec_scale);
    lmmc_vec_destroy(&shape_vec);
    lmmc_std_bool_vec_destroy(&vec_cmp);
    lmmc_mat_destroy(&matmul);
    lmmc_mat_destroy(&mat_add);
    lmmc_mat_destroy(&mat_add_scalar);
    lmmc_mat_destroy(&mat_sub);
    lmmc_mat_destroy(&mat_sub_scalar);
    lmmc_mat_destroy(&scalar_sub_mat);
    lmmc_mat_destroy(&mat_mul_elem);
    lmmc_mat_destroy(&mat_mul_scalar);
    lmmc_mat_destroy(&mat_div);
    lmmc_mat_destroy(&mat_div_scalar);
    lmmc_mat_destroy(&scalar_div_mat);
    lmmc_mat_destroy(&mat_pow_elem);
    lmmc_mat_destroy(&mat_pow_scalar);
    lmmc_mat_destroy(&mat_pow_int);
    lmmc_mat_destroy(&mat_scale);
    lmmc_std_eig_table_destroy(&eig_table);
    lmmc_std_svd_table_destroy(&svd_table);
    lmmc_std_bool_mat_destroy(&mat_cmp);
    return 0;
}
