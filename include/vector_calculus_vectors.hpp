#pragma once

#include "vector_calculus_types.hpp"

namespace LMCAS {

/**
 * @brief 计算两个向量的点积 a·b = ∑aᵢbᵢ。
 * @param[in] a 向量 a
 * @param[in] b 向量 b（维度需与 a 相同）
 * @return 点积标量表达式
 */
LMCAS_API VectorCalculusExprResult dot_product(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算两个三维向量的叉积 a×b。
 * @param[in] a 三维向量 a
 * @param[in] b 三维向量 b
 * @return 叉积向量（三分量）
 */
LMCAS_API VectorCalculusFieldResult cross_product(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算向量 a 在向量 b 上的投影向量 proj_b(a) = (a·b / b·b)·b。
 * @param[in] a 被投影向量
 * @param[in] b 投影方向向量
 * @return 投影向量；b 为零向量时返回零向量
 */
LMCAS_API VectorCalculusFieldResult vector_project(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算向量 a 在向量 b 上的标量投影 a·b / |b|。
 * @param[in] a 被投影向量
 * @param[in] b 投影方向向量
 * @return 标量投影表达式；b 为零向量时返回 nullptr
 */
LMCAS_API VectorCalculusExprResult scalar_project(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算两个向量夹角 arccos(a·b / (|a|·|b|))。
 * @param[in] a 向量 a
 * @param[in] b 向量 b
 * @return 夹角表达式（弧度）；任一向量为零时返回 nullptr
 */
LMCAS_API VectorCalculusExprResult vector_angle_symbolic(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算三个三维向量的混合积 a·(b×c)。
 * @param[in] a 向量 a
 * @param[in] b 向量 b
 * @param[in] c 向量 c
 * @return 混合积标量表达式（等于以 a,b,c 为行的行列式）
 */
LMCAS_API VectorCalculusExprResult mixed_product(
    const VectorField& a, const VectorField& b, const VectorField& c);

}
