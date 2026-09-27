#pragma once
/**
 * @file test_multivariate_factor_support.hpp
 * @brief 多元因式分解测试的共用构造函数。
 *
 * 用于 Hensel 提升前的首项系数分配（Wang's trick）集成测试、
 * 第 9 至 11 项属性测试及完整分解单元测试。
 */

#include "test_multivariate_support.hpp"
#include <rapidcheck.h>
using namespace LMCAS;

inline MultiFactorResult checked_factor_multivariate(const MultiPoly& poly)
{
    auto result = factor_multivariate_checked(poly);
    if (!result) {
        throw std::runtime_error(
            result.error().message + " while factoring " + poly.to_string());
    }
    return std::move(result.value().value);
}
