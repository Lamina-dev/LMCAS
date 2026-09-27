#include "bigint.hpp"
#include "computation_context.hpp"
#include "internal/lmmc_lifecycle.hpp"

namespace LMCAS {
namespace {

Result<unsigned int> kernel_argument(const BigInt& value, const char* operation) {
    if (value.is_negative()) {
        return Result<unsigned int>::failure(
            CasErrc::DomainError, "combinatorial arguments must be nonnegative", operation);
    }
    const auto converted = value.try_to_uint64();
    if (!converted || *converted > std::numeric_limits<unsigned int>::max()) {
        return Result<unsigned int>::failure(
            CasErrc::ResourceLimit, "argument exceeds the LMMP unsigned kernel range", operation);
    }
    return static_cast<unsigned int>(*converted);
}

Result<void> reserve_product(const BigInt& largest_factor, unsigned int factors,
                             ComputationContext& context, const char* operation) {
    auto steps = context.consume_steps(factors, operation);
    if (!steps) {
        return steps;
    }
    const auto factor_bits = largest_factor.bit_length();
    if (factors != 0 && factor_bits > std::numeric_limits<std::size_t>::max() / factors) {
        return Result<void>::failure(CasErrc::ResourceLimit,
                                     "combinatorial bit bound overflows", operation);
    }
    return context.require_integer_bits(factor_bits * factors, operation);
}

Result<BigInt> finish_combinatorial(BigInt value, ComputationContext& context,
                                   const char* operation) {
    auto access = context.consume_steps(0, operation);
    if (!access) {
        return Result<BigInt>::failure(access.error());
    }
    auto bits = context.require_integer_bits(value.bit_length(), operation);
    if (!bits) {
        return Result<BigInt>::failure(bits.error());
    }
    return Result<BigInt>::success(std::move(value));
}

BigInt require_combinatorial(Result<BigInt> result) {
    if (result) {
        return std::move(result.value());
    }
    if (result.error().code == CasErrc::DomainError) {
        throw std::domain_error(result.error().message);
    }
    if (result.error().code == CasErrc::ResourceLimit) {
        throw std::length_error(result.error().message);
    }
    throw std::runtime_error(result.error().message);
}

}

Result<BigInt> BigInt::factorial_checked(const BigInt& n,
                                       ComputationContext& context) try {
    constexpr auto operation = "BigInt::factorial_checked";
    auto access = context.consume_steps(1, operation);
    if (!access) {
        return Result<BigInt>::failure(access.error());
    }
    auto argument = kernel_argument(n, operation);
    if (!argument) {
        return Result<BigInt>::failure(argument.error());
    }
    if (argument.value() <= 1) {
        return finish_combinatorial(BigInt(1), context, operation);
    }
    auto budget = reserve_product(n, argument.value(), context, operation);
    if (!budget) {
        return Result<BigInt>::failure(budget.error());
    }
    detail::ensure_lmmc_lifecycle();
    mp_bitcnt_t binary_shift = 0;
    const auto limbs = lmmp_factorial_size_(argument.value(), &binary_shift);
    BigInt result;
    result.realloc_to(limbs);
    result.size_ = lmmp_factorial_(result.data_, binary_shift, limbs, argument.value());
    result.sign_ = POSITIVE;
    result.normalize();
    return finish_combinatorial(std::move(result), context, operation);
} catch (const std::bad_alloc&) {
    return Result<BigInt>::failure(CasErrc::ResourceLimit,
                                  "factorial allocation failed", "BigInt::factorial_checked");
} catch (const std::length_error& error) {
    return Result<BigInt>::failure(CasErrc::ResourceLimit, error.what(), "BigInt::factorial_checked");
}

Result<BigInt> BigInt::npr_checked(const BigInt& n, const BigInt& r,
                                 ComputationContext& context) try {
    constexpr auto operation = "BigInt::npr_checked";
    auto access = context.consume_steps(1, operation);
    if (!access) {
        return Result<BigInt>::failure(access.error());
    }
    auto total = kernel_argument(n, operation);
    auto rank = kernel_argument(r, operation);
    if (!total) {
        return Result<BigInt>::failure(total.error());
    }
    if (!rank) {
        return Result<BigInt>::failure(rank.error());
    }
    if (rank.value() > total.value()) {
        return finish_combinatorial(BigInt(0), context, operation);
    }
    if (rank.value() == 0) {
        return finish_combinatorial(BigInt(1), context, operation);
    }
    auto budget = reserve_product(n, rank.value(), context, operation);
    if (!budget) {
        return Result<BigInt>::failure(budget.error());
    }
    detail::ensure_lmmc_lifecycle();
    mp_bitcnt_t binary_shift = 0;
    const auto limbs = lmmp_nPr_size_(total.value(), rank.value(), &binary_shift);
    BigInt result;
    result.realloc_to(limbs);
    result.size_ = lmmp_nPr_(result.data_, binary_shift, limbs, total.value(), rank.value());
    result.sign_ = POSITIVE;
    result.normalize();
    return finish_combinatorial(std::move(result), context, operation);
} catch (const std::bad_alloc&) {
    return Result<BigInt>::failure(CasErrc::ResourceLimit,
                                  "permutation allocation failed", "BigInt::npr_checked");
} catch (const std::length_error& error) {
    return Result<BigInt>::failure(CasErrc::ResourceLimit, error.what(), "BigInt::npr_checked");
}

Result<BigInt> BigInt::ncr_checked(const BigInt& n, const BigInt& r,
                                 ComputationContext& context) try {
    constexpr auto operation = "BigInt::ncr_checked";
    auto access = context.consume_steps(1, operation);
    if (!access) {
        return Result<BigInt>::failure(access.error());
    }
    auto total = kernel_argument(n, operation);
    auto rank = kernel_argument(r, operation);
    if (!total) {
        return Result<BigInt>::failure(total.error());
    }
    if (!rank) {
        return Result<BigInt>::failure(rank.error());
    }
    if (rank.value() > total.value()) {
        return finish_combinatorial(BigInt(0), context, operation);
    }
    const unsigned int reduced_rank = std::min(rank.value(), total.value() - rank.value());
    if (reduced_rank == 0) {
        return finish_combinatorial(BigInt(1), context, operation);
    }
    auto budget = reserve_product(n, reduced_rank, context, operation);
    if (!budget) {
        return Result<BigInt>::failure(budget.error());
    }
    detail::ensure_lmmc_lifecycle();
    mp_bitcnt_t binary_shift = 0;
    const auto limbs = lmmp_nCr_size_(total.value(), reduced_rank, &binary_shift);
    BigInt result;
    result.realloc_to(limbs);
    result.size_ = lmmp_nCr_(result.data_, binary_shift, limbs, total.value(), reduced_rank);
    result.sign_ = POSITIVE;
    result.normalize();
    return finish_combinatorial(std::move(result), context, operation);
} catch (const std::bad_alloc&) {
    return Result<BigInt>::failure(CasErrc::ResourceLimit,
                                  "combination allocation failed", "BigInt::ncr_checked");
} catch (const std::length_error& error) {
    return Result<BigInt>::failure(CasErrc::ResourceLimit, error.what(), "BigInt::ncr_checked");
}

Result<BigInt> BigInt::factorial_checked(const BigInt& n) {
    ComputationContext context;
    return factorial_checked(n, context);
}

Result<BigInt> BigInt::npr_checked(const BigInt& n, const BigInt& r) {
    ComputationContext context;
    return npr_checked(n, r, context);
}

Result<BigInt> BigInt::ncr_checked(const BigInt& n, const BigInt& r) {
    ComputationContext context;
    return ncr_checked(n, r, context);
}

BigInt BigInt::factorial(const BigInt& n) {
    return require_combinatorial(factorial_checked(n));
}

BigInt BigInt::npr(const BigInt& n, const BigInt& r) {
    return require_combinatorial(npr_checked(n, r));
}

BigInt BigInt::ncr(const BigInt& n, const BigInt& r) {
    return require_combinatorial(ncr_checked(n, r));
}

}
