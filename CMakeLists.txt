cmake_minimum_required(VERSION 3.10)
project(LMCAS)

# 设置输出目录
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    message(WARNING "MSVC is deprecated, please use GCC/G++ instead.")
    add_compile_options($<$<COMPILE_LANGUAGE:C,CXX>:/utf-8 /permissive->)
endif()

if(WIN32)
    add_compile_definitions(LAMMP_WINDOWS)
endif()

# 包含 LMMC 内部的 LAMMP 子项目
add_subdirectory(LMMC/LAMMP)

# 包含 LMMC 子项目
# 禁用 LMMC 的测试和示例以节省编译时间
set(LMMC_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LMMC_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(LMMC_BUILD_SHARED OFF CACHE BOOL "" FORCE)
set(LMMC_REQUIRE_LAMMP_STATIC OFF CACHE BOOL "" FORCE)
set(LAMMP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LAMMP" CACHE PATH "" FORCE)
add_subdirectory(LMMC)

# 强制将LAMMP的输出目录设置到项目构建目录，方便调试
set_target_properties(LammpCore PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
)

# LMCAS 库
add_library(lmcas SHARED
    src/symbolic.cpp
    src/print_visitor.cpp
    src/solve_visitor.cpp
    src/matcher.cpp
    src/integration.cpp
    src/solver.cpp
    src/dllmain.c
)

target_compile_definitions(lmcas PRIVATE LAMINA_CORE_EXPORTS)

if(MSVC)
    target_compile_options(lmcas PRIVATE /utf-8 /permissive-)
endif()

# 包含路径
target_include_directories(lmcas PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LAMMP/include
    ${CMAKE_CURRENT_SOURCE_DIR}/LMMC/include
)

# 链接 LAMMP 核心库和 LMMC
target_link_libraries(lmcas PUBLIC LammpCore lmmc)

# 测试程序
add_executable(test_cas tests/test_main.cpp)
target_link_libraries(test_cas PRIVATE lmcas LammpCore)

add_executable(test_bigint_debug tests/test_bigint_debug.cpp)
target_link_libraries(test_bigint_debug PRIVATE lmcas LammpCore)

add_executable(test_normalization tests/test_normalization.cpp)
target_link_libraries(test_normalization PRIVATE lmcas LammpCore)

add_executable(test_expand_debug tests/test_expand_debug.cpp)
target_link_libraries(test_expand_debug PRIVATE lmcas LammpCore)

add_executable(test_suite tests/test_suite.cpp)
target_link_libraries(test_suite PRIVATE lmcas LammpCore)

add_executable(test_det tests/test_det.cpp)
target_link_libraries(test_det PRIVATE lmcas LammpCore)

add_executable(test_matrix_ops tests/test_matrix_ops.cpp)
target_link_libraries(test_matrix_ops PRIVATE lmcas LammpCore)

add_executable(test_solve tests/test_solve.cpp)
target_link_libraries(test_solve PRIVATE lmcas LammpCore)

add_executable(test_ast_depth tests/test_ast_depth.cpp)
target_link_libraries(test_ast_depth PRIVATE lmcas LammpCore)

add_executable(test_proof tests/test_proof.cpp)
target_link_libraries(test_proof PRIVATE lmcas LammpCore)

add_executable(test_norm tests/test_norm.cpp)
target_link_libraries(test_norm PRIVATE lmcas LammpCore)

add_executable(test_factor tests/test_factor.cpp)
target_link_libraries(test_factor PRIVATE lmcas LammpCore)

add_executable(test_arithmetic tests/test_arithmetic.cpp)
target_link_libraries(test_arithmetic PRIVATE lmcas LammpCore)

add_executable(test_fraction tests/test_fraction.cpp)
target_link_libraries(test_fraction PRIVATE lmcas LammpCore)

add_executable(test_new_features tests/test_new_features.cpp)
target_link_libraries(test_new_features PRIVATE lmcas LammpCore)

add_executable(test_complex_arithmetic tests/test_complex_arithmetic.cpp)
target_link_libraries(test_complex_arithmetic PRIVATE lmcas LammpCore)

add_executable(test_calculus tests/test_calculus.cpp)
target_link_libraries(test_calculus PRIVATE lmcas LammpCore)

add_executable(test_trig tests/test_trig.cpp)
target_link_libraries(test_trig PRIVATE lmcas LammpCore)

add_executable(test_series tests/test_series.cpp)
target_link_libraries(test_series PRIVATE lmcas LammpCore)

add_executable(test_matcher tests/test_matcher.cpp)
target_link_libraries(test_matcher PRIVATE lmcas LammpCore)

add_executable(test_poly_bridge tests/test_poly_bridge.cpp)
target_link_libraries(test_poly_bridge PRIVATE lmcas LammpCore)

add_executable(test_solve_poly tests/test_solve_visitor_poly.cpp)
target_link_libraries(test_solve_poly PRIVATE lmcas LammpCore)

add_executable(test_solver_linear tests/test_solver_linear.cpp)
target_link_libraries(test_solver_linear PRIVATE lmcas LammpCore)

add_executable(test_solver_poly tests/test_solver_poly.cpp)
target_link_libraries(test_solver_poly PRIVATE lmcas LammpCore)

add_executable(test_polynomial_gcd tests/test_polynomial_gcd.cpp)
target_link_libraries(test_polynomial_gcd PRIVATE lmcas LammpCore)

add_executable(test_number_compare tests/test_number_compare.cpp)
target_link_libraries(test_number_compare PRIVATE lmcas LammpCore)

add_executable(test_integration_fallback tests/test_integration_fallback.cpp)
target_link_libraries(test_integration_fallback PRIVATE lmcas LammpCore)

add_executable(test_bigint_add_carry tests/test_bigint_add_carry.cpp)
target_link_libraries(test_bigint_add_carry PRIVATE lmcas LammpCore)

add_executable(test_value_compare tests/test_value_compare.cpp)
target_link_libraries(test_value_compare PRIVATE lmcas LammpCore)

# add_executable(test_migration tests/test_migration.cpp)
# target_link_libraries(test_migration PRIVATE lmcas LammpCore)

# --- Benchmarks ---
add_executable(benchmark_all benchmarks/benchmark_all.cpp)
target_link_libraries(benchmark_all PRIVATE lmcas LammpCore)
target_compile_definitions(benchmark_all PRIVATE -D_SYMBOLIC_BENCHMARK)

add_executable(test_trig_rewrite tests/test_trig_rewrite.cpp)
target_link_libraries(test_trig_rewrite PRIVATE lmcas LammpCore)

add_executable(test_matcher_comm tests/test_matcher_comm.cpp)
target_link_libraries(test_matcher_comm PRIVATE lmcas LammpCore)

enable_testing()
add_test(NAME TestCAS COMMAND test_cas)
add_test(NAME TestBigIntDebug COMMAND test_bigint_debug)
add_test(NAME TestProof COMMAND test_proof)
add_test(NAME TestNorm COMMAND test_norm)
add_test(NAME TestFactor COMMAND test_factor)
add_test(NAME TestArithmetic COMMAND test_arithmetic)
add_test(NAME TestNewFeatures COMMAND test_new_features)
add_test(NAME TestComplexArithmetic COMMAND test_complex_arithmetic)
add_test(NAME TestCalculus COMMAND test_calculus)
add_test(NAME TestTrig COMMAND test_trig)
add_test(NAME TestTrigRewrite COMMAND test_trig_rewrite)
add_test(NAME TestSeries COMMAND test_series)
# add_test(NAME TestSolveExtended COMMAND test_solve_extended)
add_test(NAME TestASTDepth COMMAND test_ast_depth)
add_test(NAME TestNormalization COMMAND test_normalization)
add_test(NAME TestExpandDebug COMMAND test_expand_debug)
add_test(NAME TestSuite COMMAND test_suite)
add_test(NAME TestDet COMMAND test_det)
add_test(NAME TestMatrixOps COMMAND test_matrix_ops)
add_test(NAME TestSolve COMMAND test_solve)
add_test(NAME TestFraction COMMAND test_fraction)
add_test(NAME TestMatcher COMMAND test_matcher)
add_test(NAME TestPolyBridge COMMAND test_poly_bridge)
add_test(NAME TestSolvePoly COMMAND test_solve_poly)
add_test(NAME TestSolverLinear COMMAND test_solver_linear)
add_test(NAME TestSolverPoly COMMAND test_solver_poly)
add_test(NAME TestPolynomialGcd COMMAND test_polynomial_gcd)
add_test(NAME TestMatcherComm COMMAND test_matcher_comm)
add_test(NAME TestNumberCompare COMMAND test_number_compare)
add_test(NAME TestIntegrationFallback COMMAND test_integration_fallback)
add_test(NAME TestBigIntAddCarry COMMAND test_bigint_add_carry)
add_test(NAME TestValueCompare COMMAND test_value_compare)

# add_executable(debug_recursion debug_recursion.cpp)
# target_link_libraries(debug_recursion PRIVATE lmcas LammpCore)

add_executable(test_new_integration tests/test_new_integration.cpp)
target_link_libraries(test_new_integration PRIVATE lmcas LammpCore)
add_test(NAME TestNewIntegration COMMAND test_new_integration)

add_executable(test_cyclic_integration tests/test_cyclic_integration.cpp)
target_link_libraries(test_cyclic_integration PRIVATE lmcas LammpCore)
add_test(NAME TestCyclicIntegration COMMAND test_cyclic_integration)

add_executable(test_definite_integration tests/test_definite_integration.cpp)
target_link_libraries(test_definite_integration PRIVATE lmcas LammpCore)
add_test(NAME TestDefiniteIntegration COMMAND test_definite_integration)

add_executable(test_singularity tests/test_singularity.cpp)
target_link_libraries(test_singularity PRIVATE lmcas LammpCore)
add_test(NAME TestSingularity COMMAND test_singularity)

