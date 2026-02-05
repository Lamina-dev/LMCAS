cmake_minimum_required(VERSION 3.10)
project(LMCAS)

# 设置输出目录
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/utf-8 /permissive-)
endif()

# 包含 LAMMP 子项目
add_subdirectory(LAMMP)

# LMCAS 库
add_library(lmcas SHARED
    cas.cpp
    symbolic_core.cpp
    symbolic_simplify.cpp
    symbolic_simplify_ln.cpp
    symbolic_calculus.cpp
    symbolic_matrix.cpp
    symbolic_poly.cpp
)

target_compile_definitions(lmcas PRIVATE LAMINA_CORE_EXPORTS)

if(MSVC)
    target_compile_options(lmcas PRIVATE /utf-8 /permissive-)
endif()

# 包含路径
target_include_directories(lmcas PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/LAMMP/include
)

# 链接 LAMMP 核心库
target_link_libraries(lmcas PRIVATE LammpCore)

# 测试程序
add_executable(test_cas tests/test_main.cpp)
target_link_libraries(test_cas PRIVATE lmcas LammpCore)

add_executable(test_bigint_debug tests/test_bigint_debug.cpp)
target_link_libraries(test_bigint_debug PRIVATE lmcas LammpCore)

add_executable(test_det tests/test_det.cpp)
target_link_libraries(test_det PRIVATE lmcas LammpCore)

add_executable(test_matrix_ops tests/test_matrix_ops.cpp)
target_link_libraries(test_matrix_ops PRIVATE lmcas LammpCore)

add_executable(test_solve tests/test_solve.cpp)
target_link_libraries(test_solve PRIVATE lmcas LammpCore)

add_executable(test_proof tests/test_proof.cpp)
target_link_libraries(test_proof PRIVATE lmcas LammpCore)

add_executable(test_norm tests/test_norm.cpp)
target_link_libraries(test_norm PRIVATE lmcas LammpCore)

add_executable(test_factor tests/test_factor.cpp)
target_link_libraries(test_factor PRIVATE lmcas LammpCore)

add_executable(test_arithmetic tests/test_arithmetic.cpp)
target_link_libraries(test_arithmetic PRIVATE lmcas LammpCore)


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

# --- Benchmarks ---
add_executable(benchmark_all benchmarks/benchmark_all.cpp)
target_link_libraries(benchmark_all PRIVATE lmcas LammpCore)
target_compile_definitions(benchmark_all PRIVATE -D_SYMBOLIC_BENCHMARK)

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
add_test(NAME TestSeries COMMAND test_series)
