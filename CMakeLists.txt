cmake_minimum_required(VERSION 3.10)
project(LMCAS)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    message(WARNING
        "DO NOT USE MSVC"
    )
    add_compile_options($<$<COMPILE_LANGUAGE:C,CXX>:/utf-8>)
    add_compile_options($<$<COMPILE_LANGUAGE:C,CXX>:/permissive->)
endif()

if(WIN32)
    add_compile_definitions(LAMMP_WINDOWS)
endif()

add_subdirectory(LMMC/LAMMP)

set(LMMC_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LMMC_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(LMMC_BUILD_SHARED OFF CACHE BOOL "" FORCE)
set(LMMC_REQUIRE_LAMMP_STATIC OFF CACHE BOOL "" FORCE)
set(LAMMP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LAMMP" CACHE PATH "" FORCE)
add_subdirectory(LMMC)

set_target_properties(LammpCore PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
)

add_library(lmcas SHARED
    src/symbolic.cpp
    src/print_visitor.cpp
    src/limit_visitor.cpp
    src/solve_visitor.cpp
    src/matcher.cpp
    src/multivariate_factor.cpp
    src/multivariate_poly.cpp
    src/integration.cpp
    src/solver.cpp
    src/solve_polynomial.cpp
    src/solve_transcendental.cpp
    src/solve_mixed_transcendental.cpp
    src/newton_raphson.cpp
    src/root_of_utils.cpp
    src/interval.cpp
    src/inequality_solver.cpp
    src/parametric_solver.cpp
    src/symbolic_ode.cpp
    src/symbolic_ode_engine.cpp
    src/symbolic_geometry.cpp
    src/symbolic_vector_geometry.cpp
    src/symbolic_complex.cpp
    src/symbolic_matrix.cpp
    src/property_store.cpp
    src/relation_store.cpp
    src/inference_engine.cpp
    src/assumption_context.cpp
    src/query_interface.cpp
    src/transcendental_factor.cpp
    src/calculus_utils.cpp
    src/vector_calculus.cpp
    src/series_engine.cpp
    src/transform_engine.cpp
    src/berlekamp.cpp
    src/hensel_lift.cpp
    src/matrix_decomposition.cpp
    src/complex_analysis.cpp
    src/differential_geometry.cpp
    src/numerical_integration.cpp
    src/dllmain.c
)

target_compile_definitions(lmcas PRIVATE LAMINA_CORE_EXPORTS)

if(MSVC)
    target_compile_options(lmcas PRIVATE /utf-8 /permissive-)
endif()

target_include_directories(lmcas PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LAMMP/include
    ${CMAKE_CURRENT_SOURCE_DIR}/LMMC/include
)

target_link_libraries(lmcas PUBLIC LammpCore lmmc)

# ============================================================
# Benchmarks
# ============================================================
add_executable(benchmark_all benchmarks/benchmark_all.cpp)
target_link_libraries(benchmark_all PRIVATE lmcas LammpCore)
target_compile_definitions(benchmark_all PRIVATE -D_SYMBOLIC_BENCHMARK)

# ============================================================
# Tests (CTest)
# ============================================================
enable_testing()
add_subdirectory(tests)
