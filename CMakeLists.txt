cmake_minimum_required(VERSION 3.16)
project(LMCAS VERSION 2.0.0 LANGUAGES C CXX)

include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(LMCAS_BUILD_BENCHMARKS "Build LMCAS benchmarks" OFF)
option(LMCAS_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(LMCAS_ENABLE_TSAN "Enable ThreadSanitizer" OFF)
option(LMCAS_ENABLE_COVERAGE "Enable compiler coverage instrumentation" OFF)
option(LMCAS_WARNINGS_AS_ERRORS "Treat LMCAS compiler warnings as errors" OFF)
set(LMCAS_USE_LAMMP_ASM "AUTO" CACHE STRING "LAMMP assembly backend: AUTO | ON | OFF")
set_property(CACHE LMCAS_USE_LAMMP_ASM PROPERTY STRINGS AUTO ON OFF)

if(LMCAS_ENABLE_SANITIZERS AND LMCAS_ENABLE_TSAN)
    message(FATAL_ERROR "LMCAS_ENABLE_SANITIZERS and LMCAS_ENABLE_TSAN are mutually exclusive")
endif()

if(LMCAS_ENABLE_COVERAGE AND (LMCAS_ENABLE_SANITIZERS OR LMCAS_ENABLE_TSAN))
    message(FATAL_ERROR "LMCAS_ENABLE_COVERAGE cannot be combined with sanitizer instrumentation")
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    message(FATAL_ERROR "LMCAS 2.0 supports GCC, Clang, and MinGW; MSVC is not supported")
endif()

if(WIN32)
    add_compile_definitions(LAMMP_WINDOWS)
endif()

if(LMCAS_USE_LAMMP_ASM STREQUAL "ON")
    set(USE_ASM ON CACHE BOOL "Enable x64 assembly optimization (require NASM)" FORCE)
elseif(LMCAS_USE_LAMMP_ASM STREQUAL "OFF")
    set(USE_ASM OFF CACHE BOOL "Enable x64 assembly optimization (require NASM)" FORCE)
elseif(NOT LMCAS_USE_LAMMP_ASM STREQUAL "AUTO")
    message(FATAL_ERROR "LMCAS_USE_LAMMP_ASM must be AUTO, ON, or OFF")
endif()
add_subdirectory(LMMC/LAMMP EXCLUDE_FROM_ALL)

if(UNIX)
    target_link_libraries(LammpCore PRIVATE m)
endif()

set(LMMC_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LMMC_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(LMMC_BUILD_SHARED OFF CACHE BOOL "" FORCE)
set(LMMC_REQUIRE_LAMMP_STATIC OFF CACHE BOOL "" FORCE)
set(LAMMP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LAMMP" CACHE PATH "" FORCE)
add_subdirectory(LMMC EXCLUDE_FROM_ALL)

set_target_properties(LammpCore PROPERTIES
    PREFIX ""
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
)

add_library(lmcas SHARED
    src/symbolic.cpp
    src/numeric_evaluation.cpp
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
    src/poly_utils.cpp
    src/root_of_utils.cpp
    src/interval.cpp
    src/inequality_solver.cpp
    src/parametric_solver.cpp
    src/symbolic_ode.cpp
    src/symbolic_ode_engine.cpp
    src/symbolic_geometry.cpp
    src/symbolic_implicit_diff.cpp
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
    src/lmmc_lifecycle.cpp
    src/lsr_expr.cpp
)

if(WIN32)
    target_sources(lmcas PRIVATE src/dllmain.c)
endif()

set_target_properties(lmcas PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    CXX_EXTENSIONS OFF
)

target_compile_definitions(lmcas PRIVATE LAMINA_CORE_EXPORTS)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(lmcas PRIVATE -Wall -Wextra -Wpedantic)
    if(LMCAS_WARNINGS_AS_ERRORS)
        target_compile_options(lmcas PRIVATE -Werror)
    endif()
endif()

target_include_directories(lmcas PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LMMC/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LAMMP/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(lmcas
    PUBLIC LammpCore lmmc
)

if(LMCAS_ENABLE_SANITIZERS OR LMCAS_ENABLE_TSAN)
    if(WIN32 OR NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(FATAL_ERROR "The configured sanitizer mode requires GCC or Clang on a non-Windows platform")
    endif()

    if(LMCAS_ENABLE_SANITIZERS)
        set(_lmcas_sanitizer_flags -fsanitize=address,undefined -fno-omit-frame-pointer)
    else()
        set(_lmcas_sanitizer_flags -fsanitize=thread -fno-omit-frame-pointer)
    endif()

    target_compile_options(lmcas PUBLIC ${_lmcas_sanitizer_flags})
    target_link_options(lmcas PUBLIC ${_lmcas_sanitizer_flags})
endif()

if(LMCAS_ENABLE_COVERAGE)
    if(WIN32 OR NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(FATAL_ERROR "Coverage instrumentation requires GCC or Clang on a non-Windows platform")
    endif()

    target_compile_options(lmcas PUBLIC --coverage -O0 -g)
    target_link_options(lmcas PUBLIC --coverage)
endif()

# ============================================================
# Benchmarks
# ============================================================
if(LMCAS_BUILD_BENCHMARKS)
    add_executable(benchmark_all benchmarks/benchmark_all.cpp)
    target_link_libraries(benchmark_all PRIVATE lmcas)
    target_compile_definitions(benchmark_all PRIVATE _SYMBOLIC_BENCHMARK)
endif()

# ============================================================
# Install / package config
# ============================================================
install(TARGETS lmcas
    EXPORT LMCASTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS LammpCore
    EXPORT LMCASTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS lmmc
    EXPORT LMCASTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "symbolic_ast.hpp" EXCLUDE
        PATTERN "visitors" EXCLUDE
        PATTERN "visitors/*" EXCLUDE
)

install(DIRECTORY LMMC/LAMMP/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(DIRECTORY LMMC/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(EXPORT LMCASTargets
    NAMESPACE LMCAS::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LMCAS
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/LMCASConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/LMCASConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LMCAS
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/LMCASConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/LMCASConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/LMCASConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LMCAS
)

# ============================================================
# Tests (CTest)
# ============================================================
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
