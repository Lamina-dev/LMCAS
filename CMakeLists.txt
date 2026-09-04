cmake_minimum_required(VERSION 3.26)
project(LMCAS VERSION 2.0.0 LANGUAGES C CXX)

include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(LMCAS_BUILD_BENCHMARKS "Build LMCAS benchmarks" OFF)
option(LMCAS_BUILD_STANDALONE_TESTS "Build standalone executable tests" ${BUILD_TESTING})
option(LMCAS_BUILD_LMMC_TESTS "Build nested LMMC tests" ${BUILD_TESTING})
option(LMCAS_BUILD_GTEST_TESTS "Build GoogleTest-based tests" ${BUILD_TESTING})
option(LMCAS_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(LMCAS_ENABLE_TSAN "Enable ThreadSanitizer" OFF)
option(LMCAS_ENABLE_COVERAGE "Enable compiler coverage instrumentation" OFF)
option(LMCAS_WARNINGS_AS_ERRORS "Treat LMCAS compiler warnings as errors" OFF)
option(LMCAS_ENABLE_CLANG_TIDY "Run clang-tidy on LMCAS C++ sources" OFF)
set(LMCAS_LMMP_ASM "AUTO" CACHE STRING "LMMP assembly mode")
set_property(CACHE LMCAS_LMMP_ASM PROPERTY STRINGS GENERIC X64 ARM64 AUTO)

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
    add_compile_definitions(LMMP_WINDOWS)
endif()

if(NOT DEFINED LMMP_BUILD_TESTS)
    set(LMMP_BUILD_TESTS OFF)
endif()
if(NOT DEFINED LMMP_BUILD_EXAMPLES)
    set(LMMP_BUILD_EXAMPLES OFF)
endif()
if(NOT DEFINED LMMP_BUILD_BENCHMARKS)
    set(LMMP_BUILD_BENCHMARKS OFF)
endif()

if(NOT DEFINED LMMC_BUILD_TESTS)
    set(LMMC_BUILD_TESTS "${LMCAS_BUILD_LMMC_TESTS}")
endif()
if(NOT DEFINED LMMC_BUILD_EXAMPLES)
    set(LMMC_BUILD_EXAMPLES OFF)
endif()
if(NOT DEFINED LMMC_BUILD_SHARED)
    set(LMMC_BUILD_SHARED ON)
endif()
if(NOT DEFINED LMMP_ROOT)
    set(LMMP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LMMP")
endif()
if(NOT DEFINED LMMC_WARNINGS_AS_ERRORS)
    set(LMMC_WARNINGS_AS_ERRORS "${LMCAS_WARNINGS_AS_ERRORS}")
endif()
if(NOT DEFINED LMMC_ENABLE_SANITIZERS)
    set(LMMC_ENABLE_SANITIZERS "${LMCAS_ENABLE_SANITIZERS}")
endif()
if(NOT DEFINED LMMC_ENABLE_TSAN)
    set(LMMC_ENABLE_TSAN "${LMCAS_ENABLE_TSAN}")
endif()
if(NOT DEFINED LMMC_ENABLE_COVERAGE)
    set(LMMC_ENABLE_COVERAGE "${LMCAS_ENABLE_COVERAGE}")
endif()
if(NOT DEFINED LMMC_ENABLE_CLANG_TIDY)
    set(LMMC_ENABLE_CLANG_TIDY "${LMCAS_ENABLE_CLANG_TIDY}")
endif()
if(NOT DEFINED LMMC_LMMP_ASM)
    set(LMMC_LMMP_ASM "${LMCAS_LMMP_ASM}")
endif()
if(LMCAS_BUILD_LMMC_TESTS)
    add_subdirectory(LMMC)
else()
    add_subdirectory(LMMC EXCLUDE_FROM_ALL)
endif()

if(NOT TARGET liblmmp)
    message(FATAL_ERROR "LMMP did not define its liblmmp target")
endif()
set_target_properties(liblmmp PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
)
foreach(_lmcas_config DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set_target_properties(liblmmp PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_${_lmcas_config}
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
        LIBRARY_OUTPUT_DIRECTORY_${_lmcas_config}
            "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
        ARCHIVE_OUTPUT_DIRECTORY_${_lmcas_config}
            "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    )
endforeach()

set(LMCAS_CORE_SOURCES
    src/bigint.cpp
    src/expression_analysis.cpp
    src/exact_matrix.cpp
    src/unit.cpp
    src/residual_verification.cpp
    src/quantity.cpp
    src/symbolic_set.cpp
    src/symbolic_core.cpp
    src/normalization_arithmetic.cpp
    src/normalization_functions.cpp
    src/normalization_structural.cpp
    src/numeric_evaluation.cpp
    src/print_visitor.cpp
    src/matcher.cpp
    src/lmmc_lifecycle.cpp
)

set(LMCAS_ASSUMPTION_SOURCES
    src/property_store.cpp
    src/relation_store.cpp
    src/inference_engine_core.cpp
    src/inference_engine_arithmetic.cpp
    src/inference_engine_bounds.cpp
    src/inference_engine_monotonicity.cpp
    src/assumption_context.cpp
    src/query_interface.cpp
)

set(LMCAS_ALGEBRA_SOURCES
    src/symbolic_gcd.cpp
    src/symbolic_algebra_adapter.cpp
    src/symbolic_factor.cpp
    src/multivariate_factor.cpp
    src/multivariate_factor_checked.cpp
    src/multivariate_factor_patterns.cpp
    src/multivariate_hensel.cpp
    src/multivariate_poly.cpp
    src/solver_linear.cpp
    src/solver_groebner.cpp
    src/solver_assumptions.cpp
    src/solve_polynomial.cpp
    src/solve_transcendental.cpp
    src/solve_mixed_transcendental.cpp
    src/solve_visitor.cpp
    src/newton_raphson.cpp
    src/exact_sturm.cpp
    src/exact_algebraic.cpp
    src/exact_root.cpp
    src/complex_root_isolation.cpp
    src/poly_utils.cpp
    src/root_of_utils.cpp
    src/interval.cpp
    src/inequality_solver.cpp
    src/inequality_solver_numeric.cpp
    src/inequality_solver_polynomial.cpp
    src/parametric_solver.cpp
    src/transcendental_factor.cpp
    src/transcendental_rewrite.cpp
    src/zassenhaus_combine.cpp
    src/berlekamp.cpp
    src/hensel_lift.cpp
)

set(LMCAS_CALCULUS_SOURCES
    src/limit_visitor.cpp
    src/limit_visitor_core.cpp
    src/limit_result.cpp
    src/symbolic_calculus_adapter.cpp
    src/integration_table.cpp
    src/integration_basic_strategies.cpp
    src/integration_trigonometric.cpp
    src/integrator.cpp
    src/integration_rational.cpp
    src/integration_special_functions.cpp
    src/multiple_integral.cpp
    src/integration_weierstrass.cpp
    src/integration_trig_substitution.cpp
    src/calculus_utils.cpp
    src/numerical_integration.cpp
)

set(LMCAS_ANALYSIS_SOURCES
    src/symbolic_ode.cpp
    src/ode_classification.cpp
    src/ode_first_order.cpp
    src/ode_linear.cpp
    src/ode_advanced.cpp
    src/symbolic_geometry.cpp
    src/symbolic_implicit_diff.cpp
    src/symbolic_vector_geometry.cpp
    src/symbolic_complex.cpp
    src/vector_calculus_support.cpp
    src/vector_calculus_differential.cpp
    src/vector_calculus_integrals.cpp
    src/vector_calculus_integral_api.cpp
    src/vector_calculus_extrema.cpp
    src/vector_calculus_algebra.cpp
    src/series_engine.cpp
    src/series_convergence.cpp
    src/series_expansions.cpp
    src/transform_engine.cpp
    src/complex_analysis.cpp
    src/differential_geometry.cpp
)

set(LMCAS_LINEAR_ALGEBRA_SOURCES
    src/symbolic_matrix.cpp
    src/matrix_decomposition.cpp
)

set(LMCAS_FACADE_SOURCES
    src/lsr_expr_complex_evaluation.cpp
    src/lsr_expr_equivalence.cpp
    src/lsr_expr_errors.cpp
    src/lsr_expr_parsing.cpp
    src/lsr_expr_sets.cpp
    src/lsr_expr_solving.cpp
    src/lsr_expr_facade.cpp
    src/lsr_expr_arithmetic.cpp
    src/lsr_expr_transform.cpp
    src/lsr_expr_values.cpp
    src/lsr_value_facade.cpp
)

set(LMCAS_COMPONENT_TARGETS
    lmcas_core
    lmcas_assumptions
    lmcas_algebra
    lmcas_calculus
    lmcas_linear_algebra
    lmcas_analysis
    lmcas_facade
)

add_library(lmcas_core OBJECT ${LMCAS_CORE_SOURCES})
add_library(lmcas_assumptions OBJECT ${LMCAS_ASSUMPTION_SOURCES})
add_library(lmcas_algebra OBJECT ${LMCAS_ALGEBRA_SOURCES})
add_library(lmcas_calculus OBJECT ${LMCAS_CALCULUS_SOURCES})
add_library(lmcas_linear_algebra OBJECT ${LMCAS_LINEAR_ALGEBRA_SOURCES})
add_library(lmcas_analysis OBJECT ${LMCAS_ANALYSIS_SOURCES})
add_library(lmcas_facade OBJECT ${LMCAS_FACADE_SOURCES})

add_library(lmcas SHARED)
foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
    target_sources(lmcas PRIVATE $<TARGET_OBJECTS:${component}>)
endforeach()


set_target_properties(lmcas PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    CXX_EXTENSIONS OFF
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
)

foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
    target_compile_definitions(${component} PRIVATE LAMINA_CORE_EXPORTS)
    target_include_directories(${component} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/LMMC/include
        ${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LMMP/include
        ${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LMMP/include/lmmp
    )
    set_target_properties(${component} PROPERTIES
        CXX_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE ON
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN YES
    )
endforeach()

target_compile_definitions(lmcas PRIVATE LAMINA_CORE_EXPORTS)

set(_lmcas_warning_options
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-Wall>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-Wextra>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-Wpedantic>
)
if(LMCAS_ENABLE_CLANG_TIDY)
    find_program(LMCAS_CLANG_TIDY_EXECUTABLE NAMES clang-tidy REQUIRED)
    foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
        set_property(TARGET ${component} PROPERTY CXX_CLANG_TIDY
            "${LMCAS_CLANG_TIDY_EXECUTABLE};--config-file=${CMAKE_CURRENT_SOURCE_DIR}/.clang-tidy")
    endforeach()
endif()
if(LMCAS_WARNINGS_AS_ERRORS)
    list(APPEND _lmcas_warning_options
        $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-Werror>)
endif()
foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
    target_compile_options(${component} PRIVATE ${_lmcas_warning_options})
endforeach()
target_compile_options(lmcas PRIVATE
    $<$<COMPILE_LANG_AND_ID:C,GNU,Clang,AppleClang>:-Wall>
    $<$<COMPILE_LANG_AND_ID:C,GNU,Clang,AppleClang>:-Wextra>
    $<$<COMPILE_LANG_AND_ID:C,GNU,Clang,AppleClang>:-Wpedantic>
    $<$<AND:$<BOOL:${LMCAS_WARNINGS_AS_ERRORS}>,$<COMPILE_LANG_AND_ID:C,GNU,Clang,AppleClang>>:-Werror>
)

target_include_directories(lmcas PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LMMC/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LMMP/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LMMP/include/lmmp>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/lmmp>
)

target_link_libraries(lmcas
    PUBLIC lmmp lmmc
)

function(lmcas_apply_instrumentation target_name)
    if(LMCAS_ENABLE_SANITIZERS OR LMCAS_ENABLE_TSAN)
        if(WIN32 OR
           NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            message(FATAL_ERROR "The configured sanitizer mode requires GCC or Clang for both C and C++ on a non-Windows platform")
        endif()
        if(LMCAS_ENABLE_SANITIZERS)
            set(_lmcas_instrumentation_flags
                -fsanitize=address,undefined -fno-omit-frame-pointer)
        else()
            set(_lmcas_instrumentation_flags
                -fsanitize=thread -fno-omit-frame-pointer)
        endif()
        target_compile_options(${target_name} PRIVATE
            ${_lmcas_instrumentation_flags})
        target_link_options(${target_name} PRIVATE
            ${_lmcas_instrumentation_flags})
    elseif(LMCAS_ENABLE_COVERAGE)
        if(WIN32 OR
           NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            message(FATAL_ERROR "Coverage requires GCC or Clang for both C and C++ on a non-Windows platform")
        endif()
        target_compile_options(${target_name} PRIVATE --coverage -O0 -g)
        target_link_options(${target_name} PRIVATE --coverage)
    endif()
endfunction()

foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
    if(LMCAS_ENABLE_SANITIZERS OR LMCAS_ENABLE_TSAN)
        if(LMCAS_ENABLE_SANITIZERS)
            target_compile_options(${component} PRIVATE
                -fsanitize=address,undefined -fno-omit-frame-pointer)
        else()
            target_compile_options(${component} PRIVATE
                -fsanitize=thread -fno-omit-frame-pointer)
        endif()
    elseif(LMCAS_ENABLE_COVERAGE)
        target_compile_options(${component} PRIVATE --coverage -O0 -g)
    endif()
endforeach()
lmcas_apply_instrumentation(lmcas)

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

install(TARGETS liblmmp
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

file(GLOB _lmcas_public_headers CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
list(FILTER _lmcas_public_headers EXCLUDE REGEX "/symbolic_ast\\.hpp$")
install(FILES ${_lmcas_public_headers}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

file(GLOB _lmmp_public_headers CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LMMP/include/lmmp/*.h")
install(FILES ${_lmmp_public_headers}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lmmp)

install(DIRECTORY LMMC/include/lmmc/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lmmc
    FILES_MATCHING PATTERN "*.h")

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
    add_test(
        NAME source_structure
        COMMAND ${CMAKE_COMMAND}
            -DLMCAS_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CheckSourceStructure.cmake
    )
    add_subdirectory(tests)
endif()
