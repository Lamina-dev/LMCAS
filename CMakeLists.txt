cmake_minimum_required(VERSION 3.26)
project(LMCAS LANGUAGES C CXX)

option(BUILD_TESTING "Build development tests" OFF)
include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

configure_file(Doxyfile.in Doxyfile @ONLY)

option(LMCAS_BUILD_BENCHMARKS "Build LMCAS benchmarks" OFF)
option(LMCAS_BUILD_TESTS "Build LMCAS GoogleTest and RapidCheck tests" ${BUILD_TESTING})
option(LMCAS_BUILD_LMMC_TESTS "Build nested LMMC tests" ${BUILD_TESTING})
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
    message(FATAL_ERROR "MSVC ABI compilers (including clang-cl) are unsupported; use Windows MinGW or GNU-driver Clang")
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
# Propagate the active CAS selection when reusing the build cache.
set(LMMC_LMMP_ASM "${LMCAS_LMMP_ASM}")
if(LMCAS_BUILD_LMMC_TESTS)
    add_subdirectory(LMMC)
else()
    add_subdirectory(LMMC EXCLUDE_FROM_ALL)
endif()
if(NOT LMMC_BUILD_SHARED)
    set_property(TARGET lmmc PROPERTY POSITION_INDEPENDENT_CODE ON)
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

include(cmake/ComponentSources.cmake)

set(LMCAS_COMPONENT_TARGETS
    lmcas_exact_numeric
    lmcas_core
    lmcas_assumptions
    lmcas_matrix_kernel
    lmcas_algebra
    lmcas_calculus
    lmcas_linear_algebra
    lmcas_analysis
    lmcas_facade
)

add_library(lmcas_exact_numeric OBJECT ${LMCAS_EXACT_NUMERIC_SOURCES})
add_library(lmcas_core OBJECT ${LMCAS_CORE_SOURCES})
add_library(lmcas_assumptions OBJECT ${LMCAS_ASSUMPTION_SOURCES})
add_library(lmcas_matrix_kernel OBJECT ${LMCAS_MATRIX_KERNEL_SOURCES})
add_library(lmcas_algebra OBJECT ${LMCAS_ALGEBRA_SOURCES})
add_library(lmcas_calculus OBJECT ${LMCAS_CALCULUS_SOURCES})
add_library(lmcas_linear_algebra OBJECT ${LMCAS_LINEAR_ALGEBRA_SOURCES})
add_library(lmcas_analysis OBJECT ${LMCAS_ANALYSIS_SOURCES})
add_library(lmcas_facade OBJECT ${LMCAS_FACADE_SOURCES})


include(cmake/ComponentHeaders.cmake)

target_link_libraries(lmcas_core PRIVATE lmcas_exact_numeric)
target_link_libraries(lmcas_assumptions PRIVATE lmcas_exact_numeric lmcas_core)
target_link_libraries(lmcas_matrix_kernel PRIVATE
    lmcas_exact_numeric lmcas_core lmcas_assumptions)
target_link_libraries(lmcas_algebra PRIVATE
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_matrix_kernel)
target_link_libraries(lmcas_linear_algebra PRIVATE
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_matrix_kernel lmcas_algebra)
target_link_libraries(lmcas_calculus PRIVATE
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_algebra)
target_link_libraries(lmcas_analysis PRIVATE
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_algebra
    lmcas_linear_algebra lmcas_calculus)
target_link_libraries(lmcas_facade PRIVATE
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_algebra
    lmcas_linear_algebra lmcas_calculus lmcas_analysis)

set(_lmcas_component_sources
    ${LMCAS_EXACT_NUMERIC_SOURCES}
    ${LMCAS_CORE_SOURCES}
    ${LMCAS_ASSUMPTION_SOURCES}
    ${LMCAS_MATRIX_KERNEL_SOURCES}
    ${LMCAS_ALGEBRA_SOURCES}
    ${LMCAS_CALCULUS_SOURCES}
    ${LMCAS_LINEAR_ALGEBRA_SOURCES}
    ${LMCAS_ANALYSIS_SOURCES}
    ${LMCAS_FACADE_SOURCES}
)

add_library(lmcas SHARED)
target_compile_features(lmcas PUBLIC cxx_std_17)
foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
    target_compile_features(${component} PUBLIC cxx_std_17)
    target_sources(lmcas PRIVATE $<TARGET_OBJECTS:${component}>)
endforeach()


set_target_properties(lmcas PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
)
if(APPLE)
    set_target_properties(lmcas PROPERTIES
        MACOSX_RPATH ON
        INSTALL_NAME_DIR "@rpath"
    )
endif()

foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
    target_compile_definitions(${component} PRIVATE LMCAS_EXPORTS)
    target_include_directories(${component} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
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

if(LMCAS_BUILD_TESTS)
    # Production components have hidden visibility; tests need separately
    # compiled, exported symbols without changing the installed lmcas ABI.
    add_library(lmcas_test_internal SHARED ${_lmcas_component_sources})
    target_compile_features(lmcas_test_internal PUBLIC cxx_std_17)
    target_compile_definitions(lmcas_test_internal PUBLIC LMCAS_STATIC_DEFINE)
    target_include_directories(lmcas_test_internal PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_include_directories(lmcas_test_internal PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/LMMC/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LMMP/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/LMMC/LMMP/include/lmmp")
    target_link_libraries(lmcas_test_internal PUBLIC lmmp lmmc)
    set_target_properties(lmcas_test_internal PROPERTIES
        CXX_EXTENSIONS OFF
        CXX_VISIBILITY_PRESET default
        VISIBILITY_INLINES_HIDDEN NO
        WINDOWS_EXPORT_ALL_SYMBOLS ON)
    if(WIN32)
        # These tests observe executable-owned allocation and thread teardown.
        add_library(lmcas_test_local STATIC)
        set_target_properties(lmcas_test_local PROPERTIES LINKER_LANGUAGE CXX)
        if(MINGW)
            set_property(TARGET lmcas_test_local PROPERTY STATIC_LIBRARY_OPTIONS --thin)
        endif()
        foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
            target_sources(lmcas_test_local PRIVATE $<TARGET_OBJECTS:${component}>)
        endforeach()
        target_compile_definitions(lmcas_test_local INTERFACE LMCAS_STATIC_DEFINE)
        target_include_directories(lmcas_test_local INTERFACE
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/src")
        target_link_libraries(lmcas_test_local PUBLIC lmmp lmmc)
    endif()
endif()

target_compile_definitions(lmcas PRIVATE LMCAS_EXPORTS)

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
           NOT CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
            message(FATAL_ERROR "The configured sanitizer mode requires GCC, Clang, or AppleClang for both C and C++ on a non-Windows platform")
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
           NOT CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
            message(FATAL_ERROR "Coverage requires GCC, Clang, or AppleClang for both C and C++ on a non-Windows platform")
        endif()
        target_compile_options(${target_name} PRIVATE
            --coverage -fprofile-update=atomic -O0 -g)
        target_link_options(${target_name} PRIVATE --coverage)
    endif()
endfunction()

foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
    lmcas_apply_instrumentation(${component})
endforeach()
lmcas_apply_instrumentation(lmcas)
if(LMCAS_BUILD_TESTS)
    target_compile_options(lmcas_test_internal PRIVATE ${_lmcas_warning_options})
    lmcas_apply_instrumentation(lmcas_test_internal)
endif()

# ============================================================
# Benchmarks
# ============================================================
if(LMCAS_BUILD_BENCHMARKS)
    add_executable(benchmark_all benchmarks/benchmark_all.cpp)
    target_link_libraries(benchmark_all PRIVATE lmcas)
    target_compile_definitions(benchmark_all PRIVATE _SYMBOLIC_BENCHMARK)
    target_compile_options(benchmark_all PRIVATE ${_lmcas_warning_options})
    lmcas_apply_instrumentation(benchmark_all)
    set_target_properties(benchmark_all PROPERTIES CXX_EXTENSIONS OFF)
    if(BUILD_TESTING)
        add_test(NAME benchmark_smoke COMMAND benchmark_all)
        set_tests_properties(benchmark_smoke PROPERTIES LABELS "benchmark" TIMEOUT 300)
    endif()
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
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
list(TRANSFORM _lmcas_public_headers
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
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

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/LMCASConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LMCAS
)

# ============================================================
# Tests (CTest)
# ============================================================
if(LMCAS_BUILD_TESTS OR LMMC_BUILD_TESTS)
    enable_testing()
endif()
if(LMCAS_BUILD_TESTS)
    add_subdirectory(tests)
endif()
