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
option(LMCAS_ENABLE_QUALITY_GATES "Enable reproducible first-party quality gates" OFF)
if(LMCAS_ENABLE_QUALITY_GATES)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    include("${CMAKE_CURRENT_SOURCE_DIR}/LMMC/cmake/QualityGates.cmake")
    lmmc_quality_require_tools()
endif()
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
    message(FATAL_ERROR "LMCAS supports GCC, Clang, and MinGW; MSVC is not supported")
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
if(LMCAS_ENABLE_QUALITY_GATES)
    set(LMMC_ENABLE_QUALITY_GATES ON)
elseif(NOT DEFINED LMMC_ENABLE_QUALITY_GATES)
    set(LMMC_ENABLE_QUALITY_GATES OFF)
endif()
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

# Component dependencies are intentionally acyclic and point from higher
# layers to lower layers. Keep the final shared library as the only packaged
# C++ target; object components are build boundaries, not public libraries.
set_property(TARGET lmcas_exact_numeric PROPERTY LMCAS_LAYER_RANK 0)
set_property(TARGET lmcas_core PROPERTY LMCAS_LAYER_RANK 1)
set_property(TARGET lmcas_assumptions PROPERTY LMCAS_LAYER_RANK 2)
set_property(TARGET lmcas_matrix_kernel PROPERTY LMCAS_LAYER_RANK 3)
set_property(TARGET lmcas_algebra PROPERTY LMCAS_LAYER_RANK 4)
set_property(TARGET lmcas_linear_algebra PROPERTY LMCAS_LAYER_RANK 5)
set_property(TARGET lmcas_calculus PROPERTY LMCAS_LAYER_RANK 5)
set_property(TARGET lmcas_analysis PROPERTY LMCAS_LAYER_RANK 6)
set_property(TARGET lmcas_facade PROPERTY LMCAS_LAYER_RANK 7)

include(cmake/ComponentHeaders.cmake)

function(lmcas_component_depends component)
    get_target_property(_component_rank ${component} LMCAS_LAYER_RANK)
    set_property(TARGET ${component} APPEND PROPERTY LMCAS_COMPONENT_DEPENDENCIES ${ARGN})
    foreach(dependency IN LISTS ARGN)
        get_target_property(_dependency_rank ${dependency} LMCAS_LAYER_RANK)
        if(_dependency_rank GREATER_EQUAL _component_rank)
            message(FATAL_ERROR
                "${component} cannot depend on same/higher layer ${dependency}")
        endif()
        target_link_libraries(${component} PRIVATE ${dependency})
    endforeach()
endfunction()

lmcas_component_depends(lmcas_core lmcas_exact_numeric)
lmcas_component_depends(lmcas_assumptions lmcas_exact_numeric lmcas_core)
lmcas_component_depends(lmcas_matrix_kernel
    lmcas_exact_numeric lmcas_core lmcas_assumptions)
lmcas_component_depends(lmcas_algebra
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_matrix_kernel)
lmcas_component_depends(lmcas_linear_algebra
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_matrix_kernel lmcas_algebra)
lmcas_component_depends(lmcas_calculus
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_algebra)
lmcas_component_depends(lmcas_analysis
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_algebra
    lmcas_linear_algebra lmcas_calculus)
lmcas_component_depends(lmcas_facade
    lmcas_exact_numeric lmcas_core lmcas_assumptions lmcas_algebra
    lmcas_linear_algebra lmcas_calculus lmcas_analysis)

# Every production translation unit belongs to exactly one component.
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
list(LENGTH _lmcas_component_sources _lmcas_component_source_count)
list(REMOVE_DUPLICATES _lmcas_component_sources)
list(LENGTH _lmcas_component_sources _lmcas_unique_source_count)
if(NOT _lmcas_component_source_count EQUAL _lmcas_unique_source_count)
    message(FATAL_ERROR
        "A production source is assigned to more than one LMCAS component")
endif()
file(GLOB_RECURSE _lmcas_discovered_sources CONFIGURE_DEPENDS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(SORT _lmcas_component_sources)
list(SORT _lmcas_discovered_sources)
if(NOT _lmcas_component_sources STREQUAL _lmcas_discovered_sources)
    message(FATAL_ERROR
        "Component source lists must cover every recursive src/*.c and src/*.cpp file exactly once")
endif()

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
    add_library(lmcas_test_internal STATIC)
    foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
        target_sources(lmcas_test_internal PRIVATE $<TARGET_OBJECTS:${component}>)
    endforeach()
    target_compile_features(lmcas_test_internal PUBLIC cxx_std_17)
    target_compile_definitions(lmcas_test_internal INTERFACE LMCAS_STATIC_DEFINE)
    target_include_directories(lmcas_test_internal PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_link_libraries(lmcas_test_internal PUBLIC lmmp lmmc)
    set_target_properties(lmcas_test_internal PROPERTIES
        POSITION_INDEPENDENT_CODE ON CXX_EXTENSIONS OFF)
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
    if(LMCAS_ENABLE_QUALITY_GATES)
        target_compile_options(${component} PRIVATE -g)
    endif()
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

set(_lmcas_owned_public_headers)
foreach(component IN LISTS LMCAS_COMPONENT_TARGETS)
    get_target_property(_component_files ${component} SOURCES)
    foreach(_component_file IN LISTS _component_files)
        if(IS_ABSOLUTE "${_component_file}")
            file(RELATIVE_PATH _component_file
                "${CMAKE_CURRENT_SOURCE_DIR}" "${_component_file}")
        endif()
        cmake_path(NORMAL_PATH _component_file)
        if(_component_file MATCHES "^include/[^/]+\\.hpp$")
            list(APPEND _lmcas_owned_public_headers "${_component_file}")
        endif()
    endforeach()
endforeach()
list(LENGTH _lmcas_owned_public_headers _lmcas_owned_public_header_count)
list(REMOVE_DUPLICATES _lmcas_owned_public_headers)
list(LENGTH _lmcas_owned_public_headers _lmcas_unique_public_header_count)
if(NOT _lmcas_owned_public_header_count EQUAL _lmcas_unique_public_header_count)
    message(FATAL_ERROR
        "A public header is assigned to more than one LMCAS component")
endif()

file(GLOB _lmcas_public_headers CONFIGURE_DEPENDS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
list(SORT _lmcas_owned_public_headers)
list(SORT _lmcas_public_headers)
if(NOT _lmcas_owned_public_headers STREQUAL _lmcas_public_headers)
    message(FATAL_ERROR
        "Component header lists must cover every public include/*.hpp file exactly once")
endif()
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

if(LMCAS_ENABLE_QUALITY_GATES)
    lmmc_add_quality_gates(lmcas "${CMAKE_CURRENT_SOURCE_DIR}")
    include(cmake/Architecture.cmake)
    lmcas_add_architecture_gate()
endif()
