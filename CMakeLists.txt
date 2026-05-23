cmake_minimum_required(VERSION 3.10)
project(LMCAS)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    message(WARNING "MSVC is deprecated, please use GCC/G++ instead.")
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
    src/solve_visitor.cpp
    src/matcher.cpp
    src/integration.cpp
    src/solver.cpp
    src/solve_polynomial.cpp
    src/solve_transcendental.cpp
    src/newton_raphson.cpp
    src/root_of_utils.cpp
    src/interval.cpp
    src/inequality_solver.cpp
    src/parametric_solver.cpp
    src/symbolic_ode.cpp
    src/symbolic_geometry.cpp
    src/symbolic_vector_geometry.cpp
    src/symbolic_complex.cpp
    src/symbolic_matrix.cpp
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

add_executable(benchmark_all benchmarks/benchmark_all.cpp)
target_link_libraries(benchmark_all PRIVATE lmcas LammpCore)
target_compile_definitions(benchmark_all PRIVATE -D_SYMBOLIC_BENCHMARK)

add_executable(test_trig_rewrite tests/test_trig_rewrite.cpp)
target_link_libraries(test_trig_rewrite PRIVATE lmcas LammpCore)

add_executable(test_matcher_comm tests/test_matcher_comm.cpp)
target_link_libraries(test_matcher_comm PRIVATE lmcas LammpCore)

add_executable(test_review_fixes tests/test_review_fixes.cpp)
target_link_libraries(test_review_fixes PRIVATE lmcas LammpCore)

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
add_test(NAME TestReviewFixes COMMAND test_review_fixes)

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

add_executable(test_solve_cubic tests/test_solve_cubic.cpp)
target_link_libraries(test_solve_cubic PRIVATE lmcas LammpCore)
add_test(NAME TestSolveCubic COMMAND test_solve_cubic)

add_executable(test_solve_quartic tests/test_solve_quartic.cpp)
target_link_libraries(test_solve_quartic PRIVATE lmcas LammpCore)
add_test(NAME TestSolveQuartic COMMAND test_solve_quartic)

add_executable(test_find_rational_roots tests/test_find_rational_roots.cpp)
target_link_libraries(test_find_rational_roots PRIVATE lmcas LammpCore)
add_test(NAME TestFindRationalRoots COMMAND test_find_rational_roots)

add_executable(test_solve_factoring tests/test_solve_factoring.cpp)
target_link_libraries(test_solve_factoring PRIVATE lmcas LammpCore)
add_test(NAME TestSolveFactoring COMMAND test_solve_factoring)

add_executable(test_solve_transcendental tests/test_solve_transcendental.cpp)
target_link_libraries(test_solve_transcendental PRIVATE lmcas LammpCore)
add_test(NAME TestSolveTranscendental COMMAND test_solve_transcendental)

add_executable(test_newton_raphson tests/test_newton_raphson.cpp)
target_link_libraries(test_newton_raphson PRIVATE lmcas LammpCore)
add_test(NAME TestNewtonRaphson COMMAND test_newton_raphson)

add_executable(test_rootof tests/test_rootof.cpp)
target_link_libraries(test_rootof PRIVATE lmcas LammpCore)
add_test(NAME TestRootOf COMMAND test_rootof)

add_executable(test_sturm_isolation tests/test_sturm_isolation.cpp)
target_link_libraries(test_sturm_isolation PRIVATE lmcas LammpCore)
add_test(NAME TestSturmIsolation COMMAND test_sturm_isolation)

add_executable(test_dispatcher_routing tests/test_dispatcher_routing.cpp)
target_link_libraries(test_dispatcher_routing PRIVATE lmcas LammpCore)
add_test(NAME TestDispatcherRouting COMMAND test_dispatcher_routing)

add_executable(test_solve_extended tests/test_solve_extended.cpp)
target_link_libraries(test_solve_extended PRIVATE lmcas LammpCore)
add_test(NAME TestSolveExtended COMMAND test_solve_extended)

add_executable(test_interval tests/test_interval.cpp)
target_link_libraries(test_interval PRIVATE lmcas LammpCore)
add_test(NAME TestInterval COMMAND test_interval)

add_executable(test_interval_serialize tests/test_interval_serialize.cpp)
target_link_libraries(test_interval_serialize PRIVATE lmcas LammpCore)
add_test(NAME TestIntervalSerialize COMMAND test_interval_serialize)

add_executable(test_inequality tests/test_inequality.cpp)
target_link_libraries(test_inequality PRIVATE lmcas LammpCore)
add_test(NAME TestInequality COMMAND test_inequality)

add_executable(test_parametric_inequality tests/test_parametric_inequality.cpp)
target_link_libraries(test_parametric_inequality PRIVATE lmcas LammpCore)
add_test(NAME TestParametricInequality COMMAND test_parametric_inequality)

add_executable(test_parametric_system tests/test_parametric_system.cpp)
target_link_libraries(test_parametric_system PRIVATE lmcas LammpCore)
add_test(NAME TestParametricSystem COMMAND test_parametric_system)

add_executable(test_symbolic_ode tests/test_symbolic_ode.cpp)
target_link_libraries(test_symbolic_ode PRIVATE lmcas LammpCore)
add_test(NAME TestSymbolicODE COMMAND test_symbolic_ode)

add_executable(test_symbolic_geometry tests/test_symbolic_geometry.cpp)
target_link_libraries(test_symbolic_geometry PRIVATE lmcas LammpCore)
add_test(NAME TestSymbolicGeometry COMMAND test_symbolic_geometry)

add_executable(test_symbolic_vector_geometry tests/test_symbolic_vector_geometry.cpp)
target_link_libraries(test_symbolic_vector_geometry PRIVATE lmcas LammpCore)
add_test(NAME TestSymbolicVectorGeometry COMMAND test_symbolic_vector_geometry)

add_executable(test_symbolic_complex tests/test_symbolic_complex.cpp)
target_link_libraries(test_symbolic_complex PRIVATE lmcas LammpCore)
add_test(NAME TestSymbolicComplex COMMAND test_symbolic_complex)

add_executable(test_poly_utils tests/test_poly_utils.cpp)
target_link_libraries(test_poly_utils PRIVATE lmcas LammpCore)
add_test(NAME TestPolyUtils COMMAND test_poly_utils)

add_executable(test_symbolic_implicit_diff tests/test_symbolic_implicit_diff.cpp)
target_link_libraries(test_symbolic_implicit_diff PRIVATE lmcas LammpCore)
add_test(NAME TestSymbolicImplicitDiff COMMAND test_symbolic_implicit_diff)

add_executable(test_symbolic_matrix_advanced tests/test_symbolic_matrix_advanced.cpp)
target_link_libraries(test_symbolic_matrix_advanced PRIVATE lmcas LammpCore)
add_test(NAME TestSymbolicMatrixAdvanced COMMAND test_symbolic_matrix_advanced)
