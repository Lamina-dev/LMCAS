include(FetchContent)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

find_package(GTest 1.15.2 CONFIG QUIET)
if(NOT TARGET GTest::gtest_main)
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
        GIT_SHALLOW TRUE
        SOURCE_SUBDIR development-only)
    FetchContent_MakeAvailable(googletest)
    set(INSTALL_GTEST OFF)
    set(BUILD_GMOCK OFF)
    set(gtest_force_shared_crt ON)
    add_subdirectory("${googletest_SOURCE_DIR}" "${googletest_BINARY_DIR}" EXCLUDE_FROM_ALL SYSTEM)
endif()

find_package(rapidcheck CONFIG QUIET)
if(NOT TARGET rapidcheck)
    FetchContent_Declare(rapidcheck
        GIT_REPOSITORY https://github.com/emil-e/rapidcheck.git
        GIT_TAG 6e8dadfdafa3a74eabb52ead87f8787f72eccd0b
        GIT_SUBMODULES ""
        SOURCE_SUBDIR development-only)
    FetchContent_MakeAvailable(rapidcheck)
    set(RC_ENABLE_TESTS OFF)
    set(RC_ENABLE_EXAMPLES OFF)
    set(RC_ENABLE_GTEST OFF)
    set(RC_ENABLE_GMOCK OFF)
    set(RC_ENABLE_CATCH OFF)
    set(RC_ENABLE_DOCTEST OFF)
    set(RC_ENABLE_BOOST OFF)
    set(RC_ENABLE_BOOST_TEST OFF)
    set(RC_INSTALL_ALL_EXTRAS OFF)
    add_subdirectory("${rapidcheck_SOURCE_DIR}" "${rapidcheck_BINARY_DIR}" EXCLUDE_FROM_ALL SYSTEM)
endif()
include(GoogleTest)

function(lmcas_register_test_source source)
    cmake_path(ABSOLUTE_PATH source BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" NORMALIZE)
    get_property(_registered GLOBAL PROPERTY LMCAS_REGISTERED_TEST_SOURCES)
    if(source IN_LIST _registered)
        message(FATAL_ERROR "Test source registered more than once: ${source}")
    endif()
    set_property(GLOBAL APPEND PROPERTY LMCAS_REGISTERED_TEST_SOURCES "${source}")
endfunction()

function(lmcas_add_test TARGET_NAME SOURCE_FILE)
    lmcas_register_test_source("${SOURCE_FILE}")
    cmake_parse_arguments(_test "PUBLIC" "LABEL" "" ${ARGN})
    add_executable(${TARGET_NAME} ${SOURCE_FILE})
    set_target_properties(${TARGET_NAME} PROPERTIES EXCLUDE_FROM_ALL FALSE)
    target_link_libraries(${TARGET_NAME} PRIVATE GTest::gtest_main rapidcheck)
    if(_test_PUBLIC)
        target_link_libraries(${TARGET_NAME} PRIVATE lmcas)
        set_property(TARGET ${TARGET_NAME} PROPERTY LMCAS_TEST_LINKAGE public)
    else()
        target_link_libraries(${TARGET_NAME} PRIVATE lmcas_test_internal)
        set_property(TARGET ${TARGET_NAME} PROPERTY LMCAS_TEST_LINKAGE internal)
    endif()
    lmcas_apply_instrumentation(${TARGET_NAME})
    target_compile_options(${TARGET_NAME} PRIVATE ${_lmcas_warning_options})
    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    if(WIN32)
        get_target_property(_rapidcheck_type rapidcheck TYPE)
        if(_rapidcheck_type STREQUAL "SHARED_LIBRARY")
            add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "$<TARGET_FILE:rapidcheck>" "$<TARGET_FILE_DIR:${TARGET_NAME}>"
                VERBATIM)
        endif()
    endif()
    set(_labels gtest)
    if(_test_LABEL)
        list(APPEND _labels "${_test_LABEL}")
    endif()
    gtest_discover_tests(${TARGET_NAME}
        TEST_PREFIX "${TARGET_NAME}."
        PROPERTIES LABELS "${_labels}"
        DISCOVERY_MODE PRE_TEST
        DISCOVERY_TIMEOUT 30)
endfunction()

function(lmcas_add_test_group GROUP_NAME)
    cmake_parse_arguments(_group "PUBLIC" "" "" ${ARGN})
    set(_linkage "")
    if(_group_PUBLIC)
        set(_linkage PUBLIC)
    endif()
    foreach(TEST_SOURCE IN LISTS _group_UNPARSED_ARGUMENTS)
        get_filename_component(TEST_STEM "${TEST_SOURCE}" NAME_WE)
        lmcas_add_test("${TEST_STEM}" "${TEST_SOURCE}" ${_linkage} LABEL "${GROUP_NAME}")
    endforeach()
endfunction()

function(lmcas_add_public_header_checks TARGET_NAME)
    set(HEADER_CHECK_SOURCES "")
    foreach(HEADER_NAME IN LISTS ARGN)
        string(MAKE_C_IDENTIFIER "${HEADER_NAME}" HEADER_ID)
        set(HEADER_SOURCE
            "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}_${HEADER_ID}.cpp"
        )
        set(HEADER_CONTENT "#include <${HEADER_NAME}>\n")
        string(APPEND HEADER_CONTENT
            "#ifdef LMCAS_INTERNAL_AST_INCLUDED\n"
            "#error Public headers must not include private AST implementations\n"
            "#endif\n")
        file(GENERATE OUTPUT "${HEADER_SOURCE}" CONTENT "${HEADER_CONTENT}")
        list(APPEND HEADER_CHECK_SOURCES "${HEADER_SOURCE}")
    endforeach()
    add_library(${TARGET_NAME} OBJECT ${HEADER_CHECK_SOURCES})
    set_target_properties(${TARGET_NAME} PROPERTIES EXCLUDE_FROM_ALL FALSE)
    target_link_libraries(${TARGET_NAME} PRIVATE lmcas)
    target_compile_options(${TARGET_NAME} PRIVATE ${_lmcas_warning_options})
    lmcas_apply_instrumentation(${TARGET_NAME})
endfunction()

lmcas_register_test_source(test_lmmc_lifecycle_windows.cpp)
add_executable(test_lmmc_lifecycle_windows test_lmmc_lifecycle_windows.cpp)
target_link_libraries(test_lmmc_lifecycle_windows PRIVATE lmcas_test_internal)
set_property(TARGET test_lmmc_lifecycle_windows PROPERTY LMCAS_TEST_LINKAGE internal)
if(WIN32 AND MINGW)
    target_link_options(test_lmmc_lifecycle_windows PRIVATE -municode)
endif()

lmcas_register_test_source(lmmc_lifecycle_windows_module.cpp)
add_library(lmcas_lifecycle_windows_module MODULE lmmc_lifecycle_windows_module.cpp)
target_link_libraries(lmcas_lifecycle_windows_module PRIVATE lmcas)
set_property(TARGET lmcas_lifecycle_windows_module PROPERTY LMCAS_TEST_LINKAGE public)
foreach(_lifecycle_target test_lmmc_lifecycle_windows lmcas_lifecycle_windows_module)
    lmcas_apply_instrumentation(${_lifecycle_target})
    target_compile_options(${_lifecycle_target} PRIVATE ${_lmcas_warning_options})
    set_target_properties(${_lifecycle_target} PROPERTIES
        EXCLUDE_FROM_ALL FALSE
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
endforeach()
add_dependencies(test_lmmc_lifecycle_windows lmcas_lifecycle_windows_module)
add_test(NAME test_lmmc_lifecycle_windows
    COMMAND test_lmmc_lifecycle_windows
        "$<TARGET_FILE:lmcas_lifecycle_windows_module>" "$<TARGET_FILE:lmcas>")
set_tests_properties(test_lmmc_lifecycle_windows PROPERTIES
    LABELS "lifecycle;windows" TIMEOUT 120 SKIP_RETURN_CODE 77)


file(GLOB_RECURSE _lmcas_selfcheck_headers CONFIGURE_DEPENDS
    RELATIVE "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/include/*.hpp")
lmcas_add_public_header_checks(lmcas_public_header_checks
    ${_lmcas_selfcheck_headers})

lmcas_add_test(test_computation_context_gtest gtest/test_computation_context.cpp PUBLIC)
lmcas_add_test(test_result_contracts_gtest gtest/test_result_contracts.cpp PUBLIC)
lmcas_add_test(test_units_and_sets_gtest gtest/test_units_and_sets.cpp)


lmcas_add_test_group(
    test_core
    test_main.cpp
    test_suite.cpp
    test_normalization.cpp
    test_expansion_regressions.cpp
    test_ast_depth.cpp
    test_ast_invariants.cpp
    test_norm.cpp
    test_proof.cpp
    test_arithmetic.cpp
    test_bigint_core.cpp
    test_bigint_add_carry.cpp
    test_bigint_safety.cpp
    test_fraction.cpp
    test_rational_safety.cpp
    test_number_compare.cpp
    test_expr_serialization.cpp
    test_value_compare.cpp
    test_numeric_evaluation.cpp
    test_computation_context.cpp
    test_domain_preserving_simplify.cpp
    test_result_contracts.cpp
    test_complex_arithmetic.cpp
    test_symbolic_complex.cpp
)

lmcas_add_test_group(
    test_public_api PUBLIC
    test_symbolic_public_header.cpp
    test_expr_parsing.cpp
    test_expr_values.cpp
    test_expr_arithmetic.cpp
    test_expr_math.cpp
    test_expr_evaluation.cpp
    test_expr_matching.cpp
    test_expr_complex.cpp
    test_expr_complex_evaluation.cpp
    test_expr_sets.cpp
    test_expr_sets_domains.cpp
    test_expr_sets_solutions.cpp
    test_expr_equivalence.cpp
)

lmcas_add_test_group(
    test_calculus_group
    test_calculus.cpp
    test_series.cpp
    test_integration_strategies.cpp
    test_cyclic_integration.cpp
    test_definite_integration.cpp
    test_integration_fallback.cpp
    test_integration_table_priority.cpp
    test_integration_table_roundtrip.cpp
    test_linear_substitution_roundtrip.cpp
    test_rational_decomposition_roundtrip.cpp
    test_trig_tan_roundtrip.cpp
    test_trig_sec_roundtrip.cpp
    test_trig_sin_cos_roundtrip.cpp
    test_special_function_ast.cpp
    test_multiple_integral_separability.cpp
    test_singularity.cpp
    test_symbolic_implicit_diff.cpp
    test_trig.cpp
    test_trig_rewrite.cpp
    test_calculus_lim.cpp
    test_squeeze_theorem.cpp
    test_piecewise_limit.cpp
    test_calculus_utils.cpp
    test_calculus_utils_differentials.cpp
    test_calculus_utils_inverse.cpp
    test_calculus_utils_geometry.cpp
    test_curvature.cpp
    test_indeterminate_limits.cpp
    test_limits_at_infinity.cpp
    test_taylor_fallback.cpp
    test_exponential_limit.cpp
    test_prop_calculus_utils.cpp
    test_trig_substitution.cpp
)

lmcas_add_test_group(
    test_solving_group
    test_solve.cpp
    test_solve_visitor_poly.cpp
    test_solver_linear.cpp
    test_solver_poly.cpp
    test_solve_cubic.cpp
    test_solve_quartic.cpp
    test_find_rational_roots.cpp
    test_solve_factoring.cpp
    test_solve_transcendental.cpp
    test_solve_extended.cpp
    test_newton_raphson.cpp
    test_rootof.cpp
    test_sturm_isolation.cpp
    test_dispatcher_routing.cpp
    test_interval.cpp
    test_interval_checked.cpp
    test_interval_serialize.cpp
    test_inequality_examples.cpp
    test_inequality_soundness.cpp
    test_inequality_completeness.cpp
    test_inequality_endpoints.cpp
    test_inequality_parametric_consistency.cpp
    test_parametric_inequality.cpp
    test_parametric_system.cpp
    test_solver.cpp
    test_interval_properties.cpp
    test_interval_expressions.cpp
    test_newton_raphson_isolation.cpp
    test_newton_raphson_residual.cpp
    test_newton_raphson_deflation.cpp
    test_rootof_isolation.cpp
    test_solve_quartic_degenerate.cpp
    test_solve_quartic_properties.cpp
    test_parametric_inequality_degenerate.cpp
)

lmcas_add_test_group(
    test_algebra_group
    test_poly_bridge.cpp
    test_polynomial_gcd.cpp
    test_poly_utils.cpp
    test_factor.cpp
    test_det.cpp
    test_matrix_ops.cpp
    test_symbolic_matrix_advanced.cpp
    test_symbolic_matrix_eigen.cpp
    test_matcher.cpp
    test_matcher_comm.cpp
    test_symbolic_ode.cpp
    test_symbolic_geometry.cpp
    test_symbolic_vector_geometry_vector.cpp
    test_symbolic_vector_geometry_line_plane.cpp
    test_symbolic_vector_geometry_line_plane_construction.cpp
    test_symbolic_vector_geometry_quadric.cpp
    test_symbolic_vector_geometry_surface.cpp
    test_polynomial_resultant.cpp
    test_core_regressions.cpp
    test_equivalence_complex_allocation.cpp
    test_irrational_regressions.cpp
    test_groebner_advanced.cpp
    test_modular_arithmetic.cpp
    test_modular_fglm.cpp
    test_monomial_order.cpp
    test_print.cpp
    test_transcendental_factor_substitution.cpp
    test_transcendental_factor_substitution_polynomial.cpp
    test_transcendental_factor_finite_field.cpp
    test_transcendental_factor_zassenhaus.cpp
    test_transcendental_factor_back_substitution.cpp
    test_transcendental_factor_pipeline.cpp
    test_hensel_lift.cpp
    test_mixed_transcendental_classification.cpp
    test_mixed_transcendental_search_interval.cpp
    test_mixed_transcendental_isolation.cpp
    test_mixed_transcendental_refinement.cpp
    test_mixed_transcendental_deduplication.cpp
    test_mixed_transcendental_pipeline.cpp
    test_mixed_transcendental_pipeline_contracts.cpp
    test_cancel.cpp
    test_multipoly_construction.cpp
    test_multipoly_construction_univariate.cpp
    test_multipoly_evaluation.cpp
    test_multipoly_content.cpp
    test_multipoly_division_degree.cpp
    test_multipoly_ring.cpp
    test_multivariate_content.cpp
    test_multivariate_sqfree.cpp
    test_multivariate_diophantine.cpp
    test_hensel_lift_smoke.cpp
    test_multivariate_factor.cpp
    test_normalization_ext.cpp
    test_hensel_lift_two_factor.cpp
    test_hensel_lift_multifactor.cpp
    test_hensel_lift_symmetric.cpp
    test_multivariate_content_gcd.cpp
    test_multivariate_content_primitive.cpp
    test_multivariate_diophantine_solutions.cpp
    test_multivariate_diophantine_properties.cpp
    test_multivariate_factor_division.cpp
    test_multivariate_factor_properties.cpp
    test_multivariate_factor_complete.cpp
    test_multivariate_sqfree_multiplicities.cpp
)

# Assumption tests retain process isolation: many exercise global caches,
# randomized generators, and deliberate failure paths.
lmcas_add_test(test_assumption_relation_store test_assumption_relation_store.cpp)
lmcas_add_test(test_property_store_sign test_property_store_sign.cpp)
lmcas_add_test(test_property_store_parity_bounded test_property_store_parity_bounded.cpp)
lmcas_add_test(test_assumption_property_store test_assumption_property_store.cpp)
lmcas_add_test(test_assumption_inference_add test_assumption_inference_add.cpp)
lmcas_add_test(test_assumption_inference_mul test_assumption_inference_mul.cpp)
lmcas_add_test(test_assumption_inference_mul_sign_sets test_assumption_inference_mul_sign_sets.cpp)
lmcas_add_test(test_assumption_inference_pow test_assumption_inference_pow.cpp)
lmcas_add_test(test_assumption_inference_func test_assumption_inference_func.cpp)
lmcas_add_test(test_assumption_scope test_assumption_scope.cpp)
lmcas_add_test(test_assumption_query test_assumption_query.cpp)
lmcas_add_test(test_assumption_api test_assumption_api.cpp)
lmcas_add_test(test_assumption_monotonicity test_assumption_monotonicity.cpp)
lmcas_add_test(test_assumption_interval test_assumption_interval.cpp)
lmcas_add_test(test_assumption_interval_certified test_assumption_interval_certified.cpp)
lmcas_add_test(test_assumption_simplify test_assumption_simplify.cpp)
lmcas_add_test(test_assumption_binder_scope test_assumption_binder_scope.cpp)
lmcas_add_test(test_assumption_solver test_assumption_solver.cpp)
lmcas_add_test(test_assumption_examples test_assumption_examples.cpp)
lmcas_add_test(test_assumption_enums test_assumption_enums.cpp)
lmcas_add_test(test_assumption_declarations test_assumption_declarations.cpp)
lmcas_add_test(test_assumption_property_rules test_assumption_property_rules.cpp)
lmcas_add_test(test_assumption_relation_props test_assumption_relation_props.cpp)
lmcas_add_test(test_assumption_relation_queries test_assumption_relation_ext.cpp)
lmcas_add_test(test_assumption_cycle_depth test_assumption_cycle_depth.cpp)
lmcas_add_test(test_assumption_properties test_assumption_properties_ext.cpp)
lmcas_add_test(test_assumption_division_sign test_assumption_division_sign.cpp)
lmcas_add_test(test_assumption_domain_inference test_assumption_domain_inference.cpp)
lmcas_add_test(test_inference_engine_queries test_inference_engine_ext.cpp)
lmcas_add_test(test_assumption_query_conditions test_assumption_query_ext.cpp)
lmcas_add_test(test_assumption_cache test_assumption_cache.cpp)
lmcas_add_test(test_assumption_conditional test_assumption_conditional.cpp)
lmcas_add_test(test_assumption_with_syntax test_assumption_with_syntax.cpp)
lmcas_add_test(test_assumption_serialization test_assumption_serialization.cpp)
lmcas_add_test(test_assumption_serialization_records test_assumption_serialization_records.cpp)
lmcas_add_test(test_assumption_context_operations test_assumption_context_ext.cpp)
lmcas_add_test(test_assumption_propagation test_assumption_propagation.cpp)
lmcas_add_test(test_assumption_diagnostics test_assumption_diagnostics.cpp)
lmcas_add_test(test_assumption_propagation_unit test_assumption_propagation_unit.cpp)
lmcas_add_test(test_assumption_system_integration test_assumption_system_integration.cpp)
lmcas_add_test(test_assumption_integration test_assumption_integration.cpp)
lmcas_add_test(test_assumption_inference_trig test_assumption_inference_trig.cpp)
lmcas_add_test(test_assumption_inference_mul_domains test_assumption_inference_mul_domains.cpp)
lmcas_add_test(test_assumption_monotonicity_intervals test_assumption_monotonicity_intervals.cpp)
lmcas_add_test(test_assumption_definiteness_rules test_assumption_definiteness_rules.cpp)
lmcas_add_test(test_assumption_domain_store test_assumption_domain_store.cpp)
lmcas_add_test(test_assumption_scope_shadowing test_assumption_scope_shadowing.cpp)
lmcas_add_test(test_assumption_interval_checked test_assumption_interval_checked.cpp)
lmcas_add_test(test_assumption_inference_checked test_assumption_inference_checked.cpp)
lmcas_add_test(test_assumption_query_extended_checked test_assumption_query_extended_checked.cpp)

lmcas_add_test_group(
    test_analysis_group
    test_vector_calculus_differential.cpp
    test_vector_calculus_jacobian_hessian.cpp
    test_vector_calculus_curve_integrals.cpp
    test_vector_calculus_surface_integrals.cpp
    test_vector_calculus_extrema.cpp
    test_prop_limits.cpp
    test_prop_vector_calculus.cpp
    test_integral_theorems.cpp
    test_integral_theorems_spatial.cpp
    test_integral_theorems_checked.cpp
    test_ode_engine_first_order.cpp
    test_ode_engine_higher_order.cpp
    test_ode_engine_classification.cpp
    test_ode_engine_advanced.cpp
    test_series_engine.cpp
    test_series_sequence_limits.cpp
    test_series_summation.cpp
    test_series_convergence.cpp
    test_series_expansions.cpp
    test_prop_series_engine.cpp
    test_z_transform.cpp
    test_transform_engine.cpp
    test_transform_contracts.cpp
    test_transform_exact_decay.cpp
    test_matrix_decomposition.cpp
    test_complex_analysis.cpp
    test_differential_geometry.cpp
    test_numerical_integration.cpp
)

file(GLOB_RECURSE _lmcas_discovered_test_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB _lmcas_package_consumer_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/package_consumer/*.cpp")
set(_lmcas_test_support_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/lmmc_lifecycle_windows_module.cpp"
    ${_lmcas_package_consumer_sources})
list(REMOVE_ITEM _lmcas_discovered_test_sources
    ${_lmcas_test_support_sources})

get_property(_lmcas_registered_test_sources
    GLOBAL PROPERTY LMCAS_REGISTERED_TEST_SOURCES)
list(REMOVE_ITEM _lmcas_registered_test_sources
    ${_lmcas_test_support_sources})
list(SORT _lmcas_discovered_test_sources)
list(SORT _lmcas_registered_test_sources)
if(NOT _lmcas_discovered_test_sources STREQUAL
   _lmcas_registered_test_sources)
    set(_lmcas_unregistered_test_sources
        ${_lmcas_discovered_test_sources})
    list(REMOVE_ITEM _lmcas_unregistered_test_sources
        ${_lmcas_registered_test_sources})
    set(_lmcas_missing_test_sources
        ${_lmcas_registered_test_sources})
    list(REMOVE_ITEM _lmcas_missing_test_sources
        ${_lmcas_discovered_test_sources})
    message(FATAL_ERROR
        "LMCAS test source lists must cover every tests/*.cpp file exactly once.\n"
        "Unregistered: ${_lmcas_unregistered_test_sources}\n"
        "Missing: ${_lmcas_missing_test_sources}")
endif()
