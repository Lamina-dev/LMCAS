option(LMCAS_FETCH_GTEST "Download GoogleTest when it is not installed" OFF)

include(FetchContent)

find_package(GTest CONFIG QUIET)
if(NOT GTest_FOUND AND LMCAS_FETCH_GTEST)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
        GIT_SHALLOW TRUE
        GIT_CONFIG http.version=HTTP/1.1
    )
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

if(GTest_FOUND OR TARGET GTest::gtest_main)
    set(LMCAS_HAVE_GTEST ON)
    include(GoogleTest)
else()
    set(LMCAS_HAVE_GTEST OFF)
    message(STATUS "GoogleTest not found; skipping GTest-based tests (enable LMCAS_FETCH_GTEST to download it)")
endif()

function(lmcas_add_gtest TARGET_NAME SOURCE_FILE)
    if(NOT LMCAS_HAVE_GTEST)
        return()
    endif()
    add_executable(${TARGET_NAME} ${SOURCE_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE lmcas GTest::gtest_main)
    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    )
    gtest_discover_tests(${TARGET_NAME} DISCOVERY_TIMEOUT 30)
endfunction()

function(lmcas_add_standalone_test TARGET_NAME SOURCE_FILE)
    if(NOT LMCAS_BUILD_STANDALONE_TESTS)
        return()
    endif()
    add_executable(${TARGET_NAME} ${SOURCE_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE lmcas)
    add_test(NAME ${TARGET_NAME} COMMAND ${TARGET_NAME})
    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    )
endfunction()

function(lmcas_add_standalone_test_group TARGET_NAME)
    if(NOT LMCAS_BUILD_STANDALONE_TESTS)
        return()
    endif()

    set(TEST_SOURCES ${ARGN})
    set(RUNNER_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}_runner.cpp")
    set(RUNNER_CONTENT "#include \"test_common.hpp\"\n#include <iostream>\n")
    set(RUNNER_CALLS "")

    foreach(TEST_SOURCE IN LISTS TEST_SOURCES)
        get_filename_component(TEST_STEM "${TEST_SOURCE}" NAME_WE)
        string(MAKE_C_IDENTIFIER "${TARGET_NAME}_${TEST_STEM}_main" TEST_ENTRYPOINT)
        set_property(
            SOURCE "${TEST_SOURCE}"
            APPEND PROPERTY COMPILE_DEFINITIONS "main=${TEST_ENTRYPOINT}"
        )
        string(APPEND RUNNER_CONTENT "int ${TEST_ENTRYPOINT}();\n")
        string(APPEND RUNNER_CALLS
            "    g_failures = 0;\n"
            "    g_passes = 0;\n"
            "    std::cout << \"\\n=== ${TEST_STEM} ===\\n\";\n"
            "    if (${TEST_ENTRYPOINT}() != 0) ++failed_programs;\n"
        )
    endforeach()

    string(APPEND RUNNER_CONTENT
        "\nint main() {\n"
        "    int failed_programs = 0;\n"
        "${RUNNER_CALLS}"
        "    return failed_programs == 0 ? 0 : 1;\n"
        "}\n"
    )
    file(GENERATE OUTPUT "${RUNNER_SOURCE}" CONTENT "${RUNNER_CONTENT}")

    add_executable(${TARGET_NAME} "${RUNNER_SOURCE}" ${TEST_SOURCES})
    target_include_directories(${TARGET_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(${TARGET_NAME} PRIVATE lmcas)
    add_test(NAME ${TARGET_NAME} COMMAND ${TARGET_NAME})
    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    )
endfunction()

function(lmcas_add_public_header_checks TARGET_NAME)
    set(HEADER_CHECK_SOURCES "")
    foreach(HEADER_NAME IN LISTS ARGN)
        string(MAKE_C_IDENTIFIER "${HEADER_NAME}" HEADER_ID)
        set(HEADER_SOURCE
            "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}_${HEADER_ID}.cpp"
        )
        set(HEADER_CONTENT
            "#define LMCAS_HEADER_UNDER_TEST \"${HEADER_NAME}\"\n"
        )
        string(APPEND HEADER_CONTENT
            "#define lmcas_public_header_isolation_anchor lmcas_public_header_isolation_anchor_${HEADER_ID}\n"
            "#include \"${CMAKE_CURRENT_SOURCE_DIR}/test_public_header_isolation.cpp\"\n"
        )
        file(GENERATE OUTPUT "${HEADER_SOURCE}" CONTENT "${HEADER_CONTENT}")
        list(APPEND HEADER_CHECK_SOURCES "${HEADER_SOURCE}")
    endforeach()
    add_library(${TARGET_NAME} OBJECT ${HEADER_CHECK_SOURCES})
    target_link_libraries(${TARGET_NAME} PRIVATE lmcas)
endfunction()

lmcas_add_public_header_checks(
    lmcas_public_header_checks
    complex_analysis.hpp
    calculus_utils.hpp
    differential_geometry.hpp
    numerical_integration.hpp
    matrix_decomposition.hpp
    series_engine.hpp
    symbolic_ode_engine.hpp
    transform_engine.hpp
    vector_calculus.hpp
    relation_store.hpp
    property_store.hpp
    interval.hpp
    query_interface.hpp
    assumption_context.hpp
    inference_engine.hpp
    poly_utils.hpp
    symbolic_implicit_diff.hpp
    lsr_expr.hpp
    quantity.hpp
    symbolic_set.hpp
)

lmcas_add_gtest(test_computation_context_gtest gtest/test_computation_context.cpp)
lmcas_add_gtest(test_result_contracts_gtest gtest/test_result_contracts.cpp)
lmcas_add_gtest(test_units_and_sets_gtest gtest/test_units_and_sets.cpp)

lmcas_add_standalone_test_group(
    test_core
    test_main.cpp
    test_suite.cpp
    test_normalization.cpp
    test_expand_debug.cpp
    test_ast_depth.cpp
    test_ast_invariants.cpp
    test_symbolic_public_header.cpp
    test_norm.cpp
    test_proof.cpp
    test_arithmetic.cpp
    test_bigint_debug.cpp
    test_bigint_add_carry.cpp
    test_bigint_safety.cpp
    test_fraction.cpp
    test_rational_safety.cpp
    test_number_compare.cpp
    test_value_compare.cpp
    test_numeric_evaluation_v2.cpp
    test_computation_context.cpp
    test_domain_preserving_simplify.cpp
    test_result_contracts.cpp
    test_lsr_expr.cpp
    test_complex_arithmetic.cpp
    test_symbolic_complex.cpp
)

lmcas_add_standalone_test_group(
    test_calculus_group
    test_calculus.cpp
    test_series.cpp
    test_new_integration.cpp
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
    test_integration_enhanced.cpp
    test_singularity.cpp
    test_symbolic_implicit_diff.cpp
    test_trig.cpp
    test_trig_rewrite.cpp
    test_calculus_lim.cpp
    test_squeeze_theorem.cpp
    test_piecewise_limit.cpp
    test_calculus_utils.cpp
    test_curvature.cpp
    test_indeterminate_limits.cpp
    test_limits_at_infinity.cpp
    test_taylor_fallback.cpp
    test_debug_limit.cpp
    test_prop_calculus_utils.cpp
    test_trig_substitution.cpp
)

lmcas_add_standalone_test_group(
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
    test_inequality.cpp
    test_parametric_inequality.cpp
    test_parametric_system.cpp
    test_solver.cpp
)

lmcas_add_standalone_test_group(
    test_algebra_group
    test_poly_bridge.cpp
    test_polynomial_gcd.cpp
    test_poly_utils.cpp
    test_factor.cpp
    test_det.cpp
    test_matrix_ops.cpp
    test_symbolic_matrix_advanced.cpp
    test_matcher.cpp
    test_matcher_comm.cpp
    test_symbolic_ode.cpp
    test_symbolic_geometry.cpp
    test_symbolic_vector_geometry.cpp
    test_new_features.cpp
    test_review_fixes.cpp
    test_groebner_advanced.cpp
    test_modular_arithmetic.cpp
    test_modular_fglm.cpp
    test_monomial_order.cpp
    test_print.cpp
    test_transcendental_factor.cpp
    test_hensel_lift.cpp
    test_mixed_transcendental.cpp
    test_cancel.cpp
    test_multipoly.cpp
    test_multivariate_content.cpp
    test_multivariate_sqfree.cpp
    test_multivariate_diophantine.cpp
    test_hensel_lift_smoke.cpp
    test_multivariate_factor.cpp
    test_normalization_ext.cpp
)

# Assumption tests retain process isolation: many exercise global caches,
# randomized generators, and deliberate failure paths.
lmcas_add_standalone_test(test_assumption_relation_store test_assumption_relation_store.cpp)
lmcas_add_standalone_test(test_property_store_sign test_property_store_sign.cpp)
lmcas_add_standalone_test(test_property_store_parity_bounded test_property_store_parity_bounded.cpp)
lmcas_add_standalone_test(test_assumption_property_store test_assumption_property_store.cpp)
lmcas_add_standalone_test(test_assumption_inference_add test_assumption_inference_add.cpp)
lmcas_add_standalone_test(test_assumption_inference_mul test_assumption_inference_mul.cpp)
lmcas_add_standalone_test(test_assumption_inference_pow test_assumption_inference_pow.cpp)
lmcas_add_standalone_test(test_assumption_inference_func test_assumption_inference_func.cpp)
lmcas_add_standalone_test(test_assumption_scope test_assumption_scope.cpp)
lmcas_add_standalone_test(test_assumption_query test_assumption_query.cpp)
lmcas_add_standalone_test(test_assumption_api test_assumption_api.cpp)
lmcas_add_standalone_test(test_assumption_monotonicity test_assumption_monotonicity.cpp)
lmcas_add_standalone_test(test_assumption_interval test_assumption_interval.cpp)
lmcas_add_standalone_test(test_assumption_simplify test_assumption_simplify.cpp)
lmcas_add_standalone_test(test_assumption_solver test_assumption_solver.cpp)
lmcas_add_standalone_test(test_assumption_examples test_assumption_examples.cpp)
lmcas_add_standalone_test(test_assumption_enums test_assumption_enums.cpp)
lmcas_add_standalone_test(test_assumption_declarations test_assumption_declarations.cpp)
lmcas_add_standalone_test(test_assumption_property_rules test_assumption_property_rules.cpp)
lmcas_add_standalone_test(test_assumption_relation_props test_assumption_relation_props.cpp)
lmcas_add_standalone_test(test_assumption_relation_queries test_assumption_relation_ext.cpp)
lmcas_add_standalone_test(test_assumption_cycle_depth test_assumption_cycle_depth.cpp)
lmcas_add_standalone_test(test_assumption_properties test_assumption_properties_ext.cpp)
lmcas_add_standalone_test(test_assumption_division_sign test_assumption_division_sign.cpp)
lmcas_add_standalone_test(test_assumption_domain_inference test_assumption_domain_inference.cpp)
lmcas_add_standalone_test(test_inference_engine_queries test_inference_engine_ext.cpp)
lmcas_add_standalone_test(test_assumption_query_conditions test_assumption_query_ext.cpp)
lmcas_add_standalone_test(test_assumption_cache test_assumption_cache.cpp)
lmcas_add_standalone_test(test_assumption_conditional test_assumption_conditional.cpp)
lmcas_add_standalone_test(test_assumption_with_syntax test_assumption_with_syntax.cpp)
lmcas_add_standalone_test(test_assumption_serialization test_assumption_serialization.cpp)
lmcas_add_standalone_test(test_assumption_context_operations test_assumption_context_ext.cpp)
lmcas_add_standalone_test(test_assumption_propagation test_assumption_propagation.cpp)
lmcas_add_standalone_test(test_assumption_diagnostics test_assumption_diagnostics.cpp)
lmcas_add_standalone_test(test_assumption_propagation_unit test_assumption_propagation_unit.cpp)
lmcas_add_standalone_test(test_assumption_system_integration test_assumption_system_integration.cpp)
lmcas_add_standalone_test(test_assumption_integration test_assumption_integration.cpp)

lmcas_add_standalone_test_group(
    test_analysis_group
    test_vector_calculus.cpp
    test_prop_limits.cpp
    test_prop_vector_calculus.cpp
    test_integral_theorems.cpp
    test_ode_engine.cpp
    test_series_engine.cpp
    test_prop_series_engine.cpp
    test_z_transform.cpp
    test_transform_engine.cpp
    test_matrix_decomposition.cpp
    test_complex_analysis.cpp
    test_differential_geometry.cpp
    test_numerical_integration.cpp
)
