if(NOT DEFINED LMCAS_SOURCE_DIR)
    message(FATAL_ERROR "LMCAS_SOURCE_DIR is required")
endif()

set(failures "")
foreach(entry IN ITEMS
    "include/symbolic_ast.hpp|80"
    "include/integration.hpp|80"
    "include/visitors/limit_visitor.hpp|200"
    "include/visitors/normalization_visitor.hpp|200"
    "include/bigint.hpp|400"
    "src/bigint.cpp|900"
    "src/integration_table.cpp|900"
    "src/integration_rational.cpp|900"
    "src/inference_engine_arithmetic.cpp|1200"
    "src/lsr_expr.cpp|1650"
    "src/lsr_expr_facade.cpp|950"
    "src/multivariate_factor_patterns.cpp|400"
    "src/multivariate_factor.cpp|1500"
    "src/multivariate_hensel.cpp|650"
    "src/ode_classification.cpp|800"
    "src/ode_linear.cpp|700"
    "src/symbolic_core.cpp|1400"
    "src/symbolic_factor.cpp|950"
    "src/transcendental_factor.cpp|1000"
    "src/transcendental_rewrite.cpp|850"
    "src/zassenhaus_combine.cpp|700"
    "src/solver_groebner.cpp|1200"
    "src/vector_calculus_integrals.cpp|1500")
    string(REPLACE "|" ";" fields "${entry}")
    list(GET fields 0 relative_path)
    list(GET fields 1 maximum_lines)
    set(path "${LMCAS_SOURCE_DIR}/${relative_path}")
    file(READ "${path}" contents)
    string(REGEX REPLACE "[^\n]" "" line_breaks "${contents}")
    string(LENGTH "${line_breaks}" line_count)
    math(EXPR line_count "${line_count} + 1")
    if(line_count GREATER maximum_lines)
        string(APPEND failures
            "\n  ${relative_path}: ${line_count} lines (maximum ${maximum_lines})")
    endif()
endforeach()

file(GLOB_RECURSE implementation_files "${LMCAS_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE public_headers "${LMCAS_SOURCE_DIR}/include/*.hpp")
foreach(path IN LISTS implementation_files public_headers)
    if(path IN_LIST implementation_files)
        set(maximum_lines 1650)
    else()
        set(maximum_lines 800)
    endif()
    file(READ "${path}" contents)
    string(REGEX REPLACE "[^\n]" "" line_breaks "${contents}")
    string(LENGTH "${line_breaks}" line_count)
    math(EXPR line_count "${line_count} + 1")
    if(line_count GREATER maximum_lines)
        file(RELATIVE_PATH relative_path "${LMCAS_SOURCE_DIR}" "${path}")
        string(APPEND failures
            "\n  ${relative_path}: ${line_count} lines (maximum ${maximum_lines})")
    endif()
endforeach()

file(READ "${LMCAS_SOURCE_DIR}/include/assumption_context.hpp" assumption_context_header)
if(assumption_context_header MATCHES "query_interface\\.hpp")
    string(APPEND failures
        "\n  include/assumption_context.hpp: assumptions must not depend on the query facade")
endif()

file(READ "${LMCAS_SOURCE_DIR}/src/symbolic_core.cpp" symbolic_core_source)
foreach(forbidden_header IN ITEMS
    "integration.hpp"
    "solver.hpp"
    "solve_polynomial.hpp"
    "assumption_context.hpp"
    "query_interface.hpp"
    "series_engine.hpp"
    "transform_engine.hpp"
    "vector_calculus.hpp"
    "symbolic_ode.hpp"
    "multivariate_factor.hpp")
    if(symbolic_core_source MATCHES "${forbidden_header}")
        string(APPEND failures
            "\n  src/symbolic_core.cpp: core must not depend on ${forbidden_header}")
    endif()
endforeach()

foreach(calculus_source IN ITEMS
    integrator.cpp
    integration_table.cpp
    integration_basic_strategies.cpp
    integration_trigonometric.cpp
    integration_rational.cpp
    integration_special_functions.cpp
    integration_weierstrass.cpp
    integration_trig_substitution.cpp
    multiple_integral.cpp)
    file(READ "${LMCAS_SOURCE_DIR}/src/${calculus_source}" contents)
    foreach(forbidden_header IN ITEMS
        "solver.hpp"
        "series_engine.hpp"
        "transform_engine.hpp"
        "vector_calculus.hpp"
        "symbolic_ode.hpp")
        if(contents MATCHES "${forbidden_header}")
            string(APPEND failures
                "\n  src/${calculus_source}: calculus must not depend on ${forbidden_header}")
        endif()
    endforeach()
endforeach()

foreach(path IN LISTS implementation_files public_headers)
    file(READ "${path}" contents)
    file(RELATIVE_PATH relative_path "${LMCAS_SOURCE_DIR}" "${path}")
    if(contents MATCHES "InferenceErrorPropagation" OR
       contents MATCHES "MultipleIntegralEngine")
        string(APPEND failures
            "\n  ${relative_path}: duplicate control abstraction is forbidden")
    endif()
endforeach()

foreach(path IN LISTS implementation_files public_headers)
    file(READ "${path}" contents)
    file(RELATIVE_PATH relative_path "${LMCAS_SOURCE_DIR}" "${path}")
    foreach(forbidden_name IN ITEMS
        "ProductNode_Op"
        "NormalizationVisitorImpl"
        "differentiate_legacy"
        "solve_dispatch_checked"
        "ComplexExprResult"
        "SeriesExprResult"
        "SymbolicExprResult"
        "MatrixExprResult"
        "VectorExprResult")
        if(contents MATCHES "${forbidden_name}")
            string(APPEND failures
                "\n  ${relative_path}: retired name ${forbidden_name} is forbidden")
        endif()
    endforeach()
endforeach()

foreach(path IN LISTS implementation_files public_headers)
    file(READ "${path}" contents)
    file(RELATIVE_PATH relative_path "${LMCAS_SOURCE_DIR}" "${path}")
    if(contents MATCHES "std::(cout|cerr|clog)" OR
       contents MATCHES "(^|[^A-Za-z_])(printf|fprintf|puts)[ \\t\\r\\n]*\\(")
        string(APPEND failures
            "\n  ${relative_path}: library diagnostics must use ComputationContext")
    endif()
endforeach()

if(failures)
    message(FATAL_ERROR "Source structure limits exceeded:${failures}")
endif()
