cmake_minimum_required(VERSION 3.26)
project(lmcas_package_consumer LANGUAGES C CXX)

set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(LMCAS CONFIG REQUIRED)
include("${CMAKE_CURRENT_LIST_DIR}/../../LMMC/cmake/PackageConsumerChecks.cmake")
package_check_target(LMCAS::lmcas cxx_std_17)
package_check_target(LMCAS::lmmc c_std_11)
package_check_target(LMCAS::liblmmp "")

add_executable(lmcas_package_consumer
    main.cpp
    expr_checks.cpp
    expr_arithmetic_checks.cpp
    expr_complex_checks.cpp
    expr_sets_checks.cpp
    expr_contracts_checks.cpp
    lmmc_stdlib_checks.cpp
    lmmc_linalg_checks.cpp
    lmmc_matrix_checks.cpp
    lmmc_set_checks.cpp
)
target_link_libraries(lmcas_package_consumer PRIVATE LMCAS::lmcas)


get_target_property(_lmcas_include_dirs LMCAS::lmcas INTERFACE_INCLUDE_DIRECTORIES)
set(_lmcas_install_include "")
foreach(_include_dir IN LISTS _lmcas_include_dirs)
    if(EXISTS "${_include_dir}/symbolic.hpp")
        set(_lmcas_install_include "${_include_dir}")
        break()
    endif()
endforeach()
if(NOT _lmcas_install_include)
    message(FATAL_ERROR "LMCAS::lmcas does not expose an installed include directory")
endif()

foreach(_private_dir internal visitors lmmp/impl lmmc/internal)
    if(EXISTS "${_lmcas_install_include}/${_private_dir}")
        message(FATAL_ERROR
            "Private header directory was installed: ${_lmcas_install_include}/${_private_dir}")
    endif()
endforeach()

package_check_headers(LMCAS::lmcas "${_lmcas_install_include}" "*.hpp" CXX)
package_check_headers(LMCAS::lmmc "${_lmcas_install_include}" "lmmc/*.h" C CXX)
package_check_headers(LMCAS::liblmmp "${_lmcas_install_include}" "lmmp/*.h" C CXX)
package_check_report(lmcas_package_consumer
    LMCAS::lmcas LMCAS::lmmc LMCAS::liblmmp)
