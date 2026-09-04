cmake_minimum_required(VERSION 3.26)
project(lmcas_package_consumer LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(LMCAS CONFIG REQUIRED)

add_executable(lmcas_package_consumer
    main.cpp
    lmmc_linalg_checks.cpp
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

foreach(_private_dir internal visitors lmmp/impl)
    if(EXISTS "${_lmcas_install_include}/${_private_dir}")
        message(FATAL_ERROR
            "Private header directory was installed: ${_lmcas_install_include}/${_private_dir}")
    endif()
endforeach()

file(GLOB _lmcas_headers CONFIGURE_DEPENDS
    "${_lmcas_install_include}/*.hpp")
file(GLOB _lmmc_headers CONFIGURE_DEPENDS
    "${_lmcas_install_include}/lmmc/*.h")
file(GLOB _lmmp_headers CONFIGURE_DEPENDS
    "${_lmcas_install_include}/lmmp/*.h")

function(add_standalone_header_check header language dependency)
    get_filename_component(_name "${header}" NAME)
    string(MAKE_C_IDENTIFIER "${language}_${_name}" _target)
    if(language STREQUAL "C")
        set(_extension c)
    else()
        set(_extension cpp)
    endif()
    set(_source "${CMAKE_CURRENT_BINARY_DIR}/standalone_${_target}.${_extension}")
    file(WRITE "${_source}" "#include <${header}>\n")
    add_library("standalone_${_target}" OBJECT "${_source}")
    target_link_libraries("standalone_${_target}" PRIVATE "${dependency}")
endfunction()

foreach(_header IN LISTS _lmcas_headers)
    get_filename_component(_name "${_header}" NAME)
    add_standalone_header_check("${_name}" CXX LMCAS::lmcas)
endforeach()
foreach(_header IN LISTS _lmmc_headers)
    get_filename_component(_name "${_header}" NAME)
    add_standalone_header_check("lmmc/${_name}" C LMCAS::lmmc)
endforeach()
foreach(_header IN LISTS _lmmp_headers)
    get_filename_component(_name "${_header}" NAME)
    add_standalone_header_check("lmmp/${_name}" C LMCAS::liblmmp)
endforeach()
