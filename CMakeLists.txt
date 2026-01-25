cmake_minimum_required(VERSION 3.10)
project(LMCAS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/utf-8 /permissive-)
endif()

# 包含 LAMMP 子项目
add_subdirectory(LAMMP)

# LMCAS 库
add_library(lmcas SHARED
    cas.cpp
    symbolic.cpp
)

target_compile_definitions(lmcas PRIVATE LAMINA_CORE_EXPORTS)

if(MSVC)
    target_compile_options(lmcas PRIVATE /utf-8 /permissive-)
endif()

# 包含路径
target_include_directories(lmcas PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/LAMMP/include
)

# 链接 LAMMP 核心库
target_link_libraries(lmcas PRIVATE LammpCore)

# 测试程序
add_executable(test_cas tests/test_main.cpp)
target_link_libraries(test_cas PRIVATE lmcas LammpCore)

add_executable(test_bigint_debug tests/test_bigint_debug.cpp)
target_link_libraries(test_bigint_debug PRIVATE lmcas LammpCore)

add_executable(test_proof tests/test_proof.cpp)
target_link_libraries(test_proof PRIVATE lmcas LammpCore)

