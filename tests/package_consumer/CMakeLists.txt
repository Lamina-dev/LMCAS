cmake_minimum_required(VERSION 3.16)
project(lmcas_package_consumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(LMCAS CONFIG REQUIRED)

add_executable(lmcas_package_consumer main.cpp)
target_link_libraries(lmcas_package_consumer PRIVATE LMCAS::lmcas)
