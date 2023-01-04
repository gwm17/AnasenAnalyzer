cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
    message("Building release")
else()
    message("Building debug")
endif()

project(AnasenAnalyzer)

find_package(ROOT REQUIRED COMPONENTS GenVector)

set(AA_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(AA_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_subdirectory(vendor/catima)
add_subdirectory(vendor/fmt)
add_subdirectory(src)