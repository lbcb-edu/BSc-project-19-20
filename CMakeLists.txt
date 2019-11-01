cmake_minimum_required(VERSION 3.9)

project(orange_mapper LANGUAGES CXX)

set(orange_mapper_VERSION_MAJOR 0)
set(orange_mapper_VERSION_MINOR 1)
set(orange_mapper_VERSION_PATCH 0)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

configure_file (
    "${PROJECT_SOURCE_DIR}/src/orange_mapper_conf.h.in"
    "${PROJECT_BINARY_DIR}/conf/orange_mapper_conf.h"
)

include_directories(
    "${PROJECT_BINARY_DIR}/conf"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")

# External dependencies
find_package(Threads)

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/googletest EXCLUDE_FROM_ALL)

# Project subdirs 
add_subdirectory(src)
add_subdirectory(test)

