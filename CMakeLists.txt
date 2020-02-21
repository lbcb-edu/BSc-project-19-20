cmake_minimum_required(VERSION 3.2)

project(
    pink
        VERSION 1.0
        LANGUAGES CXX
)

set(PINK_MAPPER_VERSION_MAJOR 0)
set(PINK_MAPPER_VERSION_MINOR 1)
set(PINK_MAPPER_VERSION_PATCH 0)

# TODO: change default output paths for binaries/libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(pink_mapper_build_tests "Build pink mapper unit tests" ON)

add_executable(
    pink_mapper
    ${PROJECT_SOURCE_DIR}/src/pink_mapper.cpp
)

add_subdirectory(
    vendor/bioparser EXCLUDE_FROM_ALL
)

add_subdirectory(
    vendor/thread_pool EXCLUDE_FROM_ALL
)

add_library(
    pink_alignment src/pink_alignment.cpp
)

add_library(
    pink_minimizers src/pink_minimizers.cpp
)


target_link_libraries(
    pink_mapper
        bioparser
	pink_alignment
	pink_minimizers
	thread_pool
)

if (pink_mapper_build_tests)
    include_directories(${PROJECT_SOURCE_DIR}/src)

    add_executable(pink_mapper_test test/pink_mapper_test.cpp)

    add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)

    target_link_libraries(pink_mapper_test pink_alignment pink_minimizers thread_pool gtest_main)

endif (pink_mapper_build_tests)
