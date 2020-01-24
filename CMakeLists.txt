cmake_minimum_required(VERSION 3.1)

project(blue_mapper)

set(BLUE_MAPPER_VERSION_MAJOR 0)
set(BLUE_MAPPER_VERSION_MINOR 2)
set(BLUE_MAPPER_VERSION_PATCH 0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-Ofast -g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")

find_package(Threads)

configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/config.hpp"
)

add_subdirectory("vendor/bioparser" EXCLUDE_FROM_ALL)
add_subdirectory("vendor/thread_pool" EXCLUDE_FROM_ALL)
add_subdirectory("vendor/googletest" EXCLUDE_FROM_ALL)
add_subdirectory("blue_alignment" EXCLUDE_FROM_ALL)
add_subdirectory("blue_minimizers" EXCLUDE_FROM_ALL)
add_subdirectory("blue_common" EXCLUDE_FROM_ALL)

add_executable(blue_mapper "src/main.cpp" "src/matcher.cpp")

target_link_libraries(
  blue_mapper PUBLIC 
  blue_alignment 
  blue_common
  blue_minimizers
  bioparser
  thread_pool
)

target_include_directories(
  blue_mapper PUBLIC 
  "${PROJECT_BINARY_DIR}/blue_alignment/include"
  "${PROJECT_BINARY_DIR}/blue_common/include"
  "${PROJECT_BINARY_DIR}/blue_minimizers/include"
  "${PROJECT_BINARY_DIR}/include"
  "${PROJECT_BINARY_DIR}/vendor/bioparser/include"
  "${PROJECT_BINARY_DIR}/vendor/thread_pool/include"
)

add_executable(alignment_tests "blue_alignment/test/alignment.t.cpp")

target_link_libraries(
  alignment_tests PUBLIC
  blue_alignment
  blue_common
  gtest
  pthread
  ${CMAKE_THREAD_LIBS_INIT}
)

target_include_directories(
  alignment_tests PUBLIC
  "${PROJECT_BINARY_DIR}/blue_alignment/include"
  "${PROJECT_BINARY_DIR}/blue_common/include"
  "${PROJECT_BINARY_DIR}/vendor/googletest/googletest/include"
)

add_executable(minimizers_tests "blue_minimizers/test/minimizers.t.cpp")

target_link_libraries(
  minimizers_tests PUBLIC
  blue_common
  blue_minimizers
  gtest
  pthread
  ${CMAKE_THREAD_LIBS_INIT}
)

target_include_directories(
  minimizers_tests PUBLIC
  "${PROJECT_BINARY_DIR}/blue_common/include"
  "${PROJECT_BINARY_DIR}/blue_minimizers/include"
  "${PROJECT_BINARY_DIR}/vendor/googletest/googletest/include"
)