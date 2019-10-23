cmake_minimum_required(VERSION 3.1)

project(blue_mapper)

set(BLUE_MAPPER_VERSION_MAJOR 0)
set(BLUE_MAPPER_VERSION_MINOR 1)
set(BLUE_MAPPER_VERSION_PATCH 0)

set(CMAKE_CXX_STANDARD 17)

configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/config.hpp"
)

include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/vendor/bioparser/include")
add_subdirectory("vendor/bioparser" EXCLUDE_FROM_ALL)

add_executable(blue_mapper "src/main.cpp")
target_link_libraries(blue_mapper bioparser)