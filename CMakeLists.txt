cmake_minimum_required (VERSION 2.6)
project (Mapper)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set (Mapper_VERSION_MAJOR 0)
set (Mapper_VERSION_MINOR 1)
configure_file (
  "${PROJECT_SOURCE_DIR}/MapperConfig.h.in"
  "${PROJECT_BINARY_DIR}/src/MapperConfig.h"
  )

include_directories("${PROJECT_BINARY_DIR}/vendor/bioparser/include")
add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_executable (Mapper src/brown_mapper.cpp)
target_link_libraries (Mapper bioparser)
