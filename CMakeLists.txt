cmake_minimum_required (VERSION 2.6)
project (Mapper)

set (Mapper_VERSION_MAJOR 0)
set (Mapper_VERSION_MINOR 1)
configure_file (
  "${PROJECT_SOURCE_DIR}/MapperConfig.h.in"
  "${PROJECT_BINARY_DIR}/MapperConfig.h"
  )

include_directories("${PROJECT_BINARY_DIR}")
add_subdirectory(bioparser EXCLUDE_FROM_ALL)

add_executable (Mapper brown_mapper.cpp)
target_link_libraries (Mapper bioparser)
