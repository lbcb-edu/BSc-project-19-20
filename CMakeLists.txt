cmake_minimum_required(VERSION 3.5)

project(white_mapper)

set(CMAKE_CXX_STANDARD 11)

set(WHITE_MAPPER_VERSION_MAJOR 0)
set(WHITE_MAPPER_VERSION_MINOR 1)

configure_file("${PROJECT_SOURCE_DIR}/config.h.in"
	       "${PROJECT_BINARY_DIR}/include/config.h")


include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/vendor/bioparser/include")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_executable(white_mapper src/white_mapper.cpp)

target_link_libraries(white_mapper bioparser)
