cmake_minimum_required(VERSION 3.5)

project(white_mapper)

set(CMAKE_CXX_STANDARD 11)

set(WHITE_MAPPER_VERSION_MAJOR 0)
set(WHITE_MAPPER_VERSION_MINOR 1)

configure_file("${PROJECT_SOURCE_DIR}/config.h.in"
	       "${PROJECT_BINARY_DIR}/include/config.h")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(white_alignment EXCLUDE_FROM_ALL)
add_subdirectory(vendor/googletest EXCLUDE_FROM_ALL)

add_executable(white_mapper src/white_mapper.cpp)

target_link_libraries(
		white_mapper PUBLIC
		bioparser
		white_alignment
		)

target_include_directories(
		white_mapper PUBLIC
		"${PROJECT_SOURCE_DIR}/include"
		"${PROJECT_BINARY_DIR}/vendor/bioparser/include"
		"${PROJECT_BINARY_DIR}/white_alignment"
		"${PROJECT_BINARY_DIR}/vendor/googletest/googletest/include"
		"${PROJECT_BINARY_DIR}/vendor/googletest/googlemock/include"
		)

add_executable(alignment_test tests/alignment_test.cpp)

target_link_libraries(alignment_test white_alignment gtest)


