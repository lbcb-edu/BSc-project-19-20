cmake_minimum_required (VERSION 3.2)
project (Mapper)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (Mapper_VERSION_MAJOR 0)
set (Mapper_VERSION_MINOR 4)
configure_file (
  "${PROJECT_SOURCE_DIR}/MapperConfig.h.in"
  "${PROJECT_BINARY_DIR}/src/MapperConfig.h"
  )

include_directories ("${PROJECT_BINARY_DIR}/vendor/bioparser/include")
add_subdirectory (vendor/bioparser EXCLUDE_FROM_ALL)
#add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
add_subdirectory (brown_alignment)
add_subdirectory (brown_minimizers)
option(brown_tests "tests" ON)

add_executable (Mapper src/brown_mapper.cpp)
target_link_libraries (Mapper bioparser)
target_link_libraries (Mapper brown_alignment)
target_link_libraries (Mapper brown_minimizers)
#target_link_libraries (Mapper thread_pool)
target_include_directories(Mapper PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/brown_alignment"
                          "${PROJECT_SOURCE_DIR}/brown_minimizers"
                          )

if(brown_tests)
	include_directories(${PROJECT_SOURCE_DIR}/src)

	add_executable(brown_tests test/brown_tests.cpp)

	add_subdirectory(vendor/googletest EXCLUDE_FROM_ALL)

	target_link_libraries(brown_tests brown_alignment gtest_main brown_minimizers)
endif(brown_tests)
