cmake_minimum_required(VERSION 3.2)

project(
    pink
        VERSION 1.0
        LANGUAGES CXX
)

set(PINK_MAPPER_VERSION_MAJOR 0)
set(PINK_MAPPER_VERSION_MINOR 1)
set(PINK_MAPPER_VERSION_PATCH 0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(
    pink_mapper
    ${PROJECT_SOURCE_DIR}/src/pink_mapper.cpp   
)

add_subdirectory(
    vendor/bioparser EXCLUDE_FROM_ALL
)


target_link_libraries(
    pink_mapper
        bioparser
)
