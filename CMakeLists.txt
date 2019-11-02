cmake_minimum_required(VERSION 3.5)

project(blue_mapper)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_executable(blue_mapper src/blue_mapper.cpp)

