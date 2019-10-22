add_executable(
    orange_mapper
        orange_mapper.cc
)

add_library(orange_mapper_conf INTERFACE)
target_include_directories(orange_mapper_conf INTERFACE ./)

target_link_libraries(
    orange_mapper
        PRIVATE
            orange_mapper_conf
            bioparser
)
