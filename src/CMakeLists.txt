add_executable(
    orange_mapper
        orange_mapper.cc
)

add_library(orange_mapper_conf INTERFACE)
target_include_directories(orange_mapper_conf INTERFACE ./)

add_library(
    orange_alignment 
        orange_alignment.h
        orange_alignment.cc
)

add_library(
    orange_minimizers
        orange_minimizers.h
        orange_minimizers.cc
)

target_link_libraries(
    orange_mapper
        PRIVATE
            orange_mapper_conf
            orange_minimizers
            orange_alignment
            bioparser
)
