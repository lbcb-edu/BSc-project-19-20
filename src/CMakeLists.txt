add_executable(
    orange_mapper
        orange_mapper.cc
)

target_link_libraries(
    orange_mapper
        PRIVATE
            bioparser
)
