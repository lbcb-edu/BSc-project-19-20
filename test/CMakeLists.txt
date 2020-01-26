include(GoogleTest)

add_executable(
    orange_test
        ornage_test.cc
)

target_include_directories(
    orange_test
        PRIVATE
            ../src
)

target_link_libraries(
    orange_test
        gtest_main

        orange_alignment
        orange_minimizers
)

enable_testing()

gtest_add_tests(
    TARGET orange_test
)