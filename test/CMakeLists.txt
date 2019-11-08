include(GoogleTest)

add_executable(
    orange_alignment_test
        ornage_alignment_test.cc
)

target_include_directories(
    orange_alignment_test
        PRIVATE
            ../src
)

target_link_libraries(
    orange_alignment_test
        Threads::Threads
        gtest_main

        orange_alignment
)

enable_testing()

gtest_add_tests(
    TARGET orange_alignment_test
)