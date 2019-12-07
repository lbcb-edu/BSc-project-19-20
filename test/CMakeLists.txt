include(GoogleTest)

add_executable(
    orange_alignment_test
        ornage_alignment_test.cc
)

add_executable(
    orange_minimizers_test
        ornage_minimizers_test.cc
)

target_include_directories(
    orange_alignment_test
        PRIVATE
            ../src
)

target_include_directories(
    orange_minimizers_test
        PRIVATE
            ../src
)

target_link_libraries(
    orange_alignment_test
        Threads::Threads
        gtest_main

        orange_alignment
)

target_link_libraries(
    orange_minimizers_test
        Threads::Threads    
        gtest_main

        orange_minimizers
)

enable_testing()

gtest_add_tests(
    TARGET orange_alignment_test
    TARGET orange_minimizers_test
)