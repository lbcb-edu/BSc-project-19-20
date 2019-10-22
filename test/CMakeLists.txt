include(GoogleTest)

add_executable(
    orange_mapper_tests
        orange_mapper_tests.cc
)

target_link_libraries(
    orange_mapper_tests
        Threads::Threads
        gtest_main
)

enable_testing()

gtest_add_tests(
    TARGET orange_mapper_tests
    TEST_SUFFIX .noArgs
    TEST_LIST   noArgsTests
)