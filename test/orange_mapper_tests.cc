#include <gtest/gtest.h>

TEST(DummySuite, DummyDoubleEQ) {
    ASSERT_DOUBLE_EQ(1.0, 1.0);
}

/**
 * @brief Main testing entry point
 */
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv); 
    return RUN_ALL_TESTS();
}