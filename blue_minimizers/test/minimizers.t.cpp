#include <gtest/gtest.h>

#include <minimizers/minimizers.hpp>

namespace blue {
namespace {

TEST(MinimizerTests, InitialTest) {
  ASSERT_EQ((KMerInfo{10, 11, true}),  // lmao @ macros
            minimizers("", SequenceLength{0}, KType{10}, WindowLength{11})[0]);
}

}  // namespace
}  // namespace blue

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}