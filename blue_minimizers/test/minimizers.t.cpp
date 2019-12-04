#include <iostream>

#include <gtest/gtest.h>
#include <minimizers/minimizers.hpp>

namespace blue {
namespace {

TEST(MinimizerTests, InitialTest) {
  for (auto&& kmer :
       minimizers("AAACTG", SequenceLength{6}, KType{2}, WindowLength{2})) {
    ::std::cout << ::std::get<0>(kmer) << " " << ::std::get<1>(kmer) << " "
                << ::std::get<2>(kmer) << "\n";
  }
}

}  // namespace
}  // namespace blue

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}