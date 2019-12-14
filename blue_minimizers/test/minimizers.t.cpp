#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include <minimizers/minimizers.hpp>

namespace blue {
namespace {

TEST(MinimizerTests, InitialTest) {
  auto kmers = minimizers("AGCTGCACTGATGCTA", SequenceLength{16}, KType{5},
                          WindowLength{8});

  for (auto&& kmer : kmers) {
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