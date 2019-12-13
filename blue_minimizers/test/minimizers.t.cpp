#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include <minimizers/minimizers.hpp>

namespace blue::detail {
namespace {

TEST(MinimizersTest, StorageTest) {
  QuadField qf{19};

  qf.Set(0, 0b11);
  qf.Set(1, 0b01);
  ASSERT_EQ(0b1101, qf.At(0, 2));
  ASSERT_EQ(0b01, qf.At(1, 1));

  qf.Set(15, 0b10);
  qf.Set(16, 0b01);
  qf.Set(17, 0b11);
  ASSERT_EQ(0b100111, qf.At(15, 3));

  qf.Set(18, 0b01);
  ASSERT_EQ(0b1101, qf.At(17, 2));

  qf = QuadField{33};

  qf.Set(31, 0b11);
  ASSERT_EQ(0b11, qf.At(16, 16));

  ::std::vector<unsigned> v{
      0b11, 0b11, 0b10, 0b01, 0b00, 0b11, 0b00, 0b11, 0b10, 0b11, 0b00, 0b11,
      0b11, 0b10, 0b01, 0b00, 0b11, 0b00, 0b11, 0b10, 0b11, 0b00, 0b11, 0b00,
      0b11, 0b10, 0b00, 0b11, 0b00, 0b11, 0b10, 0b00, 0b11, 0b00, 0b10, 0b10,
      0b01, 0b00, 0b11, 0b00, 0b11, 0b10, 0b11, 0b00, 0b11, 0b00, 0b11, 0b10};

  qf = QuadField{v.size()};
  for (int i = 0; i < v.size(); ++i) {
    qf.Set(i, v[i]);
  }

  for (int k = 1; k <= 16; ++k) {
    for (int i = 0; i < v.size() - k; ++i) {
      unsigned val = 0;
      for (int j = 0; j < k; ++j) {
        val <<= 2;
        val |= v[i + j];
      }
      ASSERT_EQ(val, qf.At(i, k));
    }
  }
}

TEST(MinimizerTests, InitialTest) {
  for (auto&& kmer :
       minimizers("AAACTG", SequenceLength{6}, KType{2}, WindowLength{2})) {
    ::std::cout << ::std::get<0>(kmer) << " " << ::std::get<1>(kmer) << " "
                << ::std::get<2>(kmer) << "\n";
  }
}

}  // namespace
}  // namespace blue::detail

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}