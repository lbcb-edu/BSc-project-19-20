#include <gtest/gtest.h>

#include <iostream>
#include <string>

#include "../include/alignment/alignment.hpp"

namespace algn {
namespace {

TEST(ContiguousArrayTest, IndexingTest) {
  using namespace detail;
  Contiguous2DArray<int> int_arr({2, 3});

  ASSERT_EQ(2, int_arr.rows);
  ASSERT_EQ(3, int_arr.cols);

  for (auto i = 0; i < int_arr.rows; ++i)
    for (auto j = 0; j < int_arr.cols; ++j) int_arr[i][j] = i * j;

  for (auto i = 0; i < int_arr.rows; ++i)
    for (auto j = 0; j < int_arr.cols; ++j) ASSERT_EQ(i * j, int_arr[i][j]);

  auto const c_arr = ::std::move(int_arr);

  ASSERT_EQ(2, c_arr.rows);
  ASSERT_EQ(3, c_arr.cols);

  for (auto i = 0; i < c_arr.rows; ++i)
    for (auto j = 0; j < c_arr.cols; ++j) ASSERT_EQ(i * j, c_arr[i][j]);
}

TEST(NeedleWunschTest, GeneralTest) {
  ::std::string cigar;
  unsigned int target_begin;

  // example taken from
  // https://www.fer.unizg.hr/_download/repository/bioinformatika_skripta_v1.2.pdf#page=24
  auto score = PairwiseAlignment("TGCATAT", 7, "ATCCGAT", 7,
                                 AlignmentType::kNeedlemanWunsch, 2, -1, -1,
                                 cigar, target_begin);

  ASSERT_EQ(0, target_begin);
  ASSERT_EQ("D=X=XI==", cigar);
  ASSERT_EQ(4, score);
}

TEST(SmithWatermanTest, GeneralTest) {
  ::std::string cigar;
  unsigned int target_begin;

  auto score = PairwiseAlignment("ACCTAAGG", 8, "GGCTCAATCA", 10,
                                 AlignmentType::kSmithWaterman, 2, -1, -2,
                                 cigar, target_begin);

  ASSERT_EQ(2, target_begin);
  ASSERT_EQ(6, score);
  ASSERT_EQ("=D==", cigar);
}

}  // namespace
}  // namespace algn

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}