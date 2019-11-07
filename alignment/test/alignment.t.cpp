#include <gtest/gtest.h>

#include <iostream>

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

}  // namespace
}  // namespace algn

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}