#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <utility>

namespace algn {

enum AlignmentType { kNeedlemanWunsch, kSmithWaterman, kSuffixPrefix };

int PairwiseAlignment(const char* query, unsigned int query_length,
                      const char* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap);

int PairwiseAlignment(const char* query, unsigned int query_length,
                      const char* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap,
                      ::std::string& cigar, unsigned int& target_begin);

namespace detail {

template <typename T>
class Contiguous2DArray {
 private:
  ::std::pair<::std::size_t, ::std::size_t> dim_;
  ::std::unique_ptr<T[]> data_;

 public:
  const ::std::size_t& rows = dim_.first;
  const ::std::size_t& cols = dim_.second;

  Contiguous2DArray(::std::size_t rows, ::std::size_t cols)
      : dim_(rows, cols), data_(new T[dim_.first * dim_.second]) {}

  Contiguous2DArray(const Contiguous2DArray& other) = delete;
  Contiguous2DArray& operator=(const Contiguous2DArray& other) = delete;

  Contiguous2DArray(Contiguous2DArray&& other) = default;
  Contiguous2DArray& operator=(Contiguous2DArray&& other) = default;

  // clang-format off
  T* operator[](::std::size_t index) { 
    return &data_[index * dim_.second]; 
  }

  const T* operator[](::std::size_t index) const {
    return &data_[index * dim_.second];
  }

  // utility
  T& last() { 
    return data_[rows * cols - 1]; 
  }

  const T& last() const { 
    return data_[rows * cols - 1]; 
  }
  // clang-format on
};

}  // namespace detail
}  // namespace algn