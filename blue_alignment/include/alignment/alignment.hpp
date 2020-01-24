#pragma once

#include <algorithm>
#include <common/strong_type.hpp>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace blue {

enum AlignmentType { kNeedlemanWunsch, kSmithWaterman, kOverlap };

using Query = StrongType<const char*, struct QueryTag>;
using QueryLength = StrongType<unsigned int, struct QueryLengthTag>;

using Target = StrongType<const char*, struct TargetTag>;
using TargetLength = StrongType<unsigned int, struct TargetLengthTag>;

using Match = StrongType<int, struct MatchTag>;
using Mismatch = StrongType<int, struct MismatchTag>;
using Gap = StrongType<int, struct GapTag>;

int PairwiseAlignment(Query query, QueryLength query_length, Target target,
                      TargetLength target_length, AlignmentType type,
                      Match match, Mismatch mismatch, Gap gap);

int PairwiseAlignment(Query query, QueryLength query_length, Target target,
                      TargetLength target_length, AlignmentType type,
                      Match match, Mismatch mismatch, Gap gap,
                      ::std::string& cigar, unsigned int& target_begin);

namespace detail {

template <typename T>
class Contiguous2DArray {
 private:
  ::std::pair<::std::size_t, ::std::size_t> dim_;
  ::std::unique_ptr<T[]> data_;

 public:
  const ::std::size_t rows = dim_.first;
  const ::std::size_t cols = dim_.second;

  Contiguous2DArray(::std::size_t rows, ::std::size_t cols)
      : dim_(rows, cols), data_(new T[dim_.first * dim_.second]) {}

  Contiguous2DArray(const Contiguous2DArray& other) = delete;
  Contiguous2DArray& operator=(const Contiguous2DArray& other) = delete;

  Contiguous2DArray(Contiguous2DArray&& other) = default;
  Contiguous2DArray& operator=(Contiguous2DArray&& other) = default;

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
};

}  // namespace detail
}  // namespace blue