#pragma once

#include <algorithm>
#include <bitset>
#include <memory>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <common/strong_type.hpp>

namespace blue {

using KMerInfo = ::std::tuple<unsigned, unsigned, bool>;

using SequenceLength = StrongType<unsigned, struct SequenceLengthTag>;
using KType = StrongType<unsigned, struct KTypeTag>;
using WindowLength = StrongType<unsigned, struct WindowLengthTag>;

::std::vector<KMerInfo> minimizers(const char* sequence,
                                   SequenceLength sequence_length, KType k,
                                   WindowLength window_length);

//
//
//
//
//

namespace detail {  // declared here so the class is testable

using UnsignedIntegralType = unsigned;

// used in QuadField and quad_helpers functiong group
constexpr ::std::size_t kBits = sizeof(UnsignedIntegralType) * 8;

namespace quad_helpers {

inline UnsignedIntegralType Prefix(UnsignedIntegralType data,
                                   ::std::size_t count) noexcept {
  return data >> (kBits - count);
}

inline UnsignedIntegralType Infix(UnsignedIntegralType data,
                                  ::std::size_t index,
                                  ::std::size_t count) noexcept {
  return (data >> (kBits - index - count)) & ((1 << count) - 1);
}

};  // namespace quad_helpers

class QuadField {
 private:
  constexpr static ::std::size_t points_ = kBits / 2;

  ::std::size_t size_;
  ::std::unique_ptr<UnsignedIntegralType[]> data_;

 public:
  QuadField(::std::size_t size)
      : size_(size),
        data_(new UnsignedIntegralType[size / points_ +
                                       (size % points_ ? 1 : 0)]) {}

  QuadField(const QuadField&) = delete;
  QuadField& operator=(const QuadField&) = delete;

  QuadField(QuadField&&) = default;
  QuadField& operator=(QuadField&&) = default;

  UnsignedIntegralType At(::std::size_t index, ::std::size_t count) const {
    auto idx = index / points_;

    if (!(index % points_))
      return quad_helpers::Prefix(data_[idx], count * 2);

    if (idx == (index + count - 1) / points_)
      return quad_helpers::Infix(data_[idx], index * 2, count * 2);

    auto rem = (index + count) % points_;

    return (quad_helpers::Infix(data_[idx], index * 2, (points_ - index) * 2)
            << (rem * 2)) |
           (quad_helpers::Prefix(data_[idx + 1], rem * 2));
  }

  void Set(::std::size_t index, UnsignedIntegralType value) {
    data_[index / points_] &= ~(0b11 << (kBits - (index % points_) * 2 - 2));
    data_[index / points_] |= value << (kBits - (index % points_) * 2 - 2);
  }

  // also prints out the remainder of last integer (˵ ͡° ͜ʖ ͡°˵)
  friend ::std::ostream& operator<<(::std::ostream& out, const QuadField& qf) {
    for (int i = 0; i <= qf.size_ / qf.points_; ++i) {
      out << ::std::bitset<kBits>{qf.data_[i]};
      if (i != qf.size_ / qf.points_)
        out << ":";
    }
    return out;
  }
};

}  // namespace detail

}  // namespace blue