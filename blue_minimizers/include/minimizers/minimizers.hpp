#pragma once

#include <algorithm>
#include <bitset>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <iostream>

#include <common/strong_type.hpp>

namespace blue {

using KMerInfo = ::std::tuple<unsigned, unsigned, bool>;

using SequenceLength = StrongType<unsigned, struct SequenceLengthTag>;
using KType = StrongType<unsigned, struct KTypeTag>;
using WindowLength = StrongType<unsigned, struct WindowLengthTag>;

::std::vector<KMerInfo> minimizers(const char* sequence,
                                   SequenceLength sequence_length, KType k,
                                   WindowLength window_length);

}  // namespace blue

namespace std {

template <>
struct hash<::blue::KMerInfo> {
  ::std::size_t operator()(const ::blue::KMerInfo& ki) const noexcept {
    return ::std::apply([](auto&&... args) { return combine(0, args...); }, ki);
  }

 private:
  template <typename T>
  static ::std::size_t combine(const ::std::size_t seed, const T arg) noexcept {
    return seed ^
           (::std::hash<T>{}(arg) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
  }

  template <typename T, typename... Args>
  static ::std::size_t combine(const ::std::size_t seed, const T arg,
                               Args&&... args) noexcept {
    return combine(combine(seed, arg), args...);
  }
};

}  // namespace std