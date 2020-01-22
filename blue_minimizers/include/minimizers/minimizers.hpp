#pragma once

#include <algorithm>
#include <bitset>
#include <common/strong_type.hpp>
#include <functional>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace blue {

struct KMerInfo;

using SequenceLength = StrongType<unsigned, struct SequenceLengthTag>;
using KType = StrongType<unsigned, struct KTypeTag>;
using WindowLength = StrongType<unsigned, struct WindowLengthTag>;

::std::vector<KMerInfo> minimizers(const char* sequence,
                                   SequenceLength sequence_length, KType k,
                                   WindowLength window_length);

//

struct KMerInfo {
  unsigned kmer;
  unsigned pos;
  bool rc;

  KMerInfo() noexcept = default;
  KMerInfo(const unsigned kmer, const unsigned pos, const bool rc) noexcept
      : kmer(kmer), pos(pos), rc(rc) {}

  friend bool operator<(const KMerInfo& a, const KMerInfo& b) noexcept {
    return a.kmer < b.kmer || (a.kmer == b.kmer && a.pos < b.pos) ||
           (a.kmer == b.kmer && a.pos == b.pos && !a.rc && b.rc);
  }

  friend bool operator==(const KMerInfo& a, const KMerInfo& b) noexcept {
    return a.kmer == b.kmer && a.pos == b.pos && a.rc == b.rc;
  }

  friend void swap(KMerInfo& a, KMerInfo& b) noexcept {
    using ::std::swap;
    swap(a.kmer, b.kmer);
    swap(a.pos, b.pos);
    swap(a.rc, b.rc);
  }

  friend ::std::ostream& operator<<(::std::ostream& out,
                                    const ::blue::KMerInfo k) {
    out << k.kmer << " " << k.pos << " " << k.rc;
    return out;
  }
};

}  // namespace blue

namespace std {

template <>
struct hash<::blue::KMerInfo> {
  ::std::size_t operator()(const ::blue::KMerInfo& ki) const noexcept {
    ::std::hash<unsigned> uh;
    ::std::size_t seed{0};
    seed += uh(ki.kmer) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed += uh(ki.pos) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed += ::std::hash<bool>{}(ki.rc) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

}  // namespace std