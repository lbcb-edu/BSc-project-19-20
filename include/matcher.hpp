#pragma once

#include <alignment/alignment.hpp>
#include <cassert>
#include <memory>
#include <minimizers/minimizers.hpp>
#include <string>
#include <vector>

#include "mapper.hpp"

namespace matcher {

using ReferenceIndex =
    ::std::unordered_map<unsigned, ::std::vector<::blue::KMerInfo>>;
using FragmentIndex = ::std::vector<::blue::KMerInfo>;

struct Sequence {
  unsigned begin;
  unsigned end;
  unsigned len;

  Sequence() noexcept = default;
  Sequence(const unsigned begin, const unsigned end,
           const unsigned len) noexcept
      : begin(begin), end(end), len(len) {}
};

struct KMerMatch {
  ::blue::KMerInfo frag;
  ::blue::KMerInfo ref;

  KMerMatch() noexcept = default;
  KMerMatch(const ::blue::KMerInfo& frag, const ::blue::KMerInfo& ref) noexcept
      : frag(frag), ref(ref) {}
};

struct Region {
  KMerMatch begin;
  KMerMatch end;
};

ReferenceIndex CreateReferenceIndex(
    const ::std::unique_ptr<::mapper::Sequence>& reference, const unsigned k,
    const unsigned w, const double f);

::std::vector<::blue::KMerInfo> CreateFragmentIndex(
    const ::std::unique_ptr<::mapper::Sequence>& fragment, const unsigned k,
    const unsigned w);

Region BestMatch(const ReferenceIndex& ri, FragmentIndex&& fi);

::std::string Align(const ::mapper::Sequence& ref,
                    const ::mapper::Sequence& frag, Region&& r, const bool c,
                    const unsigned k, const ::blue::AlignmentType at);

//

template <typename T, typename C>
Sequence LIS(const ::std::vector<T>& v, const C& comp) {
  ::std::vector<Sequence> sequences{{0, 0, 1}};
  for (unsigned i = 1; i < v.size(); ++i)
    if (comp(v[i], v[sequences.front().end])) {
      sequences[0] = {i, i, 1};
    } else if (auto [p, q, l] = sequences.back(); comp(v[q], v[i])) {
      sequences.emplace_back(p, i, l + 1);
    } else {
      auto iter =
          ::std::lower_bound(sequences.begin(), sequences.end(), v[i],
                             [&v, &comp](auto const& seq, auto const& value) {
                               return comp(v[seq.end], value);
                             });

      if (iter == sequences.begin()) {
        *iter = {i, i, 1};
      } else {
        auto [u, v, k] = *(iter - 1);
        *iter = {u, i, k + 1};
      }
    }
  return sequences.back();
}

}  // namespace matcher