#include "matcher.hpp"

namespace matcher {

ReferenceIndex CreateReferenceIndex(
    const ::std::unique_ptr<::mapper::Sequence>& reference, const unsigned k,
    const unsigned w, const double f) {
  ::std::unordered_map<::blue::KMerInfo, unsigned> occurences;
  for (auto& kmer : ::blue::minimizers(
           reference->sequence.data(),
           ::blue::SequenceLength{
               static_cast<unsigned>(reference->sequence.length())},
           ::blue::KType{k}, ::blue::WindowLength{w}))
    ++occurences[kmer];

  ::std::vector<::std::pair<::blue::KMerInfo, unsigned>> index(
      ::std::make_move_iterator(occurences.begin()),
      ::std::make_move_iterator(occurences.end()));

  ::std::sort(index.begin(), index.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
  index.resize(index.size() * (1.0 - f));

  ::std::sort(index.begin(), index.end(), [](const auto& a, const auto& b) {
    return a.first.pos < b.first.pos;
  });

  ReferenceIndex ret;
  for (auto& [kmer, _] : index)
    ret[kmer.kmer].emplace_back(::std::move(kmer));

  return ret;
}

::std::vector<::blue::KMerInfo> CreateFragmentIndex(
    const ::std::unique_ptr<::mapper::Sequence>& fragment, const unsigned k,
    const unsigned w) {
  return ::blue::minimizers(fragment->sequence.data(),
                            ::blue::SequenceLength{static_cast<unsigned>(
                                fragment->sequence.length())},
                            ::blue::KType{k}, ::blue::WindowLength{w});
}

namespace detail {

bool r(const KMerMatch& a, const KMerMatch& b) noexcept {
  return a.ref.pos < b.ref.pos;
}

bool i(const KMerMatch& a, const KMerMatch& b) noexcept {
  return a.ref.pos > b.ref.pos;
}

}  // namespace detail

Region BestMatch(const ReferenceIndex& ri, FragmentIndex fi) {
  ::std::vector<::std::vector<KMerMatch>> matches(2);
  constexpr bool (*cmp[])(const KMerMatch&, const KMerMatch&) = {detail::r,
                                                                 detail::i};

  for (auto f_kmer : fi) {
    auto iter = ri.find(f_kmer.kmer);
    auto idx = f_kmer.rc & 1;
    if (iter != ri.end())
      for (auto r_kmer : iter->second)
        matches[idx ^ r_kmer.rc].emplace_back(f_kmer, r_kmer);
  }

  Sequence best_seq = {0, 0, 0};
  auto type = -1;

  for (auto i = 0; i < matches.size(); ++i)
    if (matches[i].size()) {
      auto seq = LIS(matches[i], cmp[i]);
      if (seq.len > best_seq.len)
        best_seq = seq, type = i;
    }

  if (type == -1)
    return {};

  return {matches[type][best_seq.begin], matches[type][best_seq.end]};
}

}  // namespace matcher