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

Region BestMatch(const ReferenceIndex& ri, FragmentIndex&& fi) {
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

namespace detail {

char complement(const char c) noexcept {
  constexpr char map[] = {'T',  '\0', 'G',  '\0', '\0', '\0', 'C',
                          '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                          '\0', '\0', '\0', '\0', '\0', 'A'};
  return map[c - 'A'];
}

}  // namespace detail

::std::string Align(const ::mapper::Sequence& ref,
                    const ::mapper::Sequence& frag, Region&& r, const bool c,
                    const unsigned k, const ::blue::AlignmentType at) {
  if (r.begin.frag.pos == r.end.frag.pos)
    return {};

  auto rc = r.begin.frag.rc != r.begin.ref.rc;
  auto flen = r.end.frag.pos - r.begin.frag.pos + k;
  auto rlen = ((rc ? -1 : 1) * (static_cast<int>(r.end.ref.pos) -
                                static_cast<int>(r.begin.ref.pos))) +
              k;

  rlen = ::std::min(1'000'000u / flen, rlen);

  ::std::unique_ptr<char[]> frag_rc;
  if (rc) {
    frag_rc.reset(new char[flen]);
    for (auto i = 0; i < flen; ++i)
      frag_rc[i] = detail::complement(
          frag.sequence[frag.sequence.length() - 1 - r.end.frag.pos]);
  }

  ::std::string cigar;
  unsigned target_begin;

  using namespace ::blue;

  auto score = c ? PairwiseAlignment(
                       Query{rc ? frag_rc.get()
                                : (frag.sequence.data() + r.begin.frag.pos)},
                       QueryLength{flen},
                       Target{ref.sequence.data() +
                              (rc ? r.end.ref.pos : r.begin.ref.pos)},
                       TargetLength{rlen}, at, Match{1}, Mismatch{1}, Gap{1},
                       cigar, target_begin)
                 : PairwiseAlignment(
                       Query{rc ? frag_rc.get()
                                : (frag.sequence.data() + r.begin.frag.pos)},
                       QueryLength{flen},
                       Target{ref.sequence.data() +
                              (rc ? r.end.ref.pos : r.begin.ref.pos)},
                       TargetLength{rlen}, at, Match{1}, Mismatch{1}, Gap{1});

  ::std::string ret = frag.name;
  ret += "\t" + ::std::to_string(frag.sequence.length());
  ret += "\t" + ::std::to_string(r.begin.frag.pos);
  ret += "\t" + ::std::to_string(r.end.frag.pos + k);

  ret += rc ? "\t-" : "\t+";

  ret += "\t" + ref.name;
  ret += "\t" + ::std::to_string(ref.sequence.length());
  ret += "\t" +
         ::std::to_string(rc ? (rlen - 1 - r.begin.ref.pos) : r.begin.ref.pos);
  ret += "\t" + ::std::to_string(
                    (rc ? (rlen - 1 - r.end.ref.pos) : r.end.ref.pos) + k);

  ret += "\t0";
  ret += "\t" + ::std::to_string(score);

  ret += "\t" + ref.quality;

  if (c) {
    ret += "\tcg:Z:" + cigar;
  }

  return ret;
}

}  // namespace matcher