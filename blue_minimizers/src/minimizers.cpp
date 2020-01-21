#include "minimizers/minimizers.hpp"

namespace blue {

namespace detail {

inline unsigned Complement(const unsigned c) noexcept {
  return 3 - c;
}

inline unsigned CustomOrdering(const char c) noexcept {
  switch (c) {
    case 'C':
      return 0;
    case 'A':
      return 1;
    case 'T':
      return 2;
    case 'G':
      return 3;
  }

  return 0;
}

::std::vector<KMerInfo> GenerateFrom(const char* seq, unsigned l, unsigned k) {
  unsigned reg{0}, inv{0};

  using ::std::get;

  for (auto i = 0; i < k; ++i)
    reg = (reg << 2) | CustomOrdering(seq[i]),
    inv = (inv << 2) | Complement(CustomOrdering(seq[l - i - 1]));

  ::std::vector<KMerInfo> ret;

  ret.reserve(l - k + 1);

  bool less = reg <= inv;

  ret.emplace_back(less ? reg : inv, 0, !less);

  unsigned mask = k == 16 ? ~0 : ((1 << (2 * k)) - 1);

  for (auto i = k; i < l; ++i) {
    reg = ((reg << 2) | CustomOrdering(seq[i])) & mask,
    inv = ((inv << 2) | Complement(CustomOrdering(seq[l - i - 1]))) & mask;

    less = reg <= inv;

    ret.emplace_back(less ? reg : inv, i - k + 1, !less);
  }

  return ret;
}

template <typename T>
::std::vector<int> MinByK(const ::std::vector<T>& v, unsigned k) {
  const auto sz = v.size();

  ::std::vector<int> mins(sz, -1);
  ::std::vector<int> s{0};

  for (auto i = 1; i < sz; ++i) {
    while (s.size() && v[s.back()] > v[i])
      mins[s.back()] = i - 1, s.pop_back();

    s.push_back(i);
  }

  while (s.size())
    mins[s.back()] = sz - 1, s.pop_back();

  ::std::vector<int> ret;

  s.shrink_to_fit();
  ret.reserve(sz - k + 1);

  auto j = 0;
  for (auto i = 0; i <= sz - k; ++i) {
    while (j < i || mins[j] < i + k - 1)
      ++j;

    ret.push_back(j);
  }

  return ret;
}

}  // namespace detail

::std::vector<KMerInfo> minimizers(const char* sequence,
                                   SequenceLength sequence_length, KType k,
                                   WindowLength window_length) {
  auto slen = sequence_length.get();
  auto wlen = window_length.get();
  auto klen = k.get();

  ::std::vector<KMerInfo> kmers(
      std::move(detail::GenerateFrom(sequence, slen, klen)));

  decltype(kmers) min_kmers;

  min_kmers.push_back(kmers.front());

  for (auto i = 1; i < wlen - 1; ++i)
    if (kmers[i] < min_kmers.back())
      min_kmers.push_back(kmers[i]);

  for (auto&& idx : detail::MinByK(kmers, wlen))
    min_kmers.push_back(kmers[idx]);

  auto old_end = min_kmers.size();
  min_kmers.push_back(kmers.back());

  for (auto i = 1; i < wlen - 1; ++i)
    if (kmers[kmers.size() - i - 1] < min_kmers.back())
      min_kmers.push_back(kmers[kmers.size() - i - 1]);

  auto cutoff = ::std::begin(min_kmers);
  ::std::advance(cutoff, old_end);

  ::std::reverse(cutoff, ::std::end(min_kmers));

  min_kmers.erase(::std::unique(::std::begin(min_kmers), ::std::end(min_kmers)),
                  ::std::end(min_kmers));

  return min_kmers;
}

}  // namespace blue