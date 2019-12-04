#include "minimizers/minimizers.hpp"

namespace blue {

namespace detail {

inline char Complement(const char c) noexcept {
  return '0' + '3' - c;
}

inline char CustomOrdering(const char c) noexcept {
  switch (c) {
    case 'C':
      return '0';
    case 'A':
      return '1';
    case 'T':
      return '2';
    case 'G':
      return '3';
    default:
      return '?';
  }
}

template <typename StringType>
inline ::std::string_view ViewFrom(const StringType& str, ::std::size_t pos,
                                   ::std::size_t count) {
  return ::std::string_view{str.data() + pos, count};
}

template <typename StringType>
inline KMerInfo KMerFrom(const StringType& win, const StringType& cwin,
                         unsigned k) {
  ::std::string_view candidate{"Z"};
  KMerInfo ret{0, k, false};

  for (int i = 0; i <= win.size() - k; ++i) {
    ::std::string_view kk = ViewFrom(win, i, k);
    ::std::string_view ck = ViewFrom(cwin, cwin.size() - k - i, k);

    bool complement = ck < kk;

    if (complement)
      kk = ck;

    if (kk < candidate) {
      candidate = kk;
      ret = {i, k, complement};
    }
  }

  return ret;
}

inline KMerInfo AdjustOrigin(KMerInfo ki, unsigned pos) {
  ::std::get<0>(ki) += pos - (::std::get<2>(ki) ? ::std::get<1>(ki) - 1 : 0);
  return ki;
}

}  // namespace detail

::std::vector<KMerInfo> minimizers(const char* sequence,
                                   SequenceLength sequence_length, KType k,
                                   WindowLength window_length) {
  auto len = sequence_length.get();
  auto kk = k.get();
  auto wlen = window_length.get();
  auto l = wlen + kk - 1;

  ::std::string seq(len, ' ');
  ::std::string com(len, ' ');

  ::std::transform(sequence, sequence + len, ::std::begin(seq),
                   detail::CustomOrdering);
  ::std::transform(::std::rbegin(seq), ::std::rend(seq), ::std::begin(com),
                   detail::Complement);

  ::std::vector<KMerInfo> ret;
  ret.reserve(seq.size() - window_length.get());

  for (int i = 0; i <= len - l; ++i)
    ret.push_back(
        detail::AdjustOrigin(detail::KMerFrom(detail::ViewFrom(seq, i, l),
                                              detail::ViewFrom(com, i, l), kk),
                             i));

  for (int u = 1; u < wlen; ++u) {
    auto ul = u + kk - 1;

    ret.push_back(detail::KMerFrom(detail::ViewFrom(seq, 0, ul),
                                   detail::ViewFrom(com, 0, ul), kk));

    auto origin = len - ul;

    ret.push_back(detail::AdjustOrigin(
        detail::KMerFrom(detail::ViewFrom(seq, origin, ul),
                         detail::ViewFrom(com, origin, ul), kk),
        origin));
  }

  return ret;
}

}  // namespace blue