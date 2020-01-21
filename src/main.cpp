#include <algorithm>
#include <alignment/alignment.hpp>
#include <bioparser/bioparser.hpp>
#include <cassert>
#include <cctype>
#include <iostream>
#include <limits>
#include <minimizers/minimizers.hpp>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <thread_pool/thread_pool.hpp>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "config.hpp"
#include "mapper.hpp"

namespace util {

::std::string ToLower(::std::string s) noexcept {
  for (auto&& c : s)
    c = ::std::tolower(c);
  return s;
}

bool EndsWith(const ::std::string& str, const ::std::string& ext) noexcept {
  return str.find(ext) == str.length() - ext.length();
}

::std::ostream& DisplayHelp(::std::ostream& out) noexcept {
  out << "Mapper that displays input genome statistics."
      << "\n"
      << "Arguments: [fragments] [reference] [alignmentType] [match] "
         "[mismatch] [gap] [threads] <cigar>"
      << "\n"
      << " - fragments: input file in a FASTA or FASTQ format"
      << "\n"
      << " - reference: input file in a FASTA format"
      << "\n"
      << " - alignmentType: the alignment algorithm used for aligning two "
         "random "
         "sequences from the first file options: nw (Needleman-Wunsch), sw "
         "(Smith-Waterman) or ov (overlap algorithm)"
      << "\n"
      << " - match: the cost of matching bases"
      << "\n"
      << " - mismatch: the cost of mismatching bases"
      << "\n"
      << " - gap: the cost of a gap in either source or target sequence"
      << "\n"
      << " - threads: number of threads to be used"
      << "\n"
      << " - cigar: optional argument, provide 'c' if you want cigar to be "
         "printed"
      << "\n";

  return out;
}

::std::string Version() noexcept {
  ::std::ostringstream ss;
  ss << BLUE_MAPPER_VERSION_MAJOR << "." << BLUE_MAPPER_VERSION_MINOR << "."
     << BLUE_MAPPER_VERSION_PATCH;

  return ss.str();
}

::std::ostream& operator<<(
    ::std::ostream& out,
    const ::mapper::Sequences<::mapper::Sequence>& sequences) noexcept {
  ::std::size_t min = ::std::numeric_limits<::std::size_t>::max(), max = 0,
                total = 0, temp;

  for (auto&& seq : sequences) {
    min = ::std::min(min, temp = seq->sequence.size());
    max = ::std::max(max, temp);
    total += temp;
  }

  auto avg = total / (double)sequences.size();

  out << " - number of sequences: " << sequences.size() << "\n"
      << " - average sequence length: " << avg << "\n"
      << " - max sequence length: " << max << "\n"
      << " - min sequence length: " << min;

  return out;
}

struct Terminator {
  int status;
};

::std::ostream& operator<<(::std::ostream& out, Terminator t) {
  out << "\n";
  ::std::exit(t.status);
}

int ParseInteger(const char* n) {
  try {
    return ::std::stoi(n);
  } catch (::std::invalid_argument& exc) {
    ::std::cerr << "Error while attempting to parse " << n << " as a number."
                << Terminator{EXIT_FAILURE};
  }

  return {};
}

template <typename IntegralType>
IntegralType Rnd(const IntegralType bound) {
  static ::std::mt19937 mt{::std::random_device{}()};
  static ::std::uniform_int_distribution<IntegralType> uniform{0, bound - 1};

  return uniform(mt);
}

using ReferenceIndex =
    ::std::unordered_map<unsigned, ::std::vector<::blue::KMerInfo>>;
using FragmentIndex = ::std::vector<::blue::KMerInfo>;
using Sequence = ::std::tuple<int, int, int>;

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
    return ::std::get<1>(a.first) < ::std::get<1>(b.first);
  });

  ReferenceIndex ret;
  for (auto& [kmer, _] : index)
    ret[::std::get<0>(kmer)].emplace_back(::std::move(kmer));

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

template <typename T, typename C>
Sequence LIS(const ::std::vector<T>& v, const C& comp) {
  ::std::vector<Sequence> sequences{{0, 0, 1}};
  for (int i = 1; i < v.size(); ++i)
    if (comp(v[i], v[::std::get<1>(sequences.front())])) {
      sequences[0] = {i, i, 1};
    } else if (auto [p, q, l] = sequences.back(); comp(v[q], v[i])) {
      sequences.emplace_back(p, i, l + 1);
    } else {
      auto iter =
          ::std::lower_bound(sequences.begin(), sequences.end(), v[i],
                             [&v, &comp](auto const& tuple, auto const& value) {
                               return comp(v[::std::get<1>(tuple)], value);
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

using KMerMatch = ::std::pair<::blue::KMerInfo, ::blue::KMerInfo>;

namespace detail {

bool r(const KMerMatch& a, const KMerMatch& b) noexcept {
  return ::std::get<1>(a.second) < ::std::get<1>(b.second);
}

bool i(const KMerMatch& a, const KMerMatch& b) noexcept {
  return ::std::get<1>(a.second) > ::std::get<1>(b.second);
}

bool rr(const KMerMatch& a, const KMerMatch& b) noexcept {
  return ::std::get<1>(a.first) < ::std::get<1>(b.first) &&
         ::std::get<1>(a.second) < ::std::get<1>(b.second);
}

bool ri(const KMerMatch& a, const KMerMatch& b) noexcept {
  return ::std::get<1>(a.first) < ::std::get<1>(b.first) &&
         ::std::get<1>(a.second) > ::std::get<1>(b.second);
}

bool ir(const KMerMatch& a, const KMerMatch& b) noexcept {
  return ::std::get<1>(a.first) > ::std::get<1>(b.first) &&
         ::std::get<1>(a.second) < ::std::get<1>(b.second);
}

bool ii(const KMerMatch& a, const KMerMatch& b) noexcept {
  return ::std::get<1>(a.first) > ::std::get<1>(b.first) &&
         ::std::get<1>(a.second) > ::std::get<1>(b.second);
}

}  // namespace detail

using Region = ::std::pair<KMerMatch, KMerMatch>;

Region BestMatch(const ReferenceIndex& ri, FragmentIndex fi) {
  ::std::vector<::std::vector<KMerMatch>> matches(2);
  constexpr bool (*cmp[])(const KMerMatch&, const KMerMatch&) = {detail::r,
                                                                 detail::i};

  for (auto f_kmer : fi) {
    auto iter = ri.find(::std::get<0>(f_kmer));
    auto idx = ::std::get<2>(f_kmer) & 1;
    if (iter != ri.end())
      for (auto r_kmer : iter->second)
        matches[idx ^ ::std::get<2>(r_kmer)].emplace_back(f_kmer, r_kmer);
  }

  Sequence best_seq = {0, 0, 0};
  auto type = -1;

  for (auto i = 0; i < matches.size(); ++i)
    if (matches[i].size()) {
      auto seq = LIS(matches[i], cmp[i]);
      if (::std::get<2>(seq) > ::std::get<2>(best_seq))
        best_seq = seq, type = i;
    }

  if (type == -1)
    return {};

  return {matches[type][::std::get<0>(best_seq)],
          matches[type][::std::get<1>(best_seq)]};
}

}  // namespace util

::std::ostream& operator<<(::std::ostream& out, const ::blue::KMerInfo k) {
  out << ::std::get<0>(k) << " " << ::std::get<1>(k) << " " << ::std::get<2>(k);
  return out;
}

int main(int argc, char** argv, char** env) {
  /*
   *
   * Argument validation.
   *
   */

  ::std::string arg1;

  if (argc == 2 &&
      ("-v" == (arg1 = ::util::ToLower(argv[1])) || "--version" == arg1))
    ::std::cerr << ::util::Version() << ::util::Terminator{EXIT_SUCCESS};
  else if ((argc != 3 && argc < 7) || "-h" == arg1 || "--help" == arg1)
    ::util::DisplayHelp(::std::cerr) << ::util::Terminator{EXIT_SUCCESS};

  // clang-format off
  const ::std::unordered_set<::std::string>
    fasta_ext{".fasta", ".fasta.gz", ".fa", ".fa.gz"},
    fastq_ext{".fastq", ".fastq.gz", ".fq", ".fq.gz"};
  // clang-format on

  const ::std::string src{argv[1]}, dest{argv[2]};

  int using_fasta = -1;

  for (auto&& ext : fasta_ext)
    if (::util::EndsWith(src, ext)) {
      using_fasta = 1;
      break;
    }

  if (-1 == using_fasta) {
    for (auto&& ext : fastq_ext)
      if (::util::EndsWith(src, ext)) {
        using_fasta = 0;
        break;
      }
  }

  if (-1 == using_fasta)
    ::std::cerr << "Source file should be in FASTA or FASTQ format."
                << ::util::Terminator{EXIT_FAILURE};

  bool destination_valid = false;

  for (auto&& ext : fasta_ext)
    if (::util::EndsWith(dest, ext)) {
      destination_valid = true;
      break;
    }

  if (!destination_valid)
    ::std::cerr << "Destination file should be in FASTA format."
                << ::util::Terminator{EXIT_FAILURE};

  ::blue::AlignmentType algorithm;

  if (const ::std::string a_type{::util::ToLower(argv[3])}; a_type == "nw")
    algorithm = ::blue::AlignmentType::kNeedlemanWunsch;
  else if (a_type == "sw")
    algorithm = ::blue::AlignmentType::kSmithWaterman;
  else if (a_type == "ov")
    algorithm = ::blue::AlignmentType::kOverlap;
  else
    ::std::cerr << "Invalid alignment algorithm."
                << ::util::Terminator{EXIT_FAILURE};

  int match_cost = ::util::ParseInteger(argv[4]),
      mismatch_cost = ::util::ParseInteger(argv[5]),
      gap_cost = ::util::ParseInteger(argv[6]);

  /*
   *
   * Sequence loading from input files and statistic output.
   *
   */

  ::std::cerr << "Starting with sequence loading."
              << "\n";

  auto fragments = using_fasta ? ::mapper::Parse<::mapper::Sequence>(src)
                               : ::mapper::Parse<::mapper::QSequence>(src);

  ::std::cerr << "Loaded fragments."
              << "\n";

  auto reference =
      ::std::move(::mapper::Parse<::mapper::Sequence>(dest).front());

  ::std::cerr << "Loaded reference."
              << "\n\n";

  unsigned k = argc >= 8 ? ::std::stoi(argv[7]) : 15;
  unsigned w = argc >= 9 ? ::std::stoi(argv[8]) : 5;
  auto f = argc >= 10 ? ::std::stod(argv[9]) : 0.001;

  auto c = argc >= 12 && argv[11][0] == 'c';
  auto t = ::std::stoi(argv[10]);

  ::std::cerr << "Creating minimizer index..."
              << "\n";

  auto pool = ::thread_pool::createThreadPool(t);

  auto reference_index_f = pool->submit(::util::CreateReferenceIndex,
                                        ::std::ref(reference), k, w, f);

  ::std::vector<::std::future<::std::vector<::blue::KMerInfo>>>
      fragment_indices_f;
  fragment_indices_f.reserve(fragments.size());

  ::std::transform(
      fragments.begin(), fragments.end(),
      ::std::back_inserter(fragment_indices_f), [&](const auto& f) {
        return pool->submit(::util::CreateFragmentIndex, ::std::ref(f), k, w);
      });

  const auto reference_index = reference_index_f.get();

  ::std::vector<::std::future<::util::Region>> regions_f;
  regions_f.reserve(fragment_indices_f.size());

  ::std::transform(fragment_indices_f.begin(), fragment_indices_f.end(),
                   ::std::back_inserter(regions_f),
                   [&](auto& f) -> ::std::future<::util::Region> {
                     return pool->submit(::util::BestMatch,
                                         ::std::ref(reference_index), f.get());
                   });

  ::std::vector<::util::Region> regions;
  regions.reserve(regions_f.size());

  ::std::transform(regions_f.begin(), regions_f.end(),
                   ::std::back_inserter(regions),
                   [](auto& f) { return f.get(); });

  for (auto [p, q] : regions) {
    auto [pf, pr] = p;
    auto [qf, qr] = q;

    ::std::cout << pf << "; " << pr << "\n";
    ::std::cout << qf << ", " << qr << "\n";
    ::std::cout << "\n";
  }
}