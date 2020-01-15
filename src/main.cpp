#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <bioparser/bioparser.hpp>

#include <alignment/alignment.hpp>
#include <minimizers/minimizers.hpp>

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
         "[mismatch] [gap]."
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

}  // namespace util

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

  auto reference = ::mapper::Parse<::mapper::Sequence>(dest);

  ::std::cerr << "Loaded reference."
              << "\n\n";

  ::std::cerr << "TODO..."
              << "\n";

  auto k = argc >= 8 ? ::std::stoi(argv[7]) : 15;
  auto w = argc >= 9 ? ::std::stoi(argv[8]) : 5;
  auto f = argc >= 10 ? ::std::stod(argv[9]) : 0.001;

  ::std::vector<::std::vector<::blue::KMerInfo>> kmers;
  kmers.reserve(fragments.size());

  for (auto&& frag : fragments)
    kmers.emplace_back(minimizers(
        frag->sequence.data(),
        blue::SequenceLength{static_cast<unsigned>(frag->sequence.length())},
        blue::KType{k}, blue::WindowLength{w}));

  ::std::unordered_map<::blue::KMerInfo, unsigned> occurences;

  for (auto& kmer : minimizers(reference.front()->sequence.data(),
                               blue::SequenceLength{static_cast<unsigned>(
                                   reference.front()->sequence.length())},
                               blue::KType{k}, blue::WindowLength{w}))
    ++occurences[kmer];

  ::std::vector<::std::pair<::blue::KMerInfo, unsigned>> reference_index(
      ::std::make_move_iterator(occurences.begin()),
      ::std::make_move_iterator(occurences.end()));

  occurences.clear();

  ::std::sort(reference_index.begin(), reference_index.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

  ::std::cerr << "Size before: " << reference_index.size() << "\n";
  reference_index.resize((1.0 - f) * reference_index.size());
  ::std::cerr << "Size after:  " << reference_index.size() << "\n";
}