#include <bioparser/bioparser.hpp>
#include <cctype>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

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

using namespace ::mapper;

::std::ostream& operator<<(::std::ostream& out,
                           const Sequences<Sequence>& sequences) noexcept {
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
  else if ((argc != 3 && argc != 7) || "-h" == arg1 || "--help" == arg1)
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

  using ::util::operator<<;

  ::std::cerr << "Fragment statistics: "
              << "\n"
              << fragments << "\n"
              << "Reference statistics: "
              << "\n"
              << reference << "\n";

  /*
   *
   * Sequence alignment.
   *
   */

  auto source = *fragments[::util::Rnd(fragments.size())];
  auto target = *fragments[::util::Rnd(fragments.size())];

  auto source_length = static_cast<unsigned>(source.sequence.size());
  auto target_length = static_cast<unsigned>(target.sequence.size());

  reference.clear();
  reference.shrink_to_fit();

  ::std::cerr << "\n"
              << "Aligning two random sequences, query length: "
              << source_length << ", target length: " << target_length << "\n";

  // there is much more memory that the program owns but does not use
  // it is leftover from storing sequences
  ::std::cerr << "Estimated (active) memory usage: "
              << (8 * target_length * source_length) / (1 << 20)
              << " megabytes."
              << "\n";

  ::std::string cigar;
  unsigned target_begin;

  using namespace ::blue;

  auto score = PairwiseAlignment(
      Query{source.sequence.c_str()}, QueryLength{source_length},
      Target{target.sequence.c_str()}, TargetLength{target_length}, algorithm,
      Match{match_cost}, Mismatch{mismatch_cost}, Gap{gap_cost}, cigar,
      target_begin);

  ::std::cerr << "Done."
              << "\n";
  ::std::cerr << "Alignment score: " << score << "\n";

  ::std::cerr << "Cigar length is " << cigar.size() << "."
              << " Do you want to print it (y/n)? ";

  ::std::string response;
  ::std::cin >> response;

  if (::util::ToLower(response) == "y")
    ::std::cerr << cigar << "\n";

  // return 0;

  ::std::cerr << "Calculating kMers..."
              << "\n";

  ::std::unordered_map<unsigned, unsigned> occurences;
  unsigned counter{0};

  // TODO: add command line arguments and other statistics

  ::std::ios_base::sync_with_stdio(false);
  ::std::cin.tie(nullptr);

  for (auto&& frag : fragments) {
    for (auto&& [kmer, u, v] :
         minimizers(frag->sequence.data(),
                    blue::SequenceLength{
                        static_cast<unsigned>(frag->sequence.length())},
                    blue::KType{15}, blue::WindowLength{5}))
      ++occurences[kmer];
  }

  ::std::cout << "\nNumber of distinct kmers: " << occurences.size() << ".\n";
}