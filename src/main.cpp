#include <bioparser/bioparser.hpp>
#include <cctype>
#include <future>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "config.hpp"
#include "mapper.hpp"

namespace util {

::std::string ToLower(::std::string s) noexcept {
  for (auto&& c : s) c = ::std::tolower(c);
  return s;
}

::std::ostream& DisplayHelp(::std::ostream& out) noexcept {
  out << "Mapper that displays input genome statistics." << ::std::endl
      << "Arguments: [fragments] [reference]." << ::std::endl
      << " - fragments: input file in a FASTA or FASTQ format" << ::std::endl
      << " - reference: input file in a FASTA format";

  return out;
}

::std::string Version() noexcept {
  ::std::ostringstream ss;
  ss << BLUE_MAPPER_VERSION_MAJOR << "." << BLUE_MAPPER_VERSION_MINOR << "."
     << BLUE_MAPPER_VERSION_PATCH;

  return ss.str();
}

bool EndsWith(const ::std::string& str, const ::std::string& ext) noexcept {
  return str.find(ext) == str.length() - ext.length();
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

  out << " - number of sequences: " << sequences.size() << ::std::endl
      << " - average sequence length: " << avg << ::std::endl
      << " - max sequence length: " << max << ::std::endl
      << " - min sequence length: " << min;

  return out;
}

}  // namespace util

int main(int argc, char** argv, char** env) {
  ::std::string arg1;

  if (argc == 1 || "-h" == (arg1 = ::util::ToLower(argv[1])) ||
      "--help" == arg1) {
    ::util::DisplayHelp(::std::cerr) << ::std::endl;
    ::std::exit(EXIT_SUCCESS);
  } else if ("-v" == arg1 || "--Version" == arg1) {
    ::std::cerr << ::util::Version() << ::std::endl;
    ::std::exit(EXIT_SUCCESS);
  } else if (argc != 3) {
    ::std::cerr << "This program requires 2 arguments. See --help."
                << ::std::endl;
    ::std::exit(EXIT_FAILURE);
  }

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

  if (-1 == using_fasta) {
    ::std::cerr << "Source file should be in FASTA or FASTQ format."
                << ::std::endl;
    ::std::exit(EXIT_FAILURE);
  }

  bool destination_valid = false;

  for (auto&& ext : fasta_ext)
    if (::util::EndsWith(dest, ext)) {
      destination_valid = true;
      break;
    }

  if (!destination_valid) {
    ::std::cerr << "Destination file should be in FASTA format." << ::std::endl;
    ::std::exit(EXIT_FAILURE);
  }

  auto fragments = using_fasta ? ::mapper::Parse<::mapper::Sequence>(src)
                               : ::mapper::Parse<::mapper::QSequence>(src);

  ::std::cerr << "Loaded fragments." << ::std::endl;

  auto reference = ::mapper::Parse<::mapper::Sequence>(dest);

  ::std::cerr << "Loaded reference." << ::std::endl;

  using ::util::operator<<;

  ::std::cerr << "Fragment statistics: " << ::std::endl
              << fragments << ::std::endl
              << "Reference statistics: " << ::std::endl
              << reference << ::std::endl;
}