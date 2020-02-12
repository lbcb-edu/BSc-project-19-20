#pragma once

#include <minimizers/minimizers.hpp>
#include <random>
#include <sstream>
#include <string>
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
      << "Required arguments (if not calling help or version): -s [fragments] "
         "-r [reference]"
      << "\n"
      << "Optional arguments (must come after required arguments):"
      << " -h/--help: display this help"
      << "\n"
      << " -v/--version: display mapper's current version"
      << "\n"
      << " -a/--alignment: takes an argument that specifies the alignment "
         "algorithm; nw - Needleman-Wunsch, sw - Smith-Waterman, ov - overlap"
      << "\n"
      << " -k/--kmerlen: specify the length of the kmers used for matching"
      << "\n"
      << " -w/--windowlen: specify the length of window used in kmer generation"
      << "\n"
      << " -f/--freqignored: specify the percentage of ignored kmers amongst "
         "most occuring"
      << "\n"
      << " -t/--threads: specify the number of threads used for finding "
         "alignment regions"
      << "\n"
      << " -c/--cigar: optional argument, provide 'c' if you want cigar to be "
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
