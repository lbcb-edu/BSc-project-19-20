#include <getopt.h>
#include <unistd.h>

#include <alignment/alignment.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread_pool/thread_pool.hpp>
#include <unordered_set>
#include <vector>

#include "config.hpp"
#include "mapper.hpp"
#include "matcher.hpp"
#include "util.hpp"

int main(int argc, char** argv, char** env) {
  /*
   *
   * Parsing arguments.
   *
   */

  int opt;
  constexpr ::option options[] = {
      {"help", no_argument, nullptr, 'h'},
      {"version", no_argument, nullptr, 'v'},
      {"alignment", required_argument, nullptr, 'a'},
      {"kmerlen", required_argument, 0, 'k'},
      {"windowlen", required_argument, 0, 'w'},
      {"freqignored", required_argument, 0, 'f'},
      {"threads", required_argument, 0, 't'},
      {"cigar", no_argument, 0, 'c'}};

  ::blue::AlignmentType algorithm;

  unsigned k = 15;
  unsigned w = 5;
  auto f = 0.001;

  auto c = false;
  auto t = 1;

  ::std::string src, dest;

  while ((opt = getopt_long(argc, argv, "s:r:hva:k:w:f:t:c", options,
                            nullptr)) != -1)
    switch (opt) {
      case 'h':
        ::util::DisplayHelp(::std::cerr) << ::util::Terminator{EXIT_SUCCESS};

      case 'v':
        ::std::cerr << ::util::Version() << ::util::Terminator{EXIT_SUCCESS};

      case 's':
        src = optarg;
        break;

      case 'r':
        dest = optarg;
        break;

      case 'a': {
        ::std::string_view av = optarg;
        algorithm = av == "nw"
                        ? ::blue::AlignmentType::kNeedlemanWunsch
                        : av == "sw" ? ::blue::AlignmentType::kSmithWaterman
                                     : ::blue::AlignmentType::kOverlap;
        if (algorithm == ::blue::AlignmentType::kOverlap && av != "ov")
          ::std::cerr << "Invalid alignment algorithm."
                      << ::util::Terminator{EXIT_FAILURE};
        break;
      }

      case 'k':
        try {
          k = ::std::stoi(optarg);
          if (k > sizeof(unsigned) / 2 || k <= 0)
            ::std::cerr << "KMer length must be in range [1, 16]."
                        << ::util::Terminator{EXIT_FAILURE};
        } catch (::std::invalid_argument& exc) {
          ::std::cerr
              << "KMer length should be a positive integer less than 17."
              << ::util::Terminator{EXIT_FAILURE};
        }
        break;

      case 'w':
        try {
          w = ::std::stoi(optarg);
        } catch (::std::invalid_argument& exc) {
          ::std::cerr << "Window length should be a positive integer."
                      << ::util::Terminator{EXIT_FAILURE};
        }
        break;

      case 'f':
        try {
          f = ::std::stod(optarg);
          if (f < 0.0 || f > 1.0)
            ::std::cerr
                << "Ignored percentage of most occuring minimizers can not be "
                   "negative or greater than 1."
                << ::util::Terminator{EXIT_FAILURE};
        } catch (::std::invalid_argument& exc) {
          ::std::cerr << "Ignored percentage of most occuring minimizers must "
                         "be a real number."
                      << ::util::Terminator{EXIT_FAILURE};
        }
        break;

      case 't':
        try {
          t = ::std::stoi(optarg);
          if (t < 0)
            ::std::cerr << "Number of threads must be a positive integer."
                        << ::util::Terminator{EXIT_FAILURE};
        } catch (::std::invalid_argument& exc) {
          ::std::cerr << "Number of threads must be a positive integer."
                      << ::util::Terminator{EXIT_FAILURE};
        }
        break;

      case 'c':
        c = true;
        break;

      default:
        ::std::cerr << "Unkown command line option. Use -h to display possible "
                       "arguments."
                    << ::util::Terminator{EXIT_FAILURE};
    }

  /*
   *
   * Checking input file format.
   *
   */

  if (!src.length() || !dest.length())
    ::std::cerr << "Fragment/reference genome file missing. See --help."
                << ::util::Terminator{EXIT_FAILURE};

  ::std::string arg1;

  // clang-format off
  const ::std::unordered_set<::std::string>
    fasta_ext{".fasta", ".fasta.gz", ".fa", ".fa.gz"},
    fastq_ext{".fastq", ".fastq.gz", ".fq", ".fq.gz"};
  // clang-format on

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

  /*
   *
   * Sequence loading from input files.
   *
   */

  ::std::cerr << "Starting with sequence loading."
              << "\n";

  auto fragments = using_fasta ? ::mapper::Parse<::bioparser::FastaParser>(src)
                               : ::mapper::Parse<::bioparser::FastqParser>(src);

  ::std::cerr << "Loaded fragments."
              << "\n";

  auto reference =
      ::std::move(::mapper::Parse<::bioparser::FastaParser>(dest).front());

  ::std::cerr << "Loaded reference."
              << "\n";

  /*
   *
   * Minimizer indexing.
   *
   */

  ::std::cerr << "Creating minimizer indexes and finding matches."
              << "\n";

  auto pool = ::thread_pool::createThreadPool(t);

  auto reference_index = pool->submit(::matcher::CreateReferenceIndex,
                                      ::std::ref(reference), k, w, f)
                             .get();

  ::std::vector<::std::future<::std::string>> pafs_f(fragments.size());

  for (int i = 0; i < fragments.size(); ++i)
    pafs_f[i] = pool->submit(
        [&](auto&& fragment) {
          return ::matcher::Align(
              *reference, *fragment,
              ::matcher::BestMatch(
                  reference_index,
                  ::matcher::CreateFragmentIndex(fragment, k, w)),
              c, k, algorithm);
        },
        ::std::move(fragments[i]));

  ::std::ios_base::sync_with_stdio(false);

  ::std::cerr << "Printing matches to stdout."
              << "\n\n";

  for (auto& pf : pafs_f)
    ::std::cout << pf.get() << "\n";
}