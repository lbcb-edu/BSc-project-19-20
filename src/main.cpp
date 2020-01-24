#include <alignment/alignment.hpp>
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

  auto fragments = using_fasta ? ::mapper::Parse<::bioparser::FastaParser>(src)
                               : ::mapper::Parse<::bioparser::FastqParser>(src);

  ::std::cerr << "Loaded fragments."
              << "\n";

  auto reference =
      ::std::move(::mapper::Parse<::bioparser::FastaParser>(dest).front());

  ::std::cerr << "Loaded reference."
              << "\n\n";

  /*
   *
   * Minimizer indexing.
   *
   */

  unsigned k = argc >= 8 ? ::util::ParseInteger(argv[7]) : 15;
  unsigned w = argc >= 9 ? ::util::ParseInteger(argv[8]) : 5;
  auto f = argc >= 10 ? ::util::ParseInteger(argv[9]) : 0.001;

  auto c = argc >= 12 && argv[11][0] == 'c';
  auto t = ::util::ParseInteger(argv[10]);

  ::std::cerr << "Creating minimizer index..."
              << "\n";

  auto pool = ::thread_pool::createThreadPool(t);

  auto reference_index = pool->submit(::matcher::CreateReferenceIndex,
                                      ::std::ref(reference), k, w, f)
                             .get();

  ::std::vector<::std::future<::std::vector<::blue::KMerInfo>>>
      fragment_indices_f;
  fragment_indices_f.reserve(fragments.size());

  ::std::vector<::std::future<::matcher::Region>> regions_f;
  regions_f.reserve(fragments.size());

  ::std::vector<::std::future<::std::string>> pafs_f;
  pafs_f.reserve(fragments.size());

  for (int i = 0; i < fragments.size(); i += t) {
    ::std::vector<::std::vector<::blue::KMerInfo>> f_index(t);
    for (int j = 0; j < t && i + j < fragments.size(); ++j)
      fragment_indices_f.push_back(pool->submit(
          ::matcher::CreateFragmentIndex, ::std::ref(fragments[i + j]), k, w));

    for (int j = 0; j < t && i + j < fragments.size(); ++j)
      regions_f.push_back(pool->submit([&, idx = i + j] {
        return ::matcher::BestMatch(reference_index,
                                    fragment_indices_f[idx].get());
      }));

    for (int j = 0; j < t && i + j < fragments.size(); ++j)
      pafs_f.push_back(pool->submit([&, idx = i + j] {
        return ::matcher::Align(*reference, *fragments[idx],
                                regions_f[idx].get(), c, k, algorithm);
      }));
  }

  ::std::ios_base::sync_with_stdio(false);
  ::std::cerr << "Found alignments:"
              << "\n\n";

  for (auto& f : pafs_f)
    if (auto str = f.get(); str.length())
      ::std::cerr << str << "\n";
}