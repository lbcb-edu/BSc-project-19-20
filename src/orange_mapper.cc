#include <unordered_map>
#include <functional>
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <random>

#include <bioparser/bioparser.hpp>

#include "orange_alignment.h"
#include "orange_minimizers.h"
#include "orange_mapper_conf.h"

namespace orange {
namespace mapper {

/**
 * @brief Supported genome file formats
 */
enum FileType { kFasta, kFastq, kUnsupported };

/**
 * Struct array representing long optins
 */
struct option const long_options[] = {{"help", no_argument, NULL, 'h'},
                                      {"version", no_argument, NULL, 'v'},
                                      {"type", required_argument, NULL, 't'},
                                      {"match", required_argument, 0, 'm'},
                                      {"mismatch", required_argument, 0, 's'},
                                      {"gap", required_argument, 0, 'g'},
                                      {"kmer", required_argument, NULL, 'k'},
                                      {"window", required_argument, NULL, 'w'},
                                      {NULL, 0, NULL, 0}};

/**
 * @brief Holds FASTA/FASTQ sequecnes
 */
struct Sequence {
    /**
     * @brief Construct a new FastAQ object
     *
     * @param name sequence name
     * @param name_len sequence name length
     *
     * @param seq geneom sequence
     * @param seq_len genom sequence length
     */
    Sequence(char const* name, std::uint32_t name_len, char const* seq,
             std::uint32_t seq_len)
        : name_{name}, name_len_{name_len}, seq_{seq}, seq_len_{seq_len} {}

    char const* name_;              //<<< sequence name
    std::uint32_t const name_len_;  //<<< sequence name length

    char const* seq_;              //<<< genome sequence
    std::uint32_t const seq_len_;  //<<< genome sequence
};

/**
 * @brief FastA sequence alias.
 */
// Doable since FASTA doesn't contain quality information
using FastA = Sequence;

/**
 * @brief FASTQ genome sequences.
 *
 * @details Extends @ref orange::mapper::Sequence
 *      with string quiality information
 */
struct FastQ : public Sequence {
    FastQ(char const* name, std::uint32_t name_len, char const* seq,
          std::uint32_t seq_len, char const* quality_str,
          std::uint32_t quality_str_len)
        : Sequence(name, name_len, seq, seq_len),
          quality_str_{quality_str},
          quality_str_len_{quality_str_len} {}

    char const* quality_str_;
    std::uint32_t const quality_str_len_;
};

/**
 * @brief Vector of pointers to parsed @ref orange::mapper::Sequence
 */
using VecSeqPtr = std::vector<std::unique_ptr<Sequence>>;

/**
 * @brief Vector of pointers to parsed @ref orange::mapper::FastQ
 *
 */
using VecFastQPtr = std::vector<std::unique_ptr<FastQ>>;

/**
 * @brief Simple intarface for creating @ref bioparser FASTA parser
 */
auto createFastaParser =
    std::function<std::unique_ptr<bioparser::Parser<FastA>>(
        std::string const& path)>{
        &bioparser::createParser<bioparser::FastaParser, FastA>};

/**
 * @brief Simple interface for creating @ref bioparser FASTQ parser
 */
auto createFastqParser =
    std::function<std::unique_ptr<bioparser::Parser<FastQ>>(
        std::string const& path)>{
        &bioparser::createParser<bioparser::FastqParser, FastQ>};

/**
 * @brief Prints program version to stderr,
 *      exits upon completion
 */
auto printVersion() {
    std::cerr << orange_mapper_VERSION_MAJOR << '.'
              << orange_mapper_VERSION_MINOR << '.'
              << orange_mapper_VERSION_PATCH << '\n';
    std::exit(EXIT_SUCCESS);
}

/**
 * @brief Prints program usage information to stderr,
 *      exits upon completion
 */
auto printHelp() {
    std::cerr << "Genome sequence mapper.\n"
              << "Reads are supported in FASTA and FASTQ formats while "
              << "reference is expected to be FASTA.\n\n"
              << "Usage:\n\t orange_mapper <reads> <reference> "
              << "-t <alignment_type> -m <match> -s <mismatch> -g <gap>\n\n"
              << "Options:\n\t-h\thelp\n\t-v\tverions\n";
    std::exit(EXIT_SUCCESS);
}

/**
 * @brief parses alignment type from comman line argument
 *
 * @param type command line argument value
 * @return auto @ref orange::alignment::AlignmentType
 */
auto parseAlignType(char const* type) {
    auto str_val = std::string_view{type};

    if (str_val == "global")
        return alignment::AlignmentType::kGlobal;
    if (str_val == "local")
        return alignment::AlignmentType::kLocal;
    if (str_val == "semi-global")
        return alignment::AlignmentType::kSemiGlobal;

    throw std::invalid_argument(
        "Unknow alignment type.\n"
        "Use orange_mapper -h for further information");
}

/**
 * @brief Parses command line options.
 *
 * @details After option parsing,
 *      all non option elements are placed at the end of the array.
 *
 *      </br>
 *      Associated actions with each option are executed by
 *      this function considering the option was set at launch.
 *
 * @param argc number of command line arguments
 * @param argv commane line arguments
 *
 * @param conf sequence alignemnt config passed over command line options
 *
 * @return int index of first non-option element in argv
 */
auto parseOptions(int argc, char* argv[], alignment::AlignConf& a_conf,
                  minimizers::MinimizerConf& m_conf) {
    auto opt = int{};
    while ((opt = getopt_long(argc, argv, "hvm:s:g:t:", long_options, NULL)) !=
           -1) {
        switch (opt) {
            case 'v':
                printVersion();
                break;
            case 'h':
                printHelp();
                break;
            case 't':
                a_conf.type_ = parseAlignType(optarg);
                break;
            case 'm':
                a_conf.match_ = std::atoi(optarg);
                break;
            case 's':
                a_conf.mismatch_ = std::atoi(optarg);
                break;
            case 'g':
                a_conf.gap_ = std::atoi(optarg);
                break;
            case 'k':
                m_conf.k_ = std::atoi(optarg);
                break;
            case 'w':
                m_conf.win_len_ = std::atoi(optarg);
                break;

            default:
                throw std::invalid_argument(
                    "Uknown option.\n\tUse: 'orange_mapper -h'");
                break;
        }
    }

    return optind;
}

/**
 * @brief Get the @ref orange::mapper::FileType
 *
 * @param file path to file
 * @return corresponding @ref orange::mapper::FileType
 */
auto parseFileType(std::string_view file) {
    auto format_ex = std::invalid_argument(
        "Unsuppored file format/extension!\n"
        "\t Supported: FASTA (reads and reference), FASTQ "
        "(reads)\n\nUse -h for help\n");

    // Find furthest extension
    auto dot = file.rfind('.');
    if (dot == std::string_view::npos)
        throw std::invalid_argument("Missing file extension.");
    else if (file.size() - dot < 2)
        throw format_ex;
    else if (file[dot + 1] == 'g' && file[dot + 2] == 'z') {
        // Ignoring .gz
        file.remove_suffix(file.size() - dot);
        dot = file.rfind('.');
    }

    file.remove_prefix(dot);
    if (file == ".fasta" || file == ".fa")
        return FileType::kFasta;
    else if (file == ".fastq" || file == ".fq")
        return FileType::kFastq;

    throw format_ex;
}

/**
 * @brief Loads FASTA/FASTQ files into
 *
 * @param path_to_file path to file contaning FASTA/FASTQ sequences
 *
 * @return @ref orange::mapper::VecSeqPtr
 */
auto loadFile(std::string const& path_to_file, FileType const& type) {
    auto parse = [](auto parser, auto& objects) { parser->parse(objects, -1); };

    if (type == FileType::kFasta) {
        auto objects = VecSeqPtr{};
        parse(createFastaParser(path_to_file), objects);

        return objects;
    } else {
        auto objects = VecFastQPtr{};
        parse(createFastqParser(path_to_file), objects);

        return VecSeqPtr{std::make_move_iterator(objects.begin()),
                         std::make_move_iterator(objects.end())};
    }
}

/**
 * @brief Prints stats for loaded sequences
 *
 * @details Calculated stats: <br\>
 *          <ol>
 *              <li>Number of sequences</li>
 *              <li>Average sequence length</li>
 *              <li>Minimum sequence length</li>
 *              <li>Maximum sequence length</li>
 *          <\ol>
 *
 * @param vec_seqs reference to loaded set of sequences
 */
auto printStats(std::string_view const& origin, VecSeqPtr const& vec_seq) {
    auto n_seq = vec_seq.size();

    auto avg_len = double{};
    auto min_len = std::uint32_t{std::numeric_limits<std::uint32_t>::max()};
    auto max_len = std::uint32_t{std::numeric_limits<std::uint32_t>::min()};

    std::for_each(begin(vec_seq), end(vec_seq),
                  [&avg_len, &min_len, &max_len](auto& seq_ptr) {
                      auto const& seq_len = seq_ptr->seq_len_;

                      avg_len += seq_len;
                      min_len = std::min(min_len, seq_len);
                      max_len = std::max(max_len, seq_len);
                  });

    avg_len /= n_seq;

    std::cerr << '\n'
              << origin << " Stats:\n"
              << "\tNumber or sequences: " << n_seq << '\n'
              << "\tAverage length: " << avg_len << '\n'
              << "\tMinimum length: " << min_len << '\n'
              << "\tMaximum length: " << max_len << '\n';
}

/**
 * @brief
 *
 * @param reads
 * @param ref
 * @param align_type
 * @param match
 * @param mismatch
 * @param gap
 * @return auto
 */
auto printRngAlign(Sequence const& query, Sequence const& target,
                   alignment::AlignConf conf) {
    auto cigar = std::string{};
    auto target_begin = std::uint32_t{};

    auto const align_score = alignment::pairwiseAlignment(
        query.seq_, query.seq_len_, target.seq_, target.seq_len_, conf.type_,
        conf.match_, conf.mismatch_, conf.gap_, cigar, target_begin);

    std::cerr << "Random alignment score: " << align_score << '\n'
              << "CIGAR:\n\t" << cigar << "\n\n";
}

auto printMinimizerStats(VecSeqPtr const& reads,
                         minimizers::MinimizerConf const& conf) {
    auto minims = minimizers::KMerValUMap<std::uint64_t>{};
    using KMerCnt = std::pair<minimizers::KMer, std::uint64_t>;

    for (auto const& it : reads) {
        auto& read = *it.get();
        for (auto const& minim : minimizers::minimizers(
                 read.seq_, read.seq_len_, conf.k_, conf.win_len_)) {
            ++minims[std::get<0>(minim)];
        }
    }

    std::cerr << "Number of distincs minimizers for reads: " << minims.size()
              << '\n';

    auto vec = std::vector<KMerCnt>{minims.begin(), minims.end()};
    std::sort(vec.begin(), vec.end(),
              [](auto const& l, auto const& r) { return l.second > r.second; });

    auto n_singletons_lambda = [&vec]() -> std::uint64_t {
        auto it = std::lower_bound(
            vec.begin(), vec.end(),
            std::make_pair(std::uint64_t{0}, std::uint64_t{1}),
            [](auto const& l, auto const& r) { return l.second > r.second; });

        if (it == vec.end())
            throw std::range_error("No singletons!");
        else
            return static_cast<std::uint64_t>(vec.end() - it);
    };

    auto most_freq = vec[static_cast<std::size_t>(conf.f_ * vec.size())].second;
    auto n_singletons = n_singletons_lambda();

    std::cerr << "Fraction: " << 1.0 * n_singletons / most_freq << '\n';
}

}  // namespace mapper
}  // namespace orange

/**
 * @brief Main program entry point
 *
 * @param argc number of command line arguments
 * @param argv command line arguments (options and/or sources)
 *
 * @return int execution success status
 */
int main(int argc, char* argv[]) {
    using namespace orange;

    try {
        // Alignment configuration
        auto a_conf = alignment::AlignConf{};
        auto m_conf = minimizers::MinimizerConf{};

        auto arg_index = mapper::parseOptions(argc, argv, a_conf, m_conf);

        // Check number of passed arguments
        if (argc - arg_index < 2) {
            std::cerr << "Invalid number of arguments.\n"
                         "\tRequired two input files:\n"
                         "\t\t1. Set of fragments (FASTA/FASTQ)\n"
                         "\t\t2. Reference genome (FASTA)\n\n"
                         "Followed up with alignment algorithm specification "
                         "{global, local, semi-global}\n"
                         "and match, mismatch, gap scoring values for the "
                         "coresponding type.\n"
                         "\tEg. orange_mapper escherichia_coli_r7_reads.fastq"
                         " escherichia_coli_reference.fasta -t global -m 1 -s "
                         "-1 -g -1\n";
            return EXIT_FAILURE;
        }

        // Extracting resource paths
        auto path_to_reads = std::string{argv[arg_index]};
        auto path_to_ref = std::string{argv[arg_index + 1]};

        // Determine scan formats
        auto reads_type = mapper::parseFileType(path_to_reads);
        auto ref_type = mapper::parseFileType(path_to_ref);

        // Report start of file loading
        std::cout << "Started loading files\n";

        // Load and parse data from files
        auto reads = mapper::loadFile(path_to_reads, reads_type);
        auto ref = mapper::loadFile(path_to_ref, ref_type);

        // Report end of file loading
        std::cout << "Finsihed loading files\n\ngi";

        // Prinit minimizers stats for reads
        std::cerr << "Procesing reads minimizers data...\n";
        mapper::printMinimizerStats(reads, m_conf);

    } catch (std::exception const& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// cmake -DCMAKE_BUILD_TYPE=Release ..