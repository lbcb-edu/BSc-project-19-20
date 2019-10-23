#include <unordered_map>
#include <functional>
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <string>
#include <thread>
#include <vector>
#include <memory>
#include <array>

#include <bioparser/bioparser.hpp>

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
 * @brief Doable since FASTA doesn't contain quality information
 */
using FastA = Sequence;

/**
 * @brief FASTQ genome sequences.
 *
 * @details Extends @ref orange::mapper::Sequence
 *      with cigar string quiality information
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
 * @brief Simple intarface for creating @ref bioparser FASTA parser
 */
auto createFastaParser =
    std::function<std::unique_ptr<bioparser::Parser<FastA>>(
        const std::string& path)>{
        &bioparser::createParser<bioparser::FastaParser, FastA>};

/**
 * @brief Simple interface for creating @ref bioparser FASTQ parser
 */
auto createFastqParser =
    std::function<std::unique_ptr<bioparser::Parser<FastQ>>(
        const std::string& path)>{
        &bioparser::createParser<bioparser::FastqParser, FastQ>};

/**
 * @brief Prints program version to stderr
 *
 * @return void
 */
auto printVersion() {
    std::cerr << orange_mapper_VERSION_MAJOR << '.'
              << orange_mapper_VERSION_MINOR << '.'
              << orange_mapper_VERSION_PATCH << '\n';
}

/**
 * @brief Prints program usage information to stderr
 *
 * @return void
 */
auto printHelp() {
    // TODO: Implent some serious help
    std::cerr << "HELP!\n";
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
 * @return int index of first non-option element in argv
 */
auto parseOptions(int argc, char* argv[]) {
    auto opt = int{};
    while ((opt = getopt_long(argc, argv, "hv", long_options, NULL)) != -1) {
        switch (opt) {
            case 'v':
                printVersion();
                break;
            case 'h':
                printHelp();
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
    // Ignore .gz
    auto dot = file.rfind('.');
    if (file[dot + 1] == 'g' && file[dot + 2] == 'z') {
        file.remove_suffix(file.size() - dot);
        dot = file.rfind('.');
    }

    file.remove_prefix(dot);
    if (file == ".fasta" || file == ".fa")
        return FileType::kFasta;
    else if (file == ".fastq" || file == ".fq")
        return FileType::kFastq;

    return FileType::kUnsupported;
}

/**
 * @brief Loads FASTA/FASTQ files into std::vector
 *
 * @param path_to_file path to file contaning FASTA/FASTQ sequences
 *
 * @return @ref orange::mapper::VecSeqPtr
 */
auto loadFile(std::string const& path_to_file, FileType const& type) {
    auto objects = VecSeqPtr{};
    auto parse = [&objects](auto parser) { parser->parse(objects, -1); };

    if (type == FileType::kFasta)
        parse(createFastaParser(path_to_file));
    else
        parse(createFastaParser(path_to_file));

    return objects;
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
 * @param vec_seqs
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
        auto file_index = mapper::parseOptions(argc, argv);

        // Check number of passed arguments
        if (argc - file_index < 2) {
            std::cerr << "Invalid number of arguments.\n"
                         "\tRequired two input files:\n"
                         "\t\t1. Set of fragments (FASTA/FASTQ)\n"
                         "\t\t2. Reference genom (FASTA)\n\n"
                         "\tEg. orange_mapper escherichia_coli_r7_reads.fastq "
                         "escherichia_coli_reference.fasta\n";
        }

        // Extracting resource paths
        auto path_to_reads = std::string{argv[file_index]};
        auto path_to_ref = std::string{argv[file_index + 1]};

        // Determine scan formats
        auto reads_type = mapper::parseFileType(path_to_reads);
        auto ref_type = mapper::parseFileType(path_to_ref);

        // Check if formats are supported
        if (reads_type == mapper::FileType::kUnsupported ||
            ref_type == mapper::FileType::kUnsupported) {
            std::cerr << "Unsuppored file format!\n"
                         "\t Supported: FASTA (reads and reference), FASTQ "
                         "(reference)\n";

            return EXIT_FAILURE;
        }

        // Report start of file loading
        std::cout << "Started loading files\n";

        // Load and parse data from files
        // TODO: Split this process accross multiple threads
        auto reads = mapper::loadFile(path_to_reads, reads_type);
        auto ref = mapper::loadFile(path_to_ref, ref_type);

        // Report end of file loading
        std::cout << "Finsihed loading files\n";

        // Print stats
        mapper::printStats(path_to_reads, reads);
        mapper::printStats(path_to_ref, ref);

    } catch (std::exception const& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}