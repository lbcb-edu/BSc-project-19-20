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
 * Struct array representing long optins
 */
struct option long_options[] = {{"help", no_argument, NULL, 'h'},
                                {"version", no_argument, NULL, 'v'},
                                {NULL, 0, NULL, 0}};

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

template <typename Format, typename ParserType>
auto parse_fastaq(char const* path_to_file) {
    
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

    } catch (std::exception const& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}