#include <string_view>
#include <functional>
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <thread>
#include <memory>
#include <array>

#include <bioparser/bioparser.hpp>

#include "orange_mapper_conf.h"

namespace orange {

/**
 * Struct array representing long optins
 */
struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {NULL, 0, NULL, 0}
};

/**
 * @brief Prints program version to stderr
 * 
 * @return void
 */
auto printVersion() {
    std::cerr << ORANGE_MAPPER_VERSION_MAJOR << '.' << 
        ORANGE_MAPPER_VERSION_MINOR << '.' << 
        ORANGE_MAPPER_VERSION_PATCH << '\n';
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
 * @return int index of first non-option element in argv array
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

}

/**
 * @brief Main program entry point
 * 
 * @param argc number of command line arguments
 * @param argv command line arguments (options and/or sources)
 * 
 * @return int execution success status 
 */
int main(int argc, char* argv[]) {

    try {
        auto file_index = orange::parseOptions(argc, argv);

        // Check number of passed arguments
        if (argc - file_index < 2) {
            std::cerr << "Invalid number of arguments.\n" <<
                "\tRequired two input files:\n"
                "\t\t1. Set of fragments (FASTA/FASTQ)\n" <<
                "\t\t2. Reference genom (FASTA)\n\n" <<
                "\tEg. orange_mapper escherichia_coli_r7_reads.fastq " <<
                    "escherichia_coli_reference.fasta\n";
        } 

    } catch (std::exception const& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}