#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <vector>
#include <tuple>
#include <set>
#include <map>
#include <functional>
#include <cmath>
#include <cstdlib>

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "pink_alignment.h"
#include "pink_minimizers.h"

#define VERSION "v0.1.0"

#define EXTENSIONS_A_SIZE 4
#define EXTENSIONS_Q_SIZE 4

#define MATCH 3
#define MISMATCH -1
#define GAP -2
#define K 15
#define W 5
#define F 0.001
#define CIGAR false
#define BAND_WIDTH 10000
#define DEFAULT_QUALITY 255
#define T 1

std::string extensions_a[]{".fasta", ".fa", ".fasta.gz", ".fa.gz"};
std::string extensions_q[]{".fastq", ".fq", ".fastq.gz", ".fq.gz"};

void help() {
    std::cout << "This is help option." << std::endl
              << "This program accepts two files as floating arguments - first file (which contains a set of fragments) in FASTA or FASTQ format, and second one (which contains a corresponding reference genome) in FASTA format."
              << std::endl
              << "This program has following options:\n" << std::endl
              << "    -h   --help\n       prints the help menu\n" << std::endl
              << "    -v   --version\n        prints the version number\n" << std::endl
              << "    -G   --global\n       sets Alignment type to global, which is default type.\n" << std::endl
              << "    -S   --semi_global\n       sets Alignment type to semi-global. Default type is global.\n"
              << std::endl
              << "    -L   --local\n       sets Alignment type to local. Default type is local.\n" << std::endl
              << "    -m   --match\n       sets match argument for pairwise alignment to following argument. Default value is 3.\n"
              << std::endl
              << "    -s   --mismatch\n       sets mismatch argument for pairwise alignment to following argument. Default value is -1.\n"
              << std::endl
              << "    -g   --gap\n       sets gap argument for pairwise alignment to following argument. Default value is -2.\n"
              << std::endl
              << "    -k\n       sets k for k-mers. Default value is 15.\n" << std::endl
              << "    -w\n       sets window length for finding interior minimizers. Default value is 5.\n" << std::endl
              << "    -f\n       sets percentage of minimizers which are not taken in account. Default value is 0.001\n"
              << std::endl
              << "    -c\n       includes the CIGAR strings of the alignment in the mapper printing\n" << std::endl
              << "    -t\n       sets number of threads with whom mapper is parallelized\n" << std::endl;
}

void version() {
    std::cout << "version: " << VERSION << std::endl;
}

void errorMessage() {
    std::cerr << "Your input is not ok. Please try again. You can ask for help with \"-h\" or \"--help\".\n";
    exit(1);
}

bool is_extension_ok(std::string const& argument, char type) {
    std::size_t extension_start = argument.find_last_of('.');
    std::string extension = argument.substr(extension_start);
    if (extension == ".gz") {
        std::string before_extension = argument.substr(0, extension_start - 1);
        extension_start = before_extension.find_last_of('.');
        extension = argument.substr(extension_start);
    }
    std::string *p;
    switch (type) {
        case 'a':
            p = std::find(extensions_a, extensions_a + EXTENSIONS_A_SIZE, extension);
            if (p == extensions_a + EXTENSIONS_A_SIZE)
                return false;
            else
                return true;
        case 'q':
            p = std::find(extensions_q, extensions_q + EXTENSIONS_Q_SIZE, extension);
            if (p == extensions_q + EXTENSIONS_Q_SIZE)
                return false;
            else
                return true;
        default:
            return false;
    }
}

class Fast {
public:
    std::string name;
    std::string sequence;
    std::string quality;

    Fast(
            const char *name, std::uint32_t name_length,
            const char *sequence, std::uint32_t sequence_length) :
            name{std::string(name, name_length)},
            sequence{std::string(sequence, sequence_length)} {}

    Fast(
            const char *name, std::uint32_t name_length,
            const char *sequence, std::uint32_t sequence_length,
            const char *quality, std::uint32_t quality_length) :
            name{std::string(name, name_length)},
            sequence{std::string(sequence, sequence_length)},
            quality{std::string(quality, quality_length)} {}
};

std::vector<std::unique_ptr<Fast>> parse_fasta(std::string const& fastaFile) {
    std::vector<std::unique_ptr<Fast>> fasta_objects;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(fastaFile);
    // read the whole file
    fasta_parser->parse(fasta_objects, -1);

    return fasta_objects;
}

std::vector<std::unique_ptr<Fast>> parse_fastq(std::string const& fastqFile) {
    std::vector<std::unique_ptr<Fast>> fastq_objects;
    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fast>(fastqFile);
    // read a predefined size of bytes
    std::uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
    while (true) {
        auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
        if (!status) {
            break;
        }
    }

    return fastq_objects;
}

void print_statistics(std::string const& file) {
    std::vector<std::unique_ptr<Fast>> fast_objects;
    if (is_extension_ok(file, 'a')) {
        fast_objects = parse_fasta(file);
    } else {
        fast_objects = parse_fastq(file);
    }

    unsigned int num_of_seq, sum, min_length, max_length;
    float average_length;

    num_of_seq = fast_objects.size();

    sum = 0;
    min_length = (fast_objects[0]->sequence).length();
    max_length = (fast_objects[0]->sequence).length();
    for (int i = 0; i < num_of_seq; i++) {
        sum += (fast_objects[i]->sequence).length();

        if ((fast_objects[i]->sequence).length() < min_length)
            min_length = (fast_objects[i]->sequence).length();

        if ((fast_objects[i]->sequence).length() > max_length)
            max_length = (fast_objects[i]->sequence).length();

    }
    average_length = (float) sum / num_of_seq;


    std::cerr << "File name: " << file << std::endl;
    std::cerr << "Number of sequences: " << num_of_seq << std::endl;
    std::cerr << "Average length: " << average_length << std::endl;
    std::cerr << "Minimal length: " << min_length << std::endl;
    std::cerr << "Maximal length: " << max_length << std::endl;
}

auto comparator_by_amount = [](std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>> &a,
                               std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>> &b) {
    return (a.second).size() < (b.second).size();
};

auto comparator_by_query_position_ascending = [](std::tuple<unsigned int, unsigned int, bool> &a,
                                                 std::tuple<unsigned int, unsigned int, bool> &b) {
    return std::get<0>(a) < std::get<0>(b);
};

auto comparator_by_query_position_descending = [](std::tuple<unsigned int, unsigned int, bool> &a,
                                                  std::tuple<unsigned int, unsigned int, bool> &b) {
    return std::get<0>(a) > std::get<0>(b);
};

uint32_t create_reference_genome_minimizer_index(std::vector<std::unique_ptr<Fast>> const& fast_objects,
                                             unsigned int k, unsigned int w, double f,
                                             std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &target_minimizer_index) {

    std::cerr << "Creating refernce genome minimizer index" << std::endl;

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = pink::minimizers(
            (fast_objects.front()->sequence).c_str(), (fast_objects.front()->sequence).length(), k, w);

    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator index_iterator;

    for (auto minimizer : minimizers) {
        index_iterator = target_minimizer_index.find(std::get<0>(minimizer));

        if (index_iterator == target_minimizer_index.end()) {
            std::vector<std::pair<unsigned int, bool>> locations_strands;
            locations_strands.emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
            target_minimizer_index.insert(std::make_pair(std::get<0>(minimizer), locations_strands));
        } else {
            (index_iterator->second).emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
        }
    }

    std::cerr << "Reference genome minimizer index done!" << std::endl;

    std::vector<uint32_t> occurrences;
    for(auto minimizer : target_minimizer_index) {
        occurrences.emplace_back(minimizer.second.size());
    }
    std::sort(occurrences.begin(), occurrences.end());
    unsigned int ignore = std::round(occurrences.size() * f);

    return occurrences[occurrences.size() - 1 - ignore];
}

void create_fragment_minimizer_index(std::unique_ptr<Fast> const& fast_objects_i, unsigned int k, unsigned int w,
                                std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &fragment_minimizer_index) {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = pink::minimizers(
            (fast_objects_i->sequence).c_str(), (fast_objects_i->sequence).length(), k, w);

    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator index_iterator;

    for (auto minimizer : minimizers) {
        index_iterator = fragment_minimizer_index.find(std::get<0>(minimizer));

        if (index_iterator == fragment_minimizer_index.end()) {
            std::vector<std::pair<unsigned int, bool>> locations_strands;
            locations_strands.emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
            fragment_minimizer_index.insert(std::make_pair(std::get<0>(minimizer), locations_strands));
        } else {
            (index_iterator->second).emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
        }
    }
}

void make_match_groups(std::vector<std::tuple<unsigned int, unsigned int, bool>> const& matches,
                       std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> &match_groups) {
    std::vector<std::tuple<unsigned int, unsigned int, bool>> match_group;

    match_group.emplace_back(matches[0]);
    for (int i = 1; i < matches.size(); i++) {
        if (abs(std::get<0>(matches[i]) - std::get<0>(matches[i - 1])) > BAND_WIDTH) {
            match_groups.emplace_back(match_group);
            match_group.clear();
            match_group.emplace_back(matches[i]);
        } else {
            match_group.emplace_back(matches[i]);
        }
    }
    match_groups.emplace_back(match_group);

}

void find_matches(std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> const& fragment_minimizer_index,
                  std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> const& target_minimizer_index, int ignoring,
                  std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> &match_groups) {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> matches;

    for (auto minimizer : fragment_minimizer_index) {
        auto it = target_minimizer_index.find(minimizer.first);
        if (it != target_minimizer_index.end() && (it->second).size() <= ignoring) {
            for (auto location_strand_query : minimizer.second) {
                for (auto location_strand_target : it->second) {
                    matches.emplace_back(std::make_tuple(location_strand_query.first, location_strand_target.first,
                                                         location_strand_query.second ^ location_strand_target.second));
                }
            }
        }
    }

    if (!matches.empty()) {
        std::vector<std::tuple<unsigned int, unsigned int, bool>> matches_same_strand;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> matches_opposite_strand;

        for (auto match : matches) {
            if (std::get<2>(match)) {
                matches_opposite_strand.emplace_back(match);
            } else {
                matches_same_strand.emplace_back(match);
            }
        }

        if (!matches_same_strand.empty()) {
            std::sort(matches_same_strand.begin(), matches_same_strand.end(), comparator_by_query_position_ascending);
            make_match_groups(matches_same_strand, match_groups);
        }

        if (!matches_opposite_strand.empty()) {
            std::sort(matches_opposite_strand.begin(), matches_opposite_strand.end(),
                      comparator_by_query_position_descending);
            make_match_groups(matches_opposite_strand, match_groups);
        }
    }

}

std::vector<std::tuple<unsigned int, unsigned int, bool>> longest_increasing_subsequence(std::vector<std::tuple<unsigned int, unsigned int, bool>> const& matches) {
    unsigned int n = matches.size();
    std::vector<unsigned int> tail_indexes(n, 0); // Initialized with 0
    std::vector<unsigned int> prev_indexes(n, -1); // initialized with -1

    auto comparator_by_target_position = [&matches](unsigned int const &a, unsigned int const &b) {
        return std::get<1>(matches[a]) < std::get<1>(matches[b]);
    };

    int length = 1; // it will always point to empty location

    for (unsigned int i = 1; i < n; i++) {
        auto b = tail_indexes.begin(), e = tail_indexes.begin() + length;
        auto it = lower_bound(b, e, i, comparator_by_target_position);

        if (it == e) {
            prev_indexes[i] = tail_indexes[length - 1];
            tail_indexes[length++] = i;
        } else {
            *it = i;
            int pos = std::distance(b, it);
            prev_indexes[i] = tail_indexes[pos - 1];
        }
    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> result;

    for (int i = tail_indexes[length - 1]; i >= 0; i = prev_indexes[i])
        result.emplace_back(matches[i]);

    return result;
}

std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> find_region(std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> const& match_groups) {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> candidate = longest_increasing_subsequence(match_groups.front());
    int max_length;
    std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> region;

    int query_beg, query_end, target_beg, target_end;
    query_beg = std::get<0>(candidate.back()) <= std::get<0>(candidate.front()) ? std::get<0>(candidate.back()) : std::get<0>(candidate.front());
    query_end = std::get<0>(candidate.back()) > std::get<0>(candidate.front()) ? std::get<0>(candidate.back()) : std::get<0>(candidate.front());
    target_beg = std::get<1>(candidate.back()) <= std::get<1>(candidate.front()) ? std::get<1>(candidate.back()) : std::get<1>(candidate.front());
    target_end = std::get<1>(candidate.back()) > std::get<1>(candidate.front()) ? std::get<1>(candidate.back()) : std::get<1>(candidate.front());

    region = std::make_tuple(query_beg, query_end, target_beg, target_end, std::get<2>(candidate.front()));
    max_length = std::max(std::get<1>(region) - std::get<1>(region), std::get<3>(region) - std::get<2>(region));

    std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> region_help;
    for (auto const& vec : match_groups) {
        candidate = longest_increasing_subsequence(vec);
        query_beg = std::get<0>(candidate.back()) <= std::get<0>(candidate.front()) ? std::get<0>(candidate.back()) : std::get<0>(candidate.front());
        query_end = std::get<0>(candidate.back()) > std::get<0>(candidate.front()) ? std::get<0>(candidate.back()) : std::get<0>(candidate.front());
        target_beg = std::get<1>(candidate.back()) <= std::get<1>(candidate.front()) ? std::get<1>(candidate.back()) : std::get<1>(candidate.front());
        target_end = std::get<1>(candidate.back()) > std::get<1>(candidate.front()) ? std::get<1>(candidate.back()) : std::get<1>(candidate.front());

        region_help = std::make_tuple(query_beg, query_end, target_beg, target_end, std::get<2>(candidate.front()));
        if(std::max(query_end - query_beg, target_end - target_beg) > max_length)
            region = region_help;
    }
    return region;
}

std::string paf_string(std::string const& query, const char *query_name, unsigned int query_length,
                       std::string const& target, const char *target_name, unsigned int target_length,
                       pink::AlignmentType type, int match, int mismatch, int gap, bool include_cigar,
                       std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> const& region) {

    std::string ret = std::string(query_name) + '\t' +
                      std::to_string(query_length) + '\t' +
                      std::to_string(std::get<0>(region)) + '\t' +
                      std::to_string(std::get<1>(region)) + '\t' +
                      std::string(std::get<4>(region) ? "-" : "+") + '\t' +
                      std::string(target_name) + '\t' +
                      std::to_string(target_length) + '\t' +
                      std::to_string(std::get<2>(region)) + '\t' +
                      std::to_string(std::get<3>(region)) + '\t';

    unsigned int query_substring_length = std::get<1>(region) - std::get<0>(region) + 1;
    unsigned int target_substring_length = std::get<3>(region) - std::get<2>(region) + 1;

    std::string query_substring = query.substr(std::get<0>(region), query_substring_length);
    std::string target_substring = target.substr(std::get<2>(region), target_substring_length);

    unsigned int num_of_matches = std::max(std::get<1>(region) - std::get<0>(region), std::get<3>(region) - std::get<2>(region));
    unsigned int block_length = num_of_matches;
    int mapping_quality = DEFAULT_QUALITY;

    std::string cigar;
    unsigned int target_begin;
    if (include_cigar) {
        pink::pairwise_alignment(query_substring.c_str(), query_substring_length, target_substring.c_str(),
                                 target_substring_length, type, match, mismatch, gap, cigar, target_begin);

        num_of_matches = 0;
        block_length = 0;
        for (int i = 0; i < cigar.length(); i++) {
            if (cigar.at(i) == '=')
                num_of_matches += cigar.at(i - 1) - '0';
            if(isdigit(cigar.at(i)))
                block_length += cigar.at(i) - '0';
        }
    }

    ret += std::to_string(num_of_matches) + '\t' +
           std::to_string(block_length) + '\t' +
           std::to_string(mapping_quality);

    if (include_cigar)
        ret += "\tcg:Z:" + cigar;

    return ret;
}

std::string work_with_fragment(std::vector<std::unique_ptr<Fast>> const & fast_objects1,
                                std::vector<std::unique_ptr<Fast>> const& fast_objects2,
                                std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> const& target_minimizer_index, int ignoring,
                                unsigned int k, unsigned int w, pink::AlignmentType type, int match, int mismatch, int gap,
                                bool include_cigar,
                                int &i) {
    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> fragment_minimizer_index;
    create_fragment_minimizer_index(fast_objects1[i], k, w, fragment_minimizer_index);

    std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> match_groups;
    find_matches(fragment_minimizer_index, target_minimizer_index, ignoring, match_groups);

    if (match_groups.empty()) {
        return "";
    }

    // <query_beginning, query_ending, target_beginning, target_ending, true=opposite/false=same strand>
    std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> candidate = find_region(match_groups);

    const char *query = (fast_objects1[i]->sequence).c_str();
    const char *target = (fast_objects2.front()->sequence).c_str();
    unsigned int query_length = (fast_objects1[i]->sequence).length();
    unsigned int target_length = (fast_objects2.front()->sequence).length();
    const char *query_name = (fast_objects1[i]->name).c_str();
    const char *target_name = (fast_objects2.front()->name).c_str();

    std::string paf = paf_string(query, query_name, query_length, target, target_name, target_length, type, match,
                                 mismatch, gap, include_cigar, candidate);

    return paf;
}

int main(int argc, char *argv[]) {

    int opt;
    std::string alignment_type;
    pink::AlignmentType type = pink::global;
    int match = MATCH;
    int mismatch = MISMATCH;
    int gap = GAP;
    int k = K;
    int w = W;
    double f = F;
    bool include_cigar = CIGAR;
    int t = T;
    while ((opt = getopt(argc, argv, ":hvGSLm:s:g:k:w:f:ct:")) != -1) {
        switch (opt) {
            case 'h' : help(); return 0;
            case 'v' : version(); return 0;
            case 'G' : type = pink::global; break;
            case 'S' : type = pink::semi_global; break;
            case 'L' : type = pink::local; break;
            case 'm' : match = atoi(optarg); break;
            case 's' : mismatch = atoi(optarg); break;
            case 'g' : gap = atoi(optarg); break;
            case 'k' : k = atoi(optarg); break;
            case 'w' : w = atoi(optarg); break;
            case 'f' : f = atoi(optarg); break;
            case 'c' : include_cigar = true; break;
            case 't' : t = atoi(optarg); break;
            default: errorMessage();
        }

    }

    if (argc != optind + 2) {
        errorMessage();
    }

    std::string first = argv[optind];
    std::string second = argv[optind + 1];

    if (!(is_extension_ok(first, 'a') || is_extension_ok(first, 'q')) || !is_extension_ok(second, 'a')) {
        errorMessage();
    }

    std::vector<std::unique_ptr<Fast>> fast_objects1;
    if (is_extension_ok(first, 'a')) {
        fast_objects1 = parse_fasta(first);
    } else {
        fast_objects1 = parse_fastq(first);
    }

    std::vector<std::unique_ptr<Fast>> fast_objects2;
    fast_objects2 = parse_fasta(second);

    // FIRST

//     print_statistics(first);
//     std::cerr << std::endl;
//     print_statistics(second);


    // SECOND

//    int which1 = std::rand() % fast_objects1.size();
//    int which2 = std::rand() % fast_objects1.size();
//
//    const char *query = (fast_objects1[which1]->sequence).c_str();
//    const char *target = (fast_objects1[which2]->sequence).c_str();
//    unsigned int query_length = (fast_objects1[which1]->sequence).length();
//    unsigned int target_length = (fast_objects1[which2]->sequence).length();
//
//    int cost = pink::pairwise_alignment(query, query_length, target, target_length, type, match, mismatch, gap);
//
//    std::string query_string = std::string(query, query_length);
//    std::string target_string = std::string(target, target_length);
//
//    std::cout << "Cost for pairwise alignment for " << query_string << " and " << target_string << " is: " << cost << std::endl;
//
//    std::string cigar;
//    unsigned int target_begin;
//    pink::pairwise_alignment(query, query_length, target, target_length, type, match, mismatch, gap, cigar, target_begin);
//    std::cout << "Cigar string: " << cigar << std::endl;


    // THIRD

//    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;
//
//    for (int i = 0; i < fast_objects1.size(); i++) {
//
//        std::vector<std::tuple<unsigned int, unsigned int, bool>> current = pink::minimizers(
//                (fast_objects1[i]->sequence).c_str(), (fast_objects1[i]->sequence).length(), k, w);
//
//        for (auto minimizer : current) {
//            minimizers.emplace_back(minimizer);
//        }
//
//    }
//
//    std::map<unsigned int, int> counter;
//    std::map<unsigned int, int>::iterator counter_iterator;
//
//    for (auto minimizer : minimizers) {
//        counter_iterator = counter.find(std::get<0>(minimizer));
//
//        if (counter_iterator == counter.end()) {
//            counter.insert(std::make_pair(std::get<0>(minimizer), 1));
//        } else {
//            counter_iterator->second++;
//        }
//    }
//
//    std::cout << "Number of distinct minimizers is "
//              << counter.size()
//              << std::endl;
//
//    int num_of_singletons = 0;
//    for (auto no_minimizer : counter) {
//        if (no_minimizer.second == 1)
//            num_of_singletons++;
//    }
//
//    std::cout << "Fraction of singletons: " << (double) num_of_singletons / counter.size() << std::endl;
//
//    std::vector<uint32_t> occurrences;
//    for(auto minimizer : counter) {
//        occurrences.emplace_back(minimizer.second);
//    }
//    std::sort(occurrences.begin(), occurrences.end());
//    unsigned int ignore = std::round(occurrences.size() * f);
//
//    std::cout << "Number of occurrences of the most frequent minimizer when the top " << f
//              << " frequent minimizers are not taken in account: "
//              << occurrences[occurrences.size() - 1 - ignore] << std::endl;

    // FINAL

    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> target_minimizer_index;
    uint32_t ignoring = create_reference_genome_minimizer_index(fast_objects2, k, w, f, target_minimizer_index);

    std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);
    std::vector<std::future<std::string>> thread_futures;

    unsigned int work_for_each_thread = fast_objects1.size() / t;

    unsigned int thread_begin = 0;
    unsigned int thread_end = work_for_each_thread + (fast_objects1.size() % t);

    for (int i = 0; i < fast_objects1.size(); i++) {
        thread_futures.emplace_back(
                thread_pool->submit(work_with_fragment, std::ref(fast_objects1), std::ref(fast_objects2),
                                    std::ref(target_minimizer_index), ignoring,
                                    k, w, type, match, mismatch, gap, include_cigar,
                                    i));
    }

    for (auto &it: thread_futures) {
        auto paf = it.get();
        std::cout << paf << std::endl;
    }

    return 0;
}
