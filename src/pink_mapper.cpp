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

bool is_extension_ok(std::string argument, char type) {
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

std::vector<std::unique_ptr<Fast>> parse_fasta(std::string fastaFile) {
    std::vector<std::unique_ptr<Fast>> fasta_objects;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(fastaFile);
    // read the whole file
    fasta_parser->parse(fasta_objects, -1);

    return fasta_objects;
}

std::vector<std::unique_ptr<Fast>> parse_fastq(std::string fastqFile) {
    std::vector<std::unique_ptr<Fast>> fastq_objects;
    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fast>(fastqFile);
    // read a predefined size of bytes
    std::uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
    while (true) {
        auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
        if (status == false) {
            break;
        }
    }

    return fastq_objects;
}

void print_statistics(std::string file) {
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

void create_reference_genome_minimizer_index(const std::vector<std::unique_ptr<Fast>> &fast_objects,
                                        unsigned int k, unsigned int w, double f,
                                        std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &target_minimizer_index) {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = pink::minimizers(
            (fast_objects.front()->sequence).c_str(), (fast_objects.front()->sequence).length(), k, w);

    // key is minimizer itself, value is vector of all locations and strands where that minimizer appears
    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> counter;
    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator counter_iterator;

    for (auto minimizer : minimizers) {
        counter_iterator = counter.find(std::get<0>(minimizer));

        if (counter_iterator == counter.end()) {
            std::vector<std::pair<unsigned int, bool>> locations_strands;
            locations_strands.emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
            counter.insert(std::make_pair(std::get<0>(minimizer), locations_strands));
        } else {
            (counter_iterator->second).emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
        }
    }

    std::vector<std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>>> pairs(counter.begin(),
                                                                                           counter.end());

    sort(pairs.begin(), pairs.end(), comparator_by_amount);
    unsigned int ignore = std::round(counter.size() * f);
    pairs.resize(minimizers.size() - ignore);
    pairs.shrink_to_fit();

    for (auto current : pairs) {
        target_minimizer_index.insert(current);
    }
}

void
create_fragment_minimizer_index(const std::unique_ptr<Fast> &fast_objects_i, unsigned int k, unsigned int w,
                                std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &fragment_minimizer_index) {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = pink::minimizers((fast_objects_i->sequence).c_str(), (fast_objects_i->sequence).length(), k, w);

    // key is minimizer itself, value is vector of all locations and strands where that minimizer appears
    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> counter;
    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator counter_iterator;

    for (auto minimizer : minimizers) {
        counter_iterator = counter.find(std::get<0>(minimizer));

        if (counter_iterator == counter.end()) {
            std::vector<std::pair<unsigned int, bool>> locations_strands;
            locations_strands.emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
            counter.insert(std::make_pair(std::get<0>(minimizer), locations_strands));
        } else {
            (counter_iterator->second).emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
        }
    }

    std::vector<std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>>> pairs(counter.begin(), counter.end());

    sort(pairs.begin(), pairs.end(), comparator_by_amount);
    for (auto current : pairs) {
        fragment_minimizer_index.insert(current);
    }
}

void make_match_groups(std::vector<std::tuple<unsigned int, unsigned int, bool>> matches,
                       std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> &match_groups) {
    std::vector<std::tuple<unsigned int, unsigned int, bool>> match_group;
    match_group.emplace_back(matches[0]);
    for (int i = 1; i < matches.size(); i++) {
        if (abs(std::get<0>(matches[i]) - std::get<0>(matches[i-1])) > BAND_WIDTH){
            match_groups.emplace_back(match_group);
            match_group.clear();
            match_group.emplace_back(matches[i]);
        } else {
            match_group.emplace_back(matches[i]);
        }
    }

}

void find_matches(std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> fragment_minimizer_index,
                  std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> target_minimizer_index,
                  std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> &match_groups) {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> matches;

    for (auto minimizer : fragment_minimizer_index) {
        auto it = target_minimizer_index.find(minimizer.first);
        if (it != target_minimizer_index.end()) {
            for (auto location_strand_query : minimizer.second) {
                for (auto location_strand_target : it->second) {
                        matches.emplace_back(std::make_tuple(location_strand_query.first, location_strand_target.first, location_strand_query.second ^ location_strand_target.second));
                }
            }
        }
    }

    if(!matches.empty()) {
        std::vector<std::tuple<unsigned int, unsigned int, bool>> matches_same_strand;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> matches_opposite_strand;

        for(auto match : matches){
            if(std::get<2>(match)){
                matches_opposite_strand.emplace_back(match);
            } else {
                matches_same_strand.emplace_back(match);
            }
        }

        if(!matches_same_strand.empty()) {
            std::sort(matches_same_strand.begin(), matches_same_strand.end(), comparator_by_query_position_ascending);
            make_match_groups(matches_same_strand, match_groups);
        }

        if(!matches_opposite_strand.empty()) {
            std::sort(matches_opposite_strand.begin(), matches_opposite_strand.end(), comparator_by_query_position_descending);
            make_match_groups(matches_opposite_strand, match_groups);
        }
    }

}

int CeilIndex(std::vector<std::tuple<unsigned int, unsigned int>> &v, std::vector<unsigned int> &T, int start, int end,
              std::tuple<unsigned int, unsigned int> key) {
    while (end - start > 1) {
        int middle = start + (end - start) / 2;
        // pazi da ti ostane lokacija 2
        if (std::get<1>(v[T[middle]]) >= std::get<1>(key))
            end = middle;
        else
            start = middle;
    }

    return end;
}


std::vector<std::tuple<unsigned int, unsigned int, bool>> longest_increasing_subsequence(
        std::vector<std::tuple<unsigned int, unsigned int, bool>> matches) {
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

//        if (std::get<1>(matches[i]) < std::get<1>(matches[tail_indexes[0]])) {
//            // new smallest value
//            tail_indexes[0] = i;
//        } else if (std::get<1>(matches[i]) > std::get<1>(matches[tail_indexes[length - 1]])) {
//            // arr[i] wants to extend largest subsequence
//            prev_indexes[i] = tail_indexes[length - 1];
//            tail_indexes[length++] = i;
//        } else {
//            // arr[i] wants to be a potential candidate of
//            // future subsequence
//            // It will replace ceil value in tailIndices
//
//            int pos = CeilIndex(matches, tail_indexes, -1,
//                                length - 1, matches[i]);
//
//            prev_indexes[i] = tail_indexes[pos - 1];
//            tail_indexes[pos] = i;
//        }

    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> result;

    for (int i = tail_indexes[length - 1]; i >= 0; i = prev_indexes[i])
        result.emplace_back(matches[i]);

    std::reverse(result.begin(), result.end());

    return result;
}

std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> find_region(std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> match_groups) {
    std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> candidate_group;

    for(auto vec : match_groups) {
        std::vector<std::tuple<unsigned int, unsigned int, bool>> candidate = longest_increasing_subsequence(vec);
        candidate_group.emplace_back(candidate);
    }

    int max_length = 0;
    int index_of_max_length = 0;
    for(int i = 0; i < candidate_group.size(); i++){
        if(candidate_group[i].size() > max_length) {
            max_length = candidate_group[i].size();
            index_of_max_length = i;
        }
    }

    std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> region = std::make_tuple(
            std::get<0>(candidate_group[index_of_max_length].front()), std::get<0>(candidate_group[index_of_max_length].back()),
            std::get<1>(candidate_group[index_of_max_length].front()), std::get<1>(candidate_group[index_of_max_length].back()),
            std::get<2>(candidate_group[index_of_max_length].front()));

    return region;
}

std::string paf_string(std::string query, const char *query_name, unsigned int query_length,
                       std::string target, const char *target_name, unsigned int target_length,
                       pink::AlignmentType type, int match, int mismatch, int gap, bool include_cigar,
                       std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> region) {

    std::string ret = std::string(query_name) + '\t' +
                      std::to_string(query_length) + '\t' +
                      std::to_string(std::get<0>(region)) + '\t' +
                      std::to_string(std::get<1>(region)) + '\t' +
                      std::string(std::get<4>(region) ? "+" : "-") + '\t' +
                      std::string(target_name) + '\t' +
                      std::to_string(target_length) + '\t' +
                      std::to_string(std::get<2>(region)) + '\t' +
                      std::to_string(std::get<3>(region)) + '\t';

    unsigned int query_substring_length = std::get<1>(region) - std::get<0>(region) + 1;
    unsigned int target_substring_length = std::get<3>(region) - std::get<2>(region) + 1;

    std::string query_substring = query.substr(std::get<0>(region), query_substring_length);
    std::string target_substring = target.substr(std::get<2>(region), target_substring_length);

    unsigned int num_of_matches = std::get<1>(region) - std::get<0>(region);
    unsigned int block_length = num_of_matches;
    int mapping_quality = DEFAULT_QUALITY;

    std::string cigar;
    unsigned int target_begin;
    if (include_cigar) {
        pink::pairwise_alignment(query_substring.c_str(), query_substring_length, target_substring.c_str(),
                                 target_substring_length, type, match, mismatch, gap, cigar, target_begin);

        num_of_matches = 0;
        for (char c : cigar) {
            if (c == '=')
                num_of_matches++;
        }
        block_length = cigar.size();
    }

    ret += std::to_string(num_of_matches) + '\t' +
           std::to_string(block_length) + '\t' +
           std::to_string(mapping_quality);

    if (include_cigar)
        ret += "\tcg:Z:" + cigar + '\n';
    else
        ret += '\n';

    return ret;
}

int work_with_fragments(const std::vector<std::unique_ptr<Fast>> &fast_objects1, const std::vector<std::unique_ptr<Fast>> &fast_objects2,
                    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> target_minimizer_index,
                    unsigned int k, unsigned int w, pink::AlignmentType type, int match, int mismatch, int gap,
                    bool include_cigar,
                    int thread_begin, int thread_end) {
    for (int i = thread_begin; i < thread_end; i++) {
        std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> fragment_minimizer_index;
        create_fragment_minimizer_index(fast_objects1[i], k, w, fragment_minimizer_index);

        std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> match_groups;
        find_matches(fragment_minimizer_index, target_minimizer_index, match_groups);

        if (match_groups.empty())
            continue;

        std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> candidate = find_region(match_groups);

        const char *query = (fast_objects1[i]->sequence).c_str();
        const char *target = (fast_objects2.front()->sequence).c_str();
        unsigned int query_length = (fast_objects1[i]->sequence).length();
        unsigned int target_length = (fast_objects2.front()->sequence).length();
        const char *query_name = (fast_objects1[i]->name).c_str();
        const char *target_name = (fast_objects2.front()->name).c_str();

        std::string paf = paf_string(query, query_name, query_length, target, target_name, target_length, type, match,
                                     mismatch, gap, include_cigar, candidate);

        std::cout << paf << std::endl;
    }


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
    int t = 0;
    while ((opt = getopt(argc, argv, ":hvt:m:s:g:")) != -1) {
        switch (opt) {
            case 'h':
                help();
                return 0;
            case 'v':
                version();
                return 0;
            case 'G' :
                type = pink::global;
                break;
            case 'S' :
                type = pink::semi_global;
                break;
            case 'L' :
                type = pink::local;
                break;
            case 'm':
                match = atoi(optarg);
                break;
            case 's' :
                mismatch = atoi(optarg);
                break;
            case 'g' :
                gap = atoi(optarg);
                break;
            case 'k' :
                k = atoi(optarg);
                break;
            case 'w' :
                w = atoi(optarg);
                break;
            case 'f' :
                f = atoi(optarg);
                break;
            case 'c' :
                include_cigar = true;
                break;
            case 't' :
                t = atoi(optarg);
                break;
            default:
                errorMessage();
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
//    std::string s = "";
//    for(int i = 0; i < query_length; i++){
//        s += query[i];
//    }
//    std::string t = "";
//    for(int i = 0; i < target_length; i++){
//        t += target[i];
//    }
//
//    std::cout << "Cost for pairwise alignment for " << s << " and " << t << " is: " << cost << std::endl;
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
//    // key is minimizer itself, value is vector of all locations and strands where that minimizer appears
//    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> counter;
//    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator counter_iterator;
//
//    for (auto minimizer : minimizers) {
//        counter_iterator = counter.find(std::get<0>(minimizer));
//
//        if (counter_iterator == counter.end()) {
//            std::vector<std::pair<unsigned int, bool>> locations_strands;
//            locations_strands.emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
//            counter.insert(std::make_pair(std::get<0>(minimizer), locations_strands));
//        } else {
//            (counter_iterator->second).emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
//        }
//    }
//
//    std::vector<std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>>> pairs(counter.begin(),
//                                                                                           counter.end());
//
//    sort(pairs.begin(), pairs.end(), [=](std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>> &a,
//                                         std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>> &b) {
//             return (a.second).size() < (b.second).size();
//         }
//    );
//
//    std::cout << "Number of distinct minimizers is "
//              << counter.size()
//              << std::endl;
//
//    int num_of_singletons = 0;
//    for (auto p : pairs) {
//        if ((p.second).size() == 1)
//            num_of_singletons++;
//        else
//            break;
//    }
//
//    std::cout << "Fraction of singletons: " << (double) num_of_singletons / counter.size() << std::endl;
//
//    unsigned int ignore = std::round(counter.size() * f);
//
//    std::cout << "Number of occurrences of the most frequent minimizer when the top " << f
//              << " frequent minimizers are not taken in account: "
//              << pairs[pairs.size() - 1 - ignore].second.size() << std::endl;

    // FINAL

    std::map<unsigned int, std::vector<std::pair<unsigned int, bool>>> target_minimizer_index;
    create_reference_genome_minimizer_index(fast_objects2, k, w, f, target_minimizer_index);

    std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);
    std::vector<std::future<int>> thread_futures;

    unsigned int work_for_each_thread = fast_objects1.size() / t;

    unsigned int thread_begin = 0;
    unsigned int thread_end = work_for_each_thread + (fast_objects1.size() % t);

    for (int i = 0; i < t; i++) {
        thread_futures.emplace_back(
                thread_pool->submit(work_with_fragments, std::ref(fast_objects1), std::ref(fast_objects2),
                                    std::ref(target_minimizer_index),
                                    k, w, type, match, mismatch, gap, include_cigar,
                                    thread_begin, thread_end));
        thread_begin = thread_end;
        thread_end = thread_begin + work_for_each_thread;
    }

    for (auto &it: thread_futures) {
        it.wait();
    }

    return 0;
}