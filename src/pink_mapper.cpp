#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <vector>
#include <tuple>
#include <set>
#include <map>
#include <functional>

#include "bioparser/bioparser.hpp"
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

std::string extensions_a[] {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
std::string extensions_q[] { ".fastq", ".fq", ".fastq.gz", ".fq.gz"};

void help() {
    std::cout << "This is help option." << std::endl
              << "This program accepts two files as floating arguments - first file (which contains a set of fragments) in FASTA or FASTQ format, and second one (which contains a corresponding reference genome) in FASTA format." << std::endl
              << "This program has following options:\n" << std::endl
              << "    -h   --help\n       prints the help menu\n" << std::endl
              << "    -v   --version\n        prints the version number\n" << std::endl
              << "    -t   --type\n       sets Alignment type to following argument. Default type is global.\n" << std::endl
              << "    -m   --match\n       sets match argument for pairwise alignment to following argument. Default value is 3.\n" << std::endl
              << "    -s   --mismatch\n       sets mismatch argument for pairwise alignment to following argument. Default value is -1.\n" << std::endl
              << "    -g   --gap\n       sets gap argument for pairwise alignment to following argument. Default value is -2.\n" << std::endl
              << "    -k\n       sets k for k-mers. Default value is 15.\n" << std::endl
              << "    -w\n       sets window length for finding interior minimizers. Default value is 5.\n" << std::endl
              << "    -f\n       sets percentage of minimizers which are not taken in account. Default value is 0.001\n" << std::endl;
}

void version(){
    std::cout << "version: " << VERSION << std::endl;
}

void errorMessage() {
    std::cerr << "Your input is not ok. Please try again. You can ask for help with \"-h\" or \"--help\".\n";
    exit(1);
}

bool is_extension_ok(std::string argument, char type){
    std::size_t extension_start = argument.find_last_of('.');
    std::string extension = argument.substr(extension_start);
    if(extension == ".gz"){
        std::string before_extension = argument.substr(0, extension_start-1);
        extension_start = before_extension.find_last_of('.');
        extension = argument.substr(extension_start);
    }
    std::string * p;
    switch(type){
        case 'a':
            p = std::find(extensions_a, extensions_a + EXTENSIONS_A_SIZE, extension);
            if(p == extensions_a + EXTENSIONS_A_SIZE) 
                return false;
            else
                return true;
        case 'q':
            p = std::find(extensions_q, extensions_q + EXTENSIONS_Q_SIZE, extension);
            if(p == extensions_q + EXTENSIONS_Q_SIZE)
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
        const char* name, std::uint32_t name_length,
        const char* sequence, std::uint32_t sequence_length) :
        name{std::string(name, name_length)},
        sequence{std::string(sequence, sequence_length)}
        {}

    Fast(
        const char* name, std::uint32_t name_length,
        const char* sequence, std::uint32_t sequence_length,
        const char* quality, std::uint32_t quality_length) :
        name{std::string(name, name_length)},
        sequence{std::string(sequence, sequence_length)},
        quality{std::string(quality, quality_length)}
        {}
};

std::vector<std::unique_ptr<Fast>> parse_fasta (std::string fastaFile){
    std::vector<std::unique_ptr<Fast>> fasta_objects;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(fastaFile);
    // read the whole file
    fasta_parser->parse(fasta_objects, -1);
    
    return fasta_objects;
}

std::vector<std::unique_ptr<Fast>> parse_fastq (std::string fastqFile){
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

void print_statistics (std::string file){
    std::vector<std::unique_ptr<Fast>> fast_objects;
    if(is_extension_ok(file, 'a')){
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
    for(int i = 0; i < num_of_seq; i++){
        sum += (fast_objects[i]->sequence).length();

        if((fast_objects[i]->sequence).length() < min_length)
            min_length = (fast_objects[i]->sequence).length();
        
        if((fast_objects[i]->sequence).length() > max_length)
            max_length = (fast_objects[i]->sequence).length();

    }
    average_length = (float) sum / num_of_seq;


    std::cerr << "File name: " << file << std::endl;
    std::cerr << "Number of sequences: " << num_of_seq << std::endl;
    std::cerr << "Average length: " << average_length << std::endl;
    std::cerr << "Minimal length: " << min_length << std::endl;
    std::cerr << "Maximal length: " << max_length << std::endl;
}

int main(int argc, char* argv[]) {

    int opt;
    std::string alignment_type;
    pink::AlignmentType type = pink::global;
    int match = MATCH;
    int mismatch = MISMATCH;
    int gap = GAP;
    int k = K;
    int w = W;
    double f = F;
    while((opt = getopt(argc, argv, ":hvt:m:s:g:")) != -1)  {
        switch(opt) {
            case 'h':
                help();
                return 0;
            case 'v':
                version();
                return 0;
            case 't':
                alignment_type = optarg;
                if(alignment_type == "global")
                    type = pink::global;
                else if(alignment_type == "semi_global")
                    type = pink::semi_global;
                else if(alignment_type == "local")
                    type = pink::local;
                else
                    errorMessage();
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
            default:
                errorMessage();
        }

    }

    if(argc != optind + 2) {
        errorMessage();
    }

    std::string first = argv[optind];
    std::string second = argv[optind+1];

    if(!(is_extension_ok(first, 'a') || is_extension_ok(first, 'q')) || !is_extension_ok(second, 'a')){
        errorMessage();
    }

    // FIRST

//     print_statistics(first);
//     std::cerr << std::endl;
//     print_statistics(second);


    // SECOND

//    std::vector<std::unique_ptr<Fast>> fast_objects;
//    if(is_extension_ok(first, 'a')){
//        fast_objects = parse_fasta(first);
//    } else {
//        fast_objects = parse_fastq(first);
//    }
//
//    int which1 = std::rand() % fast_objects.size();
//    int which2 = std::rand() % fast_objects.size();
//
//    const char *query = (fast_objects[which1]->sequence).c_str();
//    const char *target = (fast_objects[which2]->sequence).c_str();
//    unsigned int query_length = (fast_objects[which1]->sequence).length();
//    unsigned int target_length = (fast_objects[which2]->sequence).length();
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

    std::vector<std::unique_ptr<Fast>> fast_objects;
    if(is_extension_ok(first, 'a')){
        fast_objects = parse_fasta(first);
    } else {
        fast_objects = parse_fastq(first);
    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;

    for(int i = 0; i < fast_objects.size(); i++){
        std::vector<std::tuple<unsigned int, unsigned int, bool>> current = pink::minimizers((fast_objects[i]->sequence).c_str(), (fast_objects[i]->sequence).length(), k, w);

        for(auto minimizer : current) {
            minimizers.emplace_back(minimizer);
        }
    }

    std::set<std::tuple<unsigned int, unsigned int, bool>> unique (minimizers.begin(), minimizers.end());
    std::map<std::tuple<unsigned int, unsigned int, bool>, int> counter;

    for(auto minimizer : minimizers) {
        if(!counter.emplace(std::make_pair(minimizer, 1)).second)
            counter[minimizer] = counter[minimizer] + 1;
    }

    std::vector<std::pair<std::tuple<unsigned int, unsigned int, bool>, int>> pairs;
    for (auto itr = counter.begin(); itr != counter.end(); ++itr)
        pairs.push_back(*itr);

    sort(pairs.begin(), pairs.end(), [=](std::pair<std::tuple<unsigned int, unsigned int, bool>, int>& a, std::pair<std::tuple<unsigned int, unsigned int, bool>, int>& b)
         {
             return a.second < b.second;
         }
    );

    std::cout << "Number of distinct minimizers is "
              << unique.size()
              << std::endl;

    int num_of_singletons = 0;
    for(auto p : pairs) {
        if(p.second == 1)
            num_of_singletons++;
        else
            break;
    }

    std::cout << "Fraction of singletons: " << (double)num_of_singletons/minimizers.size() << std::endl;

    int ignore = minimizers.size()*f;

    std::cout << "Number of occurrences of the most frequent minimizer when the top " << f << " frequent minimizers: "
              << pairs[minimizers.size() - 1 - ignore].second << std::endl;


    return 0;
}