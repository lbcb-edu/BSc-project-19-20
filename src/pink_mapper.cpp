#include<iostream>
#include<algorithm>
#include<unistd.h>

#include "bioparser/bioparser.hpp"

#define VERSION "v0.1.0"

#define EXTENSIONS_A_SIZE 4
#define EXTENSIONS_Q_SIZE 4

std::string extensions_a[] {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
std::string extensions_q[] { ".fastq", ".fq", ".fastq.gz", ".fq.gz"};

void help() {
    std::cout << "This is help option." << std::endl
              << "This program has following options:\n" << std::endl
              << "    -h   --help\n       prints the help menu\n" << std::endl
              << "    -v   --version\n        prints the version number\n" << std::endl;
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
    while((opt = getopt(argc, argv, ":hv")) != -1)  {
        switch(opt) {
            case 'h':
                help();
                return 0;
            case 'v':
                version();
                return 0;
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

    print_statistics(first);
    std::cerr << std::endl;
    print_statistics(second);

    return 0;
}