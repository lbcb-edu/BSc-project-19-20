#include "MapperConfig.h"
#include "bioparser/include/bioparser/bioparser.hpp"
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>
using namespace std;

class Example1 {
public:
  Example1(const char *name, std::uint32_t name_length, const char *sequence,
           std::uint32_t sequence_length) {}
};

class Example2 {
public:
  Example2(const char *name, std::uint32_t name_length, const char *sequence,
           std::uint32_t sequence_length, const char *quality,
           std::uint32_t quality_length) {}
};

// check for given sufix
bool has_suffix(const string &str, const string &suffix) {
  return str.size() >= suffix.size() &&
         str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// covert char array to string
string convertToString(char *a, int size) {
  int i;
  string s = "";
  for (i = 0; i < size; i++) {
    s = s + a[i];
  }
  return s;
}

int main(int argc, char **argv) {

  // check for command line args
  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    fprintf(stdout, "-v (--version)  Project version\n");
    fprintf(stdout, "-h (--help)     Help\n\n");
    fprintf(stdout, "Please provide 2 files in FASTA/FASTQ format\n");
  } else if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
    fprintf(stdout, "v%d.%d\n", Mapper_VERSION_MAJOR, Mapper_VERSION_MINOR);
  } else if (argc < 3) {
    fprintf(stderr, "2 files needed\n");
    exit(EXIT_FAILURE);
    // check for file validity
  } else if ((has_suffix(convertToString(argv[1], strlen(argv[1])), ".fasta") ||
              has_suffix(convertToString(argv[1], strlen(argv[1])), ".fastq") ||
              has_suffix(convertToString(argv[1], strlen(argv[1])), ".fa") ||
              has_suffix(convertToString(argv[1], strlen(argv[1])), ".fq") ||
              has_suffix(convertToString(argv[1], strlen(argv[1])),
                         ".fasta.gz") ||
              has_suffix(convertToString(argv[1], strlen(argv[1])), ".fa.gz") ||
              has_suffix(convertToString(argv[1], strlen(argv[1])),
                         ".fastq.gz") ||
              has_suffix(convertToString(argv[1], strlen(argv[1])),
                         ".fq.gz")) &&
             (has_suffix(convertToString(argv[2], strlen(argv[2])), ".fasta") ||
              has_suffix(convertToString(argv[2], strlen(argv[2])), ".fastq") ||
              has_suffix(convertToString(argv[2], strlen(argv[2])), ".fa") ||
              has_suffix(convertToString(argv[2], strlen(argv[2])), ".fq") ||
              has_suffix(convertToString(argv[2], strlen(argv[2])),
                         ".fasta.gz") ||
              has_suffix(convertToString(argv[2], strlen(argv[2])), ".fa.gz") ||
              has_suffix(convertToString(argv[2], strlen(argv[2])),
                         ".fastq.gz") ||
              has_suffix(convertToString(argv[2], strlen(argv[2])),
                         ".fq.gz"))) {
    // parse 1st file as FASTA
    if (has_suffix(convertToString(argv[1], strlen(argv[1])), ".fasta") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fa") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fasta.gz") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fa.gz")) {

      vector<unique_ptr<Example1>> fasta_objects;
      string str = convertToString(argv[1], strlen(argv[1]));
      string path = "/home/filip/git/Project/" + str;
      auto fasta_parser =
          bioparser::createParser<bioparser::FastaParser, Example1>(path);
      fasta_parser->parse(fasta_objects, -1);

      cout << "kul1"
           << "\n";
    }
    // parse 2. file as FASTA
    if (has_suffix(convertToString(argv[2], strlen(argv[2])), ".fasta") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fa") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fasta.gz") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fa.gz")) {

      vector<unique_ptr<Example1>> fasta_objects;
      string str = convertToString(argv[1], strlen(argv[1]));
      string path = "/home/filip/git/Project/" + str;
      auto fasta_parser =
          bioparser::createParser<bioparser::FastaParser, Example1>(path);
      fasta_parser->parse(fasta_objects, -1);

      cout << "kul2"
           << "\n";
    }
    // parse 1st file as FASTA
    if (has_suffix(convertToString(argv[1], strlen(argv[1])), ".fastq") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fq") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fastq.gz") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fq.gz")) {

      std::vector<std::unique_ptr<Example2>> fastq_objects;
      string str = convertToString(argv[1], strlen(argv[1]));
      string path = "/home/filip/git/Project/" + str;
      auto fastq_parser =
          bioparser::createParser<bioparser::FastqParser, Example2>(path);
      std::uint64_t size_in_bytes = 500 * 1024 * 1024;
      while (true) {
        auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
        if (status == false) {
          break;
        }
      }

      cout << "kul3"
           << "\n";
    }
    // parse 2nd file as FASTQ
    if (has_suffix(convertToString(argv[2], strlen(argv[2])), ".fastq") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fq") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fastq.gz") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fq.gz")) {

      std::vector<std::unique_ptr<Example2>> fastq_objects;
      string str = convertToString(argv[2], strlen(argv[2]));
      string path = "/home/filip/git/Project/" + str;
      auto fastq_parser =
          bioparser::createParser<bioparser::FastqParser, Example2>(path);
      std::uint64_t size_in_bytes = 500 * 1024 * 1024;
      while (true) {
        auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
        if (status == false) {
          break;
        }
      }

      cout << "kul4"
           << "\n";
    }
  } else {
    fprintf(stderr, "2 files needed of format FASTA/FASTQ\n");
    exit(EXIT_FAILURE);
  }
}