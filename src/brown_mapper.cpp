#include <cstring>
#include <iostream>
#include <memory>
#include <vector>
#include <getopt.h>
#include "MapperConfig.h"
#include "../vendor/bioparser/include/bioparser/bioparser.hpp"
#define MAX 1000

using namespace std;

static struct option options[] = {
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

class FASTAfile {
 public:
  string name;
  string sequence;

  FASTAfile(const char *name, std::uint32_t name_length, const char *sequence,
            std::uint32_t sequence_length)
      : name(name, name_length), sequence(sequence, sequence_length) {}
};

class FASTQfile {
 public:
  string name;
  string sequence;
  string quality;

  FASTQfile(const char *name, std::uint32_t name_length, const char *sequence,
            std::uint32_t sequence_length, const char *quality,
            std::uint32_t quality_length)
      : name(name, name_length),
        sequence(sequence, sequence_length),
        quality(quality, quality_length) {}
};

void print_fasta_stats(vector<unique_ptr<FASTAfile>> &objects) {
  long min = MAX;
  long max = 0;
  long double avg = 0;
  long sum = 0;

  for (auto &a : objects) {
    auto len = (a->sequence).length();
    sum += len;

    if (len > max) max = len;

    if (len < min) min = len;
  }

  avg = sum / objects.size();

  cout << "Number of sequnces in file: " << objects.size() << "\n";
  cout << "Average sequence length: " << avg << "\n";
  cout << "Minimun sequence length: " << min << "\n";
  cout << "Maximum sequence length: " << max << "\n";
}

void print_fastq_stats(vector<unique_ptr<FASTQfile>> &objects) {
  long min = MAX;
  long max = 0;
  long double avg = 0;
  long sum = 0;

  for (auto &a : objects) {
    auto len = (a->sequence).length();
    sum += len;

    if (len > max) max = len;

    if (len < min) min = len;

    cout << "Quality: " << a->quality << "\n";
  }

  avg = sum / objects.size();

  cout << "Number of sequnces in file: " << objects.size() << "\n";
  cout << "Average sequence length: " << avg << "\n";
  cout << "Minimun sequence length: " << min << "\n";
  cout << "Maximum sequence length: " << max << "\n";
}

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

  bool flag = false;

  char opt;
  while ((opt = getopt_long(argc, argv, "hv", options, NULL)) != -1) {
        switch (opt) {
            case 'v':
                fprintf(stdout, "v%d.%d\n", Mapper_VERSION_MAJOR, Mapper_VERSION_MINOR);               
                exit(EXIT_SUCCESS);   
            case 'h':
                fprintf(stdout, "-v (--version)  Project version\n");
                fprintf(stdout, "-h (--help)     Help\n\n");
                fprintf(stdout, "Please provide 2 files in FASTA/FASTQ format\n");
                exit(EXIT_SUCCESS);   

            default:
                fprintf(stderr, "Unknown options\n");
                exit(EXIT_FAILURE);   
        }
    }

  if (argc - optind < 2) {
      fprintf(stderr, "2 files needed\n");
      exit(EXIT_FAILURE);   

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
      vector<unique_ptr<FASTAfile>> fasta_objects;
      string path = convertToString(argv[1], strlen(argv[1]));
      auto fasta_parser =
          bioparser::createParser<bioparser::FastaParser, FASTAfile>(path);
      fasta_parser->parse(fasta_objects, -1);

      print_fasta_stats(fasta_objects);

      if (strcmp(argv[1], argv[2]) == 0) flag = true;
    }
    // parse 2. file as FASTA
    if (has_suffix(convertToString(argv[2], strlen(argv[2])), ".fasta") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fa") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fasta.gz") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fa.gz")) {

      if (flag == false) {
        vector<unique_ptr<FASTAfile>> fasta_objects;
        string path = convertToString(argv[1], strlen(argv[1]));
        auto fasta_parser =
            bioparser::createParser<bioparser::FastaParser, FASTAfile>(path);
        fasta_parser->parse(fasta_objects, -1);

        print_fasta_stats(fasta_objects);
      }

    }
    // parse 1st file as FASTQ
    if (has_suffix(convertToString(argv[1], strlen(argv[1])), ".fastq") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fq") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fastq.gz") ||
        has_suffix(convertToString(argv[1], strlen(argv[1])), ".fq.gz")) {
      std::vector<std::unique_ptr<FASTQfile>> fastq_objects;
      string path = convertToString(argv[1], strlen(argv[1]));
      auto fastq_parser =
          bioparser::createParser<bioparser::FastqParser, FASTQfile>(path);
      std::uint64_t size_in_bytes = 500 * 1024 * 1024;
      while (true) {
        auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
        if (status == false) {
          break;
        }
      }
      print_fastq_stats(fastq_objects);

      if (argv[1] == argv[2]) flag = true;
    }
    // parse 2nd file as FASTQ
    if (has_suffix(convertToString(argv[2], strlen(argv[2])), ".fastq") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fq") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fastq.gz") ||
        has_suffix(convertToString(argv[2], strlen(argv[2])), ".fq.gz")) {

      if (flag == false) {
        std::vector<std::unique_ptr<FASTQfile>> fastq_objects;
        string path = convertToString(argv[2], strlen(argv[2]));
        auto fastq_parser =
            bioparser::createParser<bioparser::FastqParser, FASTQfile>(path);
        std::uint64_t size_in_bytes = 500 * 1024 * 1024;
        while (true) {
          auto status = fastq_parser->parse(fastq_objects, size_in_bytes);
          if (status == false) {
            break;
          }
        }
        print_fastq_stats(fastq_objects);
      }
    }
  } else {
    fprintf(stderr, "2 files needed of format FASTA/FASTQ\n");
    exit(EXIT_FAILURE);
  }
}