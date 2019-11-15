#include <getopt.h>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>
#include "../brown_alignment/brown_alignment.hpp"
#include "../vendor/bioparser/include/bioparser/bioparser.hpp"
#include "MapperConfig.h"
#define MAX 1000

using namespace std;

static struct option options[] = {{"version", no_argument, 0, 'v'},
                                  {"help", no_argument, 0, 'h'},
                                  {0, 0, 0, 0}};

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
        fprintf(stdout,
                "Please provide 2 files in FASTA/FASTQ format along with "
                "alignment type and match, mismatch and gap costs\nAlignment "
                "types: 0 - local, 1 - global, 2 - semi_global\n");
        exit(EXIT_SUCCESS);

      default:
        fprintf(stderr, "Unknown options\n");
        exit(EXIT_FAILURE);
    }
  }

  if (argc - optind < 7) {
    fprintf(stderr,
            "Please provide 2 files in FASTA/FASTQ format along with alignment "
            "type and match, mismatch and gap costs\nAlignment types: 0 - "
            "local, 1 - global, 2 - semi_global\n");
    exit(EXIT_FAILURE);
  }
  // parse 1st file as FASTA
  if (has_suffix(convertToString(argv[1], strlen(argv[1])), ".fasta") ||
      has_suffix(convertToString(argv[1], strlen(argv[1])), ".fa") ||
      has_suffix(convertToString(argv[1], strlen(argv[1])), ".fasta.gz") ||
      has_suffix(convertToString(argv[1], strlen(argv[1])), ".fa.gz")) {
    vector<unique_ptr<FASTAfile>> fasta_objects;
    string path = convertToString(argv[1], strlen(argv[1]));
    brown::AlignmentType type;
    if (atoi(argv[3]) == 0) type = brown::AlignmentType::local;
    if (atoi(argv[3]) == 1) type = brown::AlignmentType::global;
    if (atoi(argv[3]) == 2) type = brown::AlignmentType::semi_global;
    int match = atoi(argv[4]);
    int mismatch = atoi(argv[5]);
    int gap = atoi(argv[6]);
    auto fasta_parser =
        bioparser::createParser<bioparser::FastaParser, FASTAfile>(path);
    fasta_parser->parse(fasta_objects, -1);
    print_fasta_stats(fasta_objects);
    srand(time(NULL));
    int rand1 = rand() % fasta_objects.size();
    int rand2 = rand() % fasta_objects.size();
    string &query = fasta_objects[rand1]->sequence;
    string &target = fasta_objects[rand2]->sequence;
    string cigar;
    unsigned int target_begin;
    cout << "Alignment score: "
         << brown::pairwise_alignment(query.c_str(), query.size(),
                                      target.c_str(), target.size(), type,
                                      match, mismatch, gap, cigar, target_begin)
         << '\n';
    cout << "CIGAR string: " << cigar << '\n';

  }

  // parse 1st file as FASTQ
  else if (has_suffix(convertToString(argv[1], strlen(argv[1])), ".fastq") ||
           has_suffix(convertToString(argv[1], strlen(argv[1])), ".fq") ||
           has_suffix(convertToString(argv[1], strlen(argv[1])), ".fastq.gz") ||
           has_suffix(convertToString(argv[1], strlen(argv[1])), ".fq.gz")) {
    std::vector<std::unique_ptr<FASTQfile>> fastq_objects;
    string path = convertToString(argv[1], strlen(argv[1]));
    brown::AlignmentType type;
    if (atoi(argv[3]) == 0) type = brown::AlignmentType::local;
    if (atoi(argv[3]) == 1) type = brown::AlignmentType::global;
    if (atoi(argv[3]) == 2) type = brown::AlignmentType::semi_global;
    int match = atoi(argv[4]);
    int mismatch = atoi(argv[5]);
    int gap = atoi(argv[6]);
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
    srand(time(NULL));
    int rand1 = rand() % fastq_objects.size();
    int rand2 = rand() % fastq_objects.size();
    string &query = fastq_objects[rand1]->sequence;
    string &target = fastq_objects[rand2]->sequence;
    string cigar;
    unsigned int target_begin;
    cout << "Alignment score: "
         << brown::pairwise_alignment(query.c_str(), query.size(),
                                      target.c_str(), target.size(), type,
                                      match, mismatch, gap, cigar, target_begin)
         << '\n';
    cout << "CIGAR string: " << cigar << '\n';
  } else {
    fprintf(stderr,
            "Please provide 2 files in FASTA/FASTQ format along with alignment "
            "type and match, mismatch and gap costs\nAlignment types: 0 - "
            "local, 1 - global, 2 - semi_global\n");
    exit(EXIT_FAILURE);
  }
  // parse 2nd file as FASTA
  if (has_suffix(convertToString(argv[2], strlen(argv[2])), ".fasta") ||
      has_suffix(convertToString(argv[2], strlen(argv[2])), ".fa") ||
      has_suffix(convertToString(argv[2], strlen(argv[2])), ".fasta.gz") ||
      has_suffix(convertToString(argv[2], strlen(argv[2])), ".fa.gz")) {
    vector<unique_ptr<FASTAfile>> fasta_objects;
    string path = convertToString(argv[1], strlen(argv[1]));
    brown::AlignmentType type;
    if (atoi(argv[3]) == 0) type = brown::AlignmentType::local;
    if (atoi(argv[3]) == 1) type = brown::AlignmentType::global;
    if (atoi(argv[3]) == 2) type = brown::AlignmentType::semi_global;
    int match = atoi(argv[4]);
    int mismatch = atoi(argv[5]);
    int gap = atoi(argv[6]);
    auto fasta_parser =
        bioparser::createParser<bioparser::FastaParser, FASTAfile>(path);
    fasta_parser->parse(fasta_objects, -1);
    print_fasta_stats(fasta_objects);
    srand(time(NULL));
    int rand1 = rand() % fasta_objects.size();
    int rand2 = rand() % fasta_objects.size();
    string &query = fasta_objects[rand1]->sequence;
    string &target = fasta_objects[rand2]->sequence;
    string cigar;
    unsigned int target_begin;
    cout << "Alignment score: "
         << brown::pairwise_alignment(query.c_str(), query.size(),
                                      target.c_str(), target.size(), type,
                                      match, mismatch, gap, cigar, target_begin)
         << '\n';
    cout << "CIGAR string: " << cigar << '\n';

  }
  // parse 2nd file as FASTQ
  else if (has_suffix(convertToString(argv[2], strlen(argv[2])), ".fastq") ||
           has_suffix(convertToString(argv[2], strlen(argv[2])), ".fq") ||
           has_suffix(convertToString(argv[2], strlen(argv[2])), ".fastq.gz") ||
           has_suffix(convertToString(argv[2], strlen(argv[2])), ".fq.gz")) {
    brown::AlignmentType type;
    if (atoi(argv[3]) == 0) type = brown::AlignmentType::local;
    if (atoi(argv[3]) == 1) type = brown::AlignmentType::global;
    if (atoi(argv[3]) == 2) type = brown::AlignmentType::semi_global;
    int match = atoi(argv[4]);
    int mismatch = atoi(argv[5]);
    int gap = atoi(argv[6]);
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
    srand(time(NULL));
    int rand1 = rand() % fastq_objects.size();
    int rand2 = rand() % fastq_objects.size();
    string &query = fastq_objects[rand1]->sequence;
    string &target = fastq_objects[rand2]->sequence;
    string cigar;
    unsigned int target_begin;
    cout << "Alignment score: "
         << brown::pairwise_alignment(query.c_str(), query.size(),
                                      target.c_str(), target.size(), type,
                                      match, mismatch, gap, cigar, target_begin)
         << '\n';
    cout << "CIGAR string: " << cigar << '\n';
  } else {
    fprintf(stderr,
            "Please provide 2 files in FASTA/FASTQ format along with alignment "
            "type and match, mismatch and gap costs\nAlignment types: 0 - "
            "local, 1 - global, 2 - semi_global\n");
    exit(EXIT_FAILURE);
  }
}