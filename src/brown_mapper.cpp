#include <getopt.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <vector>
#include "../brown_alignment/brown_alignment.hpp"
#include "../brown_minimizers/brown_minimizers.hpp"
#include "../vendor/bioparser/include/bioparser/bioparser.hpp"
#include "MapperConfig.h"
#include "omp.h"
#define MAX 1000

using namespace std;

static struct option options[] = {
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {"match", required_argument, 0, 'm'},
    {"mismatch", required_argument, 0, 'x'},
    {"gap", required_argument, 0, 'g'},
    {"k", required_argument, 0, 'k'},
    {"window_length", required_argument, 0, 'w'},
    {"top minimizers not taken into account", required_argument, 0, 'f'},
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

  cout << "Number of sequences in file: " << objects.size() << "\n";
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

  cout << "Number of sequences in file: " << objects.size() << "\n";
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

// binary search
int GetCeilIndex(vector<pair<int, int>> &matches, vector<int> &T, int l, int r,
                 int key) {
  while (r - l > 1) {
    int m = l + (r - l) / 2;
    if (matches[T[m]].second >= key)
      r = m;
    else
      l = m;
  }

  return r;
}

// LIS algorithm (nlogn) but modified to suit needs
// https://www.geeksforgeeks.org/construction-of-longest-monotonically-increasing-subsequence-n-log-n/
int LongestIncreasingSubsequence(vector<pair<int, int>> &matches, int n,
                                 int &t_begin, int &t_end, int &q_begin,
                                 int &q_end) {
  vector<int> tailIndices(n, 0);
  vector<int> prevIndices(n, -1);

  int len = 1;
  for (int i = 1; i < n; i++) {
    if (matches[i].second < matches[tailIndices[0]].second) {
      // new smallest value
      tailIndices[0] = i;
    } else if (matches[i].second > matches[tailIndices[len - 1]].second) {
      // matches[i].second wants to extend largest subsequence
      prevIndices[i] = tailIndices[len - 1];
      tailIndices[len++] = i;
    } else {
      // matches[i].second wants to be a potential condidate of
      // future subsequence
      // It will replace ceil value in tailIndices
      int pos =
          GetCeilIndex(matches, tailIndices, -1, len - 1, matches[i].second);

      prevIndices[i] = tailIndices[pos - 1];
      tailIndices[pos] = i;
    }
  }
  t_begin = matches[prevIndices[1]].second;
  q_begin = matches[prevIndices[1]].first;
  t_end = matches[tailIndices[len - 1]].second;
  q_end = matches[tailIndices[len - 1]].first;

  return len;
}

int main(int argc, char **argv) {
  bool flag = false;
  int match = 0, mismatch = 0, gap = 0, k = 15, w = 5;
  float f = 0.001;
  int size1;
  int size2;
  int t_begin1, t_end1, q_begin1, q_end1;
  int t_begin2, t_end2, q_begin2, q_end2;
  string query1;
  string query2;
  string queryReference;
  vector<int> frequentRef;
  vector<int> frequent1;
  vector<int> frequent2;
  vector<tuple<unsigned int, unsigned int, bool>> minimizersListRef;
  string PAFPathQuery;
  string PAFPathTarget1;
  string PAFPathTarget2;

  char opt;
  while ((opt = getopt_long(argc, argv, "hvm:x:g:k:w:f:", options, NULL)) !=
         -1) {
    switch (opt) {
      case 0:
        break;
      case 'v':
        fprintf(stdout, "v%d.%d\n", Mapper_VERSION_MAJOR, Mapper_VERSION_MINOR);
        exit(EXIT_SUCCESS);
      case 'h':
        fprintf(stdout, "-v (--version)  Project version\n");
        fprintf(stdout, "-h (--help)     Help\n\n");
        fprintf(stdout,
                "Please provide 2 files in FASTA/FASTQ format along with "
                "alignment type and match, mismatch and gap costs, k, window "
                "width and number of top frequent minimizers not taken into "
                "acount\nAlignment "
                "types: 0 - local, 1 - global, 2 - semi_global\n");
        exit(EXIT_SUCCESS);
      case 'm':
        match = stoi(optarg);
        break;
      case 'x':
        mismatch = stoi(optarg);
        break;
      case 'g':
        gap = stoi(optarg);
        break;
      case 'k':
        k = stoi(optarg);
        break;
      case 'w':
        w = stoi(optarg);
        break;
      case 'f':
        f = stoi(optarg);
        break;
      default:
        fprintf(stderr, "Unknown options\n");
        exit(EXIT_FAILURE);
    }
  }

  if (argc != 17) {
    fprintf(
        stderr,
        "Please provide 2 files in FASTA/FASTQ format along with "
        "alignment type and match, mismatch and gap costs, k, window width and "
        "number of top frequent minimizers not taken into acount\nAlignment "
        "types: 0 - local, 1 - global, 2 - semi_global\n");
    exit(EXIT_FAILURE);
  }
  // reference genome parsing
  if (has_suffix(convertToString(argv[16], strlen(argv[16])), ".fa")) {
    vector<unique_ptr<FASTAfile>> fasta_objects;
    string pathFa = convertToString(argv[16], strlen(argv[16]));
    PAFPathQuery = pathFa;
    auto fasta_parser =
        bioparser::createParser<bioparser::FastaParser, FASTAfile>(pathFa);
    fasta_parser->parse(fasta_objects, -1);
    string &query = fasta_objects[0]->sequence;
    queryReference = query;
    minimizersListRef = brown::minimizers(query.c_str(), query.size(), k, w);
    for (int i = 0; i < minimizersListRef.size(); i++) {
      get<1>(minimizersListRef[i]) = i;
    }

    map<int, int> printListRef;
    for (int i = 0; i < minimizersListRef.size(); i++) {
      printListRef[get<0>(minimizersListRef[i])]++;
    }
    vector<pair<int, int>> kulRef;
    copy(printListRef.begin(), printListRef.end(),
         back_inserter<vector<pair<int, int>>>(kulRef));
    sort(kulRef.begin(), kulRef.end(),
         [](const pair<int, int> &l, const pair<int, int> &r) {
           if (l.second != r.second) return l.second > r.second;

           return l.first > r.first;
         });
    for (int i = 0; i < f; i++) {
      frequentRef.push_back(kulRef[i].first);
    }
  } else {
    fprintf(
        stderr,
        "Please provide 2 files in FASTA/FASTQ format along with "
        "alignment type and match, mismatch and gap costs, k, window width and "
        "number of top frequent minimizers not taken into acount\nAlignment "
        "types: 0 - local, 1 - global, 2 - semi_global\n");
    exit(EXIT_FAILURE);
  }

  // parse 1st file as FASTA
  // oWo
  if (has_suffix(convertToString(argv[13], strlen(argv[13])), ".fasta") ||
      has_suffix(convertToString(argv[13], strlen(argv[13])), ".fa") ||
      has_suffix(convertToString(argv[13], strlen(argv[13])), ".fasta.gz") ||
      has_suffix(convertToString(argv[13], strlen(argv[13])), ".fa.gz")) {
    vector<unique_ptr<FASTAfile>> fasta_objects;
    string path = convertToString(argv[13], strlen(argv[13]));
    PAFPathTarget1 = path;
    brown::AlignmentType type;
    if (atoi(argv[15]) == 0) type = brown::AlignmentType::local;
    if (atoi(argv[15]) == 1) type = brown::AlignmentType::global;
    if (atoi(argv[15]) == 2) type = brown::AlignmentType::semi_global;
    auto fasta_parser =
        bioparser::createParser<bioparser::FastaParser, FASTAfile>(path);
    fasta_parser->parse(fasta_objects, -1);
    print_fasta_stats(fasta_objects);
    srand(time(NULL));
    int rand1 = rand() % fasta_objects.size();
    int rand2 = rand() % fasta_objects.size();
    string &query = fasta_objects[rand1]->sequence;
    string &target = fasta_objects[rand2]->sequence;
    query1 = query;
    string cigar;
    unsigned int target_begin;
    cout << "Alignment score: "
         << brown::pairwise_alignment(query.c_str(), query.size(),
                                      target.c_str(), target.size(), type,
                                      match, mismatch, gap, cigar, target_begin)
         << '\n';
    cout << "CIGAR string: " << cigar << '\n';
    cout << "Beginning: " << target_begin << "\n";
    vector<tuple<unsigned int, unsigned int, bool>> minimizersList;
    minimizersList = brown::minimizers(query.c_str(), query.size(), k, w);
    for (int i = 0; i < minimizersList.size(); i++) {
      get<1>(minimizersList[i]) = i;
    }
    map<int, int> printList;
    for (int i = 0; i < minimizersList.size(); i++) {
      printList[get<0>(minimizersList[i])]++;
    }
    int cnt = 0;
    int singleton = 0;
    cout << "Number of distinct minimizers: " << printList.size() << "\n";
    for (map<int, int>::iterator itr = printList.begin();
         itr != printList.end(); ++itr) {
      if (itr->second == 1) singleton++;
    }
    cout << "Fraction of singletons: "
         << (float)singleton / (float)printList.size() << "\n";
    vector<pair<int, int>> kul;
    copy(printList.begin(), printList.end(),
         back_inserter<vector<pair<int, int>>>(kul));
    sort(kul.begin(), kul.end(),
         [](const pair<int, int> &l, const pair<int, int> &r) {
           if (l.second != r.second) return l.second > r.second;

           return l.first > r.first;
         });

    cout << "Number of occurences of the most frequent minimizer without f "
            "most frequent: "
         << kul[f].second << "\n";
    vector<pair<int, int>> matches;
    bool flag = false;

    for (int i = 0; i < minimizersListRef.size(); i++) {
      flag = false;
      for (int k = 0; k < frequentRef.size(); k++) {
        if (get<0>(minimizersListRef[i]) == frequentRef[k]) flag = true;
      }
      if (flag == true) continue;
      for (int j = 0; j < minimizersList.size(); j++) {
        if (get<0>(minimizersListRef[i]) == get<0>(minimizersList[j]))
          matches.push_back(make_pair(i, j));
      }
    }

    size1 = LongestIncreasingSubsequence(matches, matches.size(), t_begin1,
                                         t_end1, q_begin1, q_end1);
    cout << "Fragment start and end index: " << t_begin1 << " - " << t_end1
         << "\n";
    cout << "Reference start and end index: " << q_begin1 << " - " << q_end1
         << "\n";
    cout << "size: " << size1 << "\n";
    cout << "\n";
  }

  // parse 1st file as FASTQ
  // :3
  // <3
  else if (has_suffix(convertToString(argv[13], strlen(argv[13])), ".fastq") ||
           has_suffix(convertToString(argv[13], strlen(argv[13])), ".fq") ||
           has_suffix(convertToString(argv[13], strlen(argv[13])),
                      ".fastq.gz") ||
           has_suffix(convertToString(argv[13], strlen(argv[13])), ".fq.gz")) {
    std::vector<std::unique_ptr<FASTQfile>> fastq_objects;
    string path = convertToString(argv[13], strlen(argv[13]));
    brown::AlignmentType type;
    if (atoi(argv[15]) == 0) type = brown::AlignmentType::local;
    if (atoi(argv[15]) == 1) type = brown::AlignmentType::global;
    if (atoi(argv[15]) == 2) type = brown::AlignmentType::semi_global;
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
    cout << "Beginning: " << target_begin << "\n";
    vector<tuple<unsigned int, unsigned int, bool>> minimizersList;
    minimizersList = brown::minimizers(query.c_str(), query.size(), k, w);
    for (int i = 0; i < minimizersList.size(); i++) {
      get<1>(minimizersList[i]) = i;
    }
    map<int, int> printList;
    for (int i = 0; i < minimizersList.size(); i++) {
      printList[get<0>(minimizersList[i])]++;
    }
    int cnt = 0;
    int singleton = 0;
    cout << "Number of distinct minimizers: " << printList.size() << "\n";
    for (map<int, int>::iterator itr = printList.begin();
         itr != printList.end(); ++itr) {
      if (itr->second == 1) singleton++;
    }
    cout << "Fraction of singletons: "
         << (float)singleton / (float)printList.size() << "\n";
    vector<pair<int, int>> kul;
    copy(printList.begin(), printList.end(),
         back_inserter<vector<pair<int, int>>>(kul));
    sort(kul.begin(), kul.end(),
         [](const pair<int, int> &l, const pair<int, int> &r) {
           if (l.second != r.second) return l.second > r.second;

           return l.first > r.first;
         });
    for (int i = 0; i < f; i++) {
      frequent1.push_back(kul[i].first);
    }

    cout << "Number of occurences of the most frequent minimizer without f "
            "most frequent: "
         << kul[f].second << "\n";
    cout << "\n";
  } else {
    fprintf(
        stderr,
        "Please provide 2 files in FASTA/FASTQ format along with "
        "alignment type and match, mismatch and gap costs, k, window width and "
        "number of top frequent minimizers not taken into acount\nAlignment "
        "types: 0 - local, 1 - global, 2 - semi_global\n");
    exit(EXIT_FAILURE);
  }
  // parse 2nd file as FASTA
  if (has_suffix(convertToString(argv[14], strlen(argv[14])), ".fasta") ||
      has_suffix(convertToString(argv[14], strlen(argv[14])), ".fa") ||
      has_suffix(convertToString(argv[14], strlen(argv[14])), ".fasta.gz") ||
      has_suffix(convertToString(argv[14], strlen(argv[14])), ".fa.gz")) {
    vector<unique_ptr<FASTAfile>> fasta_objects;
    string path = convertToString(argv[14], strlen(argv[14]));
    PAFPathTarget2 = path;
    brown::AlignmentType type;
    if (atoi(argv[15]) == 0) type = brown::AlignmentType::local;
    if (atoi(argv[15]) == 1) type = brown::AlignmentType::global;
    if (atoi(argv[15]) == 2) type = brown::AlignmentType::semi_global;
    auto fasta_parser =
        bioparser::createParser<bioparser::FastaParser, FASTAfile>(path);
    fasta_parser->parse(fasta_objects, -1);
    print_fasta_stats(fasta_objects);
    srand(time(NULL));
    int rand1 = rand() % fasta_objects.size();
    int rand2 = rand() % fasta_objects.size();
    string &query = fasta_objects[rand1]->sequence;
    string &target = fasta_objects[rand2]->sequence;
    query2 = query;
    string cigar;
    unsigned int target_begin;
    cout << "Alignment score: "
         << brown::pairwise_alignment(query.c_str(), query.size(),
                                      target.c_str(), target.size(), type,
                                      match, mismatch, gap, cigar, target_begin)
         << '\n';
    cout << "CIGAR string: " << cigar << '\n';
    cout << "Beginning: " << target_begin << "\n";
    vector<tuple<unsigned int, unsigned int, bool>> minimizersList;
    minimizersList = brown::minimizers(query.c_str(), query.size(), k, w);
    for (int i = 0; i < minimizersList.size(); i++) {
      get<1>(minimizersList[i]) = i;
    }
    map<int, int> printList;
    for (int i = 0; i < minimizersList.size(); i++) {
      printList[get<0>(minimizersList[i])]++;
    }
    int cnt = 0;
    int singleton = 0;
    cout << "Number of distinct minimizers: " << printList.size() << "\n";
    for (map<int, int>::iterator itr = printList.begin();
         itr != printList.end(); ++itr) {
      if (itr->second == 1) singleton++;
    }
    cout << "Fraction of singletons: "
         << (float)singleton / (float)printList.size() << "\n";
    vector<pair<int, int>> kul;
    copy(printList.begin(), printList.end(),
         back_inserter<vector<pair<int, int>>>(kul));
    sort(kul.begin(), kul.end(),
         [](const pair<int, int> &l, const pair<int, int> &r) {
           if (l.second != r.second) return l.second > r.second;

           return l.first > r.first;
         });

    cout << "Number of occurences of the most frequent minimizer without f "
            "most frequent: "
         << kul[f].second << "\n";
    vector<pair<int, int>> matches;
    bool flag = false;

    for (int i = 0; i < minimizersListRef.size(); i++) {
      flag = false;
      for (int k = 0; k < frequentRef.size(); k++) {
        if (get<0>(minimizersListRef[i]) == frequentRef[k]) flag = true;
      }
      if (flag == true) continue;
      for (int j = 0; j < minimizersList.size(); j++) {
        if (get<0>(minimizersListRef[i]) == get<0>(minimizersList[j]))
          matches.push_back(make_pair(i, j));
      }
    }
    size2 = LongestIncreasingSubsequence(matches, matches.size(), t_begin2,
                                         t_end2, q_begin2, q_end2);
    cout << "Fragment start and end index: " << t_begin2 << " - " << t_end2
         << "\n";
    cout << "Reference start and end index: " << q_begin2 << " - " << q_end2
         << "\n";
    cout << "size: " << size2 << "\n";
    cout << "\n";
  }
  // parse 2nd file as FASTQ
  else if (has_suffix(convertToString(argv[14], strlen(argv[14])), ".fastq") ||
           has_suffix(convertToString(argv[14], strlen(argv[14])), ".fq") ||
           has_suffix(convertToString(argv[14], strlen(argv[14])),
                      ".fastq.gz") ||
           has_suffix(convertToString(argv[14], strlen(argv[14])), ".fq.gz")) {
    brown::AlignmentType type;
    if (atoi(argv[15]) == 0) type = brown::AlignmentType::local;
    if (atoi(argv[15]) == 1) type = brown::AlignmentType::global;
    if (atoi(argv[15]) == 2) type = brown::AlignmentType::semi_global;
    std::vector<std::unique_ptr<FASTQfile>> fastq_objects;
    string path = convertToString(argv[14], strlen(argv[14]));
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
    cout << "Beginning: " << target_begin << "\n";
  } else {
    fprintf(
        stderr,
        "Please provide 2 files in FASTA/FASTQ format along with "
        "alignment type and match, mismatch and gap costs, k, window width and "
        "number of top frequent minimizers not taken into acount\nAlignment "
        "types: 0 - local, 1 - global, 2 - semi_global\n");
    exit(EXIT_FAILURE);
  }
  int size;
  int t_begin, t_end, q_begin, q_end;
  string cigar;
  unsigned int target_begin;

  if (size1 > size2) {
    size = size1;
    t_begin = t_begin1;
    t_end = t_end1;
    q_begin = q_begin1;
    q_end = q_end1;
    queryReference = queryReference.substr(q_begin, size);
    query1 = query1.substr(t_begin, size);
    cout << PAFPathQuery << "\n";
    cout << queryReference.size() << "\n";
    cout << q_begin << "\n";
    cout << q_end << "\n";
    cout << "+"
         << "\n";
    cout << PAFPathTarget1 << "\n";
    cout << query1.size() << "\n";
    cout << t_begin << "\n";
    cout << t_end << "\n";
    brown::pairwise_alignment(queryReference.c_str(), queryReference.size(),
                              query1.c_str(), query1.size(),
                              brown::AlignmentType::local, match, mismatch, gap,
                              cigar, target_begin);
    string number;
    int noOfMatches = 0;
    int blockLength = 0;

    for (char c : cigar) {
      if (isdigit(c)) {
        number += c;
      } else if (c == '=') {
        noOfMatches += stoi(number);
        blockLength += stoi(number);
        number = "";
      } else {
        blockLength += stoi(number);
        number = "";
      }
    }
    cout << noOfMatches << "\n";
    cout << blockLength << "\n";
    cout << "empty"
         << "\n";
    cout << "cg:Z:" << cigar << "\n";
  } else {
    size = size2;
    t_begin = t_begin2;
    t_end = t_end2;
    q_begin = q_begin2;
    q_end = q_end2;
    queryReference = queryReference.substr(q_begin, size);
    query2 = query2.substr(t_begin, size);
    cout << PAFPathQuery << "\n";
    cout << queryReference.size() << "\n";
    cout << q_begin << "\n";
    cout << q_end << "\n";
    cout << "+"
         << "\n";
    cout << PAFPathTarget2 << "\n";
    cout << query2.size() << "\n";
    cout << t_begin << "\n";
    cout << t_end << "\n";
    brown::pairwise_alignment(queryReference.c_str(), queryReference.size(),
                              query2.c_str(), query2.size(),
                              brown::AlignmentType::local, match, mismatch, gap,
                              cigar, target_begin);
    string number;
    int noOfMatches = 0;
    int blockLength = 0;

    for (char c : cigar) {
      if (isdigit(c)) {
        number += c;
      } else if (c == '=') {
        noOfMatches += stoi(number);
        blockLength += stoi(number);
        number = "";
      } else {
        blockLength += stoi(number);
        number = "";
      }
    }
    cout << noOfMatches << "\n";
    cout << blockLength << "\n";
    cout << "empty"
         << "\n";
    cout << "cg:Z:" << cigar << "\n";
  }
}