#include "brown_alignment.hpp"
#include <iostream>
#include <vector>

using namespace std;

namespace brown {

enum class valueType { deletion, insertion, match, mismatch, none };

struct matrixCell {
  int value;
  valueType cost;
};

void do_cigar(vector<vector<matrixCell>>& matrix, string& cigar, int& row,
              int& column, int& cntD, int& cntI, int& cntM, int& cntX) {
  if (matrix[row][column].cost == valueType::deletion) {
    if (cntI != 0) {
      string a = to_string(cntI);
      cigar.append(a);
      cigar.push_back('I');
      cntI = 0;
    }
    if (cntM != 0) {
      string a = to_string(cntM);
      cigar.append(a);
      cigar.push_back('=');
      cntM = 0;
    }
    if (cntX != 0) {
      string a = to_string(cntX);
      cigar.append(a);
      cigar.push_back('X');
      cntX = 0;
    }
    cntD++;
    row--;
  } else if (matrix[row][column].cost == valueType::insertion) {
    if (cntD != 0) {
      string a = to_string(cntD);
      cigar.append(a);
      cigar.push_back('D');
      cntD = 0;
    }
    if (cntM != 0) {
      string a = to_string(cntM);
      cigar.append(a);
      cigar.push_back('=');
      cntM = 0;
    }
    if (cntX != 0) {
      string a = to_string(cntX);
      cigar.append(a);
      cigar.push_back('X');
      cntX = 0;
    }
    cntI++;
    column--;
  } else if (matrix[row][column].cost == valueType::match) {
    if (cntI != 0) {
      cout << "bla";
      string a = to_string(cntI);
      cigar.append(a);
      cigar.push_back('I');
      cntI = 0;
    }
    if (cntD != 0) {
      string a = to_string(cntD);
      cigar.append(a);
      cigar.push_back('D');
      cntD = 0;
    }
    if (cntX != 0) {
      string a = to_string(cntX);
      cigar.append(a);
      cigar.push_back('X');
      cntX = 0;
    }
    cntM++;
    row--;
    column--;
  } else if (matrix[row][column].cost == valueType::mismatch) {
    if (cntI != 0) {
      cout << "bla";
      string a = to_string(cntI);
      cigar.append(a);
      cigar.push_back('I');
      cntI = 0;
    }
    if (cntM != 0) {
      string a = to_string(cntM);
      cigar.append(a);
      cigar.push_back('=');
      cntM = 0;
    }
    if (cntD != 0) {
      string a = to_string(cntD);
      cigar.append(a);
      cigar.push_back('D');
      cntD = 0;
    }
    cntX++;
    row--;
    column--;
  }
}

void create_cigar(vector<vector<matrixCell>>& matrix, AlignmentType type,
                  string& cigar, unsigned int& target_begin, int row,
                  int column, const char* query, unsigned int query_length,
                  const char* target) {
  if (row < query_length) {
    string a = to_string(query_length - row);
    cigar.append(a);
    cigar.push_back('S');
  }
  int cntD = 0, cntI = 0, cntM = 0, cntX = 0;
  if (type == AlignmentType::local) {
    while (matrix[row][column].value != 0) {
      do_cigar(matrix, cigar, row, column, cntD, cntI, cntM, cntX);
    }
  } else {
    while (matrix[row][column].cost != valueType::none) {
      do_cigar(matrix, cigar, row, column, cntD, cntI, cntM, cntX);
    }
  }
  target_begin = column;
}

void align(const char* query, unsigned int query_length, const char* target,
           unsigned int target_length, AlignmentType type, int match,
           int mismatch, int gap, vector<vector<matrixCell>>& matrix) {
  matrix[0][0].cost = valueType::none;
  matrix[0][0].value = 0;

  if (type == AlignmentType::global) {
    for (int i = 1; i < query_length + 1; i++) {
      matrix[i][0].value = gap * i;
      matrix[i][0].cost = valueType::deletion;
    }
    for (int i = 1; i < target_length + 1; i++) {
      matrix[0][i].value = gap * i;
      matrix[0][i].cost = valueType::insertion;
    }
  } else if (type == AlignmentType::semi_global) {
    for (int i = 1; i < query_length + 1; i++) {
      matrix[i][0].value = gap * i;
      matrix[i][0].cost = valueType::deletion;
    }
    for (int i = 1; i < target_length + 1; i++) {
      matrix[0][i].value = 0;
      matrix[0][i].cost = valueType::none;
    }
  } else if (type == AlignmentType::local) {
    for (int i = 1; i < target_length + 1; i++) {
      matrix[0][i].value = 0;
      matrix[0][i].cost = valueType::none;
    }
    for (int i = 1; i < query_length + 1; i++) {
      matrix[i][0].value = 0;
      matrix[i][0].cost = valueType::none;
    }
  }
  for (int i = 1; i < query_length + 1; i++) {
    for (int j = 1; j < target_length + 1; j++) {
      int matchPrev;
      bool matchCheck;
      if (query[i - 1] == target[i - 1]) {
        matchPrev = matrix[i - 1][j - 1].value + match;
        matchCheck = true;
        // cout << "kul";
      } else {
        matchPrev = matrix[i - 1][j - 1].value + mismatch;
        matchCheck = false;
      }
      int insertion = matrix[i][j - 1].value + gap;
      int deletion = matrix[i - 1][j].value + gap;

      if (matchPrev > insertion) {
        if (deletion > matchPrev) {
          matrix[i][j].value = deletion;
          matrix[i][j].cost = valueType::deletion;
        } else {
          // cout << "dostadobro";
          matrix[i][j].value = matchPrev;
          if (matchCheck == true) {
            matrix[i][j].cost = valueType::match;
            // cout << "nekaj";
          } else {
            matrix[i][j].cost = valueType::mismatch;
            // cout << "damn";
          }
        }
      } else {
        if (deletion > insertion) {
          matrix[i][j].value = deletion;
          matrix[i][j].cost = valueType::deletion;
        } else {
          matrix[i][j].value = insertion;
          matrix[i][j].cost = valueType::insertion;
        }
      }
    }
  }
}

int global_alignment(const char* query, unsigned int query_length,
                     const char* target, unsigned int target_length,
                     AlignmentType type, int match, int mismatch, int gap,
                     string& cigar, unsigned int& target_begin) {
  vector<vector<matrixCell>> matrix(query_length + 1,
                                    vector<matrixCell>(target_length + 1));

  align(query, query_length, target, target_length, type, match, mismatch, gap,
        matrix);
  create_cigar(matrix, type, cigar, target_begin, query_length, target_length,
               query, query_length, target);
  return matrix[query_length][target_length].value;
}

int local_alignment(const char* query, unsigned int query_length,
                    const char* target, unsigned int target_length,
                    AlignmentType type, int match, int mismatch, int gap,
                    string& cigar, unsigned int& target_begin) {
  vector<vector<matrixCell>> matrix(query_length + 1,
                                    vector<matrixCell>(target_length + 1));
  align(query, query_length, target, target_length, type, match, mismatch, gap,
        matrix);
  int maximum = 0;
  int max1 = 0;
  int max2 = 0;
  for (int i = 1; i < query_length + 1; i++) {
    for (int j = 1; j < target_length + 1; j++) {
      if (matrix[i][j].value > maximum) {
        maximum = matrix[i][j].value;
        max1 = i;
        max2 = j;
      }
    }
  }
  create_cigar(matrix, type, cigar, target_begin, max1, max2, query,
               query_length, target);
  return maximum;
}
// sufix of the first string overlaps with the prefix of the second
int semi_global_alignment(const char* query, unsigned int query_length,
                          const char* target, unsigned int target_length,
                          AlignmentType type, int match, int mismatch, int gap,
                          string& cigar, unsigned int& target_begin) {
  vector<vector<matrixCell>> matrix(query_length + 1,
                                    vector<matrixCell>(target_length + 1));
  align(query, query_length, target, target_length, type, match, mismatch, gap,
        matrix);
  int max = 0;
  int max1 = 0;
  for (int i = 1; i < query_length + 1; i++) {
    if (matrix[query_length][i].value > max) {
      max = matrix[query_length][i].value;
      max1 = i;
    }
  }
  create_cigar(matrix, type, cigar, target_begin, max1, target_length, query,
               query_length, target);
  return max;
}

int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type, int match, int mismatch, int gap,
                       string& cigar, unsigned int& target_begin) {
  switch (type) {
    case AlignmentType::global:
      return global_alignment(query, query_length, target, target_length, type,
                              match, mismatch, gap, cigar, target_begin);
    case AlignmentType::local:
      return local_alignment(query, query_length, target, target_length, type,
                             match, mismatch, gap, cigar, target_begin);
    case AlignmentType::semi_global:
      return semi_global_alignment(query, query_length, target, target_length,
                                   type, match, mismatch, gap, cigar,
                                   target_begin);
    default:
      fprintf(stderr, "Unknown alignment type\n");
      break;
  }
}
}  // namespace brown