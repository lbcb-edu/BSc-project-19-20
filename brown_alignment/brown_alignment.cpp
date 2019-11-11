#include "brown_alignment.hpp"
#include <iostream>
#include <vector>

using namespace std;

namespace brown {

enum class valueType { deletion, left, match, mismatch, none };

struct matrixCell {
  int value;
  valueType cost;
};

void create_cigar(vector<vector<matrixCell>> &matrix, AlignmentType type, string& cigar, unsigned int& target_begin, int row, int column) {
  if(type == AlignmentType::local) {
    while (matrix[row][column].value != 0) {
      if (matrix[row][column].cost == valueType::deletion){
        cigar.push_back('D');
        row--;
      } else if (matrix[row][column].cost == valueType::left) {
        cigar.push_back('I');
        column--;
      } else if (matrix[row][column].cost == valueType::match) {
        cigar.push_back('=');
        row--;
        column--;
      } else if (matrix[row][column].cost == valueType::mismatch) {
        cigar.push_back('X');
        row--;
        column--;
      }
    }
  } else {
    while (matrix[row][column].cost != valueType::none) {
      if (matrix[row][column].cost == valueType::deletion){
        cigar.push_back('D');
        row--;
      } else if (matrix[row][column].cost == valueType::left) {
        cigar.push_back('I');
        column--;
      } else if (matrix[row][column].cost == valueType::match) {
        cigar.push_back('=');
        row--;
        column--;
      } else if (matrix[row][column].cost == valueType::mismatch) {
        cigar.push_back('X');
        row--;
        column--;
      }
    }
  }
}

void align(const char* query, unsigned int query_length,
                 const char* target, unsigned int target_length,
                 AlignmentType type, int match, int mismatch, int gap,
                 vector<vector<matrixCell>>& matrix) {
  matrix[0][0].cost = valueType::none;
  matrix[0][0].value = 0;

  if (type == AlignmentType::global) {
    for (int i = 1; i < query_length + 1; i++) {
      matrix[i][0].value = gap * i;
      matrix[i][0].cost = valueType::deletion;
    }
    for (int i = 1; i < target_length + 1; i++) {
      matrix[0][i].value = gap * i;
      matrix[0][i].cost = valueType::left;
    }
  } else if (type == AlignmentType::semi_global) {
    for (int i = 1; i < query_length + 1; i++) {
      matrix[i][0].value = gap * i;
      matrix[i][0].cost = valueType::deletion;
    }
    for (int i = 1; i < target_length + 1; i++) {
      matrix[0][i].value = 0;
      matrix[0][i].cost = valueType::left;
    }
  } else if (type == AlignmentType::local) {
    for (int i = 1; i < target_length + 1; i++) {
      matrix[0][i].value = 0;
      matrix[0][i].cost = valueType::left;
    }
    for (int i = 1; i < query_length + 1; i++) {
      matrix[i][0].value = 0;
      matrix[i][0].cost = valueType::deletion;
    }
  }
  for (int i = 1; i < query_length + 1; i++) {
    for (int j = 1; j < target_length + 1; j++) {
      int matchPrev;
      bool matchCheck;
      if (query[i - 1] == target[i - 1]) {
        matchPrev = matrix[i - 1][j - 1].value + match;
        matchCheck = true;
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
          matrix[i][j].value = matchPrev;
          if (matchCheck == true)
            matrix[i][j].cost = valueType::match;
          else 
            matrix[i][j].cost = valueType::mismatch;
        }
      } else {
        if (deletion > insertion) {
          matrix[i][j].value = deletion;
          matrix[i][j].cost = valueType::deletion;
        } else {
          matrix[i][j].value = insertion;
          matrix[i][j].cost = valueType::left;
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

  align(query, query_length, target, target_length, type, match, mismatch,
              gap, matrix);
  create_cigar(matrix, type, cigar, target_begin, query_length, target_length);
  return matrix[query_length][target_length].value;
}

int local_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type, int match, int mismatch, int gap,
                       string& cigar, unsigned int& target_begin) {
  vector<vector<matrixCell>> matrix(query_length + 1,
                                    vector<matrixCell>(target_length + 1));
  align(query, query_length, target, target_length, type, match, mismatch,
              gap, matrix);
  int maximum = 0;
  for (int i = 1; i < query_length + 1; i++) {
    for (int j = 1; j < target_length + 1; j++) {
      if (matrix[i][j].value > maximum) maximum = matrix[i][j].value;
    }
  }
  create_cigar(matrix, type, cigar, target_begin, query_length, target_length);
  return maximum;
}
//sufix of the first string overlaps with the prefix of the second
int semi_global_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type, int match, int mismatch, int gap,
                       string& cigar, unsigned int& target_begin) {
  vector<vector<matrixCell>> matrix(query_length + 1,
                                    vector<matrixCell>(target_length + 1));
  align(query, query_length, target, target_length, type, match, mismatch,
              gap, matrix);
  int max = 0;
  for (int i = 1; i < query_length + 1; i++) {
    if (matrix[query_length][i].value > max)
      max = matrix[query_length][i].value;
  }
  create_cigar(matrix, type, cigar, target_begin, query_length, 1);
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
                                   type, match, mismatch, gap, cigar, target_begin);
    default:
      fprintf(stderr, "Unknown alignment type\n");
      break;
  }
}
}  // namespace brown