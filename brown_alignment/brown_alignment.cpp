#include "brown_alignment.hpp"
#include <iostream>
#include <vector>

using namespace std;

namespace brown {

enum class ParentType { up, left, upLeft, none };

struct matrixCell {
  int value;
  ParentType parent;
};

void init_matrix(const char* query, unsigned int query_length,
                 const char* target, unsigned int target_length,
                 AlignmentType type, int match, int mismatch, int gap,
                 vector<vector<matrixCell>>& matrix) {
  matrix[0][0].parent = ParentType::none;
  matrix[0][0].value = 0;

  if (type == AlignmentType::global) {
    for (int i = 1; i < query_length + 1; i++) {
      matrix[i][0].value = gap * i;
      matrix[i][0].parent = ParentType::up;
    }
    for (int i = 1; i < target_length + 1; i++) {
      matrix[0][i].value = gap * i;
      matrix[0][i].parent = ParentType::left;
    }
  } else if (type == AlignmentType::semi_global) {
    for (int i = 1; i < query_length + 1; i++) {
      matrix[i][0].value = gap * i;
      matrix[i][0].parent = ParentType::up;
    }
    for (int i = 1; i < target_length + 1; i++) {
      matrix[0][i].value = 0;
      matrix[0][i].parent = ParentType::left;
    }
  }
  for (int i = 1; i < query_length + 1; i++) {
    for (int j = 1; j < target_length + 1; j++) {
      int matchPrev;
      if (query[i - 1] == target[i - 1]) {
        matchPrev = matrix[i - 1][j - 1].value + match;
      } else {
        matchPrev = matrix[i - 1][j - 1].value + mismatch;
      }
      int insertion = matrix[i][j - 1].value + gap;
      int deletion = matrix[i - 1][j].value + gap;

      if (matchPrev > insertion) {
        if (deletion > matchPrev) {
          matrix[i][j].value = deletion;
          matrix[i][j].parent = ParentType::up;
        } else {
          matrix[i][j].value = matchPrev;
          matrix[i][j].parent = ParentType::upLeft;
        }
      } else {
        if (deletion > insertion) {
          matrix[i][j].value = deletion;
          matrix[i][j].parent = ParentType::up;
        } else {
          matrix[i][j].value = insertion;
          matrix[i][j].parent = ParentType::left;
        }
      }
    }
  }
}

int global_alignment(const char* query, unsigned int query_length,
                     const char* target, unsigned int target_length,
                     AlignmentType type, int match, int mismatch, int gap) {
  vector<vector<matrixCell>> matrix(query_length + 1,
                                    vector<matrixCell>(target_length + 1));

  init_matrix(query, query_length, target, target_length, type, match, mismatch,
              gap, matrix);
  return matrix[query_length][target_length].value;
}

int local_alignment(const char* query, unsigned int query_length,
                    const char* target, unsigned int target_length,
                    AlignmentType type, int match, int mismatch, int gap) {
  vector<vector<matrixCell>> matrix(query_length + 1,
                                    vector<matrixCell>(target_length + 1));
  int max = 0;
  for (int i = 1; i < query_length + 1; i++) {
    if (matrix[query_length + 1][i].value > max)
      max = matrix[query_length + 1][i].value;
  }
  return max;
}

int semi_global_alignment(const char* query, unsigned int query_length,
                          const char* target, unsigned int target_length,
                          AlignmentType type, int match, int mismatch,
                          int gap) {
  vector<vector<matrixCell>> matrix(query_length + 1,
                                    vector<matrixCell>(target_length + 1));
}

int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type, int match, int mismatch, int gap) {
  switch (type) {
    case AlignmentType::global:
      return global_alignment(query, query_length, target, target_length, type,
                              match, mismatch, gap);
    case AlignmentType::local:
      return local_alignment(query, query_length, target, target_length, type,
                             match, mismatch, gap);
    case AlignmentType::semi_global:
      return semi_global_alignment(query, query_length, target, target_length,
                                   type, match, mismatch, gap);
    default:
      fprintf(stderr, "Unknown alignment type\n");
      break;
  }
}
}  // namespace brown