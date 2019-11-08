#include "alignment/alignment.hpp"

// TODO: delet
#include <iostream>

namespace algn {
namespace detail {

enum Action { kMatch, kMismatch, kInsert, kDelete, kNone };

struct Cell {
  int score;
  Action action = kNone;
};

Contiguous2DArray<Cell> NeedlemanWunsch(const char* query,
                                        unsigned int query_length,
                                        const char* target,
                                        unsigned int target_length, int match,
                                        int mismatch, int gap) {
  Contiguous2DArray<Cell> alignment(query_length + 1, target_length + 1);

  for (int j = 0; j < alignment.cols; ++j)
    alignment[0][j] = {j * gap, Action::kInsert};

  for (int i = 1; i < alignment.rows; ++i) {
    alignment[i][0] = {i * gap, Action::kDelete};

    for (int j = 1; j < alignment.cols; ++j) {
      bool is_match = query[i - 1] == target[j - 1];
      int aligned =
              alignment[i - 1][j - 1].score + (is_match ? match : mismatch),
          inserted = alignment[i][j - 1].score + gap,
          deleted = alignment[i - 1][j].score + gap;

      // clang-format off
      if (aligned > inserted)
        if (aligned > deleted)
          alignment[i][j] = {aligned,
                             is_match ? Action::kMatch : Action::kMismatch};
        else
          alignment[i][j] = {deleted,  Action::kDelete};
      else 
        if (inserted > deleted)
          alignment[i][j] = {inserted, Action::kInsert};
      else
          alignment[i][j] = {deleted,  Action::kDelete};
      // clang-format on
    }
  }

  return alignment;
}

void ConstructCigar(const detail::Contiguous2DArray<detail::Cell>& alignment,
                    ::std::string& cigar, unsigned int& target_begin) {
  detail::Cell current = alignment.last();
  ::std::pair<::std::size_t, ::std::size_t> pos{alignment.rows - 1,
                                                alignment.cols - 1};

  bool begin_set = false;

  while (true) {
    switch (current.action) {
      case Action::kMatch:
        cigar.append("=");
        --pos.first;
        --pos.second;
        break;

      case Action::kMismatch:
        cigar.append("X");
        --pos.first;
        --pos.second;
        break;

      case Action::kInsert:
        cigar.append("D");
        --pos.second;
        break;

      case Action::kDelete:
        cigar.append("I");
        --pos.first;
        break;
    }

    if (!pos.first && !begin_set) target_begin = pos.second;
    if (pos == decltype(pos){0, 0}) break;

    current = alignment[pos.first][pos.second];
  }

  ::std::reverse(::std::begin(cigar), ::std::end(cigar));
}

}  // namespace detail

int PairwiseAlignment(const char* query, unsigned int query_length,
                      const char* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap) {
  if (type == AlignmentType::kNeedlemanWunsch)
    return detail::NeedlemanWunsch(query, query_length, target, target_length,
                                   match, mismatch, gap)
        .last()
        .score;

  return -1;
}

int PairwiseAlignment(const char* query, unsigned int query_length,
                      const char* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap,
                      ::std::string& cigar, unsigned int& target_begin) {
  if (type == AlignmentType::kNeedlemanWunsch) {
    auto alignment = detail::NeedlemanWunsch(
        query, query_length, target, target_length, match, mismatch, gap);

    detail::ConstructCigar(alignment, cigar, target_begin);
    return alignment.last().score;
  }

  return -1;
}

}  // namespace algn