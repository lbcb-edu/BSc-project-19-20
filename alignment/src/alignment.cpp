#include "alignment/alignment.hpp"

namespace algn {
namespace detail {

// clang-format off
enum Action {
  kMatch    = '=',
  kMismatch = 'X',
  kInsert   = 'D',
  kDelete   = 'I',
  kNone     = '?'
};
// clang-format on

struct Cell {
  int score = 0;
  Action action = kNone;
};

inline Cell Max(Cell a, Cell b) { return a.score > b.score ? a : b; }

inline Cell Max(Cell a, Cell b, Cell c) {
  return a.score > b.score ? Max(a, c) : Max(b, c);
}

Contiguous2DArray<Cell> NeedlemanWunsch(const char* query,
                                        unsigned int query_length,
                                        const char* target,
                                        unsigned int target_length, int match,
                                        int mismatch, int gap) {
  Contiguous2DArray<Cell> alignment{query_length + 1, target_length + 1};

  for (int j = 1; j < alignment.cols; ++j)
    alignment[0][j] = {j * gap, Action::kInsert};

  for (int i = 1; i < alignment.rows; ++i) {
    alignment[i][0] = {i * gap, Action::kDelete};

    for (int j = 1; j < alignment.cols; ++j) {
      bool is_match = query[i - 1] == target[j - 1];
      alignment[i][j] =
          Max({alignment[i - 1][j - 1].score + (is_match ? match : mismatch),
               is_match ? Action::kMatch : Action::kMismatch},
              {alignment[i][j - 1].score + gap, Action::kInsert},
              {alignment[i - 1][j].score + gap, Action::kDelete});
    }
  }

  return alignment;
}

using Position = ::std::pair<::std::size_t, ::std::size_t>;
using MaxCell = ::std::pair<Contiguous2DArray<Cell>, Position>;

MaxCell SmithWaterman(const char* query, unsigned int query_length,
                      const char* target, unsigned int target_length, int match,
                      int mismatch, int gap) {
  Contiguous2DArray<Cell> alignment{query_length + 1, target_length + 1};
  ::std::pair<int, Position> best_cell{0, {0, 0}};

  for (int i = 1; i < alignment.rows; ++i)
    for (int j = 1; j < alignment.cols; ++j) {
      bool is_match = query[i - 1] == target[j - 1];
      alignment[i][j] = Max(
          Cell{0, Action::kNone},
          Max({alignment[i - 1][j - 1].score + (is_match ? match : mismatch),
               is_match ? Action::kMatch : Action::kMismatch},
              {alignment[i][j - 1].score + gap, Action::kInsert},
              {alignment[i - 1][j].score + gap, Action::kDelete}));

      if (alignment[i][j].score > best_cell.first) {
        best_cell.first = alignment[i][j].score;
        best_cell.second = {i, j};
      }
    }

  return {::std::move(alignment), best_cell.second};
}

inline void UpdatePosition(Action action, Position& pos) {
  switch (action) {
    case Action::kMatch:
      --pos.first;
      --pos.second;
      break;

    case Action::kMismatch:
      --pos.first;
      --pos.second;
      break;

    case Action::kInsert:
      --pos.second;
      break;

    case Action::kDelete:
      --pos.first;
      break;
  }
}

void ConstructCigar(const detail::Contiguous2DArray<detail::Cell>& alignment,
                    Position pos, ::std::string& cigar,
                    unsigned int& target_begin) {
  Cell current = alignment[pos.first][pos.second];

  while (current.action != kNone) {
    cigar += current.action;
    target_begin = pos.second - 1;
    UpdatePosition(current.action, pos);
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

  if (type == AlignmentType::kSmithWaterman) {
    auto ret = detail::SmithWaterman(query, query_length, target, target_length,
                                     match, mismatch, gap);

    return ret.first[ret.second.first][ret.second.first].score;
  }

  return -1;
}

int PairwiseAlignment(const char* query, unsigned int query_length,
                      const char* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap,
                      ::std::string& cigar, unsigned int& target_begin) {
  if (type == AlignmentType::kNeedlemanWunsch) {
    auto alignment = detail::NeedlemanWunsch(
        query, query_length, target, target_length, match, mismatch, gap);

    detail::ConstructCigar(alignment, {alignment.rows - 1, alignment.cols - 1},
                           cigar, target_begin);
    return alignment.last().score;
  }

  if (type == AlignmentType::kSmithWaterman) {
    auto ret = detail::SmithWaterman(query, query_length, target, target_length,
                                     match, mismatch, gap);
    auto start = ret.first[ret.second.first][ret.second.second];

    detail::ConstructCigar(ret.first, ret.second, cigar, target_begin);
    return start.score;
  }

  return -1;
}

}  // namespace algn