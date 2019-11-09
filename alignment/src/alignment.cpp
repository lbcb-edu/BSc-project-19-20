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

inline Cell max(Cell a, Cell b, Cell c) {
  return a.score > b.score ? a.score > c.score ? a : c
                           : b.score > c.score ? b : c;
}

Contiguous2DArray<Cell> NeedlemanWunsch(const char* query,
                                        unsigned int query_length,
                                        const char* target,
                                        unsigned int target_length, int match,
                                        int mismatch, int gap) {
  Contiguous2DArray<Cell> alignment(query_length + 1, target_length + 1);

  for (int j = 1; j < alignment.cols; ++j)
    alignment[0][j] = {j * gap, Action::kInsert};

  for (int i = 1; i < alignment.rows; ++i) {
    alignment[i][0] = {i * gap, Action::kDelete};

    for (int j = 1; j < alignment.cols; ++j) {
      bool is_match = query[i - 1] == target[j - 1];
      alignment[i][j] =
          max({alignment[i - 1][j - 1].score + (is_match ? match : mismatch),
               is_match ? Action::kMatch : Action::kMismatch},
              {alignment[i][j - 1].score + gap, Action::kInsert},
              {alignment[i - 1][j].score + gap, Action::kDelete});
    }
  }

  return alignment;
}

using Position = ::std::pair<::std::size_t, ::std::size_t>;

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
                    ::std::string& cigar, unsigned int& target_begin) {
  detail::Cell current = alignment.last();
  Position pos{alignment.rows - 1, alignment.cols - 1};

  bool begin_set = false;

  while (current.action != kNone) {
    cigar += current.action;
    UpdatePosition(current.action, pos);

    if (!pos.first && !begin_set) target_begin = pos.second;

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