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

struct AlignmentContext {
  const char* query;
  unsigned int query_length;

  const char* target;
  unsigned int target_length;

  int match;
  int mismatch;
  int gap;
};

struct Cell {
  int score = 0;
  Action action = kNone;
};

inline Cell Max(Cell a, Cell b) {
  return a.score > b.score ? a : b;
}

inline Cell Max(Cell a, Cell b, Cell c) {
  return a.score > b.score ? Max(a, c) : Max(b, c);
}

using Position = ::std::pair<::std::size_t, ::std::size_t>;
using MaxCell = ::std::pair<Contiguous2DArray<Cell>, Position>;

template <typename R>
R Align(const AlignmentContext& ai,
        ::std::function<void(Contiguous2DArray<Cell>&)>
            init,  // allowed to be emtpy
        ::std::function<Cell(Cell, Position)> update,
        ::std::function<R(Contiguous2DArray<Cell>&&)> compute_return) {
  Contiguous2DArray<Cell> alignment{ai.query_length + 1, ai.target_length + 1};

  if (init)
    init(alignment);

  for (int i = 1; i < alignment.rows; ++i)
    for (int j = 1; j < alignment.cols; ++j) {
      bool is_match = ai.query[i - 1] == ai.target[j - 1];
      Cell candidate = Max(
          {alignment[i - 1][j - 1].score + (is_match ? ai.match : ai.mismatch),
           is_match ? Action::kMatch : Action::kMismatch},
          {alignment[i][j - 1].score + ai.gap, Action::kInsert},
          {alignment[i - 1][j].score + ai.gap, Action::kDelete});

      alignment[i][j] = update(candidate, {i, j});
    }

  return compute_return(::std::forward<decltype(alignment)>(alignment));
}

Contiguous2DArray<Cell> NeedlemanWunsch(const AlignmentContext& ai) {
  return Align<Contiguous2DArray<Cell>>(
      ai,
      [&ai](auto& alignment) {
        for (int j = 1; j < alignment.cols; ++j)
          alignment[0][j] = {j * ai.gap, Action::kInsert};
        for (int i = 1; i < alignment.rows; ++i)
          alignment[i][0] = {i * ai.gap, Action::kDelete};
      },
      [](auto cell, auto) { return cell; },
      [](auto&& alignment) {
        return ::std::forward<decltype(alignment)>(alignment);
      });
}  // namespace detail

MaxCell SmithWaterman(const AlignmentContext& ai) {
  ::std::pair<int, Position> best_cell{0, {0, 0}};

  return Align<MaxCell>(
      ai, [](auto&) {},
      [&best_cell](auto cell, auto pos) {
        cell = Max(Cell{}, cell);
        if (cell.score > best_cell.first)
          best_cell = {cell.score, pos};
        return cell;
      },
      [&best_cell](auto&& alignment) {
        return MaxCell{::std::forward<decltype(alignment)>(alignment),
                       best_cell.second};
      });
}

MaxCell Overlap(const AlignmentContext& ai) {
  return Align<MaxCell>(
      ai,
      [&ai](auto& alignment) {
        for (int i = 1; i < alignment.rows; ++i)
          alignment[i][0] = {i * ai.gap, Action::kDelete};
      },
      [](auto cell, auto) { return cell; },
      [](auto&& alignment) {
        auto r = alignment.rows - 1;
        ::std::pair<int, Position> best_cell{0, {r, 0}};

        for (int j = 1; j < alignment.cols; ++j)
          if (alignment[r][j].score > best_cell.first) {
            best_cell = {alignment[r][j].score, {r, j}};
          }

        return MaxCell{::std::move(alignment), best_cell.second};
      });
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
                    unsigned int& target_begin,
                    ::std::function<void(::std::string&)> cigar_modifier) {
  Cell current = alignment[pos.first][pos.second];

  while (current.action != kNone) {
    cigar += current.action;
    target_begin = pos.second - 1;
    UpdatePosition(current.action, pos);
    current = alignment[pos.first][pos.second];
  }

  if (cigar_modifier)
    cigar_modifier(cigar);

  ::std::reverse(::std::begin(cigar), ::std::end(cigar));
}

}  // namespace detail

int PairwiseAlignment(Query query, QueryLength query_length, Target target,
                      TargetLength target_length, AlignmentType type,
                      Match match, Mismatch mismatch, Gap gap) {
  const detail::AlignmentContext ai{
      query.get(), query_length.get(), target.get(), target_length.get(),
      match.get(), mismatch.get(),     gap.get()};

  if (type == AlignmentType::kNeedlemanWunsch)
    return detail::NeedlemanWunsch(ai).last().score;

  auto ret = type == AlignmentType::kSmithWaterman ? detail::SmithWaterman(ai)
                                                   : detail::Overlap(ai);

  return ret.first[ret.second.first][ret.second.first].score;
}

int PairwiseAlignment(Query query, QueryLength query_length, Target target,
                      TargetLength target_length, AlignmentType type,
                      Match match, Mismatch mismatch, Gap gap,
                      ::std::string& cigar, unsigned int& target_begin) {
  const detail::AlignmentContext ai{
      query.get(), query_length.get(), target.get(), target_length.get(),
      match.get(), mismatch.get(),     gap.get()};

  if (type == AlignmentType::kNeedlemanWunsch) {
    auto alignment = detail::NeedlemanWunsch(ai);
    detail::ConstructCigar(alignment, {alignment.rows - 1, alignment.cols - 1},
                           cigar, target_begin, {});

    return alignment.last().score;
  }

  bool isSW = type == AlignmentType::kSmithWaterman;

  auto ret = isSW ? detail::SmithWaterman(ai) : detail::Overlap(ai);
  auto start = ret.first[ret.second.first][ret.second.second];

  detail::ConstructCigar(ret.first, ret.second, cigar, target_begin,
                         [isSW](::std::string& cigar) {
                           if (isSW)
                             cigar.resize(cigar.length() - 1);
                         });

  return start.score;
}

}  // namespace algn