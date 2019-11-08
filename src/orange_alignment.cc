#include <initializer_list>
#include <functional>
#include <algorithm>
#include <limits>
#include <vector>
#include <string>
#include <array>

#include "orange_alignment.h"

namespace orange {
namespace alignment {

/**
 * @brief Describes the location of the parent
 *      cell relative to queried cell position
 */
enum ParentDir {
    kDiagonal,  //<<< Marks a diagonal move in the matrix
    kLeft,      //<<< Marks a move from the left in the matrix
    kUp,        //<<< Marks a move from above in the matrix

    kNo  //<<< The cell is located at the left or upper edge of the matrix
};

/**
 * @brief Represents score in the scoring matrix
 *      for sequence alignment algorithms without backtracking for CIGAR
 */
using Score = int;

/**
 * @brief Cell in the scoring matrix
 *      for sequence aligment mamorizing parent cell relative position
 */
struct MemCell {
    Score score_;
    ParentDir parent_;  //<<< marks the origin of the scoring
};

/**
 * @brief References surrounding cells
 *      relevant for updating the score
 */
struct CellSurr {
    MemCell& diag_;
    MemCell& left_;
    MemCell& up_;
};

template <typename Cell>
using TRow = std::vector<Cell>;

using Matrix = std::vector<TRow<MemCell>>;

template <typename Row>
auto gapRowFactory(int gap, int sz) {
    auto row = Row(sz + 1);
    for (auto i = int{1}; i <= sz; ++i) row[i] = 1 * gap;
    return row;
}

int pairwiseAlignment(char const* query, unsigned int query_length,
                      char const* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap) {
    using Row = TRow<Score>;
    using ScoreFn_t = std::function<int(Row&, Row&, const int&, int)>;

    auto mx_val = std::numeric_limits<int>::min();
    auto rows = [&type, &gap, &target_length] {
        if (type == AlignmentType::kGlobal)
            return std::array<Row, 2>{gapRowFactory<Row>(gap, target_length),
                                      Row(target_length + 1, {0})};
        return std::array<Row, 2>{Row(target_length + 1, {0}),
                                  Row(target_length + 1, {0})};
    }();

    auto scoreFn = [&type, &mx_val, &gap]() -> ScoreFn_t {
        if (type == AlignmentType::kLocal)
            return [&mx_val, &gap](Row& prev, Row& curr, const int& pos,
                                   int matched) {
                auto ret = std::max({prev[pos - 1] + matched,
                                     curr[pos - 1] + gap, prev[pos] + gap, 0});
                mx_val = std::max(mx_val, ret);

                return ret;
            };

        return [&gap](Row& prev, Row& curr, const int& pos, int matched) {
            return std::max({prev[pos - 1] + matched, curr[pos - 1] + gap,
                             prev[pos] + gap});
        };
    }();

    for (auto i = int{1}; i <= query_length; ++i) {
        auto& prev = rows[(i - 1) & 1];
        auto& curr = rows[i & 1];

        if (type != AlignmentType::kLocal)
            curr[0] = i * gap;

        for (auto j = int{1}; j <= target_length; ++j) {
            curr[j] = scoreFn(prev, curr, j,
                              (query[i] == target[j]) ? match : mismatch);
        }
    }

    if (type == AlignmentType::kLocal)
        return mx_val;
    else {
        auto& last_row = rows[query_length & 1];

        if (type == AlignmentType::kGlobal)
            return last_row.back();

        // Semi-global
        return *std::max_element(begin(last_row), end(last_row));
    }
}

int pairwiseAlignment(char const* query, unsigned int query_length,
                      char const* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap,
                      std::string& cigar, unsigned int& target_begin) {
    using Row = TRow<MemCell>;

    return 0;
}

}  // namespace alignment
}  // namespace orange