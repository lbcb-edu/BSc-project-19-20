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
template <typename Cell>
struct CellSurr {
    Cell& diag_;
    Cell& left_;
    Cell& up_;
};

template <typename Cell>
using TRow = std::vector<Cell>;

template <typename Row>
using Matrix = std::vector<Row>;

int pairwiseAlignment(char const* query, unsigned int query_length,
                      char const* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap) {
    using Row = TRow<Score>;
    using Surr = CellSurr<Score>;

    auto scoreFn = [&gap](Surr&& surr, int matched) {
        return std::max(
            {surr.diag_ + matched, surr.left_ + gap, surr.up_ + gap});
    };

    auto gapRowFactory = [&gap, &sz = target_length]() {
        auto row = Row(sz + 1);
        for (auto i = int{1}; i <= sz; ++i) row[i] = i * gap;
        return row;
    };

    auto rows = [&type, &gapRowFactory, &target_length] {
        if (type == AlignmentType::kGlobal)
            return std::array<Row, 2>{gapRowFactory(),
                                      Row(target_length + 1, {0})};
        return std::array<Row, 2>{Row(target_length + 1, {0}),
                                  Row(target_length + 1, {0})};
    }();

    auto mx_val = std::numeric_limits<int>::min();

    for (auto i = int{1}; i <= query_length; ++i) {
        auto& prev = rows[(i - 1) & 1];
        auto& curr = rows[i & 1];

        if (type != AlignmentType::kLocal)
            curr[0] = i * gap;

        for (auto j = int{1}; j <= target_length; ++j) {
            curr[j] = scoreFn({prev[j - 1], curr[j - 1], prev[j]},
                              (query[i] == target[j]) ? match : mismatch);

            if (type == AlignmentType::kLocal) {
                if (curr[j] < 0)
                    curr[j] = 0;
                mx_val = std::max(mx_val, curr[j]);
            }
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
    return 0;
}

}  // namespace alignment
}  // namespace orange