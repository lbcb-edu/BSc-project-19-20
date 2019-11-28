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
 * @brief Estimates CIGAR string lenght based on query length
 */
auto constexpr kAproxLenFactor = 0.8;

/**
 * @brief Describes the location of the parent
 *      cell relative to queried cell position
 */
enum ParentDir {
    kStop,  //<<< Marks an undefined score origin

    kDiagonal,  //<<< Marks a diagonal move in the matrix
    kLeft,      //<<< Marks a move from the left in the matrix
    kUp,        //<<< Marks a move from above in the matrix
};

/**
 * @brief Cell in the scoring matrix
 *      for sequence aligment mamorizing parent cell relative position
 */
struct Cell {
    int score_;
    ParentDir parent_;  //<<< marks the origin of the scoring

    bool operator<(const Cell& other) const { return score_ < other.score_; }
};

/**
 * @brief Holds coordinates of a cell in a matrix
 */
struct Pos {
    unsigned int r_;  //<<< Row coordinate
    unsigned int c_;  //<<< Column coordinate
};

int pairwiseAlignment(char const* query, unsigned int query_length,
                      char const* target, unsigned int target_length,
                      AlignmentType type, int match, int mismatch, int gap) {
    using Row = std::vector<int>;
    using ScoreFn = std::function<int(Row&, Row&, const int&, const int)>;

    auto mx_cell = std::numeric_limits<int>::min();
    auto rows = [&type, &gap, &target_length] {
        auto rows = std::array<Row, 2>{Row(target_length + 1, {0}),
                                       Row(target_length + 1, {0})};

        if (type == AlignmentType::kGlobal)
            for (auto i = int{1}; i <= target_length; ++i) rows[0][i] = i * gap;

        return rows;
    }();

    /**
     * @brief determines scoring function based on alignment type
     */
    auto scoreFn = [&type, &mx_cell, &gap]() -> ScoreFn {
        if (type == AlignmentType::kLocal)
            return [&mx_cell, &gap](Row& prev, Row& curr, const int& pos,
                                    const int matched) {
                auto ret = std::max({prev[pos - 1] + matched,
                                     curr[pos - 1] + gap, prev[pos] + gap, 0});
                mx_cell = std::max(mx_cell, ret);

                return ret;
            };

        return [&gap](Row& prev, Row& curr, const int& pos, const int matched) {
            return std::max({prev[pos - 1] + matched, curr[pos - 1] + gap,
                             prev[pos] + gap});
        };
    }();

    // Matrix traversal
    for (auto i = int{1}; i <= query_length; ++i) {
        auto& prev = rows[(i - 1) & 1];
        auto& curr = rows[i & 1];

        if (type != AlignmentType::kLocal)
            curr[0] = i * gap;

        for (auto j = int{1}; j <= target_length; ++j) {
            curr[j] = scoreFn(prev, curr, j,
                              (query[i - 1] == target[j - 1]) ? match : mismatch);
        }
    }

    const auto& last_row = rows[query_length & 1];
    if (type == AlignmentType::kLocal)
        return mx_cell;
    else {
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
    using Row = std::vector<Cell>;
    using Matrix = std::vector<Row>;
    using ScoreFn = std::function<Cell(const int&, const int&, const int)>;

    using namespace std::string_literals;

    auto matrix = [&type, &query_length, &target_length, &gap]() -> Matrix {
        auto matrix = Matrix(query_length + 1,
                             Row(target_length + 1, {0, ParentDir::kStop}));

        if (type == AlignmentType::kGlobal)
            for (auto i = int{1}; i <= target_length; ++i)
                matrix[0][i].score_ = i * gap;

        if (type != AlignmentType::kLocal)
            for (auto i = int{1}; i <= query_length; ++i)
                matrix[i][0].score_ = i * gap;

        return matrix;
    }();

    auto mx_cell = Cell{std::numeric_limits<int>::min(), ParentDir::kStop};
    auto mx_pos = Pos{0, 0};

    auto scoreGlobFn = [&matrix, &gap](const int& r, const int& c,
                                       const int matched) {
        return std::max({
            Cell{matrix[r - 1][c - 1].score_ + matched, ParentDir::kDiagonal},
            Cell{matrix[r][c - 1].score_ + gap, ParentDir::kLeft},
            Cell{matrix[r - 1][c].score_ + gap, ParentDir::kUp},
        });
    };

    auto scoreFn = [&type, &mx_cell, &mx_pos, &scoreGlobFn]() -> ScoreFn {
        if (type == AlignmentType::kLocal) {
            return [&mx_cell, &mx_pos, &scoreGlobFn](const unsigned int& r,
                                                     const unsigned int& c,
                                                     const int matched) {
                auto ret = std::max(scoreGlobFn(r, c, matched),
                                    Cell{0, ParentDir::kStop});

                if (mx_cell < ret) {
                    mx_cell = ret;
                    mx_pos = {r, c};
                }

                return ret;
            };
        } else
            return scoreGlobFn;
    }();

    for (auto i = std::uint32_t{1}; i <= query_length; ++i) {
        for (auto j = std::uint32_t{1}; j <= target_length; ++j) {
            matrix[i][j] =
                scoreFn(i, j, (query[i - 1] == target[j - 1]) ? match : mismatch);
        }
    }

    auto backtrack = [&matrix, &cigar, &target_begin](Pos pos) {
        auto nextPos = [](const Pos& curr, const ParentDir& dir) {
            switch (dir) {
                case ParentDir::kDiagonal:
                    return Pos{curr.r_ - 1, curr.c_ - 1};
                case ParentDir::kLeft:
                    return Pos{curr.r_, curr.c_ - 1};
                case ParentDir::kUp:
                    return Pos{curr.r_ - 1, curr.c_};
                default:
                    throw std::logic_error("Ups, something went wrong. :(");
            }
        };

        auto cigarFrag = [](ParentDir dir, int len) {
            auto n = std::to_string(len);
            switch (dir) {
                case ParentDir::kDiagonal:
                    return "M" + n;
                case ParentDir::kLeft:
                    return "D" + n;
                case ParentDir::kUp:
                    return "I" + n;

                default:
                    throw std::logic_error("Ups, something went wrong. :(");
            }

            return ""s;
        };
        cigar.reserve(matrix.size() * kAproxLenFactor);

        auto curr_dir = matrix[pos.r_][pos.c_].parent_;
        auto curr_pos = pos;
        int cnt = 1;

        while (curr_dir != kStop) {
            auto nxt_pos = nextPos(curr_pos, curr_dir);
            const auto& nxt_dir = matrix[nxt_pos.r_][nxt_pos.c_].parent_;

            if (nxt_dir != curr_dir) {
                cigar += cigarFrag(curr_dir, cnt);

                curr_dir = nxt_dir;
                cnt = 1;
            } else
                ++cnt;

            curr_pos = nxt_pos;
        }

        cigar.shrink_to_fit();

        std::reverse(begin(cigar), end(cigar));
        target_begin = curr_pos.c_;
    };

    switch (type) {
        case AlignmentType::kLocal:
            backtrack({mx_pos.r_, mx_pos.c_});
            return mx_cell.score_;

        case AlignmentType::kGlobal:
            backtrack({query_length, target_length});
            return matrix[query_length][target_length].score_;

        case AlignmentType::kSemiGlobal:
            auto& last_row = matrix[query_length];

            auto pos = static_cast<unsigned int>(
                std::max_element(begin(last_row), end(last_row)) -
                last_row.begin());

            backtrack({query_length, pos});
            return last_row[pos].score_;
    }
}

}  // namespace alignment
}  // namespace orange