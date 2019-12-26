#include <algorithm>
#include <tuple>

#include <gtest/gtest.h>

#include "orange_alignment.h"
#include "orange_minimizers.h"

TEST(PairwiseAlignmentRow, global) {
    EXPECT_EQ(orange::alignment::pairwiseAlignment(
                  "GCATGCU", 7, "GATTACA", 7,
                  orange::alignment::AlignmentType::kGlobal, 1, -1, -1),
              0);
}

TEST(PairwiseAlignmentRow, local) {
    EXPECT_EQ(orange::alignment::pairwiseAlignment(
                  "ACCTAAGG", 8, "GGCTCAATCA", 10,
                  orange::alignment::AlignmentType::kLocal, 2, -1, -2),
              6);
}

TEST(PairwiseAlignmentRow, semi_global) {
    EXPECT_EQ(orange::alignment::pairwiseAlignment(
                  "TCCG", 4, "ACTCCGAT", 8,
                  orange::alignment::AlignmentType::kSemiGlobal, 4, -1, -2),
              16);
}

TEST(PairwiseAlignmentMatrix, global) {
    std::string cigar;
    unsigned int target_start = 0;

    EXPECT_EQ(orange::alignment::pairwiseAlignment(
                  "GCATGCU", 7, "GATTACA", 7,
                  orange::alignment::AlignmentType::kGlobal, 1, -1, -1, cigar,
                  target_start),
              0);
    EXPECT_EQ(target_start, 0);
}

TEST(PairwiseAlignmentMatrix, local) {
    std::string cigar;
    unsigned int target_start = 0;

    EXPECT_EQ(orange::alignment::pairwiseAlignment(
                  "ACCTAAGG", 8, "GGCTCAATCA", 10,
                  orange::alignment::AlignmentType::kLocal, 2, -1, -2, cigar,
                  target_start),
              6);
    EXPECT_EQ(target_start, 2);
}

TEST(Minimizers, s15k3win5) {
    using namespace orange::minimizers;

    auto sequence = std::string{"GACATCATCGCCACA"};
    auto correct_ans = KMers{{1, 10, 0}, {4, 11, 0},  {6, 2, 0}, {6, 5, 0},
                             {17, 1, 0}, {17, 12, 0}, {52, 0, 0}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 3, 5);
    std::sort(algo_ans.begin(), algo_ans.end());

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s15k1win1) {
    using namespace orange::minimizers;

    auto sequence = std::string{"TACTAGCTTCCATGT"};
    auto correct_ans =
        KMers{{0, 2, 0},  {0, 5, 1}, {0, 6, 0},  {0, 9, 0},  {0, 10, 0},
              {0, 13, 1}, {1, 0, 1}, {1, 1, 0},  {1, 3, 1},  {1, 4, 0},
              {1, 7, 1},  {1, 8, 1}, {1, 11, 0}, {1, 12, 1}, {1, 14, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 1, 1);
    std::sort(algo_ans.begin(), algo_ans.end());

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s64k16win32) {
    using namespace orange::minimizers;

    auto sequence = std::string{
        "ATCTGCGGCAGTTGTATCCGATGTTTGAGACCTAGTCAGTGTCTTGATGGACGTTCCATCAACG"};
    auto correct_ans =
        KMers{{102188212, 34, 1}, {165801515, 30, 0}, {230340418, 17, 0},
              {353957445, 11, 1}, {418792748, 7, 1},  {518174452, 36, 0},
              {728714881, 42, 0}, {754051992, 2, 0},  {987485408, 48, 1},
              {1657740985, 0, 0}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 16, 32);
    std::sort(algo_ans.begin(), algo_ans.end());

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s15k7win5) {
    using namespace orange::minimizers;

    auto sequence = std::string{"AAAACAGATGGGGTT"};
    auto correct_ans =
        KMers{{24, 6, 1},   {98, 5, 1},   {395, 4, 1}, {1582, 3, 1},
              {2986, 0, 1}, {4102, 7, 1}, {5121, 8, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 7, 5);

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s64k16win8) {
    using namespace orange::minimizers;

    auto sequence = std::string{
        "GAACTCTGCTAGGAAACCTTGGCCTGCATAATCAAAAGTAATGTCGGCATTGGACACTCGCTGA"};
    auto correct_ans =
        KMers{{95474084, 37, 1},  {100303473, 6, 1},   {183546262, 16, 0},
              {186212695, 22, 0}, {219826861, 31, 1},  {401213895, 5, 1},
              {425809814, 26, 0}, {452217035, 47, 0},  {520485506, 10, 1},
              {583660866, 3, 0},  {934811331, 44, 1},  {1219657040, 2, 0},
              {1378656084, 1, 0}, {1808868141, 48, 0}, {2821135848, 0, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 16, 8);
    std::sort(algo_ans.begin(), algo_ans.end());

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s16k5w8) {
    using namespace orange::minimizers;

    auto sequence = std::string{"AGCTGCACTGATGCTA"};
    auto correct_ans =
        KMers{{75, 5, 0}, {97, 8, 1}, {114, 0, 1}, {434, 10, 0}, {625, 11, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 5, 8);
    std::sort(algo_ans.begin(), algo_ans.end());

    EXPECT_EQ(correct_ans, algo_ans);
}