#include <algorithm>
#include <string>
#include <vector>
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

    auto sequence = std::string{"TGGTCTCACTGACCC"};
    auto correct_ans = KMers{{5, 0, 1}, {5, 11, 0}, {5, 11, 0}, {7, 7, 0},
                             {7, 7, 0}, {7, 7, 0},  {7, 9, 1},  {8, 3, 1},
                             {8, 3, 1}, {8, 3, 1},  {11, 5, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 3, 5);
    std::sort(begin(algo_ans), end(algo_ans));

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s15k3win3) {
    using namespace orange::minimizers;

    auto sequence = std::string{"GTCGGAGGTCTTTGC"};
    auto correct_ans =
        KMers{{0, 10, 1}, {1, 11, 1}, {6, 12, 1}, {8, 8, 1},  {9, 1, 1},
              {10, 5, 0}, {18, 0, 1}, {18, 7, 1}, {20, 6, 1}, {23, 3, 1},
              {26, 2, 0}, {29, 4, 1}, {31, 9, 0}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 3, 3);
    std::sort(begin(algo_ans), end(algo_ans));

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s15k5win8) {
    using namespace orange::minimizers;

    auto sequence = std::string{"TGGTGGGCTCCAAAT"};
    auto correct_ans =
        KMers{{81, 0, 1},  {86, 3, 1},  {86, 3, 1},  {86, 3, 1},
              {175, 8, 1}, {175, 8, 1}, {175, 8, 1}, {344, 4, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 5, 8);
    std::sort(begin(algo_ans), end(algo_ans));

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s22k7win13) {
    using namespace orange::minimizers;

    auto sequence = std::string{"ATAGATTTCTCCTAGCTTGTCT"};
    auto correct_ans = KMers{
        {138, 5, 1}, {138, 5, 1}, {138, 5, 1},  {138, 5, 1},   {138, 5, 1},
        {138, 5, 1}, {552, 6, 1}, {2211, 7, 1}, {2555, 13, 0}, {2555, 13, 0}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 7, 13);
    std::sort(begin(algo_ans), end(algo_ans));

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s33k16win20) {
    using namespace orange::minimizers;

    auto sequence = std::string{"GGACGAGCTATCAAGGCAATGTTCAAAATATGG"};
    auto correct_ans =
        KMers{{172211456, 12, 0}, {172211456, 12, 0}, {172211456, 12, 0},
              {172211456, 12, 0}, {172211456, 12, 0}, {250872634, 17, 0},
              {412930212, 2, 0},  {412930212, 2, 0},  {412930212, 2, 0},
              {657729806, 5, 0},  {657729806, 5, 0},  {657729806, 5, 0},
              {794543151, 10, 1}, {794543151, 10, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 16, 20);
    std::sort(begin(algo_ans), end(algo_ans));

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s10k1win1) {
    using namespace orange::minimizers;

    auto sequence = std::string{"GACGACGTTG"};
    auto correct_ans =
        KMers{{0, 1, 0}, {0, 4, 0}, {0, 7, 1}, {0, 8, 1}, {1, 0, 1},
              {1, 2, 0}, {1, 3, 1}, {1, 5, 0}, {1, 6, 1}, {1, 9, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 1, 1);
    std::sort(begin(algo_ans), end(algo_ans));

    EXPECT_EQ(correct_ans, algo_ans);
}