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
    auto correct_ans = KMers{{4, 1, 0},  {4, 12, 0},  {13, 3, 0},
                             {13, 6, 0}, {17, 11, 0}, {33, 0, 0}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 3, 5);

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s15k1win1) {
    using namespace orange::minimizers;

    auto sequence = std::string{"TACTAGCTTCCATGT"};
    auto correct_ans =
        KMers{{0, 0, 1}, {0, 1, 0},  {0, 3, 1},  {0, 4, 0},  {0, 7, 1},
              {0, 8, 1}, {0, 11, 0}, {0, 12, 1}, {0, 14, 1}, {1, 2, 0},
              {1, 5, 1}, {1, 6, 0},  {1, 9, 0},  {1, 10, 0}, {1, 13, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 1, 1);

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s64k16win32) {
    using namespace orange::minimizers;

    auto sequence = std::string{
        "ATCTGCGGCAGTTGTATCCGATGTTTGAGACCTAGTCAGTGTCTTGATGGACGTTCCATCAACG"};
    auto correct_ans =
        KMers{{20352065, 10, 1},  {81408263, 9, 1},   {114610308, 39, 1},
              {325633054, 8, 1},  {822580808, 1, 1},  {932858860, 0, 0},
              {973526224, 43, 1}, {974902480, 46, 0}, {1871610293, 48, 1}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 16, 32);

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s15k7win5) {
    using namespace orange::minimizers;

    auto sequence = std::string{"AAAACAGATGGGGTT"};
    auto correct_ans =
        KMers{{18, 0, 0}, {72, 1, 0}, {291, 2, 0}, {340, 8, 1}, {1166, 3, 0}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 7, 5);

    EXPECT_EQ(correct_ans, algo_ans);
}

TEST(Minimizers, s64k16win8) {
    using namespace orange::minimizers;

    auto sequence = std::string{
        "GAACTCTGCTAGGAAACCTTGGCCTGCATAATCAAAAGTAATGTCGGCATTGGACACTCGCTGA"};
    auto correct_ans =
        KMers{{11595172, 33, 0},  {25073555, 13, 0},  {46380691, 34, 0},
              {100294220, 14, 0}, {125422081, 1, 0},  {184382024, 4, 1},
              {185522767, 35, 0}, {218149101, 29, 0}, {249188257, 39, 0},
              {331873852, 28, 1}, {401176880, 15, 0}, {536408978, 23, 1},
              {640603193, 46, 1}, {671480414, 10, 0}, {1048868472, 48, 0},
              {2178839168, 0, 0}};

    auto algo_ans = minimizers(sequence.c_str(), sequence.size(), 16, 8);

    EXPECT_EQ(correct_ans, algo_ans);
}