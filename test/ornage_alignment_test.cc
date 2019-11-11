#include <gtest/gtest.h>

#include "orange_alignment.h"

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