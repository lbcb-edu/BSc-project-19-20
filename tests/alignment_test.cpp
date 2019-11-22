//
// Created by nera on 22. 11. 2019..
//
#include <gtest/gtest.h>
#include "../white_alignment/white_alignment.cpp"


TEST(PairwiseAlignment, globalType) {
    EXPECT_EQ(white::PairwiseAlignment("GCATGCU", 7, "GATTACA", 7, white::AlignmentType::kGlobal, 1, -1, -1), 0);
}

TEST(PairwiseAlignment, globalType2) {
    EXPECT_EQ(white::PairwiseAlignment("GCAT", 4, "GCT", 3, white::AlignmentType::kGlobal, 1, -1, -1), 2);
}

TEST(PairwiseAlignment, localType) {
    EXPECT_EQ(white::PairwiseAlignment("ACCTAAGG", 8, "GGCTCAATCA", 10, white::AlignmentType::kLocal, 2, -1, -2), 6);
}

TEST(PairwiseAlignment, semiGlobalType) {
    EXPECT_EQ (white::PairwiseAlignment("TCCG", 4, "ACTCCGAT", 8, white::AlignmentType::kSemiGlobal, 4, -1, -2), 12);
}

int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}