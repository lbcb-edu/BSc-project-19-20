//
// Created by nera on 22. 11. 2019..
//
#include <gtest/gtest.h>
#include "../white_alignment/white_alignment.cpp"


TEST(PairwiseAlignment, globalType) {
    EXPECT_EQ(white::PairwiseAlignment("GCATGCU", 7, "GATTACA", 7, white::AlignmentType::kGlobal, 1, -1, -1), 0);
}

TEST(PairwiseAlignment, globalType2) {
    ::std::string cig2;
    unsigned int tb2;
    EXPECT_EQ(white::PairwiseAlignment("GCAT", 4, "GCT", 3, white::AlignmentType::kGlobal, 1, -1, -1, cig2, tb2), 2);
    EXPECT_EQ(tb2, 0);
    EXPECT_EQ(cig2, "2M1I1M");
}

TEST(PairwiseAlignment, localType) {
    EXPECT_EQ(white::PairwiseAlignment("ACCTAAGG", 8, "GGCTCAATCA", 10, white::AlignmentType::kLocal, 2, -1, -2), 6);
}

TEST(PairwiseAlignment, semiGlobalType) {
    std::string cig;
    unsigned int tb;
    EXPECT_EQ (white::PairwiseAlignment("TCCG", 4, "ACTCCGAT", 8, white::AlignmentType::kSemiGlobal, 4, -1, -2,
            cig, tb), 16);
    EXPECT_EQ (cig, "4M");
    EXPECT_EQ (tb, 2);
}

TEST(PairwiseAlignment, globalTypeCigar) {
    ::std::string cigarStr;
    unsigned int tBegin;
    int value = PairwiseAlignment("CAGCACTTGGATTCTCGG", 18,
                          "CAGCGTGG", 8,
                          white::AlignmentType::kGlobal,
                          1,
                          -1,
                          -1,
                          cigarStr,
                          tBegin);
    EXPECT_EQ (cigarStr, "3M2I1M3I1M4I1M1I2M");
    EXPECT_EQ (tBegin, 0);
}


int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}