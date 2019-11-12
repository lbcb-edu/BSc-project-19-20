#include "gtest/gtest.h"
#include "pink_alignment.h"

TEST(PairwiseAlignment, global) {
    std::string cigar;
    int target_begin;
    EXPECT_EQ(pink::pairwise_alignment("TGCATAT", 7, "ATCCGAT", 7, pink::global, 1, -1, -1, cigar, target_begin), 0);

}

TEST(PairwiseAlignment, semi_global) {
    std::string cigar;
    int target_begin;
    EXPECT_EQ (pink::pairwise_alignment ("ACCTAAGG", 8, "GGCTCAATCA", 10, orange::AlignmentType::local, 2, -1, -2, cigar, target_start), 6);
    EXPECT_EQ (target_start, 2);
}

TEST(PairwiseAlignment, local) {
    std::string cigar;
    unsigned int target_begin = 0;
    EXPECT_EQ(6, blue::pairwise_alignment("ACCTAAGG", 8, "GGCTCAATCA", 10, blue::local, 2, -1, -2, cigar, target_begin));
    EXPECT_EQ(target_begin, 2);
}