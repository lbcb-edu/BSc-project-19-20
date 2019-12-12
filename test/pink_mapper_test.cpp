#include "gtest/gtest.h"
#include "pink_alignment.h"
#include "pink_minimizers.h"

#include <algorithm>

TEST(PairwiseAlignment, global) {
    std::string cigar;
    unsigned int target_begin;
    EXPECT_EQ(pink::pairwise_alignment("TGCATAT", 7, "ATCCGAT", 7, pink::global, 1, -1, -1, cigar, target_begin), 0);

}

TEST(PairwiseAlignment, semi_global) {
    std::string cigar;
    unsigned int target_begin;
    EXPECT_EQ (pink::pairwise_alignment ("ACCCAAGGG", 9, "GGCTCCATTA", 10, pink::semi_global, 1, -1, -1, cigar, target_begin), 2);
    EXPECT_EQ (target_begin, 0);
}

TEST(PairwiseAlignment, local) {
    std::string cigar;
    unsigned int target_begin = 0;
    EXPECT_EQ(6, pink::pairwise_alignment("ACCTAAGG", 8, "GGCTCAATCA", 10, pink::local, 2, -1, -2, cigar, target_begin));
    EXPECT_EQ(target_begin, 2);
}

TEST (Minimizers, TestTGACGTACATGGACA) {
std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = pink::minimizers("TGACGTACATGGACA", 15, 3, 4);

std::sort(minimizers.begin(), minimizers.end());

EXPECT_EQ(std::get<0>(minimizers.at(0)), 1);
EXPECT_EQ(std::get<1>(minimizers.at(0)), 11);
EXPECT_EQ(std::get<2>(minimizers.at(0)), 0);
EXPECT_EQ(std::get<0>(minimizers.at(1)), 6);
EXPECT_EQ(std::get<1>(minimizers.at(1)), 7);
EXPECT_EQ(std::get<2>(minimizers.at(1)), 1);
EXPECT_EQ(std::get<0>(minimizers.at(2)), 14);
EXPECT_EQ(std::get<1>(minimizers.at(2)), 4);
EXPECT_EQ(std::get<2>(minimizers.at(2)), 0);
EXPECT_EQ(std::get<0>(minimizers.at(3)), 17);
EXPECT_EQ(std::get<1>(minimizers.at(3)), 12);
EXPECT_EQ(std::get<2>(minimizers.at(3)), 1);
EXPECT_EQ(std::get<0>(minimizers.at(4)), 19);
EXPECT_EQ(std::get<1>(minimizers.at(4)), 2);
EXPECT_EQ(std::get<2>(minimizers.at(4)), 1);
EXPECT_EQ(std::get<0>(minimizers.at(5)), 45);
EXPECT_EQ(std::get<1>(minimizers.at(5)), 0);
EXPECT_EQ(std::get<2>(minimizers.at(5)), 1);

}