#include "../vendor/googletest/googletest/include/gtest/gtest.h"
#include "../brown_alignment/brown_alignment.hpp"
#include "../brown_minimizers/brown_minimizers.hpp"
#include <string>
#include <vector>
#include <tuple>

namespace {
    TEST(BrownAlignmentTest, LocalAlignmentTest) {
        char * query = "ACCTAAGG";
        char * target = "GGCTCAATCA";

        std::string cigar;
        unsigned int target_begin = 0;

        EXPECT_EQ(6, brown::pairwise_alignment(query, 8, target, 10, brown::AlignmentType::local, 2, -1, -2, cigar, target_begin));
        EXPECT_EQ(target_begin, 2);
    }

    TEST(BrownAlignmentTest, GlobalAlignmentTest){
        char * query = "CTCTGTTCG";
        char * target = "CGTATCTTGA";
        std::string cigar;
        unsigned int target_begin = 0;

        EXPECT_EQ(-5, brown::pairwise_alignment(query, 9, target, 10, brown::AlignmentType::global, 0, -1, -1, cigar, target_begin));
    }

    TEST(BrownAlignmentTest, SemiGlobalAlignmentTest){
        char * query = "TCCG";
        char * target = "ACTCCGAT";

        std::string cigar;
        unsigned int target_begin = 0;

        EXPECT_EQ(4, brown::pairwise_alignment(query, 4, target, 8, brown::AlignmentType::semi_global, 1, -2, -1, cigar, target_begin));
        EXPECT_EQ(target_begin, 2);
    }

    TEST(BrownMinimizersTest, minimizersTest) {
        char * sequence = "TGACGTACAT";

        std::vector<std::tuple<unsigned int, unsigned int, bool>> result = brown::minimizers(sequence, 10, 3, 4);
        // for (int i = 0; i < result.size(); i++)
        // {
        //     std::cout << std::get<0>(result[i]) << "\n";
        // }
        
        EXPECT_EQ(2, result.size());

    }
}