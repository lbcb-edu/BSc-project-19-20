//
// Created by nera on 16. 01. 2020..
//
#include <gtest/gtest.h>
#include "../white_minimizers/white_minimizers.cpp"


TEST(WhiteMinimizers, TestWhiteMinimizer1) {
    std::string s1 = "TCAGGAAGAAGCAGA";

    int k = 3;
    int window_length = 4;

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;

    minimizers.push_back(std::make_tuple(10, 2, true));
    minimizers.push_back(std::make_tuple(2, 5, true));
    minimizers.push_back(std::make_tuple(2, 8, true));
    minimizers.push_back(std::make_tuple(8, 12, true));

    EXPECT_EQ (minimizers, white::minimizers(s1.c_str(), s1.size(), k, window_length));
}

TEST(WhiteMinimizers, TestWhiteMinimizer2) {
    std::string s2 = "GTCATGCACGTTCAC";

    int k = 3;
    int window_length = 4;
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;

    minimizers.push_back(std::make_tuple(1, 9, false));
    minimizers.push_back(std::make_tuple(6, 7, true));
    minimizers.push_back(std::make_tuple(14, 2, false));

    EXPECT_EQ (minimizers, white::minimizers(s2.c_str(), s2.size(), k, window_length));
}

TEST(WhiteMinimizers, TestWhiteMinimizer3) {
    std::string s3 = "GTCATGCACGTTCAC"; 

    int k = 6;
    int window_length = 5;

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;


    minimizers.push_back(std::make_tuple(110, 6, false));
    minimizers.push_back(std::make_tuple(441, 5, false));
    minimizers.push_back(std::make_tuple(445, 7, true));

    EXPECT_EQ (minimizers, white::minimizers(s3.c_str(), s3.size(), k, window_length));
}

int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}