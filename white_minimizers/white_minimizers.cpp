#include "white_minimizers.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <set>
#include <map>

namespace white{

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length,
                                                                         unsigned int k, unsigned int window_length) {
        
        std::vector<std::tuple<unsigned int, unsigned int, bool>> ret;

        std::set<std::tuple<unsigned int, unsigned int, bool>> minimizer_set;
        std::map<char, int> val = {{'A', 1}, {'T', 2}, {'G', 3}, {'C', 0}};
        std::map<char, int> rev_val = {{'A', 2}, {'T', 1}, {'G', 0}, {'C', 3}};

        unsigned int len = k + window_length - 1;
        std::tuple<unsigned int, unsigned int, bool> kmer, rev_kmer;

        for (int i = 0, i_rev = sequence_length - 1; i <= sequence_length - len && i_rev  >= len - 1; i++, i_rev--) {

            std::tuple<unsigned int, unsigned int, bool> minimal;
            for (int j = i, j_rev = i_rev; j <= i + len - k && j_rev >= i_rev - k; j++, j_rev--) {
                unsigned int k_mer = 0;
                unsigned int k_mer_rev = 0;
                for (unsigned int m = 0; m < k; m++) {
                    k_mer |= val[sequence[j + m]];
                    k_mer_rev |= rev_val[sequence[j_rev - m]];
                    if (m + 1 < k) {
                        k_mer <<= 2;
                        k_mer_rev <<= 2;
                    }
                }
                kmer = std::make_tuple(k_mer, j, true);
                rev_kmer = std::make_tuple(k_mer_rev, j_rev, false);
                minimal = std::min(kmer, rev_kmer);
            }
            minimizer_set.emplace(minimal);
        }

        for (auto minimizer : minimizer_set) {
            ret.emplace_back(minimizer);
        }
        
        return ret;
    }
}

