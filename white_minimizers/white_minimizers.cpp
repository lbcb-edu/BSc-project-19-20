#include "white_minimizers.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>

namespace white{
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length,
                                                                         unsigned int k,
                                                                         unsigned int window_length){
        std::vector<std::tuple<unsigned int, unsigned int, bool>> ret;

        std::string radni = sequence;
        int l = window_length + k - 1;
        
        // pravi(numV) i reverse(numR)
        std::vector<int> numV, numR;

        for (auto it : radni) {
            switch (it) {
                case 'C':
                    numV.push_back(0);
                    numR.push_back(3);
                case 'A':
                    numV.push_back(1);
                    numR.push_back(2);
                case 'T':
                    numV.push_back(2);
                    numR.push_back(1);
                case 'G':
                    numV.push_back(3);
                    numR.push_back(0);
            }
        }

        for (int i = 0; i < sequence_length - window_length; i++) {
            for (int j = i; j < i + window_length; j++) {
                unsigned int tmp_low = 4294967295U, index_pos = 0;
                std::vector<std::tuple<unsigned int, unsigned int, bool>> window_calc;
                unsigned int kmer = 0, kmer_rev = 0;
                for (int m = j; m < j + k; m++) {
                    kmer |= numV[m];
                    kmer_rev |= numR[m];
                    kmer <<= 2;
                    kmer_rev <<= 2;
                }
                kmer |= numV[numV.size() - 1];
                kmer_rev |= numR[numR.size() - 1];
                if (kmer < tmp_low) {
                    tmp_low = kmer;
                    index_pos = j;
                    window_calc.push_back(std::make_tuple(kmer, index_pos, 1));
                } else if (kmer_rev < tmp_low) {
                    tmp_low = kmer_rev;
                    index_pos = j;
                    window_calc.push_back(std::make_tuple(kmer_rev, index_pos, 0));
                }
                std::sort(window_calc.begin(), window_calc.end());
                ret.push_back(window_calc[0]);
            }
        }
        return ret;
    }
}

