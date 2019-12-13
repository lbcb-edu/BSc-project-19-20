#include <vector>
#include <tuple>
#include <algorithm>
#include <unordered_set>

#include <iostream>

#include "pink_minimizers.h"

namespace std {
        template<>
        struct hash<tuple<unsigned int, unsigned int, bool >> {
            inline size_t operator()(const tuple<unsigned int, unsigned int, bool> &value) const {
                return hash<int>()(get<0>(value) ^ get<1>(value) ^ get<2>(value));
            }
        };
}

namespace pink {

    int get_nucleotide_value(char nucleotide) {
        switch (nucleotide) {
            case 'C' :
                return 0;
            case 'A' :
                return 1;
            case 'T' :
                return 2;
            case 'G' :
                return 3;
            default:
                return 1; // A and T often occur more frequently - A chosen for default
        }
    }

    int get_reverse_nucleotide_value(char nucleotide) {
        switch (nucleotide) {
            case 'C' :
                return 3;
            case 'A' :
                return 2;
            case 'T' :
                return 1;
            case 'G' :
                return 0;
            default:
                return 1; // A and T often occur more frequently - A chosen for default
        }
    }

    void make_kmer_triples(std::tuple<unsigned int, unsigned int, bool>& k_mer_triple, std::tuple<unsigned int, unsigned int, bool>& k_mer_triple_reverse,
                            int j, int j_reverse, const char* sequence, unsigned int k) {
        unsigned  int k_mer = 0b0;
        unsigned  int k_mer_reverse = 0b0;
        for(unsigned int m = 0; m < k; m++) {
            k_mer = k_mer | get_nucleotide_value(sequence[j + m]);
            if(m + 1 < k)
                k_mer = k_mer << 2;
        }
        for(unsigned int m = 0; m < k; m++) {
            k_mer_reverse = k_mer_reverse | get_reverse_nucleotide_value(sequence[j_reverse - m]);
            if(m + 1 < k)
                k_mer_reverse = k_mer_reverse << 2;
        }

        k_mer_triple = std::make_tuple(k_mer, j, true);
        k_mer_triple_reverse = std::make_tuple(k_mer_reverse, j_reverse, false);
    }

    void insert(std::vector<std::tuple<unsigned int, unsigned int, bool>>& k_mers, std::vector<std::tuple<unsigned int, unsigned int, bool>>& k_mers_reverse,
                     std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers) {
        std::tuple<unsigned int, unsigned int, bool> minimal;
        if(*std::min_element(k_mers.begin(), k_mers.end()) <= *std::min_element(k_mers_reverse.begin(), k_mers_reverse.end()))
            minimal = *std::min_element(k_mers.begin(), k_mers.end());
        else
            minimal = *std::min_element(k_mers_reverse.begin(), k_mers_reverse.end());
        minimizers.emplace(minimal);
    }

    void interior_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers) {
        unsigned int l = k + window_length - 1;
        std::tuple<unsigned int, unsigned int, bool> k_mer_triple;
        std::tuple<unsigned int, unsigned int, bool> k_mer_triple_reverse;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers_reverse;
        for(int i = 0, i_reverse = sequence_length - 1; i <= sequence_length - l && i_reverse >= l - 1; i++, i_reverse--) {
            k_mers.clear();
            k_mers_reverse.clear();
            for(int j = i, j_reverse = i_reverse; j <= i + l - k && j_reverse >= i_reverse - k; j++, j_reverse--){
                make_kmer_triples(k_mer_triple, k_mer_triple_reverse, j, j_reverse, sequence, k);

                k_mers.push_back(k_mer_triple);
                k_mers_reverse.push_back(k_mer_triple_reverse);
            }

            insert(k_mers, k_mers_reverse, minimizers);
        }
    }

    void begin_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers) {
        for(unsigned int u = 1; u <= sequence_length - k + 1; u++){
            std::tuple<unsigned int, unsigned int, bool> k_mer_triple;
            std::tuple<unsigned int, unsigned int, bool> k_mer_triple_reverse;
            std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers;
            std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers_reverse;
            k_mers.clear();
            k_mers_reverse.clear();
            for(int j = 0, j_reverse = sequence_length - 1; j <= u - 1 && j_reverse >= sequence_length - u; j++, j_reverse--){
                make_kmer_triples(k_mer_triple, k_mer_triple_reverse, j, j_reverse, sequence, k);

                k_mers.push_back(k_mer_triple);
                k_mers_reverse.push_back(k_mer_triple_reverse);
            }

            insert(k_mers, k_mers_reverse, minimizers);
        }
    }

    void end_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers) {
        for(unsigned int u = 1; u <= sequence_length - k + 1; u++){
            unsigned int l = k + u - 1;
            std::tuple<unsigned int, unsigned int, bool> k_mer_triple;
            std::tuple<unsigned int, unsigned int, bool> k_mer_triple_reverse;
            std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers;
            std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers_reverse;
            k_mers.clear();
            k_mers_reverse.clear();
            for(int j = sequence_length - k, j_reverse = k - 1; j >= sequence_length - k - u + 1 && j_reverse <= u + k - 2; j--, j_reverse++){
                make_kmer_triples(k_mer_triple, k_mer_triple_reverse, j, j_reverse, sequence, k);

                k_mers.push_back(k_mer_triple);
                k_mers_reverse.push_back(k_mer_triple_reverse);
            }

            insert(k_mers, k_mers_reverse, minimizers);
        }
    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length,
                                                                         unsigned int k,
                                                                         unsigned int window_length){

        std::vector<std::tuple<unsigned int, unsigned int, bool>> ret;

        std::unordered_set<std::tuple<unsigned int, unsigned int, bool>> minimizers_set;


        interior_minimizers(sequence, sequence_length, k, window_length, minimizers_set);
        begin_minimizers(sequence, sequence_length, k, minimizers_set);
        end_minimizers(sequence, sequence_length, k, minimizers_set);

        for(auto minimizer : minimizers_set){
            ret.emplace_back(minimizer);
        }

        return ret;
    }

}

