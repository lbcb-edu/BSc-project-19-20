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

    unsigned int make_kmer(unsigned int last_kmer, char next, unsigned int mask) {
        unsigned int kmer = last_kmer << 2;
        kmer = kmer | get_nucleotide_value(next);
        kmer = kmer & mask;
        return kmer;
    }

    unsigned int make_kmer_reverse(unsigned int last_kmer_reverse, char next, unsigned int mask) {
        unsigned int kmer_reverse = last_kmer_reverse >> 2;
        kmer_reverse = kmer_reverse | get_reverse_nucleotide_value(next);
        kmer_reverse = kmer_reverse & mask;
        return kmer_reverse;
    }

    unsigned int make_first_kmer(int j, const char * sequence, unsigned int k) {
        unsigned  int k_mer = 0b0;
        for(unsigned int m = 0; m < k; m++) {
            k_mer = k_mer | get_nucleotide_value(sequence[j + m]);
            if(m + 1 < k)
                k_mer = k_mer << 2;
        }
        return k_mer;
    }

    unsigned int make_first_kmer_reverse(int j_reverse, const char * sequence, unsigned int k) {
        unsigned  int k_mer_reverse = 0b0;
        for(unsigned int m = 0; m < k; m++) {
            k_mer_reverse = k_mer_reverse | get_reverse_nucleotide_value(sequence[j_reverse - m]);
            if(m + 1 < k)
                k_mer_reverse = k_mer_reverse << 2;
        }
        return k_mer_reverse;
    }

    std::tuple<unsigned int, unsigned int, bool> make_kmer_triplet(int j, int j_reverse, const char* sequence, unsigned int k, bool first,
                                                                   unsigned int &last_kmer, unsigned int &last_kmer_reverse,
                                                                   unsigned int mask) {
        unsigned  int k_mer, k_mer_reverse;
        if(first) {
            k_mer = make_first_kmer(j, sequence, k);
            k_mer_reverse = make_first_kmer_reverse(j_reverse, sequence, k);
        } else {
            k_mer = make_kmer(last_kmer, sequence[j + k - 1], mask);
            k_mer_reverse = make_kmer_reverse(last_kmer_reverse, sequence[j_reverse - k + 1], mask);
        }
        last_kmer = k_mer;
        last_kmer_reverse = k_mer_reverse;

        if(k_mer < k_mer_reverse)
            return std::make_tuple(k_mer, j, true);
        else
            return std::make_tuple(k_mer_reverse, j_reverse, false);
    }

    void interior_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers) {
        unsigned int l = k + window_length - 1;
        std::tuple<unsigned int, unsigned int, bool> k_mer_triplet;
        std::tuple<unsigned int, unsigned int, bool> k_mer_triplet_minimal;
        unsigned int last_kmer, last_kmer_reverse;
        unsigned int mask = (1 << (2 * k)) - 1;
        for(int i = 0, i_reverse = sequence_length - 1; i <= sequence_length - l && i_reverse >= l - 1; i++, i_reverse--) {
            for(int j = i, j_reverse = i_reverse; j <= i + l - k && j_reverse >= i_reverse - k; j++, j_reverse--){
                if(i == j) {
                    k_mer_triplet_minimal = make_kmer_triplet(j, j_reverse, sequence, k, true, last_kmer, last_kmer_reverse, mask);
                }
                k_mer_triplet = make_kmer_triplet(j, j_reverse, sequence, k, false, last_kmer, last_kmer_reverse, mask);
                if(std::get<0>(k_mer_triplet) < std::get<0>(k_mer_triplet_minimal))
                    k_mer_triplet_minimal = k_mer_triplet;
            }
            minimizers.emplace(k_mer_triplet_minimal);
        }
    }

    void begin_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers) {
//        for(unsigned int u = 1; u <= sequence_length - k + 1; u++){
//            std::tuple<unsigned int, unsigned int, bool> k_mer_triplet;
//            std::tuple<unsigned int, unsigned int, bool> k_mer_triplet_minimal;
//            std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers;
//            std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers_reverse;
//            unsigned int last_kmer, last_kmer_reverse;
//            unsigned int mask = (1 << (2 * k)) - 1;
//            for(int j = 0, j_reverse = sequence_length - 1; j <= u - 1 && j_reverse >= sequence_length - u; j++, j_reverse--){
//                if(j == 0)
//                    k_mer_triplet_minimal = make_kmer_triplet(j, j_reverse, sequence, k, true, last_kmer, last_kmer_reverse, mask);
//                k_mer_triplet = make_kmer_triplet(j, j_reverse, sequence, k, false, last_kmer, last_kmer_reverse, mask);
//
//                if(std::get<0>(k_mer_triplet) < std::get<0>(k_mer_triplet_minimal))
//                    k_mer_triplet_minimal = k_mer_triplet;
//            }
//            minimizers.emplace(k_mer_triplet_minimal);
//        }
        unsigned int last_kmer, last_kmer_reverse;
        unsigned int mask = (1 << (2 * k)) - 1;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> beg_mers;
        std::tuple<unsigned int, unsigned int, bool> k_mer_triplet_minimal = make_kmer_triplet(0, sequence_length - 1, sequence, k,
                                                                                       true, last_kmer, last_kmer_reverse, mask);
        beg_mers.emplace_back(k_mer_triplet_minimal);
        for(int i = 1; i < sequence_length - k + 1; i ++){
            std::tuple<unsigned int, unsigned int, bool> k_mer_triplet = make_kmer_triplet(i, sequence_length - i - 1, sequence, k,
                                                                                           true, last_kmer, last_kmer_reverse, mask);
            if(std::get<0>(k_mer_triplet) < std::get<0>(beg_mers.back()))
                beg_mers.emplace_back(k_mer_triplet);
        }

        for(auto const& beg_mer : beg_mers)
            minimizers.emplace(beg_mer);
    }

    void end_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers) {
//        for(unsigned int u = 1; u <= sequence_length - k + 1; u++){
//            unsigned int l = k + u - 1;
//            std::tuple<unsigned int, unsigned int, bool> k_mer_triplet;
//            std::tuple<unsigned int, unsigned int, bool> k_mer_triplet_minimal;
//            std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers;
//            std::vector<std::tuple<unsigned int, unsigned int, bool>> k_mers_reverse;
//            unsigned int last_kmer, last_kmer_reverse;
//            unsigned int mask = (1 << (2 * k)) - 1;
//            for(int j = sequence_length - k, j_reverse = k - 1; j >= sequence_length - k - u + 1 && j_reverse <= u + k - 2; j--, j_reverse++){
//                if(j == sequence_length - k)
//                    k_mer_triplet_minimal = make_kmer_triplet(j, j_reverse, sequence, k, true, last_kmer, last_kmer_reverse, mask);
//
//                k_mer_triplet = make_kmer_triplet(j, j_reverse, sequence, k, false, last_kmer, last_kmer_reverse, mask);
//
//                if(std::get<0>(k_mer_triplet) < std::get<0>(k_mer_triplet_minimal))
//                    k_mer_triplet_minimal = k_mer_triplet;
//            }
//
//            minimizers.emplace(k_mer_triplet_minimal);
//        }

        unsigned int last_kmer, last_kmer_reverse;
        unsigned int mask = (1 << (2 * k)) - 1;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> beg_mers;
        std::tuple<unsigned int, unsigned int, bool> k_mer_triplet_minimal = make_kmer_triplet(sequence_length - k, k - 1, sequence, k,
                                                                                               true, last_kmer, last_kmer_reverse, mask);
        beg_mers.emplace_back(k_mer_triplet_minimal);
        for(int i = sequence_length - k - 1; i >= 0; i--){
            std::tuple<unsigned int, unsigned int, bool> k_mer_triplet = make_kmer_triplet(i, sequence_length - i - 1, sequence, k,
                                                                                           true, last_kmer, last_kmer_reverse, mask);
            if(std::get<0>(k_mer_triplet) < std::get<0>(beg_mers.back()))
                beg_mers.emplace_back(k_mer_triplet);
        }

        for(auto const& beg_mer : beg_mers)
            minimizers.emplace(beg_mer);
    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length,
                                                                         unsigned int k,
                                                                         unsigned int window_length){

        std::vector<std::tuple<unsigned int, unsigned int, bool>> ret;

        std::unordered_set<std::tuple<unsigned int, unsigned int, bool>> minimizers_set;


        begin_minimizers(sequence, sequence_length, k, minimizers_set);
        end_minimizers(sequence, sequence_length, k, minimizers_set);
        interior_minimizers(sequence, sequence_length, k, window_length, minimizers_set);

        for(auto minimizer : minimizers_set){
            ret.emplace_back(minimizer);
        }

        return ret;
    }

}