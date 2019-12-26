#include <string_view>
#include <algorithm>

#include "orange_minimizers.h"

namespace orange {
namespace minimizers {

/**
 * @brief Represents a KMer in sequence
 *      flagged as minimizer in at least
 *      one search window.
 */
using KMerFlagged = std::pair<KMer, bool>;
/**
 * @brief std::vector of @ref orange::minimizers::KMerFlagged
 */
using KMersFlagged = std::vector<KMerFlagged>;

/**
 * @brief Converts bases to their masks
 *
 * @param base
 * @return mask
 */
std::uint32_t getMask(char const base) {
    switch (base) {
        case 'A':
            return 1;
        case 'C':
            return 0;
        case 'G':
            return 3;
        case 'T':
            return 2;
        default:
            throw std::invalid_argument("Invalid sequence!");
    }
}

/**
 * @brief Returns a complement of a base
 *          in a number format.
 *
 * @param mask
 * @return complement
 */
std::uint32_t complementMask(std::uint32_t const mask) { return 3 - mask; }

/**
 * @param sequence scanned sequence
 * @param k k-mer length
 *
 * @return KMersFlagged all k-mers from the sequence
 *      with a negative flag value (they are not included by default)
 */
KMersFlagged generateKMersFlagged(std::string_view sequence, std::uint32_t k) {
    auto kmers = KMersFlagged{};
    auto curr_kmer = std::uint32_t{0};
    auto comp_kmer = std::uint32_t{0};
    auto mask = (k < 16 ? ((1 << (2 * k)) - 1) : KMerValMax);

    auto shift_add_curr = [&curr_kmer, &mask](char const& base) {
        curr_kmer <<= 2, curr_kmer &= mask;
        curr_kmer |= getMask(base);
    };

    auto shift_add_comp = [&comp_kmer, &k](char const& base) {
        auto base_mask = complementMask(getMask(base));
        comp_kmer >>= 2, comp_kmer |= base_mask << (2 * k - 2);
    };

    kmers.reserve(sequence.size() - k + 1);
    auto append_kmers = [&kmers, &curr_kmer, &comp_kmer](std::uint32_t pos) {
        kmers.push_back({{curr_kmer, pos, false}, false});
        kmers.push_back({{comp_kmer, pos, false}, false});
    };

    for (auto i = std::uint32_t{0}; i < k; ++i) {
        shift_add_curr(sequence[i]);
        shift_add_comp(sequence[i]);
    }

    append_kmers(0);
    for (auto i = k; i < sequence.size(); ++i) {
        shift_add_curr(sequence[i]);
        shift_add_comp(sequence[i]);

        append_kmers(i - k + 1);
    }

    return kmers;
}

KMers minimizers(char const* sequence, std::uint32_t sequence_length,
                 std::uint32_t k, std::uint32_t window_length) {

    auto sequence_view = std::string_view{sequence, sequence_length};
    auto kmers = generateKMersFlagged(sequence_view, k);



}
}  // namespace minimizers
}  // namespace orange