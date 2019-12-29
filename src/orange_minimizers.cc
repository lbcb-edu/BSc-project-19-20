#include <string_view>
#include <functional>
#include <queue>
#include <deque>

#include "orange_minimizers.h"

namespace orange {
namespace minimizers {

using MinimsPos = std::vector<std::size_t>;

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
 * @return KMersMarked all k-mers from the sequence
 *      with a negative flag value (they are not included by default)
 */
KMers generateKMers(char const* sequence, std::uint32_t sequence_length,
                    std::uint32_t k) {
    auto kmers = KMers{};
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

    kmers.reserve(sequence_length - k + 1);
    auto append_kmers = [&kmers, &curr_kmer, &comp_kmer](std::uint32_t pos) {
        if (curr_kmer <= comp_kmer)
            kmers.emplace_back(curr_kmer, pos, false);
        else
            kmers.emplace_back(comp_kmer, pos, true);
    };

    for (auto i = std::uint32_t{0}; i < k; ++i) {
        shift_add_curr(sequence[i]);
        shift_add_comp(sequence[i]);
    }

    append_kmers(0);
    for (auto i = k; i < sequence_length; ++i) {
        shift_add_curr(sequence[i]);
        shift_add_comp(sequence[i]);

        append_kmers(i - k + 1);
    }

    return kmers;
}

KMers minimizers(char const* sequence, std::uint32_t sequence_length,
                 std::uint32_t k, std::uint32_t window_length) {
    auto kmers = generateKMers(sequence, sequence_length, k);
    auto minims = KMers{};
    auto i = std::size_t{1};

    // Fill imaginary min_queue
    minims.push_back(kmers.front());
    for (auto j = std::size_t{1}; j < window_length; ++j)
        if (std::get<0>(minims.back()) > std::get<0>(kmers[i])) {
            minims.push_back(std::move(kmers[i]));
            i = j + 1;
        }

    /* clang-format off */
    // Simulate min_queue on minims vector
    auto prev_peek_pos = minims.size() - 1;
    while (i < kmers.size()) {
        while (minims.back() != minims[prev_peek_pos] &&
               std::get<0>(minims.back()) > std::get<0>(kmers[i]))
            minims.pop_back();
        minims.push_back(std::move(kmers[i++]));

        if (std::get<0>(minims.back()) <
            std::get<0>(minims[prev_peek_pos]))
                prev_peek_pos = minims.size() - 1;
        if (std::get<1>(minims.back()) -
            std::get<1>(minims[prev_peek_pos]) == window_length)
                ++prev_peek_pos;
    }
    /* clang-format on */

    return minims;
}
}  // namespace minimizers
}  // namespace orange