#include <string_view>
#include <functional>
#include <queue>
#include <deque>
#include <set>

#include "orange_minimizers.h"

namespace orange {
namespace minimizers {

using MinimsPos = std::set<std::size_t>;

class MinKPosQueue {
public:
    MinKPosQueue(KMers const& kmers) : kmers_{kmers} {}

    void push(std::size_t pos) {
        raw_queue_.push(pos);
        while (!min_queue_.empty() && kmers_[min_queue_.back()] > kmers_[pos])
            min_queue_.pop_back();
        min_queue_.push_back(pos);
    }

    void pop() {
        if (raw_queue_.front() == min_queue_.front())
            min_queue_.pop_front();
        raw_queue_.pop();
    }

    bool empty() const { return raw_queue_.empty(); }

    std::size_t pos_min() const { return min_queue_.front(); }

private:
    KMers const& kmers_;
    std::queue<std::size_t> raw_queue_;
    std::deque<std::size_t> min_queue_;
};

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
KMers generateKMersFlagged(std::string_view sequence, std::uint32_t k) {
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

    kmers.reserve(2 * (sequence.size() - k + 1));
    auto append_kmers = [&kmers, &curr_kmer, &comp_kmer](std::uint32_t pos) {
        kmers.emplace_back(curr_kmer, pos, false);
        kmers.emplace_back(comp_kmer, pos, true);
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

/**
 *
 * @param begin KMer range beginning
 * @param end KMer range ending
 * @param scan_start
 * @param win_len
 */
void markMinimizers(KMers const& kmers, MinimsPos& minims_pos,
                    std::size_t const begin_pos, std::size_t const end_pos,
                    std::uint32_t scan_start, std::uint32_t scan_end,
                    std::uint32_t win_len) {
    auto queue = MinKPosQueue{kmers};
    auto win_end = scan_start + win_len;

    auto mark_and_pop = [&queue, &minims_pos] {
        minims_pos.insert(queue.pos_min());
        queue.pop(), queue.pop();  // pop original, pop complement
    };

    for (auto pos = begin_pos; pos != end_pos; ++pos) {
        if (std::get<1>(kmers[pos]) >= win_end) {
            mark_and_pop();
            ++win_end;
        }

        queue.push(pos);
    }

    /* clang-format off */

    while (!queue.empty() && win_end < scan_end)
        mark_and_pop(), ++win_end;
    if (!queue.empty())
        mark_and_pop();

    /* clang-format on */
}

KMers minimizers(char const* sequence, std::uint32_t sequence_length,
                 std::uint32_t k, std::uint32_t window_length) {
    auto sequence_view = std::string_view{sequence, sequence_length};
    auto all_kmers = generateKMersFlagged(sequence_view, k);
    auto minims_pos = MinimsPos{};
    auto minims = KMers{};

    markMinimizers(all_kmers, minims_pos, 0, all_kmers.size(), 0,
                   sequence_length, window_length);

    // end minimizers, front
    for (auto u = std::uint32_t{1}; u < window_length; ++u)
        markMinimizers(all_kmers, minims_pos, 0, 2 * u, 0, u, u);

    // end minimizers, back
    if (k < window_length) {
        auto rev_start = all_kmers.size() - 2;
        for (auto u = k; u < window_length; ++u) {
            markMinimizers(all_kmers, minims_pos, rev_start, all_kmers.size(),
                           sequence_length - u, sequence_length, window_length);
            rev_start -= 2;
        }
    }

    /* clang-format off */
    for (auto const& pos : minims_pos)
        minims.push_back(std::move(all_kmers[pos]));
    /* clang-format on */
    return minims;
}
}  // namespace minimizers
}  // namespace orange