#include <string_view>
#include <functional>
#include <algorithm>
#include <queue>
#include <deque>

#include "orange_minimizers.h"

namespace orange {
namespace minimizers {

/**
 * @brief Represents a KMer in sequence
 *      flagged as minimizer in at least
 *      one search window.
 */
struct KMerMarked {
    KMer const raw_;
    bool marked_;
};

/**
 * @brief std::vector of @ref orange::minimizers::KMerMarked
 */
using KMersMarked = std::vector<KMerMarked>;

using KMMIterator = KMersMarked::iterator;

/**
 * @brief K-Mare queue with O(1) minimum element access
 */
class MinKQueue {
public:
    /**
     * @brief Adds a Kmer to the back of the queue
     *
     * @param kmare
     */
    void push(KMerMarked& kmer) {
        raw_queue_.push(kmer);
        while (!min_queue_.empty() && min_queue_.back().get().raw_ > kmer.raw_)
            min_queue_.pop_back();

        min_queue_.push_back(raw_queue_.back());
    }

    /**
     * @brief Removes an element from the front of the queue
     */
    void pop() {
        if (!min_queue_.empty() &&
            raw_queue_.front().get().raw_ == min_queue_.front().get().raw_)
            min_queue_.pop_front();
        raw_queue_.pop();
    }

    /**
     * @brief Checks if MinKQueue is empty
     *
     * @return true if it's empty
     */
    bool empty() const { return raw_queue_.empty(); }

    /**
     * @brief returns the KMerMakred reference with the smallest value,
     *      aka. minimizer of the sliding window
     *
     * @return KMer
     */
    auto min() { return min_queue_.front(); }

private:
    std::queue<std::reference_wrapper<KMerMarked>> raw_queue_;
    std::deque<std::reference_wrapper<KMerMarked>> min_queue_;
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
KMersMarked generateKMersFlagged(std::string_view sequence, std::uint32_t k) {
    auto kmers = KMersMarked{};
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

/**
 *
 * @param begin KMer range beginning
 * @param end KMer range ending
 * @param scan_start
 * @param win_len
 */
void makrMinimizers(KMMIterator const begin, KMMIterator const end,
                    std::uint32_t scan_start, std::uint32_t win_len) {
    auto queue = MinKQueue{};
    auto win_end = scan_start + win_len;

    auto mark_and_pop = [&queue] {
        queue.min().get().marked_ = true;
        queue.pop(), queue.pop();  // pop original, pop complement
    };

    for (auto curr = begin; curr != end; ++curr) {
        auto& kmer = *curr;
        if (std::get<1>(kmer.raw_) >= win_end) {
            mark_and_pop();
            ++win_end;
        }

        queue.push(kmer);
    }

    mark_and_pop();
}

KMers minimizers(char const* sequence, std::uint32_t sequence_length,
                 std::uint32_t k, std::uint32_t window_length) {
    auto sequence_view = std::string_view{sequence, sequence_length};
    auto kmers = generateKMersFlagged(sequence_view, k);
}
}  // namespace minimizers
}  // namespace orange