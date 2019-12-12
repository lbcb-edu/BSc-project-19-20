#include "orange_minimizers.h"

#include <string_view>
#include <functional>
#include <queue>
#include <deque>
#include <set>

namespace orange {
namespace minimizers {

using KMerVal = std::uint32_t;

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
    void push(KMer kmer) {
        raw_queue_.push(kmer);
        while (!min_queue_.empty() &&
               std::get<0>(min_queue_.back().get()) > std::get<0>(kmer))
            min_queue_.pop_back();

        min_queue_.push_back(raw_queue_.back());
    }

    /**
     * @brief Removes an element from the front of the queue
     */
    void pop() {
        if (!min_queue_.empty() && std::get<1>(raw_queue_.front()) >=
                                       std::get<1>(min_queue_.front().get()))
            min_queue_.pop_front();
        raw_queue_.pop();
    }

    /**
     * @brief Checks if MinKQueue is empty
     *
     * @return true if it's empty
     */
    bool empyt() const { return raw_queue_.empty(); }

    /**
     * @brief Returnt the Kmer with the smallest value,
     *      aka. minimizer of the sliding window
     *
     * @return Kmer
     */
    KMer min() const { return min_queue_.front().get(); }

private:
    std::queue<KMer> raw_queue_;  //<<< Contains all the kmers
    std::deque<std::reference_wrapper<KMer>>
        min_queue_;  //<<< Contains smallest kmers
};

/**
 * @brief Get the Complement object
 *
 * @param base Nucleoid base (A, T, C, G)
 * @return char a complement (T, A, G, C)
 */
char getComplement(char const base) {
    switch (base) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        default:
            throw std::invalid_argument("Invalid sequence!");
    }
}

/**
 * @brief Get the Mask object
 *
 * @param base  (A, C, G, T)
 * @return std::uint32_t (0, 1, 2, 3)
 */
std::uint32_t getMask(char const base) {
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            throw std::invalid_argument("Invalid sequence!");
    }
}

/**
 * @brief Streams consecutive KMers over a sequence
 */
class KMerStream {
public:
    KMerStream(std::string_view seq, std::uint32_t k, std::uint32_t begin,
               std::uint32_t end)
        : seq_{seq},
          k_{k},
          pos_{begin},
          begin_{begin},
          end_{end},
          kmer_{0},
          comp_kmer_{0} {

        auto shift = int{2} * k - 2;  // Each base holds up two bits
        for (auto i = std::uint32_t{begin_}; i < begin_ + k_; ++i) {
            kmer_ |= getMask(seq_[i]) << shift;
            comp_kmer_ |= getMask(getComplement(seq_[i])) << shift;

            shift -= 2;
        }

        // First 2 * k bits are set to 1
        if (k < 16)
            mask_ = (1 << (2 * k)) - 1;
        else
            mask_ = std::numeric_limits<std::uint32_t>::max();
    }

    void shiftKMer() {
        if (pos_ + k_ >= end_)
            return;

        // Shift and ignore bits
        kmer_ <<= 2, comp_kmer_ <<= 2;
        kmer_ &= mask_, comp_kmer_ &= mask_;

        kmer_ |= getMask(seq_[pos_ + k_]);
        comp_kmer_ |= getMask(getComplement(seq_[pos_ + k_]));

        ++pos_;
    }

    std::uint32_t pos() const { return pos_; }

    bool hasShift() const { return pos_ + k_ < end_; }

    KMer getCurrKMer() const { return KMer{kmer_, pos_, false}; }

    KMer getCurrCompKMer() const { return KMer{comp_kmer_, pos_, true}; }

private:
    std::string_view seq_;  //<<< Sequence

    std::uint32_t k_;  //<<< KMer length

    std::uint32_t pos_;  //<<< Refers to the current symbol in a sequence

    std::uint32_t begin_;
    std::uint32_t end_;

    std::uint32_t mask_;  //<<< Used for ignoring leftover bits in a shift

    KMerVal kmer_;       //<<< KMer on the original strand
    KMerVal comp_kmer_;  //<<< KMer on the complement strand
};

KMers findMinimizers(std::string_view seq, std::uint32_t k,
                     std::uint32_t win_len, std::uint32_t begin,
                     std::uint32_t end) {
    auto kmer_stream = KMerStream{seq, k, begin, end};
    auto queue = MinKQueue{};

    auto minimizers = KMers{};

    auto win_end = begin +  win_len - 1;

    auto push = [&queue, &kmer_stream] {
        queue.push(kmer_stream.getCurrKMer());
        queue.push(kmer_stream.getCurrCompKMer());
    };

    auto shift_push = [&push, &kmer_stream] {
        kmer_stream.shiftKMer();
        push();
    };

    auto min_and_pop = [&queue, &minimizers] {
        minimizers.push_back(queue.min());
        // double pop because
        // we store original and complement
        queue.pop();
        queue.pop();
    };

    // clang-format off
    push();
    while (kmer_stream.hasShift()) {
        if (kmer_stream.pos() == win_end)
            ++win_end, min_and_pop();
        shift_push();
    }

    while (win_end < end && !queue.empyt())
        ++win_end, min_and_pop();
    // clang-format on

    return minimizers;
}

/**
 * @brief Finds internal minimizers of a given sequene.
 *
 * @param seq target sequence
 * @param k size of a K-mer
 * @param win_len sliding window lenght
 * @return KMers set of minimizers
 */
KMers findMinimizers(std::string_view seq, std::uint32_t k,
                     std::uint32_t win_len) {
    return findMinimizers(seq, k, win_len, 0, seq.length());
}

KMers minimizers(char const* sequence, std::uint32_t sequence_length,
                 std::uint32_t k, std::uint32_t window_length) {
    if (k < 0 || k > 16)
        throw std::invalid_argument(
            "Unsupported Kmer size! "
            "\'k\' has to be > 0 and <= 16");

    if (window_length > sequence_length)
        throw std::invalid_argument("Window size is too large!");

    auto seq = std::string_view(sequence, sequence_length);
    auto minimizers = findMinimizers(seq, k, window_length);

    auto unique_sorted = [](auto& vec) {
        auto set = std::set<KMer>{};
        for (auto& it: vec)
            set.insert(it);
        vec.assign(set.begin(), set.end());

        return vec;
    };

    auto extend_minimizers = [&minimizers](auto& vec) {
      auto nw_sz = minimizers.size() + vec.size();

      minimizers.reserve(nw_sz);
      minimizers.insert(minimizers.end(), vec.begin(), vec.end());
    };

    for (auto u = std::uint32_t{1}; u < window_length; ++u) {
        auto minims_front = findMinimizers(seq, k, u, 0, u);
        extend_minimizers(minims_front);
    }

    if (k < window_length) {
        for (auto u = k; u < window_length ; ++u) {
            auto minims_back = findMinimizers(seq, k, u, sequence_length - u, sequence_length);
            extend_minimizers(minims_back);
        }
    }

    return unique_sorted(minimizers);
}

}  // namespace minimizers
}  // namespace orange