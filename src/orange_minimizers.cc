#include "orange_minimizers.h"

#include <string_view>
#include <functional>
#include <queue>
#include <deque>

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
        raw_queue_.push(std::move(kmer));
    }

    /**
     * @brief Removes an element from the front of the queue
     */
    void pop() {
        if (raw_queue_.front() == min_queue_.front().get())
            min_queue_.pop_front();
        raw_queue_.pop();
    }

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
    KMerStream(std::string_view seq, std::uint32_t k)
        : seq_{seq}, k_{k}, pos_{0}, kmer_{0}, comp_kmer_{0} {
        auto shift = int{2} * k - 2;  // Each base holds up two bits
        for (; pos_ < k_; ++pos_) {
            kmer_ |= getMask(seq_[pos_]) << shift;
            comp_kmer_ |= getMask(getComplement(seq_[pos_])) << shift;

            shift -= 2;
        }

        // First 2 * k bits are set to 1
        mask_ = (1 << (2 * k)) - 1;
    }

    void shiftKMer() {
        if (pos_ >= seq_.size())
            return;

        // Shift and ignore bits
        kmer_ <<= 2, comp_kmer_ <<= 2;
        kmer_ &= mask_, comp_kmer_ &= mask_;

        kmer_ |= getMask(seq_[pos_]);
        comp_kmer_ |= getMask(getComplement(seq_[pos_]));

        ++pos_;
    }

    std::uint32_t pos() const { return pos_; }

    bool hasKMer() const { return pos_ - 1 < seq_.size(); }

    KMer getCurrKMer() const { return KMer{kmer_, pos_ - k_, false}; }

    KMer getCurrCompKMer() const { return KMer{comp_kmer_, pos_ - k_, true}; }

private:
    std::string_view seq_; //<<< Sequence

    std::uint32_t k_; //<<< KMer length

    std::uint32_t pos_; //<<< Refers to the current symbol in a sequence
    std::uint32_t mask_;  //<<< Used for ignoring leftover bits in a shift

    KMerVal kmer_; //<<< KMer on the original strand
    KMerVal comp_kmer_; //<<< KMer on the complement strand
};

/**
 * @brief Finds internal minimizers of a given sequene.
 *
 * @param seq target sequence
 * @param k size of a K-mer
 * @param win_len sliding window lenght
 * @return KMers set of minimizers
 */
KMers internalMinimizers(std::string_view seq, std::uint32_t k,
                         std::uint32_t win_len) {
    auto kmer_stream = KMerStream{seq, k};
    auto queue = MinKQueue{};

    auto minimizers = KMers{};

    auto push_shift = [&queue, &kmer_stream] {
        queue.push(kmer_stream.getCurrKMer());
        queue.push(kmer_stream.getCurrCompKMer());

        kmer_stream.shiftKMer();
    };

    auto min_and_pop = [&queue, &minimizers] {
        minimizers.push_back(queue.min());
        queue.pop();
    };

    while (kmer_stream.pos() < win_len) push_shift();
    min_and_pop();

    while (kmer_stream.hasKMer()) push_shift(), min_and_pop();

    return minimizers;
}

KMers minimizers(char const* sequence, unsigned int sequence_length,
                 unsigned int k, unsigned int window_length) {
    if (k < 0 || k > 16)
        throw std::invalid_argument(
            "Unsupported Kmer size! "
            "\'k\' has to be > 0 and <= 16");

    if (window_length > sequence_length)
        throw std::invalid_argument("Window size is too large!");

    auto minimizers = internalMinimizers(
        std::string_view(sequence, sequence_length), k, window_length);

    return minimizers;
}

}  // namespace minimizers
}  // namespace orange