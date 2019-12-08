#include "orange_minimizers.h"

#include <string_view>
#include <functional>
#include <queue>
#include <deque>

namespace orange {
namespace minimizers {

using Window = std::uint32_t;

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
    std::deque<std::reference_wrapper<KMer>> min_queue_;  //<<< Contains smallest kmers
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
 * @brief Streams consecutive windows over a sequence
 */
class WindowStream {
public:
    WindowStream(std::string_view seq, std::uint32_t k, std::uint32_t win_sz)
        : seq_{seq}, k_{k}, win_sz_{win_sz}, pos_{0}, win_{0}, comp_win_{0} {
        auto shift = int{2} * k - 2;  // Each base holds up two bits
        for (; pos_ < win_sz_; ++pos_) {
            win_ |= getMask(seq_[pos_]) << shift;
            comp_win_ |= getMask(getComplement(seq_[pos_])) << shift;

            shift -= 2;
        }
    }

    void shiftWins() {
        if (pos_ >= win_sz_)
            return;

        win_ <<= 2, comp_win_ <<= 2;

        win_ |= getMask(seq_[pos_]);
        comp_win_ |= getMask(getComplement(seq_[pos_]));

        ++pos_;
    }

    bool hasWindow() const { return pos_ - 1 < win_sz_; }

    KMer getCurrWin() const { return KMer{win_, pos_ - k_, false}; }

    KMer getCurrCompWin() const { return KMer{comp_win_, pos_ - k_, true}; }

private:
    std::string_view seq_;

    std::uint32_t k_;
    std::uint32_t win_sz_;

    std::uint32_t pos_;

    Window win_;
    Window comp_win_;
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
    auto win_stream = WindowStream{seq, k, win_len};
    auto queue = MinKQueue{};

    auto minimizers = KMers{};

    while (win_stream.hasWindow()) {
        queue.push(win_stream.getCurrWin());
        queue.push(win_stream.getCurrCompWin());

        minimizers.push_back(queue.min());
    }

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
}

}  // namespace minimizers
}  // namespace orange