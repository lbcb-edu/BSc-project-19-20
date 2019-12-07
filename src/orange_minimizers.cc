#include "orange_minimizers.h"

#include <string_view>
#include <queue>
#include <deque>

namespace orange {
namespace minimizers {

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
               std::get<0>(min_queue_.back()) > std::get<0>(kmer))
            min_queue_.pop_back();
        min_queue_.push_back(std::move(kmer));
    }

    /**
     * @brief Removes an element from the front of the queue
     */
    void pop() {
        if (raw_queue_.front() == min_queue_.front())
            min_queue_.pop_front();
        raw_queue_.pop();
    }

    /**
     * @brief Returnt the Kmer with the smallest value,
     *      aka. minimizer of the sliding window
     *
     * @return Kmer
     */
    KMer min() const { return min_queue_.front(); }

private:
    std::queue<KMer> raw_queue_;  //<<< Contains all the kmers
    std::deque<KMer> min_queue_;  //<<< Contains smallest kmers
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
 * @brief Returns a bitmask representation of 
 *      a @ref orange::minimizers::KMer
 * 
 * @param kmer KMer to compress 
 * @return std::unit32_t bitmask representation of KMer
 */
std::uint32_t compressKMerVal(std::string_view kmer) {
    auto compressed = std::uint32_t{};
    auto k = kmer.size(); 

    // Every base requires two bits
    // That's the reason why we double KMer size
    auto shift = 2 * k;
    for (auto base: kmer) {
        compressed |= getMask(base) << shift;
        shift -= 2;
    }

    return compressed;
}

/**
 * @brief Parameters for sliding window algorithm
 */
struct SlideConf {
    char const* seq_;
    std::uint32_t const seq_len_;

    std::uint32_t const start_;
    bool const anchored_;
    
    std::uint32_t const k_;
    std::uint32_t const win_len_;
};

KMers slidingWindow(SlideConf const& conf) {

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