#include "orange_minimizers.h"

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
    void push(Kmer const kmer) {
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
    Kmer min() const {
        return min_queue_.front();
    }

private:
    std::queue<Kmer> raw_queue_;
    std::deque<Kmer> min_queue_;
};


std::vector<Kmer> minimizers(const char* sequence,
                                  unsigned int sequence_length, unsigned int k,
                                  unsigned int window_length);

}  // namespace minimizers
}  // namespace orange