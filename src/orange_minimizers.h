#ifndef ORANGE_MINIMIZERS_H_
#define ORANGE_MINIMIZERS_H_

#include <vector>
#include <tuple>

namespace orange {
namespace minimizers {

/**
 * @brief std::tuple describing a minimizer.
 *
 * @detalis <minimizer bit mask, position in sequence, original or complement>
 */
using Kmer = std::tuple<std::uint32_t, std::uint32_t, bool>;

/**
 * @brief Returns a std::vector of @ref orange::minimizers::minimizer
 *      for a given sequence
 *
 * @details <a
 * href="https://academic.oup.com/bioinformatics/article/20/18/3363/202143"> See
 * </a>
 *
 * @param sequence target sequence
 * @param sequence_length target sequence length
 * @param k lenght of a minimizer
 * @param window_length window length
 *
 * @return std::vector<minimizer>
 */
std::vector<Kmer> minimizers(const char* sequence,
                                  unsigned int sequence_length, unsigned int k,
                                  unsigned int window_length);

}  // namespace minimizers
}  // namespace orange

#endif /* ORANGE_MINIMIZERS_H_ */