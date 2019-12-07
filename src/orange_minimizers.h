#ifndef ORANGE_MINIMIZERS_H_
#define ORANGE_MINIMIZERS_H_

#include <vector>
#include <tuple>

namespace orange {
namespace minimizers {

/**
 * @brief std::tuple describing a K-Mer.
 *
 * @detalis <K-Mer bit mask, position in sequence, original or complement>
 */
using KMer = std::tuple<std::uint32_t, std::uint32_t, bool>;

/**
 * @brief std::vector of @ref orange::minimizers::KMer
 *      representing a sequence of K-Mers
 */
using KMers = std::vector<KMer>;

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
KMers minimizers(char const* sequence,
                                  std::uint32_t sequence_length, std::uint32_t k,
                                  std::uint32_t window_length);

}  // namespace minimizers
}  // namespace orange

#endif /* ORANGE_MINIMIZERS_H_ */