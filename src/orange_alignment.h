#ifndef ORANGE_ALIGNMENT_H_
#define ORANGE_ALIGNMENT_H_

#include <string>

namespace orange {
namespace alignment {

/**
 * @brief Determines type of sequence alignment done by
 *          @ref orange::alignment::pairwise_alignment
 */
enum AlignmentType {
    kGlobal,     ///< Global alignment using Needleman-Wunsch
    kLocal,      ///< Local alignment using Smith-Waterman
    kSemiGlobal  ///< Semi-global alignment
};

/**
 * @brief Holds alignment type and scoring scheme 
 */
struct AlignConf {
    AlignmentType type_;
    int match_;
    int mismatch_;
    int gap_;
};

/**
 * @brief Calculates alignment score between query and target sequences
 *
 * @param query qurey seuqence
 * @param query_length query sequence length
 * @param target target sequence
 * @param target_length target sequence length
 *
 * @param type @ref orange::alignment::AlignmentType
 *
 * @param match match transition cost used in an algirthm of choice
 * @param mismatch mismatch transition cost used in an algirthm of choice
 * @param gap gap transition cost used in an algirthm of choice
 *
 * @return int alignment score
 */
int pairwiseAlignment(char const* query, unsigned int query_length,
                       char const* target, unsigned int target_length,
                       AlignmentType type, int match, int mismatch, int gap);

/**
 * @brief
 *
 * @param query qurey seuqence
 * @param query_length query sequence length
 * @param target target sequence
 * @param target_length target sequence length
 *
 * @param type @ref orange::alignment::AlignmentType
 *
 * @param match match transition cost used in an algirthm of choice
 * @param mismatch mismatch transition cost used in an algirthm of choice
 * @param gap gap transition cost used in an algirthm of choic
 *
 * @param cigar return through reference, alignment cigar std::string
 * @param int retrun through reference, alignment starting postion in target
 *              sequence
 *
 * @return int alignment score
 */
int pairwiseAlignment(char const* query, unsigned int query_length,
                       char const* target, unsigned int target_length,
                       AlignmentType type, int match, int mismatch, int gap,
                       std::string& cigar, unsigned int& target_begin);

}  // namespace alignment
}  // namespace orange

#endif /* ORANGE_ALIGNMENT_H_ */