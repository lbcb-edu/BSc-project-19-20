#include "orange_alignment.h"

namespace orange {
namespace alignment {

int pairwise_alignment(char const* query, unsigned int query_length,
                       char const* target, unsigned int target_length,
                       AlignmentType type, int match, int mismatch, int gap) {
    return 0;
}

int pairwise_alignment(char const* query, unsigned int query_length,
                       char const* target, unsigned int target_length,
                       AlignmentType type, int match, int mismatch, int gap,
                       std::string& cigar, unsigned int& target_begin) {
    return 0;
}

}  // namespace alignment
}  // namespace orange