#include <iostream>
#include "brown_alignment.hpp"

namespace brown {
    
    int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap);    
}