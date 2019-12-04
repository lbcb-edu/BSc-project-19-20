#include <string>
namespace white {

    enum AlignmentType {
        kGlobal,
        kLocal,
        kSemiGlobal};

    int PairwiseAlignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match,
                           int mismatch,
                           int gap);

    int PairwiseAlignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match,
                           int mismatch,
                           int gap,
                           std::string& cigar,
                           unsigned int& target_begin);
}

