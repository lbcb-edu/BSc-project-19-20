#include <iostream>
#include <vector>
#include "brown_alignment.hpp"

using namespace std;

namespace brown {

    struct matrixCell {
        int cost;
        char parent;
    };

    int global_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       int match,
                       int mismatch,
                       int gap) {

        vector<vector<int>> matrix(query_length + 1, vector<int>(target_length + 1, 0));

    }

    int local_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       int match,
                       int mismatch,
                       int gap) {

        vector<vector<int>> matrix(query_length + 1, vector<int>(target_length + 1, 0));

    }

    int semi_global_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       int match,
                       int mismatch,
                       int gap) {

        vector<vector<int>> matrix(query_length + 1, vector<int>(target_length + 1, 0));
    }    

    int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap) {

        switch (type)
        {
        case AlignmentType::global:
            return global_alignment(query, query_length, target, target_length, match, mismatch, gap);
        case AlignmentType::local:
            return local_alignment(query, query_length, target, target_length, match, mismatch, gap);
        case  AlignmentType::semi_global:
            return semi_global_alignment(query, query_length, target, target_length, match, mismatch, gap);
        default:
            fprintf(stderr, "Unknown alignment type\n");
            break;
        }      
    }
}