#include <iostream>
#include "pink_alignment.h"

namespace pink {

    enum ParentEnum { no, up, left, diagonal };

    typedef struct {
        int value;
        ParentEnum parent;
    } cell;

    void init_matrix(cell** matrix, unsigned int m, unsigned int n, AlignmentType type, int gap){
        matrix[0][0].value = 0;
        matrix[0][0].parent = no;
        if(type == global) {
            for (unsigned int i = 0; i < m; i++) {
                matrix[i][0].value = i * gap;
                matrix[i][0].parent = up;
            }
            for (unsigned int j = 0; j < n; j++) {
                matrix[0][j].value = j * gap;
                matrix[0][j].parent = left;
            }
        } else {
            for (unsigned int i = 0; i < m; i++) {
                matrix[i][0].value = 0;
                matrix[i][0].parent = no;
            }
            for (unsigned int j = 0; j < n; j++) {
                matrix[0][j].value = 0;
                matrix[0][j].parent = no;
            }
        }
    }

    int max(int a, int b, int c){
        if(a >= b && a >= c)
            return a;
        if(b >= a && b >= c)
            return b;
        if(c >= a && c >= b)
            return c;
    }

    void make_cigar(std::string& cigar, unsigned int& target_begin, const char* query, const char* target, cell** matrix, int no_of_rows, int no_of_columns, int index_i, int index_j) {
        cigar = "";
        std::string str = "";
        int m = 0;
        int mm = 0;
        int d = 0;
        int i = 0;
        if(index_i < no_of_rows - 1){
            str = std::to_string(no_of_rows - 1 - index_i) + 'S';
            cigar = str.append(cigar);
        }
        while (matrix[index_i][index_j].parent != no) {
            switch (matrix[index_i][index_j].parent) {
                case up :
                    if (mm != 0) {
                        str = std::to_string(mm) + 'X';
                        cigar = str.append(cigar);
                        mm = 0;
                    } else if (m != 0) {
                        str = std::to_string(m) + '=';
                        cigar = str.append(cigar);
                        m = 0;
                    } else if (i != 0) {
                        str = std::to_string(i) + 'I';
                        cigar = str.append(cigar);
                        i = 0;
                    }
                    d++;
                    index_i--;
                    break;
                case left :
                    if (mm != 0) {
                        str = std::to_string(mm) + 'X';
                        cigar = str.append(cigar);
                        mm = 0;
                    } else if (m != 0) {
                        str = std::to_string(m) + '=';
                        cigar = str.append(cigar);
                        m = 0;
                    } else if (i != 0) {
                        str = std::to_string(i) + 'I';
                        cigar = str.append(cigar);
                        i = 0;
                    }
                    i++;
                    index_j--;

                    break;
                case diagonal :
                    if (d != 0) {
                        str = std::to_string(d) + 'D';
                        cigar = str.append(cigar);
                        d = 0;
                    } else if (i != 0) {
                        str = std::to_string(i) + 'I';
                        cigar = str.append(cigar);
                        i = 0;
                    }

                    if(query[index_i-1] == target[index_j-1]){
                        if(mm != 0) {
                            str = std::to_string(mm) + 'X';
                            cigar = str.append(cigar);
                            mm = 0;
                        }
                        m++;
                    } else {
                        if (m != 0) {
                            str = std::to_string(m) + '=';
                            cigar = str.append(cigar);
                            m = 0;
                        }
                        mm++;
                    }

                    index_i--;
                    index_j--;
                    break;
                default:
                    index_i--;
                    index_j--;
                    break;
            }
        }
        if(index_i != 0){
            str = std::to_string(index_i) + 'S';
            cigar = str.append(cigar);
        }
        target_begin = index_j;


    }





    int pairwise_alignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match,
                           int mismatch,
                           int gap,
                           std::string& cigar,
                           unsigned int& target_begin) {

        unsigned int no_of_rows = query_length + 1;
        unsigned int no_of_columns = target_length + 1;

        cell** matrix = new cell * [no_of_rows];
        for(unsigned int i = 0; i < no_of_rows; i++)
            matrix[i] = new cell[no_of_columns];

        if(type == global)
            init_matrix(matrix, no_of_rows, no_of_columns, type, gap);
        else // using suffix-prefix alignment algorithm for semi-global alignment
            init_matrix(matrix, no_of_rows, no_of_columns, type,0);


        for(unsigned int i = 1; i < no_of_rows; i++){
            for(unsigned int j = 1; j < no_of_columns; j++){
                int insertion = matrix[i][j-1].value + gap;
                int deletion = matrix[i-1][j].value + gap;
                int w = query[i] == target[j] ? match : mismatch;
                int mmatch = matrix[i-1][j-1].value + w;

                int val = max(insertion, deletion, mmatch);
                if(type == local)
                    if(val < 0)
                        val = 0;

                matrix[i][j].value = val;

                if(val == mmatch)
                    matrix[i][j].parent = diagonal;
                else if(val == deletion)
                    matrix[i][j].parent = up;
                else if(val == insertion)
                    matrix[i][j].parent = left;
                else
                    matrix[i][j].parent = no;
            }
        }

        int alignment_score;
        int index_i, index_j;
        if(type == global) {
            alignment_score = matrix[no_of_rows - 1][no_of_columns - 1].value;
            index_i = no_of_rows - 1;
            index_j = no_of_columns - 1;
        } else if(type == semi_global) {
            alignment_score = matrix[0][0].value;
            for(int i = 1; i < no_of_rows; i++){
                if(matrix[i][no_of_columns - 1].value > alignment_score) {
                    alignment_score = matrix[i][no_of_columns - 1].value;
                    index_i = i;
                    index_j = no_of_columns - 1;
                }
            }
            for(int j = 1; j < no_of_columns; j++){
                if(matrix[no_of_rows - 1][j].value > alignment_score) {
                    alignment_score = matrix[no_of_rows - 1][j].value;
                    index_i = no_of_rows - 1;
                    index_j = j;
                }
            }
        } else {
            alignment_score = matrix[0][0].value;
            for(int i = 1; i < no_of_rows; i++){
                for(int j = 1; j < no_of_columns; j++){
                    if(matrix[i][j].value > alignment_score) {
                        alignment_score = matrix[i][j].value;
                        index_i = i;
                        index_j = j;
                    }
                }
            }
        }

        make_cigar(cigar, target_begin, query, target, matrix, no_of_rows, no_of_columns, index_i, index_j);

        for (int i = 0; i < no_of_rows; i++) {
            delete[] matrix[i];
        }
        delete[] matrix;

        return alignment_score;
    }


    int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap) {
        std::string str = "";
        unsigned a;
        return pairwise_alignment(query, query_length, target, target_length, type, match, mismatch, gap, str, a);
    }
}