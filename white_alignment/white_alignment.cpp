#include "white_alignment.hpp"
#include <algorithm>
#include <initializer_list>
#include <iostream>

namespace white {

    enum Direction {
        kUp,
        kLeft,
        kDiagonal,
        kNone
    };

    struct Cell {
        int value;
        Direction direction;

        bool operator<(const Cell &other) const {
            return value < other.value;
        };

        bool operator>(const Cell &other) const {
            return value > other.value;
        };

        bool operator==(const Cell &other) const {
            return value == other.value;
        };
    };

    void matrixInit(int gap, int match, int mismatch,
                    const char *query, const char *target,
                    Cell **matrix, int rowCnt, int columCnt,
                    AlignmentType alignmentType) {
        int rowIndex, columnIndex;

        matrix[0][0] = {0, Direction::kNone};

        if (alignmentType == AlignmentType::kGlobal) {
            for (rowIndex = 0; rowIndex < rowCnt; rowIndex++) {
                matrix[rowIndex][0] = {gap * rowIndex, Direction::kUp};
            }

            for (columnIndex = 0; columnIndex < columCnt; columnIndex++) {
                matrix[rowIndex][0] = {gap * columnIndex, Direction::kLeft};
            }
        } else { // both local and semi-global
            for (rowIndex = 0; rowIndex < rowCnt; rowIndex++) {
                matrix[rowIndex][0] = {0, Direction::kNone};
            }

            for (columnIndex = 0; columnIndex < columCnt; columnIndex++) {
                matrix[0][columnIndex] = {0, Direction::kNone};
            }
        }
        for (columnIndex = 1; columnIndex < columCnt; columnIndex++) {
            int replacement = matrix[rowIndex - 1][columnIndex - 1].value +
                              query[rowIndex] == target[columnIndex] ? match : mismatch;
            int insertion = matrix[rowIndex][columnIndex - 1].value + gap;
            int deletion = matrix[rowIndex - 1][columnIndex].value + gap;
            Cell maxValue = std::max({
                                             Cell{replacement, Direction::kDiagonal},
                                             Cell{insertion, Direction::kLeft},
                                             Cell{deletion, Direction::kUp}
                                     });
            matrix[rowIndex][columnIndex] = maxValue;
            if (alignmentType == AlignmentType::kLocal) {
                if (maxValue.value < 0) {
                    matrix[rowIndex][columnIndex].value = 0;
                }
            }
        }
    };

    int PairwiseAlignment(const char *query, unsigned int query_length,
                          const char *target, unsigned int target_length,
                          AlignmentType type,
                          int match,
                          int mismatch,
                          int gap) {

        int rowCount = query_length + 1;
        int colCount = target_length + 1;
        Cell **matrix = new Cell * [rowCount];
        for (int i = 0; i < rowCount; ++i)
            matrix[i] = new Cell[colCount];
        matrixInit(gap, match, mismatch, query, target, matrix, rowCount, colCount, type);

        //    osl memorije
        //    for(int i = 0; i < sizeY; ++i) {
        //    delete [] ary[i];
        //    }
        //    delete [] ary;

    }
}


