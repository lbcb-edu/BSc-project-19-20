#include "white_alignment.hpp"
#include <algorithm>

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
        int rowInd = 0, columnInd = 0;

        matrix[0][0] = {0, Direction::kNone};

        if (alignmentType == AlignmentType::kGlobal) {
            for (rowInd = 1; rowInd < rowCnt; rowInd++) {
                matrix[rowInd][0] = {gap * rowInd, Direction::kUp};
            }

            for (columnInd = 1; columnInd < columCnt; columnInd++) {
                matrix[0][columnInd] = {gap * columnInd, Direction::kLeft};
            }
        } else { // both local and semi-global
            for (rowInd = 1; rowInd < rowCnt; rowInd++) {
                matrix[rowInd][0] = {0, Direction::kNone};
            }

            for (columnInd = 1; columnInd < columCnt; columnInd++) {
                matrix[0][columnInd] = {0, Direction::kNone};
            }
        }
        for (rowInd = 1; rowInd < rowCnt; rowInd++) {
            for (columnInd = 1; columnInd < columCnt; columnInd++) {
                int w = query[rowInd] == target[columnInd] ? match : mismatch;
                int replacement = matrix[rowInd - 1][columnInd - 1].value + w;

                int insertion = matrix[rowInd][columnInd - 1].value + gap;
                int deletion = matrix[rowInd - 1][columnInd].value + gap;
                Cell maxValue = std::max({
                                                 Cell{replacement, Direction::kDiagonal},
                                                 Cell{insertion, Direction::kLeft},
                                                 Cell{deletion, Direction::kUp}
                                         });
                matrix[rowInd][columnInd] = maxValue;
                if (alignmentType == AlignmentType::kLocal) {
                    if (maxValue.value < 0) {
                        matrix[rowInd][columnInd].value = 0;
                    }
                }
            }
        }

    };

    int AlignmentEnd (int &rowIndex, int &colIndex, Cell** matrix,
            int rowCount, int colCount, AlignmentType alignmentType) {
        int maxScore = 0;
        int row = 0, col = 0;
        if (alignmentType == AlignmentType::kGlobal) {
            row = rowCount - 1;
            col = colCount - 1;
            maxScore = matrix[row][col].value;
        } else if (alignmentType == AlignmentType::kSemiGlobal) {
            maxScore = matrix[1][colCount - 1].value;
            row = 0, col = colCount - 1;
            for (int i = 2; i < rowCount; i++) {
                if (matrix[i][colCount - 1].value > maxScore) {
                    maxScore = matrix[i][colCount - 1].value;
                    row = i;
                }
            }
            for (int j = 1; j < colCount; j++) {
                if (matrix[rowCount - 1][j].value > maxScore) {
                    maxScore = matrix[rowCount - 1][j].value;
                    row = rowCount - 1;
                    col = j;
                }
            }
        } else { //local
            row = 0, col = 0;
            for (int i = 1; i < rowCount; i++) {
                for (int j = 1; j < colCount; j++) {
                    if (matrix[i][j].value > maxScore) {
                        maxScore = matrix[i][j].value;
                        row = i;
                        col = j;
                    }
                }
            }
        }
        rowIndex = row;
        colIndex = col;
        return maxScore;
    }

    void next (int &rowIndex, int &colIndex, Direction dir) {
        if (dir == Direction::kLeft) {
            colIndex--;
        } else if (dir == Direction::kUp) {
            rowIndex--;
        } else if (dir == Direction::kDiagonal) {
            rowIndex--;
            colIndex--;
        } else {
            //........
        }
        return;
    }

    std::string cigarSegment (int noOfReps, Direction dir) {
        std::string Op = "", noRep = std::to_string(noOfReps);

        //to do S, =, X
        switch (dir) {
            case Direction::kDiagonal:
                Op.append("M");
                return Op + noRep;
            case Direction::kLeft:
                Op.append("D");
                return Op + noRep;
            case Direction::kUp:
                Op.append("I");
                return Op + noRep;
        }
        return "";

    }

    void cigarCreate (const char *query, const char *target,
            std::string& cigar, unsigned int& target_begin,
            Cell** matrix, int rowCount, int colCount,
            int rowIndex, int colIndex) {
        //S, =, X
        std::string seq = "";
        int currentRow = rowIndex, currentCol = colIndex, noOfRep = 1;
        int nextRow = currentRow, nextCol = currentCol;
        while (matrix[currentRow][currentCol].direction != Direction::kNone) {
            next(nextRow, nextCol, matrix[currentRow][currentCol].direction);
            if (matrix[currentRow][currentCol].direction != matrix[nextRow][nextCol].direction) {
                seq.append(cigarSegment(noOfRep, matrix[currentRow][currentCol].direction));
                noOfRep = 1;
            } else noOfRep++;


            currentRow = nextRow;
            currentCol = nextCol;
        }
        cigar = seq;

        target_begin = currentCol;



    }

    int PairwiseAlignment(const char *query, unsigned int query_length,
                          const char *target, unsigned int target_length,
                          AlignmentType type,
                          int match,
                          int mismatch,
                          int gap) {
        std::string Cigar = "";
        unsigned int t;
        int rowCount = query_length + 1;
        int colCount = target_length + 1;
        Cell **matrix = new Cell * [rowCount];
        for (int i = 0; i < rowCount; ++i)
            matrix[i] = new Cell[colCount];
        matrixInit(gap, match, mismatch, query, target, matrix, rowCount, colCount, type);
        int alignmentRow, alignmentCol;
        int optimalAlignment = AlignmentEnd(alignmentRow, alignmentCol, matrix, rowCount, colCount, type);
        cigarCreate(query, target, Cigar, t, matrix, rowCount, colCount, alignmentRow, alignmentCol);

        for (int i = 0; i < rowCount; i++) {
            delete[] matrix[i];
        }
        delete[] matrix;

        return optimalAlignment;

    }

    int PairwiseAlignment(const char* query, unsigned int query_length,
                          const char* target, unsigned int target_length,
                          AlignmentType type,
                          int match,
                          int mismatch,
                          int gap,
                          std::string& cigar,
                          unsigned int& target_begin) {

        int rowCount = query_length + 1;
        int colCount = target_length + 1;
        Cell **matrix = new Cell * [rowCount];
        for (int i = 0; i < rowCount; ++i)
            matrix[i] = new Cell[colCount];
        matrixInit(gap, match, mismatch, query, target, matrix, rowCount, colCount, type);
        int alignmentRow, alignmentCol;
        int optimalAlignment = AlignmentEnd(alignmentRow, alignmentCol, matrix, rowCount, colCount, type);
        cigarCreate(query, target, cigar, target_begin, matrix, rowCount, colCount, alignmentRow, alignmentCol);

        for (int i = 0; i < rowCount; i++) {
            delete[] matrix[i];
        }
        delete[] matrix;

        return optimalAlignment;
    }
}

