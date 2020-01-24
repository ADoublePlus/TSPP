#include "Matrix.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>

int RoundUptoPowerOf2(int num)
{
    num--;

    num |= num >> 1;
    num |= num >> 2;
    num |= num >> 4;
    num |= num >> 8;
    num |= num >> 16;

    num++;

    return num;
}

void InitializeMatrix(Matrix* mat, int num_rows, int num_cols)
{
    int rcol_align = RoundUptoPowerOf2(num_cols);
    int size = num_rows * rcol_align;

    mat->_rCol_align = rcol_align;
    mat->_size = size;
    mat->data = (float*)_aligned_malloc(size * sizeof(float), 32);
    mat->num_cols = num_cols;
    mat->num_rows = num_rows;
}

void PrintMatrix(const Matrix& mat)
{
    for (int row = 0; row < mat.num_rows; row++)
    {
        for (int col = 0; col < mat.num_cols; col++)
        {
            printf("%f ", mat.get(row, col));
        }

        printf("\n");
    }
}