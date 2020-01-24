#pragma once

#include <cstdlib>
#include <cassert>

struct Matrix
{
    float* data;
    int num_rows;
    int num_cols;
    int _size;
    int _rCol_align;

    float* operator[](int row) const
    {
        assert(row < num_rows);
        return &data[row * _rCol_align];
    }

    float& get(int row, int col)
    {
        assert(col < num_cols);
        return (&data[row * _rCol_align])[col];
    }

    float get(int row, int col) const
    {
        assert(col < num_cols);
        return (&data[row * _rCol_align])[col];
    }
};

void InitializeMatrix(Matrix* mat, int num_rows, int num_cols);
void PrintMatrix(const Matrix& mat);