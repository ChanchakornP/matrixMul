#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>
#include <pthread.h>
#include <sys/resource.h>
#include <math.h>
#define UnrollingFactor 4

typedef struct
{
    double **data;
    int rows;
    int columns;
} matrix;

typedef struct
{
    matrix A;
    matrix B;
    matrix Result;
    int start_row;
    int end_row;
    int start_col;
    int end_col;
} t_arg;

int NUM_THREADS, M_BLOCK, N_BLOCK; // Assume NUM_THREAD >= M_BLOCK * N_BLOCK

void matrixMul(matrix A, matrix B, matrix Result)
{
    double **A_data = A.data;
    double **B_data = B.data;
    double **Result_data = Result.data;
    int row_A = A.rows, col_A = A.columns, col_B = B.columns;
    for (int i = 0; i < row_A; i++)
    {
        for (int j = 0; j < col_A; j++)
        {
            for (int k = 0; k < col_B; k++)
            {
                Result_data[i][k] += A_data[i][j] * B_data[j][k];
            }
        }
    }
}

void printMatrix(matrix A, int row, int col)
{
    double **A_data = ((double **)A.data);
    for (int i = 0; i < A.rows && i < row; i++)
    {
        for (int j = 0; j < A.columns && j < col; j++)
            printf("%.2f  ", A_data[i][j]);
        printf("\n");
    }
    return;
}

void matrixMulLoopUnrolling(matrix A, matrix B, matrix Result)
{
    double **A_data = A.data;
    double **B_data = B.data;
    double **Result_data = Result.data;
    double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
    int row_A = A.rows, col_A = A.columns, col_B = B.columns, i, j, k;
    int row_A0 = A.rows / 4 * 4;
    int col_A0 = A.columns / 4 * 4;
    int col_B0 = B.columns / 4 * 4;

    for (i = 0; i < row_A0; i += UnrollingFactor)
    {
        for (j = 0; j < col_A0; j += UnrollingFactor)
        {
            a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
            a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
            a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
            a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

            for (k = 0; k < col_B0; k += UnrollingFactor)
            {
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                Result_data[i][k + 1] += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
                Result_data[i][k + 2] += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
                Result_data[i][k + 3] += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;

                Result_data[i + 1][k] += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                Result_data[i + 1][k + 1] += a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
                Result_data[i + 1][k + 2] += a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
                Result_data[i + 1][k + 3] += a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;

                Result_data[i + 2][k] += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                Result_data[i + 2][k + 1] += a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
                Result_data[i + 2][k + 2] += a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
                Result_data[i + 2][k + 3] += a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;

                Result_data[i + 3][k] += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
                Result_data[i + 3][k + 1] += a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
                Result_data[i + 3][k + 2] += a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
                Result_data[i + 3][k + 3] += a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;
            }
            // If column B is indivisible by unrolling factor. This do the rest.
            for (k = col_B0; k < col_B; k++)
            {
                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];
                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                Result_data[i + 1][k] += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                Result_data[i + 2][k] += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                Result_data[i + 3][k] += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
            }
        }
        // If column A is indivisible by unrolling factor. This do the rest, while keep k unrolling.
        for (j = col_A0; j < col_A; j++)
        {
            a00 = A_data[i][j];
            a10 = A_data[i + 1][j];
            a20 = A_data[i + 2][j];
            a30 = A_data[i + 3][j];
            for (k = 0; k < col_B0; k += UnrollingFactor)
            {
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                Result_data[i][k] += a00 * b00;
                Result_data[i][k + 1] += a00 * b01;
                Result_data[i][k + 2] += a00 * b02;
                Result_data[i][k + 3] += a00 * b03;

                Result_data[i + 1][k] += a10 * b00;
                Result_data[i + 1][k + 1] += a10 * b01;
                Result_data[i + 1][k + 2] += a10 * b02;
                Result_data[i + 1][k + 3] += a10 * b03;

                Result_data[i + 2][k] += a20 * b00;
                Result_data[i + 2][k + 1] += a20 * b01;
                Result_data[i + 2][k + 2] += a20 * b02;
                Result_data[i + 2][k + 3] += a20 * b03;

                Result_data[i + 3][k] += a30 * b00;
                Result_data[i + 3][k + 1] += a30 * b01;
                Result_data[i + 3][k + 2] += a30 * b02;
                Result_data[i + 3][k + 3] += a30 * b03;
            }
            for (k = col_B0; k < col_B; k++)
            {
                b00 = B_data[j][k];
                Result_data[i][k] += a00 * b00;
                Result_data[i + 1][k] += a10 * b00;
                Result_data[i + 2][k] += a20 * b00;
                Result_data[i + 3][k] += a30 * b00;
            }
        }
    }
    // If row A is indivisible by unrolling factor. This do the rest, while keep j, k unrolling.
    for (i = row_A0; i < row_A; i++)
    {
        for (j = 0; j < col_A0; j += UnrollingFactor)
        {
            a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
            for (k = 0; k < col_B0; k += UnrollingFactor)
            {
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                Result_data[i][k + 1] += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
                Result_data[i][k + 2] += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
                Result_data[i][k + 3] += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
            }
            for (k = col_B0; k < col_B; k++)
            {
                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];
                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
            }
        }
        for (j = col_A0; j < col_A; j++)
        {
            a00 = A_data[i][j];
            for (k = 0; k < col_B0; k += UnrollingFactor)
            {
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                Result_data[i][k] += a00 * b00;
                Result_data[i][k + 1] += a00 * b01;
                Result_data[i][k + 2] += a00 * b02;
                Result_data[i][k + 3] += a00 * b03;
            }
            for (k = col_B0; k < col_B; k++)
            {
                b00 = B_data[j][k];
                Result_data[i][k] += a00 * b00;
            }
        }
    }
}

void *_matrixMulLoopUnrollingParallel(void *arg)
{
    t_arg *args = (t_arg *)arg;
    double **A_data = args->A.data;
    double **B_data = args->B.data;
    double **Result_data = args->Result.data;
    double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;

    int i, j, k;
    int start_row = args->start_row;
    int start_col = args->start_col;
    int end_row = args->end_row;
    int end_col = args->end_col;
    int last_row_A_to_unroll = (end_row - start_row) / UnrollingFactor * UnrollingFactor + start_row;
    int last_col_A_to_unroll = args->A.columns / UnrollingFactor * UnrollingFactor;
    int last_col_B_to_unroll = (end_col - start_col) / UnrollingFactor * UnrollingFactor + start_col;
    for (i = start_row; i < last_row_A_to_unroll; i += UnrollingFactor)
    {
        for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
        {
            a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
            a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
            a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
            a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

            for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
            {
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                Result_data[i][k + 1] += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
                Result_data[i][k + 2] += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
                Result_data[i][k + 3] += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;

                Result_data[i + 1][k] += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                Result_data[i + 1][k + 1] += a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
                Result_data[i + 1][k + 2] += a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
                Result_data[i + 1][k + 3] += a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;

                Result_data[i + 2][k] += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                Result_data[i + 2][k + 1] += a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
                Result_data[i + 2][k + 2] += a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
                Result_data[i + 2][k + 3] += a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;

                Result_data[i + 3][k] += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
                Result_data[i + 3][k + 1] += a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
                Result_data[i + 3][k + 2] += a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
                Result_data[i + 3][k + 3] += a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;
            }
            // If column B is indivisible by unrolling factor. This do the rest.
            for (k = last_col_B_to_unroll; k < end_col; k++)
            {
                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];
                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                Result_data[i + 1][k] += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                Result_data[i + 2][k] += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                Result_data[i + 3][k] += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
            }
        }
        // If column A is indivisible by unrolling factor. This do the rest, while keep k unrolling.
        for (j = last_col_A_to_unroll; j < args->A.columns; j++)
        {
            a00 = A_data[i][j];
            a10 = A_data[i + 1][j];
            a20 = A_data[i + 2][j];
            a30 = A_data[i + 3][j];
            for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
            {
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                Result_data[i][k] += a00 * b00;
                Result_data[i][k + 1] += a00 * b01;
                Result_data[i][k + 2] += a00 * b02;
                Result_data[i][k + 3] += a00 * b03;

                Result_data[i + 1][k] += a10 * b00;
                Result_data[i + 1][k + 1] += a10 * b01;
                Result_data[i + 1][k + 2] += a10 * b02;
                Result_data[i + 1][k + 3] += a10 * b03;

                Result_data[i + 2][k] += a20 * b00;
                Result_data[i + 2][k + 1] += a20 * b01;
                Result_data[i + 2][k + 2] += a20 * b02;
                Result_data[i + 2][k + 3] += a20 * b03;

                Result_data[i + 3][k] += a30 * b00;
                Result_data[i + 3][k + 1] += a30 * b01;
                Result_data[i + 3][k + 2] += a30 * b02;
                Result_data[i + 3][k + 3] += a30 * b03;
            }
            for (k = last_col_B_to_unroll; k < end_col; k++)
            {
                b00 = B_data[j][k];
                Result_data[i][k] += a00 * b00;
                Result_data[i + 1][k] += a10 * b00;
                Result_data[i + 2][k] += a20 * b00;
                Result_data[i + 3][k] += a30 * b00;
            }
        }
    }
    for (i = last_row_A_to_unroll; i < end_row; i++)
    {
        for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
        {
            a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
            for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
            {
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                Result_data[i][k + 1] += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
                Result_data[i][k + 2] += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
                Result_data[i][k + 3] += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
            }
            for (k = last_col_B_to_unroll; k < end_col; k++)
            {
                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];
                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
            }
        }
        for (j = last_col_A_to_unroll; j < args->A.columns; j++)
        {
            a00 = A_data[i][j];
            for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
            {
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                Result_data[i][k] += a00 * b00;
                Result_data[i][k + 1] += a00 * b01;
                Result_data[i][k + 2] += a00 * b02;
                Result_data[i][k + 3] += a00 * b03;
            }
            for (k = last_col_B_to_unroll; k < end_col; k++)
            {
                b00 = B_data[j][k];
                Result_data[i][k] += a00 * b00;
            }
        }
    }
    return NULL;
}

void compareMatrix(matrix A, matrix B)
{
    if (A.rows != B.rows || A.columns != B.columns)
    {
        printf("Dimension of two matrices are different");
        return;
    }
    int cnt = 0;
    double **A_data = ((double **)A.data);
    double **B_data = ((double **)B.data);
    for (int i = 0; i < A.rows; i++)
        for (int j = 0; j < A.columns; j++)
            if ((A_data[i][j] - B_data[i][j]) * (A_data[i][j] - B_data[i][j]) > 1.0E-10)
            {
                printf("%.6f\n", (A_data[i][j] - B_data[i][j]) * (A_data[i][j] - B_data[i][j]));
                cnt++;
            }
}

void initializeMatrixWithRandom(matrix *mat, int rows, int columns)
{
    double *Buff0 = malloc(rows * columns * sizeof(double));
    double **Buff = malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        Buff[i] = Buff0 + i * columns;
    }
    mat->data = Buff;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            mat->data[i][j] = (double)rand() / RAND_MAX;
        }
    }
    mat->rows = rows;
    mat->columns = columns;
}

void initializeMatrixWithZero(matrix *mat, int rows, int columns)
{
    double *Buff0 = malloc(rows * columns * sizeof(double));
    double **Buff = malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        Buff[i] = Buff0 + i * columns;
    }
    mat->data = Buff;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            mat->data[i][j] = 0.0;
        }
    }
    mat->rows = rows;
    mat->columns = columns;
}

void matrixMulSeq(matrix A, matrix B, matrix C, matrix D)
{
    int cost1 = A.rows * A.columns * B.columns + A.columns * C.rows * C.columns; // cost of (A*B)*C
    int cost2 = B.rows * B.columns * C.columns + A.rows * A.columns * C.columns; // cost of A*(B*C)
    matrix TMP, TMP2;

    initializeMatrixWithZero(&TMP, A.rows, B.columns);
    matrixMul(A, B, TMP);
    matrixMul(TMP, C, D);
}

void matrixMulParallel(matrix A, matrix B, matrix C, matrix D)
{
    // In the case dimension is indivisible by num_threads, last thread will handle the rest
    int num_threads_col = sqrt(NUM_THREADS);
    int num_threads_row = NUM_THREADS / num_threads_col;
    pthread_t threads[num_threads_col][num_threads_row];
    t_arg arg[num_threads_row][num_threads_col];
    matrix TMP;
    int cost1 = A.rows * A.columns * B.columns + A.columns * C.rows * C.columns; // cost of (A*B)*C
    int cost2 = B.rows * B.columns * C.columns + A.rows * A.columns * C.columns; // cost of A*(B*C)

    // Assume cost1 < cost2
    initializeMatrixWithZero(&TMP, A.rows, B.columns);
    for (int i = 0; i < num_threads_row; i++)
    {
        for (int j = 0; j < num_threads_col; j++)
        {
            int block_size_row = TMP.rows / num_threads_row;
            int block_size_col = TMP.columns / num_threads_col;
            int start_row = block_size_row * i;
            int start_col = block_size_col * j;
            // If it is the last index -> C.rows
            // else block_size_ * (_ + 1)
            int end_row = i == num_threads_row - 1 ? TMP.rows : block_size_row * (i + 1);
            int end_col = j == num_threads_col - 1 ? TMP.columns : block_size_col * (j + 1);
            arg[i][j].A = A;
            arg[i][j].B = B;
            arg[i][j].Result = TMP;
            arg[i][j].start_row = start_row;
            arg[i][j].start_col = start_col;
            arg[i][j].end_row = end_row;
            arg[i][j].end_col = end_col;
            // printf("Running Thread[%d][%d]\n", i, j);
            // printf("start_row = %d start_col = %d\n", start_row, start_col);
            // printf("end_row = %d  end_col = %d\n", end_row, end_col);
            pthread_create(&threads[i][j], NULL, _matrixMulLoopUnrollingParallel, (void *)&arg[i][j]);
        }
    }

    for (int i = 0; i < num_threads_row; i++)
        for (int j = 0; j < num_threads_col; j++)
        {
            pthread_join(threads[i][j], NULL);
        }

    // printf("Finish first matrix multiplication\n");
    for (int i = 0; i < num_threads_row; i++)
    {
        for (int j = 0; j < num_threads_col; j++)
        {
            int block_size_row = D.rows / num_threads_row;
            int block_size_col = D.columns / num_threads_col;
            int start_row = block_size_row * i;
            int start_col = block_size_col * j;
            // If it is the last index -> C.rows
            // else block_size_ * (_ + 1)
            int end_row = i == num_threads_row - 1 ? D.rows : block_size_row * (i + 1);
            int end_col = j == num_threads_col - 1 ? D.columns : block_size_col * (j + 1);
            arg[i][j].A = TMP;
            arg[i][j].B = C;
            arg[i][j].Result = D;
            arg[i][j].start_row = start_row;
            arg[i][j].start_col = start_col;
            arg[i][j].end_row = end_row;
            arg[i][j].end_col = end_col;
            // printf("Running Thread[%d][%d]\n", i, j);
            // printf("start_row = %d start_col = %d\n", start_row, start_col);
            // printf("end_row = %d  end_col = %d\n", end_row, end_col);
            pthread_create(&threads[i][j], NULL, _matrixMulLoopUnrollingParallel, (void *)&arg[i][j]);
        }
    }
    for (int i = 0; i < num_threads_row; i++)
        for (int j = 0; j < num_threads_col; j++)
        {
            pthread_join(threads[i][j], NULL);
        }

    return;
}

int main(int argc, char *argv[])
{
    srand(6048); // my SID = 540906048

    // check a number of inputs 1 + 6 = 7
    if (argc != 7)
    {
        return 1;
    }

    // declare variables
    int m, k, l, n;
    char *endptr;
    char *mode;
    // Get the Input arguments
    m = strtol(argv[1], &endptr, 10);
    if (endptr == argv[1] || *endptr != '\0')
    {
        fprintf(stderr, "The first argument is incorrect!");
        return -1;
    }
    k = strtol(argv[2], &endptr, 10);
    if (endptr == argv[2] || *endptr != '\0')
    {
        fprintf(stderr, "The second argument is incorrect!");
        return -1;
    }
    l = strtol(argv[3], &endptr, 10);
    if (endptr == argv[3] || *endptr != '\0')
    {
        fprintf(stderr, "The third argument is incorrect!");
        return -1;
    }
    n = strtol(argv[4], &endptr, 10);
    if (endptr == argv[4] || *endptr != '\0')
    {
        fprintf(stderr, "The fourth argument is incorrect!");
        return -1;
    }
    mode = argv[5];
    if (strcmp(mode, "float") != 0 && strcmp(mode, "double") != 0)
    {
        fprintf(stderr, "The fifth argument is incorrect!");
        return -1;
    }

    NUM_THREADS = strtol(argv[6], &endptr, 10);
    if (endptr == argv[6] || *endptr != '\0' || NUM_THREADS > 64 || NUM_THREADS < 0)
    {
        fprintf(stderr, "The sixth argument is incorrect!");
        return -1;
    }

    matrix A, B, C, D, D2;
    initializeMatrixWithRandom(&A, m, k); // Init with random numbers (0-1).
    initializeMatrixWithRandom(&B, k, l);
    initializeMatrixWithRandom(&C, l, n);
    initializeMatrixWithZero(&D, m, n);  // Init with all zeros
    initializeMatrixWithZero(&D2, m, n); // Init with all zeros

    matrixMulSeq(A, B, C, D);
    matrixMulParallel(A, B, C, D2);
    compareMatrix(D, D2);
    // printMatrix(D, 4, 4);
    // printMatrix(D2, 4, 4);
}