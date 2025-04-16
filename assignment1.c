#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <sys/resource.h>
#include <math.h>

#define UnrollingFactor 4
#define ALIGNMENT 1024
int IS_DOUBLE;
typedef struct
{
    void **data;
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

int NUM_THREADS;

void print_memory_usage()
{
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    printf("Memory usage: %ld KB\n", usage.ru_maxrss);
}

void freeMatrix(matrix *A)
{
    free(A->data[0]);
    free(A->data);
}

float matrixMulDouble(matrix A, matrix B, matrix Result)
{
    float seconds, microseconds, elapsed;
    struct timeval start_time, end_time;
    double **A_data = (double **)A.data;
    double **B_data = (double **)B.data;
    double **Result_data = (double **)Result.data;
    int row_A = A.rows, col_A = A.columns, col_B = B.columns;
    gettimeofday(&start_time, 0);
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
    gettimeofday(&end_time, 0);
    // timer end here
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    return elapsed;
}

float matrixMulFloat(matrix A, matrix B, matrix Result)
{
    float seconds, microseconds, elapsed;
    struct timeval start_time, end_time;

    float **A_data = (float **)A.data;
    float **B_data = (float **)B.data;
    float **Result_data = (float **)Result.data;
    int row_A = A.rows, col_A = A.columns, col_B = B.columns;

    gettimeofday(&start_time, 0);
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
    gettimeofday(&end_time, 0);
    // timer end here
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    return elapsed;
}

void printMatrix(matrix A, int row, int col)
{
    if (IS_DOUBLE)
    {
        double **A_data = ((double **)A.data);
        for (int i = 0; i < A.rows && i < row; i++)
        {
            for (int j = 0; j < A.columns && j < col; j++)
                printf("%.2f  ", A_data[i][j]);
            printf("\n");
        }
    }
    else
    {
        float **A_data = ((float **)A.data);
        for (int i = 0; i < A.rows && i < row; i++)
        {
            for (int j = 0; j < A.columns && j < col; j++)
                printf("%.2f  ", A_data[i][j]);
            printf("\n");
        }
    }
    printf("\n");

    return;
}

float matrixMulLoopUnrollingDouble(matrix A, matrix B, matrix Result)
{
    float seconds, microseconds, elapsed;
    struct timeval start_time, end_time;

    double **A_data = (double **)A.data;
    double **B_data = (double **)B.data;
    double **Result_data = (double **)Result.data;
    register double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    register double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
    register int i, j, k;
    int row_A = A.rows, col_A = A.columns, col_B = B.columns;
    int row_A0 = A.rows / 4 * 4;
    int col_A0 = A.columns / 4 * 4;
    int col_B0 = B.columns / 4 * 4;

    gettimeofday(&start_time, 0);
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
    gettimeofday(&end_time, 0);
    // timer end here
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    return elapsed;
}

float matrixMulLoopUnrollingDoubleIKJ(matrix A, matrix B, matrix Result)
{
    // This function is not being used.
    float seconds, microseconds, elapsed;
    struct timeval start_time, end_time;

    double **A_data = (double **)A.data;
    double **B_data = (double **)B.data;
    double **Result_data = (double **)Result.data;
    double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
    double c00, c01, c02, c03, c10, c11, c12, c13, c20, c21, c22, c23, c30, c31, c32, c33;
    int row_A = A.rows, col_A = A.columns, col_B = B.columns, i, j, k;
    int row_A0 = A.rows / 4 * 4;
    int col_A0 = A.columns / 4 * 4;
    int col_B0 = B.columns / 4 * 4;

    gettimeofday(&start_time, 0);

    for (i = 0; i < row_A0; i += UnrollingFactor)
    {
        for (k = 0; k < col_B0; k += UnrollingFactor)
        {
            c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
            c10 = 0.0, c11 = 0.0, c12 = 0.0, c13 = 0.0;
            c20 = 0.0, c21 = 0.0, c22 = 0.0, c23 = 0.0;
            c30 = 0.0, c31 = 0.0, c32 = 0.0, c33 = 0.0;
            for (j = 0; j < col_A0; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30, c01 += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31, c02 += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32, c03 += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
                c10 += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30, c11 += a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31, c12 += a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32, c13 += a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;
                c20 += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30, c21 += a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31, c22 += a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32, c23 += a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;
                c30 += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30, c31 += a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31, c32 += a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32, c33 += a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;
            }
            // Handle the rest of j
            for (j = col_A0; j < A.columns; j++)
            {
                a00 = A_data[i][j];
                a10 = A_data[i + 1][j];
                a20 = A_data[i + 2][j];
                a30 = A_data[i + 3][j];

                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                c00 += a00 * b00, c01 += a00 * b01, c02 += a00 * b02, c03 += a00 * b03;
                c10 += a10 * b00, c11 += a10 * b01, c12 += a10 * b02, c13 += a10 * b03;
                c20 += a20 * b00, c21 += a20 * b01, c22 += a20 * b02, c23 += a20 * b03;
                c30 += a30 * b00, c31 += a30 * b01, c32 += a30 * b02, c33 += a30 * b03;
            }

            Result_data[i][k] = c00, Result_data[i][k + 1] = c01, Result_data[i][k + 2] = c02, Result_data[i][k + 3] = c03;
            Result_data[i + 1][k] = c10, Result_data[i + 1][k + 1] = c11, Result_data[i + 1][k + 2] = c12, Result_data[i + 1][k + 3] = c13;
            Result_data[i + 2][k] = c20, Result_data[i + 2][k + 1] = c21, Result_data[i + 2][k + 2] = c22, Result_data[i + 2][k + 3] = c23;
            Result_data[i + 3][k] = c30, Result_data[i + 3][k + 1] = c31, Result_data[i + 3][k + 2] = c32, Result_data[i + 3][k + 3] = c33;
        }
    }
    gettimeofday(&end_time, 0);
    // timer end here
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    return elapsed;
}

float matrixMulLoopUnrollingFloat(matrix A, matrix B, matrix Result)
{
    float seconds, microseconds, elapsed;
    struct timeval start_time, end_time;

    float **A_data = (float **)A.data;
    float **B_data = (float **)B.data;
    float **Result_data = (float **)Result.data;
    register float a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    register float b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
    register int i, j, k;
    int row_A = A.rows, col_A = A.columns, col_B = B.columns;
    int row_A0 = A.rows / 4 * 4;
    int col_A0 = A.columns / 4 * 4;
    int col_B0 = B.columns / 4 * 4;

    gettimeofday(&start_time, 0);
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
    gettimeofday(&end_time, 0);
    // timer end here
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    return elapsed;
}

void *_matrixMulLoopUnrollingParallelDouble(void *arg)
{
    t_arg *args = (t_arg *)arg;

    double **A_data = (double **)args->A.data;
    double **B_data = (double **)args->B.data;
    double **Result_data = (double **)args->Result.data;
    register double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    register double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;

    register int i, j, k;
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

void *_matrixMulLoopUnrollingParallelDoubleIKJV1(void *arg)
{
    // This function is not being used.
    t_arg *args = (t_arg *)arg;

    double **A_data = (double **)args->A.data;
    double **B_data = (double **)args->B.data;
    double **Result_data = (double **)args->Result.data;
    register double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    register double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;

    register int i, j, k;
    int start_row = args->start_row;
    int start_col = args->start_col;
    int end_row = args->end_row;
    int end_col = args->end_col;
    int last_row_A_to_unroll = (end_row - start_row) / UnrollingFactor * UnrollingFactor + start_row;
    int last_col_A_to_unroll = args->A.columns / UnrollingFactor * UnrollingFactor;
    int last_col_B_to_unroll = (end_col - start_col) / UnrollingFactor * UnrollingFactor + start_col;

    for (i = start_row; i < last_row_A_to_unroll; i += UnrollingFactor)
    {
        for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
        {
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

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

            // Handle the rest of j
            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                a10 = A_data[i + 1][j];
                a20 = A_data[i + 2][j];
                a30 = A_data[i + 3][j];

                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                Result_data[i][k] += a00 * b00, Result_data[i][k + 1] += a00 * b01, Result_data[i][k + 2] += a00 * b02, Result_data[i][k + 3] += a00 * b03;
                Result_data[i + 1][k] += a10 * b00, Result_data[i + 2][k + 1] += a10 * b01, Result_data[i + 1][k + 2] += a10 * b02, Result_data[i + 1][k + 3] += a10 * b03;
                Result_data[i + 2][k] += a20 * b00, Result_data[i + 2][k + 2] += a20 * b01, Result_data[i + 2][k + 2] += a20 * b02, Result_data[i + 2][k + 3] += a20 * b03;
                Result_data[i + 3][k] += a30 * b00, Result_data[i + 2][k + 3] += a30 * b01, Result_data[i + 3][k + 2] += a30 * b02, Result_data[i + 3][k + 3] += a30 * b03;
            }
        }
        // Handle the rest of k
        for (k = last_col_B_to_unroll; k < end_col; k++)
        {
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];

                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                Result_data[i + 1][k] += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                Result_data[i + 2][k] += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                Result_data[i + 3][k] += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
            }

            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                a10 = A_data[i + 1][j];
                a20 = A_data[i + 2][j];
                a30 = A_data[i + 3][j];
                b00 = B_data[j][k];

                Result_data[i][k] += a00 * b00;
                Result_data[i + 1][k] += a10 * b00;
                Result_data[i + 2][k] += a20 * b00;
                Result_data[i + 3][k] += a30 * b00;
            }
            Result_data[i][k] = Result_data[i][k];
            Result_data[i + 1][k] = Result_data[i + 1][k];
            Result_data[i + 2][k] = Result_data[i + 2][k];
            Result_data[i + 3][k] = Result_data[i + 3][k];
        }
    }
    // Handle the rest of I
    for (i = last_row_A_to_unroll; i < end_row; i++)
    {
        for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
        {
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30, Result_data[i][k + 1] += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31, Result_data[i][k + 2] += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32, Result_data[i][k + 3] += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
            }

            // Handle the rest of j
            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                Result_data[i][k] += a00 * b00, Result_data[i][k + 1] += a00 * b01, Result_data[i][k + 2] += a00 * b02, Result_data[i][k + 3] += a00 * b03;
            }
        }
        // Handle the rest of k
        for (k = last_col_B_to_unroll; k < end_col; k++)
        {
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];

                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];

                Result_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
            }

            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                b00 = B_data[j][k];

                Result_data[i][k] += a00 * b00;
            }
            Result_data[i][k] = Result_data[i][k];
        }
    }

    return NULL;
}

void *_matrixMulLoopUnrollingParallelDoubleIKJV2(void *arg)
{
    // This function is not being used.
    t_arg *args = (t_arg *)arg;

    double **restrict A_data = (double **)args->A.data;
    double **restrict B_data = (double **)args->B.data;
    double **restrict Result_data = (double **)args->Result.data;
    double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
    register double c00, c01, c02, c03, c10, c11, c12, c13, c20, c21, c22, c23, c30, c31, c32, c33;

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
        for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
        {
            c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
            c10 = 0.0, c11 = 0.0, c12 = 0.0, c13 = 0.0;
            c20 = 0.0, c21 = 0.0, c22 = 0.0, c23 = 0.0;
            c30 = 0.0, c31 = 0.0, c32 = 0.0, c33 = 0.0;
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30, c01 += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31, c02 += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32, c03 += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
                c10 += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30, c11 += a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31, c12 += a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32, c13 += a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;
                c20 += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30, c21 += a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31, c22 += a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32, c23 += a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;
                c30 += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30, c31 += a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31, c32 += a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32, c33 += a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;
            }

            // Handle the rest of j
            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                a10 = A_data[i + 1][j];
                a20 = A_data[i + 2][j];
                a30 = A_data[i + 3][j];

                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                c00 += a00 * b00, c01 += a00 * b01, c02 += a00 * b02, c03 += a00 * b03;
                c10 += a10 * b00, c11 += a10 * b01, c12 += a10 * b02, c13 += a10 * b03;
                c20 += a20 * b00, c21 += a20 * b01, c22 += a20 * b02, c23 += a20 * b03;
                c30 += a30 * b00, c31 += a30 * b01, c32 += a30 * b02, c33 += a30 * b03;
            }
            Result_data[i][k] = c00, Result_data[i][k + 1] = c01, Result_data[i][k + 2] = c02, Result_data[i][k + 3] = c03;
            Result_data[i + 1][k] = c10, Result_data[i + 1][k + 1] = c11, Result_data[i + 1][k + 2] = c12, Result_data[i + 1][k + 3] = c13;
            Result_data[i + 2][k] = c20, Result_data[i + 2][k + 1] = c21, Result_data[i + 2][k + 2] = c22, Result_data[i + 2][k + 3] = c23;
            Result_data[i + 3][k] = c30, Result_data[i + 3][k + 1] = c31, Result_data[i + 3][k + 2] = c32, Result_data[i + 3][k + 3] = c33;
        }
        // Handle the rest of k
        for (k = last_col_B_to_unroll; k < end_col; k++)
        {
            c00 = 0.0;
            c10 = 0.0;
            c20 = 0.0;
            c30 = 0.0;
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                c10 += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                c20 += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                c30 += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
            }

            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                a10 = A_data[i + 1][j];
                a20 = A_data[i + 2][j];
                a30 = A_data[i + 3][j];
                b00 = B_data[j][k];

                c00 += a00 * b00;
                c10 += a10 * b00;
                c20 += a20 * b00;
                c30 += a30 * b00;
            }
            Result_data[i][k] = c00;
            Result_data[i + 1][k] = c10;
            Result_data[i + 2][k] = c20;
            Result_data[i + 3][k] = c30;
        }
    }
    // Handle the rest of I
    for (i = last_row_A_to_unroll; i < end_row; i++)
    {
        for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
        {
            c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30, c01 += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31, c02 += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32, c03 += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
            }

            // Handle the rest of j
            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                c00 += a00 * b00, c01 += a00 * b01, c02 += a00 * b02, c03 += a00 * b03;
            }
            Result_data[i][k] = c00, Result_data[i][k + 1] = c01, Result_data[i][k + 2] = c02, Result_data[i][k + 3] = c03;
        }
        // Handle the rest of k
        for (k = last_col_B_to_unroll; k < end_col; k++)
        {
            c00 = 0.0;
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];

                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
            }

            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                b00 = B_data[j][k];

                c00 += a00 * b00;
            }
            Result_data[i][k] = c00;
        }
    }

    return NULL;
}

void *_matrixMulLoopUnrollingParallelFloatIKJV2(void *arg)
{
    // This function is not being used.
    t_arg *args = (t_arg *)arg;

    float **restrict A_data = (float **)args->A.data;
    float **restrict B_data = (float **)args->B.data;
    float **restrict Result_data = (float **)args->Result.data;
    float a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    float b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
    register float c00, c01, c02, c03, c10, c11, c12, c13, c20, c21, c22, c23, c30, c31, c32, c33;

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
        for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
        {
            c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
            c10 = 0.0, c11 = 0.0, c12 = 0.0, c13 = 0.0;
            c20 = 0.0, c21 = 0.0, c22 = 0.0, c23 = 0.0;
            c30 = 0.0, c31 = 0.0, c32 = 0.0, c33 = 0.0;
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30, c01 += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31, c02 += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32, c03 += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
                c10 += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30, c11 += a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31, c12 += a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32, c13 += a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;
                c20 += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30, c21 += a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31, c22 += a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32, c23 += a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;
                c30 += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30, c31 += a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31, c32 += a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32, c33 += a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;
            }

            // Handle the rest of j
            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                a10 = A_data[i + 1][j];
                a20 = A_data[i + 2][j];
                a30 = A_data[i + 3][j];

                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                c00 += a00 * b00, c01 += a00 * b01, c02 += a00 * b02, c03 += a00 * b03;
                c10 += a10 * b00, c11 += a10 * b01, c12 += a10 * b02, c13 += a10 * b03;
                c20 += a20 * b00, c21 += a20 * b01, c22 += a20 * b02, c23 += a20 * b03;
                c30 += a30 * b00, c31 += a30 * b01, c32 += a30 * b02, c33 += a30 * b03;
            }
            Result_data[i][k] = c00, Result_data[i][k + 1] = c01, Result_data[i][k + 2] = c02, Result_data[i][k + 3] = c03;
            Result_data[i + 1][k] = c10, Result_data[i + 1][k + 1] = c11, Result_data[i + 1][k + 2] = c12, Result_data[i + 1][k + 3] = c13;
            Result_data[i + 2][k] = c20, Result_data[i + 2][k + 1] = c21, Result_data[i + 2][k + 2] = c22, Result_data[i + 2][k + 3] = c23;
            Result_data[i + 3][k] = c30, Result_data[i + 3][k + 1] = c31, Result_data[i + 3][k + 2] = c32, Result_data[i + 3][k + 3] = c33;
        }
        // Handle the rest of k
        for (k = last_col_B_to_unroll; k < end_col; k++)
        {
            c00 = 0.0;
            c10 = 0.0;
            c20 = 0.0;
            c30 = 0.0;
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];

                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                c10 += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                c20 += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                c30 += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
            }

            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                a10 = A_data[i + 1][j];
                a20 = A_data[i + 2][j];
                a30 = A_data[i + 3][j];
                b00 = B_data[j][k];

                c00 += a00 * b00;
                c10 += a10 * b00;
                c20 += a20 * b00;
                c30 += a30 * b00;
            }
            Result_data[i][k] = c00;
            Result_data[i + 1][k] = c10;
            Result_data[i + 2][k] = c20;
            Result_data[i + 3][k] = c30;
        }
    }
    // Handle the rest of I
    for (i = last_row_A_to_unroll; i < end_row; i++)
    {
        for (k = start_col; k < last_col_B_to_unroll; k += UnrollingFactor)
        {
            c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30, c01 += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31, c02 += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32, c03 += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
            }

            // Handle the rest of j
            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];

                c00 += a00 * b00, c01 += a00 * b01, c02 += a00 * b02, c03 += a00 * b03;
            }
            Result_data[i][k] = c00, Result_data[i][k + 1] = c01, Result_data[i][k + 2] = c02, Result_data[i][k + 3] = c03;
        }
        // Handle the rest of k
        for (k = last_col_B_to_unroll; k < end_col; k++)
        {
            c00 = 0.0;
            for (j = 0; j < last_col_A_to_unroll; j += UnrollingFactor)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];

                b00 = B_data[j][k];
                b10 = B_data[j + 1][k];
                b20 = B_data[j + 2][k];
                b30 = B_data[j + 3][k];

                c00 += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
            }

            for (j = last_col_A_to_unroll; j < args->A.columns; j++)
            {
                a00 = A_data[i][j];
                b00 = B_data[j][k];

                c00 += a00 * b00;
            }
            Result_data[i][k] = c00;
        }
    }

    return NULL;
}

void *_matrixMulLoopUnrollingParallelFloat(void *arg)
{
    t_arg *args = (t_arg *)arg;
    float **A_data = (float **)args->A.data;
    float **B_data = (float **)args->B.data;
    float **Result_data = (float **)args->Result.data;
    register float a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
    register float b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;

    register int i, j, k;
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
    if (IS_DOUBLE)
    {
        double **A_data = ((double **)A.data);
        double **B_data = ((double **)B.data);
        for (int i = 0; i < A.rows; i++)
            for (int j = 0; j < A.columns; j++)
                // if the abs diff larger than 1.0E-10, it is not equal
                if ((A_data[i][j] - B_data[i][j]) * (A_data[i][j] - B_data[i][j]) > 1.0E-10)
                {
                    printf("%.6f\n", (A_data[i][j] - B_data[i][j]) * (A_data[i][j] - B_data[i][j]));
                    cnt++;
                }
    }
    else
    {
        float **A_data = ((float **)A.data);
        float **B_data = ((float **)B.data);
        for (int i = 0; i < A.rows; i++)
            for (int j = 0; j < A.columns; j++)
            {
                // https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
                float diff = fabsf(A_data[i][j] - B_data[i][j]);
                float maxVal = fmaxf(fabsf(A_data[i][j]), fabsf(B_data[i][j]));
                // roughtly compare -> the first one -> diff is greater than 1e-3
                // precise compare -> the second one -> considers the magnitude
                if (diff > 1e-3 && diff / (maxVal + 1e-6) > 1e-3)
                {
                    printf("%.3f\n", fabs(A_data[i][j] - B_data[i][j]));
                    cnt++;
                }
            }
    }
    if (cnt == 0)
        printf("Compare Matrices Done. There are no differences!\n");
    else
        printf("Results are incorrect! The number of different elements is %d\n", cnt);
}

void initializeMatrixWithRandom(matrix *mat, int rows, int columns)
{
    if (IS_DOUBLE)
    {
        double *Buff0;
        double **Buff;
        if (rows * columns * sizeof(double) % ALIGNMENT == 0)
        {
            Buff0 = aligned_alloc(ALIGNMENT, rows * columns * sizeof(double));
            Buff = aligned_alloc(ALIGNMENT, rows * sizeof(double *));
        }
        else
        {
            Buff0 = malloc(rows * columns * sizeof(double));
            Buff = malloc(rows * sizeof(double *));
        }

        for (int i = 0; i < rows; i++)
        {
            Buff[i] = Buff0 + i * columns;
        }
        mat->data = (void **)Buff;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                ((double **)mat->data)[i][j] = (double)rand() / RAND_MAX;
            }
        }
    }
    else
    {
        float *Buff0;
        float **Buff;
        if (rows * columns * sizeof(float) % ALIGNMENT == 0)
        {
            Buff0 = aligned_alloc(ALIGNMENT, rows * columns * sizeof(float));
            Buff = aligned_alloc(ALIGNMENT, rows * sizeof(float *));
        }
        else
        {
            Buff0 = malloc(rows * columns * sizeof(float));
            Buff = malloc(rows * sizeof(float *));
        }

        for (int i = 0; i < rows; i++)
        {
            Buff[i] = Buff0 + i * columns;
        }
        mat->data = (void **)Buff;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                ((float **)mat->data)[i][j] = (float)rand() / RAND_MAX;
            }
        }
    }
    mat->rows = rows;
    mat->columns = columns;
}

void initializeMatrixWithZero(matrix *mat, int rows, int columns)
{
    if (IS_DOUBLE)
    {
        double *Buff0;
        double **Buff;
        if (rows * columns * sizeof(double) % ALIGNMENT == 0)
        {
            Buff0 = aligned_alloc(ALIGNMENT, rows * columns * sizeof(double));
            Buff = aligned_alloc(ALIGNMENT, rows * sizeof(double *));
        }
        else
        {
            Buff0 = malloc(rows * columns * sizeof(double));
            Buff = malloc(rows * sizeof(double *));
        }

        for (int i = 0; i < rows; i++)
        {
            Buff[i] = Buff0 + i * columns;
        }
        mat->data = (void **)Buff;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                ((double **)mat->data)[i][j] = 0.0;
            }
        }
    }
    else
    {
        float *Buff0;
        float **Buff;
        if (rows * columns * sizeof(float) % ALIGNMENT == 0)
        {
            Buff0 = aligned_alloc(ALIGNMENT, rows * columns * sizeof(float));
            Buff = aligned_alloc(ALIGNMENT, rows * sizeof(float *));
        }
        else
        {
            Buff0 = malloc(rows * columns * sizeof(float));
            Buff = malloc(rows * sizeof(float *));
        }

        for (int i = 0; i < rows; i++)
        {
            Buff[i] = Buff0 + i * columns;
        }
        mat->data = (void **)Buff;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                ((float **)mat->data)[i][j] = 0.0;
            }
        }
    }
    mat->rows = rows;
    mat->columns = columns;
}

float matrixMulSeq(matrix A, matrix B, matrix C, matrix D)
{
    int cost1 = A.rows * A.columns * B.columns + A.rows * B.columns * C.columns; // cost of (A*B)*C
    int cost2 = B.rows * B.columns * C.columns + A.rows * A.columns * C.columns; // cost of A*(B*C)
    matrix TMP;
    float elapsed = 0.0;

    // struct timeval start_time, end_time;
    // float seconds, microseconds, elapsed;
    float (*matrixMul)(matrix A, matrix B, matrix Result) = IS_DOUBLE ? matrixMulLoopUnrollingDouble : matrixMulLoopUnrollingFloat;
    if (cost1 <= cost2)
    {
        // Compute (AxB)xC
        initializeMatrixWithZero(&TMP, A.rows, B.columns);
        elapsed += matrixMul(A, B, TMP);
        elapsed += matrixMul(TMP, C, D);
    }
    else
    {
        // Compute Ax(BxC)
        initializeMatrixWithZero(&TMP, B.rows, C.columns);
        elapsed += matrixMul(B, C, TMP);
        elapsed += matrixMul(A, TMP, D);
    }

    printf("Sequential computing takes %f seconds to finish the computation.\n", elapsed);

    freeMatrix(&TMP);
    return elapsed;
}

float runParallelMultiplication(matrix A, matrix B, matrix Result, int num_threads_row, int num_threads_col, void *(*worker)(void *), t_arg args[][num_threads_col], pthread_t threads[][num_threads_col])
{
    int block_size_row = Result.rows / num_threads_row;
    int block_size_col = Result.columns / num_threads_col;
    float seconds, microseconds, elapsed;
    struct timeval start_time, end_time;

    for (int i = 0; i < num_threads_row; i++)
    {
        for (int j = 0; j < num_threads_col; j++)
        {
            int start_row = block_size_row * i;
            int start_col = block_size_col * j;
            int end_row = (i == num_threads_row - 1) ? Result.rows : block_size_row * (i + 1);
            int end_col = (j == num_threads_col - 1) ? Result.columns : block_size_col * (j + 1);

            args[i][j].A = A;
            args[i][j].B = B;
            args[i][j].Result = Result;
            args[i][j].start_row = start_row;
            args[i][j].start_col = start_col;
            args[i][j].end_row = end_row;
            args[i][j].end_col = end_col;
        }
    }

    // try to the time spend avoid declaring variable and setting up thread arguments
    // timer start from here
    gettimeofday(&start_time, 0);
    for (int i = 0; i < num_threads_row; i++)
    {
        for (int j = 0; j < num_threads_col; j++)
        {
            if (pthread_create(&threads[i][j], NULL, worker, (void *)&args[i][j]) != 0)
            {
                fprintf(stderr, "Thread create failed\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    for (int i = 0; i < num_threads_row; i++)
    {
        for (int j = 0; j < num_threads_col; j++)
        {
            pthread_join(threads[i][j], NULL);
        }
    }
    gettimeofday(&end_time, 0);
    // timer end here
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    return elapsed;
}

float matrixMulParallel(matrix A, matrix B, matrix C, matrix D)
{
    int num_threads_row = NUM_THREADS / sqrt(NUM_THREADS);
    int num_threads_col = NUM_THREADS / num_threads_row;
    // printf("Number thread col = %d\nNumber thread row = %d\n", num_threads_col, num_threads_row);

    pthread_t threads[num_threads_row][num_threads_col];
    t_arg arg[num_threads_row][num_threads_col];

    matrix TMP;
    float elapsed = 0.0;

    int cost1 = A.rows * A.columns * B.columns + A.rows * B.columns * C.columns; // cost of (A*B)*C
    int cost2 = B.rows * B.columns * C.columns + A.rows * A.columns * C.columns; // cost of A*(B*C)

    void *(*_matrixMulLoopUnrollingParallel)(void *arg) = IS_DOUBLE ? _matrixMulLoopUnrollingParallelDouble : _matrixMulLoopUnrollingParallelFloat;

    if (cost1 <= cost2)
    {
        // Compute (A*B)*C
        initializeMatrixWithZero(&TMP, A.rows, B.columns);
        elapsed += runParallelMultiplication(A, B, TMP, num_threads_row, num_threads_col, _matrixMulLoopUnrollingParallel, arg, threads);
        elapsed += runParallelMultiplication(TMP, C, D, num_threads_row, num_threads_col, _matrixMulLoopUnrollingParallel, arg, threads);
    }
    else
    {
        // Compute A*(B*C)
        initializeMatrixWithZero(&TMP, B.rows, C.columns);
        elapsed += runParallelMultiplication(B, C, TMP, num_threads_row, num_threads_col, _matrixMulLoopUnrollingParallel, arg, threads);
        elapsed += runParallelMultiplication(A, TMP, D, num_threads_row, num_threads_col, _matrixMulLoopUnrollingParallel, arg, threads);
    }

    printf("Parallel computing takes %f seconds to finish the computation.\n", elapsed);

    freeMatrix(&TMP);
    return elapsed;
}

int main(int argc, char *argv[])
{
    srand(6048); // my SID = 540906048

    // check a number of inputs 1 + 6 = 7
    if (argc != 7)
    {
        fprintf(stderr, "The input is incorrect!");
        return 1;
    }

    // declare variables
    int m, k, l, n;
    char *endptr;
    char *mode;
    // Get the Input arguments
    m = strtol(argv[1], &endptr, 10);
    if (endptr == argv[1] || *endptr != '\0' || m <= 0)
    {
        fprintf(stderr, "The first argument is incorrect!");
        return -1;
    }
    k = strtol(argv[2], &endptr, 10);
    if (endptr == argv[2] || *endptr != '\0' || k <= 0)
    {
        fprintf(stderr, "The second argument is incorrect!");
        return -1;
    }
    l = strtol(argv[3], &endptr, 10);
    if (endptr == argv[3] || *endptr != '\0' || l <= 0)
    {
        fprintf(stderr, "The third argument is incorrect!");
        return -1;
    }
    n = strtol(argv[4], &endptr, 10);
    if (endptr == argv[4] || *endptr != '\0' || n <= 0)
    {
        fprintf(stderr, "The fourth argument is incorrect!");
        return -1;
    }
    mode = argv[5];
    if (strcmp(mode, "float") == 0)
    {
        IS_DOUBLE = 0;
    }
    else if (strcmp(mode, "double") == 0)
    {
        IS_DOUBLE = 1;
    }
    else
    {
        fprintf(stderr, "The fifth argument is incorrect!");
        return -1;
    }

    NUM_THREADS = strtol(argv[6], &endptr, 10);
    if (endptr == argv[6] || *endptr != '\0' || NUM_THREADS > 64 || NUM_THREADS <= 0)
    {
        fprintf(stderr, "The sixth argument is incorrect!");
        return -1;
    }

    printf("User Inputs:\n\tSize %d %d %d %d\n\tmode=%s\n\tnum_threads=%d\n", m, k, l, n, mode, NUM_THREADS);
    matrix A, B, C, D, D2;
    initializeMatrixWithRandom(&A, m, k); // Init with random numbers (0-1).
    initializeMatrixWithRandom(&B, k, l);
    initializeMatrixWithRandom(&C, l, n);
    initializeMatrixWithZero(&D, m, n);  // Init with all zeros
    initializeMatrixWithZero(&D2, m, n); // Init with all zeros

    float seqTime, parTime, speedUp;
    parTime = matrixMulParallel(A, B, C, D2);
    seqTime = matrixMulSeq(A, B, C, D);
    compareMatrix(D, D2);
    speedUp = seqTime / parTime;
    printf("Speed-up is %.4f\n", speedUp);
    // printMatrix(D, 4, 4);
    // printMatrix(D2, 4, 4);
    print_memory_usage();
}