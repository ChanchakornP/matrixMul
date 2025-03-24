#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>

typedef struct
{
    void **data;
    int rows;
    int columns;
} Matrix;

bool IS_DOUBLE;
int NUM_THREADS;

void initializeMatrix(Matrix *matrix, int rows, int columns);
void freeMatrix(Matrix *matrix);
void matrixMul(Matrix A, Matrix B, Matrix Buffer);
void matrixMulSequencial(Matrix A, Matrix B, Matrix C, Matrix D);
void matrixMulUnrolledSequencial(Matrix A, Matrix B, Matrix C, Matrix D);

void initializeMatrix(Matrix *matrix, int rows, int columns)
{
    if (IS_DOUBLE)
    {
        double *Buff0 = malloc(rows * columns * sizeof(double));
        double **Buff = malloc(rows * sizeof(double *));
        for (int i = 0; i < rows; i++)
        {
            Buff[i] = Buff0 + i * columns;
        }
        matrix->data = (void **)Buff;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                ((double **)matrix->data)[i][j] = 0.0;
            }
        }
    }
    else
    {
        float *Buff0 = malloc(rows * columns * sizeof(float));
        float **Buff = malloc(rows * sizeof(float *));
        for (int i = 0; i < rows; i++)
        {
            Buff[i] = Buff0 + i * columns;
        }
        matrix->data = (void **)Buff;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                ((float **)matrix->data)[i][j] = 0.0;
            }
        }
    }
    matrix->rows = rows;
    matrix->columns = columns;
}

void initializeMatrixRandomNumber(Matrix *matrix, int rows, int columns)
{
    if (IS_DOUBLE)
    {
        double *Buff0 = malloc(rows * columns * sizeof(double));
        double **Buff = malloc(rows * sizeof(double *));
        for (int i = 0; i < rows; i++)
        {
            Buff[i] = Buff0 + i * columns;
        }
        matrix->data = (void **)Buff;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                ((double **)matrix->data)[i][j] = (double)rand() / RAND_MAX;
            }
        }
    }
    else
    {
        float *Buff0 = malloc(rows * columns * sizeof(float));
        float **Buff = malloc(rows * sizeof(float *));
        for (int i = 0; i < rows; i++)
        {
            Buff[i] = Buff0 + i * columns;
        }
        matrix->data = (void **)Buff;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                ((float **)matrix->data)[i][j] = (float)rand() / RAND_MAX;
            }
        }
    }
    matrix->rows = rows;
    matrix->columns = columns;
}

void freeMatrix(Matrix *matrix)
{
    if (IS_DOUBLE)
    {
        free(((double **)matrix->data)[0]);
        free(matrix->data);
    }
    else
    {
        free(((float **)matrix->data)[0]);
        free(matrix->data);
    }
}

void matrixMul(Matrix A, Matrix B, Matrix Buffer)
{
    struct timeval start_time, end_time;
    float seconds, microseconds, elapsed;
    if (IS_DOUBLE)
    {
        gettimeofday(&start_time, 0);
        double **buffer_data = ((double **)Buffer.data);
        double **A_data = ((double **)A.data);
        double **B_data = ((double **)B.data);
        for (int i = 0; i < A.rows; i++)
        {
            for (int k = 0; k < A.columns; k++)
            {
                for (int j = 0; j < B.columns; j++)
                {
                    {
                        buffer_data[i][j] += A_data[i][k] * B_data[k][j];
                    }
                }
            }
        }
        gettimeofday(&end_time, 0);
        seconds = end_time.tv_sec - start_time.tv_sec;
        microseconds = end_time.tv_usec - start_time.tv_usec;
        elapsed = seconds + 1e-6 * microseconds;
        printf("matmul double takes %f seconds to finish the computation.\n\n", elapsed);
    }
    else
    {

        gettimeofday(&start_time, 0);
        float **buffer_data = ((float **)Buffer.data);
        float **A_data = ((float **)A.data);
        float **B_data = ((float **)B.data);
        for (int i = 0; i < A.rows; i++)
        {
            for (int k = 0; k < A.columns; k++)
            {
                for (int j = 0; j < B.columns; j++)
                {
                    {
                        buffer_data[i][j] += A_data[i][k] * B_data[k][j];
                    }
                }
            }
        }
        gettimeofday(&end_time, 0);
        seconds = end_time.tv_sec - start_time.tv_sec;
        microseconds = end_time.tv_usec - start_time.tv_usec;
        elapsed = seconds + 1e-6 * microseconds;
        printf("matmul float takes %f seconds to finish the computation.\n\n", elapsed);
    }
}

void matrixMulUnrolled(Matrix A, Matrix B, Matrix Buffer)
{
    struct timeval start_time, end_time;
    float seconds, microseconds, elapsed;
    if (IS_DOUBLE)
    {
        gettimeofday(&start_time, 0);
        double ai00, ai01, ai02, ai03, ai10, ai11, ai12, ai13, ai20, ai21, ai22, ai23, ai30, ai31, ai32, ai33;
        double bj00, bj01, bj02, bj03, bj10, bj11, bj12, bj13, bj20, bj21, bj22, bj23, bj30, bj31, bj32, bj33;

        double **buffer_data = ((double **)Buffer.data);
        double **A_data = (double **)A.data;
        double **B_data = (double **)B.data;

        int i, j, k;
        int NRA0 = A.rows / 4 * 4;
        int NCA0 = A.columns / 4 * 4;
        int NCB0 = B.columns / 4 * 4;

        for (i = 0; i < NRA0; i += 4)
        {
            for (k = 0; k < NCA0; k += 4)
            {
                ai00 = A_data[i][k];
                ai01 = A_data[i][k + 1];
                ai02 = A_data[i][k + 2];
                ai03 = A_data[i][k + 3];
                ai10 = A_data[i + 1][k];
                ai11 = A_data[i + 1][k + 1];
                ai12 = A_data[i + 1][k + 2];
                ai13 = A_data[i + 1][k + 3];
                ai20 = A_data[i + 2][k];
                ai21 = A_data[i + 2][k + 1];
                ai22 = A_data[i + 2][k + 2];
                ai23 = A_data[i + 2][k + 3];
                ai30 = A_data[i + 3][k];
                ai31 = A_data[i + 3][k + 1];
                ai32 = A_data[i + 3][k + 2];
                ai33 = A_data[i + 3][k + 3];

                for (j = 0; j < NCB0; j += 4)
                {
                    bj00 = B_data[k][j];
                    bj01 = B_data[k][j + 1];
                    bj02 = B_data[k][j + 2];
                    bj03 = B_data[k][j + 3];
                    bj10 = B_data[k + 1][j];
                    bj11 = B_data[k + 1][j + 1];
                    bj12 = B_data[k + 1][j + 2];
                    bj13 = B_data[k + 1][j + 3];
                    bj20 = B_data[k + 2][j];
                    bj21 = B_data[k + 2][j + 1];
                    bj22 = B_data[k + 2][j + 2];
                    bj23 = B_data[k + 2][j + 3];
                    bj30 = B_data[k + 3][j];
                    bj31 = B_data[k + 3][j + 1];
                    bj32 = B_data[k + 3][j + 2];
                    bj33 = B_data[k + 3][j + 3];

                    buffer_data[i][j] = buffer_data[i][j] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                    buffer_data[i][j + 1] = buffer_data[i][j + 1] + ai00 * bj01 + ai01 * bj11 + ai02 * bj21 + ai03 * bj31;
                    buffer_data[i][j + 2] = buffer_data[i][j + 2] + ai00 * bj02 + ai01 * bj12 + ai02 * bj22 + ai03 * bj32;
                    buffer_data[i][j + 3] = buffer_data[i][j + 3] + ai00 * bj03 + ai01 * bj13 + ai02 * bj23 + ai03 * bj33;

                    buffer_data[i + 1][j] = buffer_data[i + 1][j] + ai10 * bj00 + ai11 * bj10 + ai12 * bj20 + ai13 * bj30;
                    buffer_data[i + 1][j + 1] = buffer_data[i + 1][j + 1] + ai10 * bj01 + ai11 * bj11 + ai12 * bj21 + ai13 * bj31;
                    buffer_data[i + 1][j + 2] = buffer_data[i + 1][j + 2] + ai10 * bj02 + ai11 * bj12 + ai12 * bj22 + ai13 * bj32;
                    buffer_data[i + 1][j + 3] = buffer_data[i + 1][j + 3] + ai10 * bj03 + ai11 * bj13 + ai12 * bj23 + ai13 * bj33;

                    buffer_data[i + 2][j] = buffer_data[i + 2][j] + ai20 * bj00 + ai21 * bj10 + ai22 * bj20 + ai23 * bj30;
                    buffer_data[i + 2][j + 1] = buffer_data[i + 2][j + 1] + ai20 * bj01 + ai21 * bj11 + ai22 * bj21 + ai23 * bj31;
                    buffer_data[i + 2][j + 2] = buffer_data[i + 2][j + 2] + ai20 * bj02 + ai21 * bj12 + ai22 * bj22 + ai23 * bj32;
                    buffer_data[i + 2][j + 3] = buffer_data[i + 2][j + 3] + ai20 * bj03 + ai21 * bj13 + ai22 * bj23 + ai23 * bj33;

                    buffer_data[i + 3][j] = buffer_data[i + 3][j] + ai30 * bj00 + ai31 * bj10 + ai32 * bj20 + ai33 * bj30;
                    buffer_data[i + 3][j + 1] = buffer_data[i + 3][j + 1] + ai30 * bj01 + ai31 * bj11 + ai32 * bj21 + ai33 * bj31;
                    buffer_data[i + 3][j + 2] = buffer_data[i + 3][j + 2] + ai30 * bj02 + ai31 * bj12 + ai32 * bj22 + ai33 * bj32;
                    buffer_data[i + 3][j + 3] = buffer_data[i + 3][j + 3] + ai30 * bj03 + ai31 * bj13 + ai32 * bj23 + ai33 * bj33;
                }

                // For elelments in remaining j columns
                for (j = NCB0; j < B.columns; j++)
                {
                    bj00 = B_data[k][j];
                    bj10 = B_data[k + 1][j];
                    bj20 = B_data[k + 2][j];
                    bj30 = B_data[k + 3][j];
                    buffer_data[i][j] = buffer_data[i][j] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                    buffer_data[i + 1][j] = buffer_data[i + 1][j] + ai10 * bj00 + ai11 * bj10 + ai12 * bj20 + ai13 * bj30;
                    buffer_data[i + 2][j] = buffer_data[i + 2][j] + ai20 * bj00 + ai21 * bj10 + ai22 * bj20 + ai23 * bj30;
                    buffer_data[i + 3][j] = buffer_data[i + 3][j] + ai30 * bj00 + ai31 * bj10 + ai32 * bj20 + ai33 * bj30;
                }
            }

            // for the remaining k
            for (k = NCA0; k < A.columns; k++)
            {
                ai00 = A_data[i][k];
                ai10 = A_data[i + 1][k];
                ai20 = A_data[i + 2][k];
                ai30 = A_data[i + 3][k];

                for (j = 0; j < NCB0; j += 4)
                {
                    bj00 = B_data[k][j];
                    bj01 = B_data[k][j + 1];
                    bj02 = B_data[k][j + 2];
                    bj03 = B_data[k][j + 3];

                    buffer_data[i][j] += ai00 * bj00;
                    buffer_data[i][j + 1] += ai00 * bj01;
                    buffer_data[i][j + 2] += ai00 * bj02;
                    buffer_data[i][j + 3] += ai00 * bj03;

                    buffer_data[i + 1][j] += ai10 * bj00;
                    buffer_data[i + 1][j + 1] += ai10 * bj01;
                    buffer_data[i + 1][j + 2] += ai10 * bj02;
                    buffer_data[i + 1][j + 3] += ai10 * bj03;

                    buffer_data[i + 2][j] += ai20 * bj00;
                    buffer_data[i + 2][j + 1] += ai20 * bj01;
                    buffer_data[i + 2][j + 2] += ai20 * bj02;
                    buffer_data[i + 2][j + 3] += ai20 * bj03;

                    buffer_data[i + 3][j] += ai30 * bj00;
                    buffer_data[i + 3][j + 1] += ai30 * bj01;
                    buffer_data[i + 3][j + 2] += ai30 * bj02;
                    buffer_data[i + 3][j + 3] += ai30 * bj03;
                }

                // For elelments in remaining j columns
                for (j = NCB0; j < B.columns; j++)
                {
                    bj00 = B_data[k][j];
                    buffer_data[i][j] += ai00 * bj00;
                    buffer_data[i + 1][j] += ai10 * bj00;
                    buffer_data[i + 2][j] += ai20 * bj00;
                    buffer_data[i + 3][j] += ai30 * bj00;
                }
            }
        }

        // For elements in remaining i rows
        for (i = NRA0; i < A.rows; i++)
        {
            for (k = 0; k < NCA0; k += 4)
            {
                ai00 = A_data[i][k];
                ai01 = A_data[i][k + 1];
                ai02 = A_data[i][k + 2];
                ai03 = A_data[i][k + 3];
                for (j = 0; j < NCB0; j += 4)
                {
                    bj00 = B_data[k][j];
                    bj01 = B_data[k][j + 1];
                    bj02 = B_data[k][j + 2];
                    bj03 = B_data[k][j + 3];
                    bj10 = B_data[k + 1][j];
                    bj11 = B_data[k + 1][j + 1];
                    bj12 = B_data[k + 1][j + 2];
                    bj13 = B_data[k + 1][j + 3];
                    bj20 = B_data[k + 2][j];
                    bj21 = B_data[k + 2][j + 1];
                    bj22 = B_data[k + 2][j + 2];
                    bj23 = B_data[k + 2][j + 3];
                    bj30 = B_data[k + 3][j];
                    bj31 = B_data[k + 3][j + 1];
                    bj32 = B_data[k + 3][j + 2];
                    bj33 = B_data[k + 3][j + 3];

                    buffer_data[i][j] = buffer_data[i][j] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                    buffer_data[i][j + 1] = buffer_data[i][j + 1] + ai00 * bj01 + ai01 * bj11 + ai02 * bj21 + ai03 * bj31;
                    buffer_data[i][j + 2] = buffer_data[i][j + 2] + ai00 * bj02 + ai01 * bj12 + ai02 * bj22 + ai03 * bj32;
                    buffer_data[i][j + 3] = buffer_data[i][j + 3] + ai00 * bj03 + ai01 * bj13 + ai02 * bj23 + ai03 * bj33;
                }

                // For elelments in remaining j columns
                for (j = NCB0; j < B.columns; j++)
                {
                    bj00 = B_data[k][j];
                    bj10 = B_data[k + 1][j];
                    bj20 = B_data[k + 2][j];
                    bj30 = B_data[k + 3][j];
                    buffer_data[i][j] = buffer_data[i][j] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                }
            }

            // for the remaining k
            for (k = NCA0; k < A.columns; k++)
            {
                ai00 = A_data[i][k];
                for (j = 0; j < B.columns; j++)
                    buffer_data[i][j] += ai00 * B_data[k][j];
            }
        }
        gettimeofday(&end_time, 0);
        seconds = end_time.tv_sec - start_time.tv_sec;
        microseconds = end_time.tv_usec - start_time.tv_usec;
        elapsed = seconds + 1e-6 * microseconds;
        printf("matmul with unrolled by the factor of 4 double takes %f seconds to finish the computation.\n\n", elapsed);
    }
    else
    {
        gettimeofday(&start_time, 0);
        float ai00, ai01, ai02, ai03, ai10, ai11, ai12, ai13, ai20, ai21, ai22, ai23, ai30, ai31, ai32, ai33;
        float bj00, bj01, bj02, bj03, bj10, bj11, bj12, bj13, bj20, bj21, bj22, bj23, bj30, bj31, bj32, bj33;

        float **buffer_data = ((float **)Buffer.data);
        float **A_data = (float **)A.data;
        float **B_data = (float **)B.data;

        int i, j, k;
        int NRA0 = A.rows / 4 * 4;
        int NCA0 = A.columns / 4 * 4;
        int NCB0 = B.columns / 4 * 4;

        for (i = 0; i < NRA0; i += 4)
        {
            for (k = 0; k < NCA0; k += 4)
            {
                ai00 = A_data[i][k];
                ai01 = A_data[i][k + 1];
                ai02 = A_data[i][k + 2];
                ai03 = A_data[i][k + 3];
                ai10 = A_data[i + 1][k];
                ai11 = A_data[i + 1][k + 1];
                ai12 = A_data[i + 1][k + 2];
                ai13 = A_data[i + 1][k + 3];
                ai20 = A_data[i + 2][k];
                ai21 = A_data[i + 2][k + 1];
                ai22 = A_data[i + 2][k + 2];
                ai23 = A_data[i + 2][k + 3];
                ai30 = A_data[i + 3][k];
                ai31 = A_data[i + 3][k + 1];
                ai32 = A_data[i + 3][k + 2];
                ai33 = A_data[i + 3][k + 3];

                for (j = 0; j < NCB0; j += 4)
                {
                    bj00 = B_data[k][j];
                    bj01 = B_data[k][j + 1];
                    bj02 = B_data[k][j + 2];
                    bj03 = B_data[k][j + 3];
                    bj10 = B_data[k + 1][j];
                    bj11 = B_data[k + 1][j + 1];
                    bj12 = B_data[k + 1][j + 2];
                    bj13 = B_data[k + 1][j + 3];
                    bj20 = B_data[k + 2][j];
                    bj21 = B_data[k + 2][j + 1];
                    bj22 = B_data[k + 2][j + 2];
                    bj23 = B_data[k + 2][j + 3];
                    bj30 = B_data[k + 3][j];
                    bj31 = B_data[k + 3][j + 1];
                    bj32 = B_data[k + 3][j + 2];
                    bj33 = B_data[k + 3][j + 3];

                    buffer_data[i][j] = buffer_data[i][j] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                    buffer_data[i][j + 1] = buffer_data[i][j + 1] + ai00 * bj01 + ai01 * bj11 + ai02 * bj21 + ai03 * bj31;
                    buffer_data[i][j + 2] = buffer_data[i][j + 2] + ai00 * bj02 + ai01 * bj12 + ai02 * bj22 + ai03 * bj32;
                    buffer_data[i][j + 3] = buffer_data[i][j + 3] + ai00 * bj03 + ai01 * bj13 + ai02 * bj23 + ai03 * bj33;

                    buffer_data[i + 1][j] = buffer_data[i + 1][j] + ai10 * bj00 + ai11 * bj10 + ai12 * bj20 + ai13 * bj30;
                    buffer_data[i + 1][j + 1] = buffer_data[i + 1][j + 1] + ai10 * bj01 + ai11 * bj11 + ai12 * bj21 + ai13 * bj31;
                    buffer_data[i + 1][j + 2] = buffer_data[i + 1][j + 2] + ai10 * bj02 + ai11 * bj12 + ai12 * bj22 + ai13 * bj32;
                    buffer_data[i + 1][j + 3] = buffer_data[i + 1][j + 3] + ai10 * bj03 + ai11 * bj13 + ai12 * bj23 + ai13 * bj33;

                    buffer_data[i + 2][j] = buffer_data[i + 2][j] + ai20 * bj00 + ai21 * bj10 + ai22 * bj20 + ai23 * bj30;
                    buffer_data[i + 2][j + 1] = buffer_data[i + 2][j + 1] + ai20 * bj01 + ai21 * bj11 + ai22 * bj21 + ai23 * bj31;
                    buffer_data[i + 2][j + 2] = buffer_data[i + 2][j + 2] + ai20 * bj02 + ai21 * bj12 + ai22 * bj22 + ai23 * bj32;
                    buffer_data[i + 2][j + 3] = buffer_data[i + 2][j + 3] + ai20 * bj03 + ai21 * bj13 + ai22 * bj23 + ai23 * bj33;

                    buffer_data[i + 3][j] = buffer_data[i + 3][j] + ai30 * bj00 + ai31 * bj10 + ai32 * bj20 + ai33 * bj30;
                    buffer_data[i + 3][j + 1] = buffer_data[i + 3][j + 1] + ai30 * bj01 + ai31 * bj11 + ai32 * bj21 + ai33 * bj31;
                    buffer_data[i + 3][j + 2] = buffer_data[i + 3][j + 2] + ai30 * bj02 + ai31 * bj12 + ai32 * bj22 + ai33 * bj32;
                    buffer_data[i + 3][j + 3] = buffer_data[i + 3][j + 3] + ai30 * bj03 + ai31 * bj13 + ai32 * bj23 + ai33 * bj33;
                }

                // For elelments in remaining j columns
                for (j = NCB0; j < B.columns; j++)
                {
                    bj00 = B_data[k][j];
                    bj10 = B_data[k + 1][j];
                    bj20 = B_data[k + 2][j];
                    bj30 = B_data[k + 3][j];
                    buffer_data[i][j] = buffer_data[i][j] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                    buffer_data[i + 1][j] = buffer_data[i + 1][j] + ai10 * bj00 + ai11 * bj10 + ai12 * bj20 + ai13 * bj30;
                    buffer_data[i + 2][j] = buffer_data[i + 2][j] + ai20 * bj00 + ai21 * bj10 + ai22 * bj20 + ai23 * bj30;
                    buffer_data[i + 3][j] = buffer_data[i + 3][j] + ai30 * bj00 + ai31 * bj10 + ai32 * bj20 + ai33 * bj30;
                }
            }

            // for the remaining k
            for (k = NCA0; k < A.columns; k++)
            {
                ai00 = A_data[i][k];
                ai10 = A_data[i + 1][k];
                ai20 = A_data[i + 2][k];
                ai30 = A_data[i + 3][k];

                for (j = 0; j < NCB0; j += 4)
                {
                    bj00 = B_data[k][j];
                    bj01 = B_data[k][j + 1];
                    bj02 = B_data[k][j + 2];
                    bj03 = B_data[k][j + 3];

                    buffer_data[i][j] += ai00 * bj00;
                    buffer_data[i][j + 1] += ai00 * bj01;
                    buffer_data[i][j + 2] += ai00 * bj02;
                    buffer_data[i][j + 3] += ai00 * bj03;

                    buffer_data[i + 1][j] += ai10 * bj00;
                    buffer_data[i + 1][j + 1] += ai10 * bj01;
                    buffer_data[i + 1][j + 2] += ai10 * bj02;
                    buffer_data[i + 1][j + 3] += ai10 * bj03;

                    buffer_data[i + 2][j] += ai20 * bj00;
                    buffer_data[i + 2][j + 1] += ai20 * bj01;
                    buffer_data[i + 2][j + 2] += ai20 * bj02;
                    buffer_data[i + 2][j + 3] += ai20 * bj03;

                    buffer_data[i + 3][j] += ai30 * bj00;
                    buffer_data[i + 3][j + 1] += ai30 * bj01;
                    buffer_data[i + 3][j + 2] += ai30 * bj02;
                    buffer_data[i + 3][j + 3] += ai30 * bj03;
                }

                // For elelments in remaining j columns
                for (j = NCB0; j < B.columns; j++)
                {
                    bj00 = B_data[k][j];
                    buffer_data[i][j] += ai00 * bj00;
                    buffer_data[i + 1][j] += ai10 * bj00;
                    buffer_data[i + 2][j] += ai20 * bj00;
                    buffer_data[i + 3][j] += ai30 * bj00;
                }
            }
        }

        // For elements in remaining i rows
        for (i = NRA0; i < A.rows; i++)
        {
            for (k = 0; k < NCA0; k += 4)
            {
                ai00 = A_data[i][k];
                ai01 = A_data[i][k + 1];
                ai02 = A_data[i][k + 2];
                ai03 = A_data[i][k + 3];
                for (j = 0; j < NCB0; j += 4)
                {
                    bj00 = B_data[k][j];
                    bj01 = B_data[k][j + 1];
                    bj02 = B_data[k][j + 2];
                    bj03 = B_data[k][j + 3];
                    bj10 = B_data[k + 1][j];
                    bj11 = B_data[k + 1][j + 1];
                    bj12 = B_data[k + 1][j + 2];
                    bj13 = B_data[k + 1][j + 3];
                    bj20 = B_data[k + 2][j];
                    bj21 = B_data[k + 2][j + 1];
                    bj22 = B_data[k + 2][j + 2];
                    bj23 = B_data[k + 2][j + 3];
                    bj30 = B_data[k + 3][j];
                    bj31 = B_data[k + 3][j + 1];
                    bj32 = B_data[k + 3][j + 2];
                    bj33 = B_data[k + 3][j + 3];

                    buffer_data[i][j] = buffer_data[i][j] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                    buffer_data[i][j + 1] = buffer_data[i][j + 1] + ai00 * bj01 + ai01 * bj11 + ai02 * bj21 + ai03 * bj31;
                    buffer_data[i][j + 2] = buffer_data[i][j + 2] + ai00 * bj02 + ai01 * bj12 + ai02 * bj22 + ai03 * bj32;
                    buffer_data[i][j + 3] = buffer_data[i][j + 3] + ai00 * bj03 + ai01 * bj13 + ai02 * bj23 + ai03 * bj33;
                }

                // For elelments in remaining j columns
                for (j = NCB0; j < B.columns; j++)
                {
                    bj00 = B_data[k][j];
                    bj10 = B_data[k + 1][j];
                    bj20 = B_data[k + 2][j];
                    bj30 = B_data[k + 3][j];
                    buffer_data[i][j] = buffer_data[i][j] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                }
            }

            // for the remaining k
            for (k = NCA0; k < A.columns; k++)
            {
                ai00 = A_data[i][k];
                for (j = 0; j < B.columns; j++)
                    buffer_data[i][j] += ai00 * B_data[k][j];
            }
        }
        gettimeofday(&end_time, 0);
        seconds = end_time.tv_sec - start_time.tv_sec;
        microseconds = end_time.tv_usec - start_time.tv_usec;
        elapsed = seconds + 1e-6 * microseconds;
        printf("matmul with unrolled by the factor of 4 float takes %f seconds to finish the computation.\n\n", elapsed);
    }
}

void matrixMulSequencial(Matrix A, Matrix B, Matrix C, Matrix D)
{
    int cost1 = A.rows * A.columns * B.columns + A.columns * C.rows * C.columns; // cost of (A*B)*C
    int cost2 = B.rows * B.columns * C.columns + A.rows * A.columns * C.columns; // cost of A*(B*C)
    Matrix T;

    if (cost1 < cost2)
    {
        initializeMatrix(&T, A.rows, B.columns);
        matrixMul(A, B, T);
        matrixMul(T, C, D);
    }
    else
    {
        initializeMatrix(&T, B.rows, C.columns);
        matrixMul(B, C, T);
        matrixMul(A, T, D);
    }

    freeMatrix(&T);
}
void matrixMulUnrolledSequencial(Matrix A, Matrix B, Matrix C, Matrix D)
{
    int cost1 = A.rows * A.columns * B.columns + A.columns * C.rows * C.columns; // cost of (A*B)*C
    int cost2 = B.rows * B.columns * C.columns + A.rows * A.columns * C.columns; // cost of A*(B*C)
    Matrix T;

    if (cost1 < cost2)
    {
        initializeMatrix(&T, A.rows, B.columns);
        matrixMulUnrolled(A, B, T);
        matrixMulUnrolled(T, C, D);
    }
    else
    {
        initializeMatrix(&T, B.rows, C.columns);
        matrixMulUnrolled(B, C, T);
        matrixMulUnrolled(A, T, D);
    }

    freeMatrix(&T);
}

void printMatrix(Matrix A, int row, int col)
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

void compareMatrix(Matrix A, Matrix B)
{
    if (A.rows != B.rows || A.columns != B.columns)
    {
        printf("Dimension of two matrices are different");
        return;
    }

    printf("Starting comparison\n\n");
    int cnt = 0;
    if (IS_DOUBLE)
    {
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
    else
    {
        float **A_data = ((float **)A.data);
        float **B_data = ((float **)B.data);
        for (int i = 0; i < A.rows; i++)
            for (int j = 0; j < A.columns; j++)
                if ((A_data[i][j] - B_data[i][j]) * (A_data[i][j] - B_data[i][j]) > 1.0E-10)
                {
                    printf("%.6f\n", (A_data[i][j] - B_data[i][j]) * (A_data[i][j] - B_data[i][j]));
                    cnt++;
                }
    }
    if (cnt == 0)
        printf("Done. There are no differences!\n");
    else
        printf("Results are incorrect! The number of different elements is %d\n", cnt);
}

void resetMatrix(Matrix A)
{
    void **A_data;
    A_data = A.data;

    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.columns; j++)
        {
            if (IS_DOUBLE)
            {
                ((double **)A_data)[i][j] = 0;
            }
            else
            {
                ((float **)A_data)[i][j] = 0;
            }
        }
    }
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
    Matrix MatA, MatB, MatC, MatD, MatD2;
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
    if (strcmp(mode, "double") == 0)
    {
        IS_DOUBLE = 1;
    }
    else
    {
        IS_DOUBLE = 0;
    }

    NUM_THREADS = strtol(argv[6], &endptr, 10);
    if (endptr == argv[6] || *endptr != '\0' || NUM_THREADS > 64 || NUM_THREADS < 0)
    {
        fprintf(stderr, "The sixth argument is incorrect!");
        return -1;
    }
    // end of getting the inputs

    initializeMatrixRandomNumber(&MatA, m, k); // Init with random numbers (0-1).
    initializeMatrixRandomNumber(&MatB, k, l);
    initializeMatrixRandomNumber(&MatC, l, n);
    initializeMatrix(&MatD, m, n); // Init with all zeros
    initializeMatrix(&MatD2, m, n);

    printf("Is double : %d\n", IS_DOUBLE);
    // Display first 4x4 Matrix
    printf("Printing Matrix A\n");
    printMatrix(MatA, 4, 4);
    printf("Printing Matrix B\n");
    printMatrix(MatB, 4, 4);
    printf("Printing Matrix C\n");
    printMatrix(MatC, 4, 4);
    printf("Printing Matrix D\n");
    printMatrix(MatD, 4, 4);

    matrixMulSequencial(MatA, MatB, MatC, MatD);
    matrixMulUnrolledSequencial(MatA, MatB, MatC, MatD2);

    compareMatrix(MatD, MatD2);
    return 0;
}