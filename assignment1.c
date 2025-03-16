#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <stdbool.h>

typedef struct
{
    void **data;
    int rows;
    int columns;
} Matrix;

void matrixMulSequencial(Matrix A, Matrix B, Matrix C, Matrix D);
void matrixMulParallel(Matrix A, Matrix B, Matrix C, Matrix D);
void matrixMul(Matrix A, Matrix B, Matrix Buffer);
void matrixMulUnrolledFour(Matrix A, Matrix B, Matrix Buffer);

bool IS_DOUBLE;
int NUM_THREADS;

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
    void *A0, *B0, *C0, *D0, *TMP0;
    void **A, **B, **C, **D, **TMP;
    Matrix MatA, MatB, MatC, MatD;

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
    if (endptr == argv[6] || *endptr != '\0')
    {
        fprintf(stderr, "The sixth argument is incorrect!");
        return -1;
    }
    // end of getting the inputs

    // Initialize Arrays
    MatA.rows = m;
    MatA.columns = k;
    MatB.rows = k;
    MatB.columns = l;
    MatC.rows = l;
    MatC.columns = n;
    MatD.rows = m;
    MatD.columns = n;

    // allocate memory by the mode
    IS_DOUBLE = strcmp(mode, "double") == 0;
    if (IS_DOUBLE)
    {
        A0 = malloc(m * k * sizeof(double));
        B0 = malloc(k * l * sizeof(double));
        C0 = malloc(l * n * sizeof(double));
        D0 = malloc(m * n * sizeof(double));
        A = malloc(m * sizeof(double *));
        B = malloc(k * sizeof(double *));
        C = malloc(l * sizeof(double *));
        D = malloc(m * sizeof(double *));
        MatA.data = (double **)A;
        MatB.data = (double **)B;
        MatC.data = (double **)C;
        MatD.data = (double **)D;
    }
    else
    {
        A0 = malloc(m * k * sizeof(float));
        B0 = malloc(k * l * sizeof(float));
        C0 = malloc(l * n * sizeof(float));
        D0 = malloc(m * n * sizeof(float));
        A = malloc(m * sizeof(float *));
        B = malloc(k * sizeof(float *));
        C = malloc(l * sizeof(float *));
        D = malloc(m * sizeof(float *));
        MatA.data = (float **)A;
        MatB.data = (float **)B;
        MatC.data = (float **)C;
        MatD.data = (float **)D;
    }

    for (int i = 0; i < m; i++)
    {
        A[i] = A0 + i * k;
        D[i] = D0 + i * n;
    }
    for (int i = 0; i < k; i++)
    {
        B[i] = B0 + i * l;
    }
    for (int i = 0; i < l; i++)
    {
        C[i] = C0 + i * n;
    }

    // If the mode is double
    if (IS_DOUBLE)
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < k; j++)
            {
                ((double **)A)[i][j] = (float)rand() / RAND_MAX;
            }
        }

        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < l; j++)
            {
                ((double **)B)[i][j] = (float)rand() / RAND_MAX;
            }
        }

        for (int i = 0; i < l; i++)
        {
            for (int j = 0; j < n; j++)
            {
                ((double **)C)[i][j] = (float)rand() / RAND_MAX;
            }
        }
        for (int i = 0; i < l; i++)
        {
            for (int j = 0; j < n; j++)
            {
                ((double **)D)[i][j] = 0;
            }
        }
    }
    else
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < k; j++)
            {
                ((float **)A)[i][j] = (float)rand() / RAND_MAX;
            }
        }

        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < l; j++)
            {
                ((float **)B)[i][j] = (float)rand() / RAND_MAX;
            }
        }

        for (int i = 0; i < l; i++)
        {
            for (int j = 0; j < n; j++)
            {
                ((float **)C)[i][j] = (float)rand() / RAND_MAX;
            }
        }
        for (int i = 0; i < l; i++)
        {
            for (int j = 0; j < n; j++)
            {
                ((float **)D)[i][j] = 0;
            }
        }
    }
    // end of initializing arrays

    // find the optimal multiplication order
    // (A * B) * C costs = m*k*l + m*l*n
    // A * (B * C) costs = m*k*n + k*l*n
    if ((m * k * l + m * l * n) < (k * l * n + m * k * n))
    {
        matrixMulSequencial(MatA, MatB, MatC, MatD); // (A * B) * C is faster
        // matrixMulParallel(A, B, C, D, numThreads);
    }
    else
    {
        matrixMulSequencial(MatB, MatC, MatA, MatD); // (B * C) * A is faster
        // matrixMulParallel(B, C, A, D, numThreads);
    }

    return 0;
}

void matrixMul(Matrix A, Matrix B, Matrix Buffer)
{
    if (IS_DOUBLE)
    {
        for (int i = 0; i < A.rows; i++)
        {
            for (int j = 0; j < B.columns; j++)
            {
                ((double **)Buffer.data)[i][j] = 0;
                for (int k = 0; k < A.columns; k++)
                {
                    ((double **)Buffer.data)[i][j] += ((double **)A.data)[i][k] * ((double **)B.data)[k][j];
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < A.rows; i++)
        {
            for (int j = 0; j < B.columns; j++)
            {
                ((float **)Buffer.data)[i][j] = 0;
                for (int k = 0; k < A.columns; k++)
                {
                    ((float **)Buffer.data)[i][j] += ((float **)A.data)[i][k] * ((float **)B.data)[k][j];
                }
            }
        }
    }
}

void matrixMulUnrolledFour(Matrix A, Matrix B, Matrix Buffer)
{
    if (IS_DOUBLE)
    {
        double ai00, ai01, ai02, ai03;
        double ai10, ai11, ai12, ai13;
        double ai20, ai21, ai22, ai23;
        double ai30, ai31, ai32, ai33;
        double bj00, bj01, bj02, bj03;
        double bj10, bj11, bj12, bj13;
        double bj20, bj21, bj22, bj23;
        double bj30, bj31, bj32, bj33;

        for (int i = 0; i < A.rows; i += 4)
        {

            for (int j = 0; j < B.columns; j += 4)
            {
                ai00 = ((double **)A.data)[i][j], ai01 = ((double **)A.data)[i][j + 1], ai02 = ((double **)A.data)[i][j + 2], ai03 = ((double **)A.data)[i][j + 3];
                ai10 = ((double **)A.data)[i + 1][j], ai11 = ((double **)A.data)[i + 1][j + 1], ai12 = ((double **)A.data)[i + 1][j + 2], ai13 = ((double **)A.data)[i + 1][j + 3];
                ai20 = ((double **)A.data)[i + 2][j], ai21 = ((double **)A.data)[i + 2][j + 1], ai22 = ((double **)A.data)[i + 2][j + 2], ai23 = ((double **)A.data)[i + 2][j + 3];
                ai30 = ((double **)A.data)[i + 3][j], ai31 = ((double **)A.data)[i + 3][j + 1], ai32 = ((double **)A.data)[i + 3][j + 2], ai33 = ((double **)A.data)[i + 3][j + 3];
                for (int k = 0; k < A.columns; k += 4)
                {
                    bj00 = ((double **)B.data)[j][k], bj01 = ((double **)B.data)[j][k + 1], bj02 = ((double **)B.data)[j][k + 2], bj03 = ((double **)B.data)[j][k + 3];
                    bj10 = ((double **)B.data)[j + 1][k], bj11 = ((double **)B.data)[j + 1][k + 1], bj12 = ((double **)B.data)[j + 1][k + 2], bj13 = ((double **)B.data)[j + 1][k + 3];
                    bj20 = ((double **)B.data)[j + 2][k], bj21 = ((double **)B.data)[j + 2][k + 1], bj22 = ((double **)B.data)[j + 2][k + 2], bj23 = ((double **)B.data)[j + 2][k + 3];
                    bj30 = ((double **)B.data)[j + 3][k], bj31 = ((double **)B.data)[j + 3][k + 1], bj32 = ((double **)B.data)[j + 3][k + 2], bj33 = ((double **)B.data)[j + 3][k + 3];

                    ((double **)Buffer.data)[i][k] = ((double **)Buffer.data)[i][k] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                    ((double **)Buffer.data)[i][k + 1] = ((double **)Buffer.data)[i][k + 1] + ai00 * bj01 + ai01 * bj11 + ai02 * bj21 + ai03 * bj31;
                    ((double **)Buffer.data)[i][k + 2] = ((double **)Buffer.data)[i][k + 2] + ai00 * bj02 + ai01 * bj12 + ai02 * bj22 + ai03 * bj32;
                    ((double **)Buffer.data)[i][k + 3] = ((double **)Buffer.data)[i][k + 3] + ai00 * bj03 + ai01 * bj13 + ai02 * bj23 + ai03 * bj33;

                    ((double **)Buffer.data)[i + 1][k] = ((double **)Buffer.data)[i + 1][k] + ai10 * bj00 + ai11 * bj10 + ai12 * bj20 + ai13 * bj30;
                    ((double **)Buffer.data)[i + 1][k + 1] = ((double **)Buffer.data)[i + 1][k + 1] + ai10 * bj01 + ai11 * bj11 + ai12 * bj21 + ai13 * bj31;
                    ((double **)Buffer.data)[i + 1][k + 2] = ((double **)Buffer.data)[i + 1][k + 2] + ai10 * bj02 + ai11 * bj12 + ai12 * bj22 + ai13 * bj32;
                    ((double **)Buffer.data)[i + 1][k + 3] = ((double **)Buffer.data)[i + 1][k + 3] + ai10 * bj03 + ai11 * bj13 + ai12 * bj23 + ai13 * bj33;

                    ((double **)Buffer.data)[i + 2][k] = ((double **)Buffer.data)[i + 2][k] + ai20 * bj00 + ai21 * bj10 + ai22 * bj20 + ai23 * bj30;
                    ((double **)Buffer.data)[i + 2][k + 1] = ((double **)Buffer.data)[i + 2][k + 1] + ai20 * bj01 + ai21 * bj11 + ai22 * bj21 + ai23 * bj31;
                    ((double **)Buffer.data)[i + 2][k + 2] = ((double **)Buffer.data)[i + 2][k + 2] + ai20 * bj02 + ai21 * bj12 + ai22 * bj22 + ai23 * bj32;
                    ((double **)Buffer.data)[i + 2][k + 3] = ((double **)Buffer.data)[i + 2][k + 3] + ai20 * bj03 + ai21 * bj13 + ai22 * bj23 + ai23 * bj33;

                    ((double **)Buffer.data)[i + 3][k] = ((double **)Buffer.data)[i + 3][k] + ai30 * bj00 + ai31 * bj10 + ai32 * bj20 + ai33 * bj30;
                    ((double **)Buffer.data)[i + 3][k + 1] = ((double **)Buffer.data)[i + 3][k + 1] + ai30 * bj01 + ai31 * bj11 + ai32 * bj21 + ai33 * bj31;
                    ((double **)Buffer.data)[i + 3][k + 2] = ((double **)Buffer.data)[i + 3][k + 2] + ai30 * bj02 + ai31 * bj12 + ai32 * bj22 + ai33 * bj32;
                    ((double **)Buffer.data)[i + 3][k + 3] = ((double **)Buffer.data)[i + 3][k + 3] + ai30 * bj03 + ai31 * bj13 + ai32 * bj23 + ai33 * bj33;
                }
            }
        }
    }
    else
    {
        float ai00, ai01, ai02, ai03;
        float ai10, ai11, ai12, ai13;
        float ai20, ai21, ai22, ai23;
        float ai30, ai31, ai32, ai33;
        float bj00, bj01, bj02, bj03;
        float bj10, bj11, bj12, bj13;
        float bj20, bj21, bj22, bj23;
        float bj30, bj31, bj32, bj33;

        for (int i = 0; i < A.rows; i += 4)
        {

            for (int j = 0; j < B.columns; j += 4)
            {
                ai00 = ((float **)A.data)[i][j], ai01 = ((float **)A.data)[i][j + 1], ai02 = ((float **)A.data)[i][j + 2], ai03 = ((float **)A.data)[i][j + 3];
                ai10 = ((float **)A.data)[i + 1][j], ai11 = ((float **)A.data)[i + 1][j + 1], ai12 = ((float **)A.data)[i + 1][j + 2], ai13 = ((float **)A.data)[i + 1][j + 3];
                ai20 = ((float **)A.data)[i + 2][j], ai21 = ((float **)A.data)[i + 2][j + 1], ai22 = ((float **)A.data)[i + 2][j + 2], ai23 = ((float **)A.data)[i + 2][j + 3];
                ai30 = ((float **)A.data)[i + 3][j], ai31 = ((float **)A.data)[i + 3][j + 1], ai32 = ((float **)A.data)[i + 3][j + 2], ai33 = ((float **)A.data)[i + 3][j + 3];
                for (int k = 0; k < A.columns; k += 4)
                {
                    bj00 = ((float **)B.data)[j][k], bj01 = ((float **)B.data)[j][k + 1], bj02 = ((float **)B.data)[j][k + 2], bj03 = ((float **)B.data)[j][k + 3];
                    bj10 = ((float **)B.data)[j + 1][k], bj11 = ((float **)B.data)[j + 1][k + 1], bj12 = ((float **)B.data)[j + 1][k + 2], bj13 = ((float **)B.data)[j + 1][k + 3];
                    bj20 = ((float **)B.data)[j + 2][k], bj21 = ((float **)B.data)[j + 2][k + 1], bj22 = ((float **)B.data)[j + 2][k + 2], bj23 = ((float **)B.data)[j + 2][k + 3];
                    bj30 = ((float **)B.data)[j + 3][k], bj31 = ((float **)B.data)[j + 3][k + 1], bj32 = ((float **)B.data)[j + 3][k + 2], bj33 = ((float **)B.data)[j + 3][k + 3];

                    ((float **)Buffer.data)[i][k] = ((float **)Buffer.data)[i][k] + ai00 * bj00 + ai01 * bj10 + ai02 * bj20 + ai03 * bj30;
                    ((float **)Buffer.data)[i][k + 1] = ((float **)Buffer.data)[i][k + 1] + ai00 * bj01 + ai01 * bj11 + ai02 * bj21 + ai03 * bj31;
                    ((float **)Buffer.data)[i][k + 2] = ((float **)Buffer.data)[i][k + 2] + ai00 * bj02 + ai01 * bj12 + ai02 * bj22 + ai03 * bj32;
                    ((float **)Buffer.data)[i][k + 3] = ((float **)Buffer.data)[i][k + 3] + ai00 * bj03 + ai01 * bj13 + ai02 * bj23 + ai03 * bj33;

                    ((float **)Buffer.data)[i + 1][k] = ((float **)Buffer.data)[i + 1][k] + ai10 * bj00 + ai11 * bj10 + ai12 * bj20 + ai13 * bj30;
                    ((float **)Buffer.data)[i + 1][k + 1] = ((float **)Buffer.data)[i + 1][k + 1] + ai10 * bj01 + ai11 * bj11 + ai12 * bj21 + ai13 * bj31;
                    ((float **)Buffer.data)[i + 1][k + 2] = ((float **)Buffer.data)[i + 1][k + 2] + ai10 * bj02 + ai11 * bj12 + ai12 * bj22 + ai13 * bj32;
                    ((float **)Buffer.data)[i + 1][k + 3] = ((float **)Buffer.data)[i + 1][k + 3] + ai10 * bj03 + ai11 * bj13 + ai12 * bj23 + ai13 * bj33;

                    ((float **)Buffer.data)[i + 2][k] = ((float **)Buffer.data)[i + 2][k] + ai20 * bj00 + ai21 * bj10 + ai22 * bj20 + ai23 * bj30;
                    ((float **)Buffer.data)[i + 2][k + 1] = ((float **)Buffer.data)[i + 2][k + 1] + ai20 * bj01 + ai21 * bj11 + ai22 * bj21 + ai23 * bj31;
                    ((float **)Buffer.data)[i + 2][k + 2] = ((float **)Buffer.data)[i + 2][k + 2] + ai20 * bj02 + ai21 * bj12 + ai22 * bj22 + ai23 * bj32;
                    ((float **)Buffer.data)[i + 2][k + 3] = ((float **)Buffer.data)[i + 2][k + 3] + ai20 * bj03 + ai21 * bj13 + ai22 * bj23 + ai23 * bj33;

                    ((float **)Buffer.data)[i + 3][k] = ((float **)Buffer.data)[i + 3][k] + ai30 * bj00 + ai31 * bj10 + ai32 * bj20 + ai33 * bj30;
                    ((float **)Buffer.data)[i + 3][k + 1] = ((float **)Buffer.data)[i + 3][k + 1] + ai30 * bj01 + ai31 * bj11 + ai32 * bj21 + ai33 * bj31;
                    ((float **)Buffer.data)[i + 3][k + 2] = ((float **)Buffer.data)[i + 3][k + 2] + ai30 * bj02 + ai31 * bj12 + ai32 * bj22 + ai33 * bj32;
                    ((float **)Buffer.data)[i + 3][k + 3] = ((float **)Buffer.data)[i + 3][k + 3] + ai30 * bj03 + ai31 * bj13 + ai32 * bj23 + ai33 * bj33;
                }
            }
        }
    }
}

void matrixMulSequencial(Matrix A, Matrix B, Matrix C, Matrix D)
{
    Matrix T;
    if (IS_DOUBLE)
    {
        double *Buff0 = malloc(A.rows * B.columns * sizeof(double));
        double **Buff = malloc(A.rows * sizeof(double *));
        for (int i = 0; i < A.rows; i++)
        {
            Buff[i] = Buff0 + i * B.columns;
        }
        T.data = Buff;
        T.rows = A.rows;
        T.columns = B.columns;
        for (int i = 0; i < T.rows; i++)
        {
            for (int j = 0; j < T.columns; j++)
            {
                ((double **)T.data)[i][j] = 0;
            }
        }
        matrixMul(A, B, T);
        matrixMul(T, C, D);
        free(Buff0);
        free(Buff);
    }
    else
    {
        float *Buff0 = malloc(B.columns * C.rows * sizeof(float));
        float **Buff = malloc(B.columns * sizeof(float *));
        for (int i = 0; i < B.columns; i++)
        {
            Buff[i] = Buff0 + i * C.rows;
        }
        T.data = Buff;
        T.rows = A.rows;
        T.columns = B.columns;
        for (int i = 0; i < T.rows; i++)
        {
            for (int j = 0; j < T.columns; j++)
            {
                ((float **)T.data)[i][j] = 0;
            }
        }
        matrixMul(B, C, T);
        matrixMul(A, T, D);
        free(Buff0);
        free(Buff);
    }
}

void matrixMulParallel(Matrix A, Matrix B, Matrix C, Matrix D)
{
    // pthread_t threads[numThreads];
    // for (int i = 0; i < numThreads; i++)
    // {
    //     pthread_create(&threads[i], NULL, (void *(*)(void *))matrixMultiplication, NULL);
    // }
}