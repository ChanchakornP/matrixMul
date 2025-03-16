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
        for (int i = 0; i < A.rows; i++)
        {
            for (int j = 0; j < B.columns; j++)
            {
                double **buffer_data = ((double **)Buffer.data); // dereference first
                buffer_data[i][j] = 0;
                for (int k = 0; k < A.columns; k++)
                {
                    buffer_data[i][j] += ((double **)A.data)[i][k] * ((double **)B.data)[k][j];
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
        for (int i = 0; i < A.rows; i++)
        {
            for (int j = 0; j < B.columns; j++)
            {
                float **buffer_data = ((float **)Buffer.data); // dereference first
                buffer_data[i][j] = 0;
                for (int k = 0; k < A.columns; k++)
                {
                    buffer_data[i][j] += ((float **)A.data)[i][k] * ((float **)B.data)[k][j];
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
        double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
        double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;

        double **buffer_data = ((double **)Buffer.data); // dereference first
        double **A_data = (double **)A.data;
        double **B_data = (double **)B.data;

        gettimeofday(&start_time, 0);
        for (int i = 0; i < A.rows; i += 4)
        {
            for (int j = 0; j < A.columns; j += 4)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];
                for (int k = 0; k < B.columns; k += 4)
                {
                    b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                    b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                    b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                    b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                    buffer_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                    buffer_data[i][k + 1] += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
                    buffer_data[i][k + 2] += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
                    buffer_data[i][k + 3] += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;

                    buffer_data[i + 1][k] += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                    buffer_data[i + 1][k + 1] += a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
                    buffer_data[i + 1][k + 2] += a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
                    buffer_data[i + 1][k + 3] += a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;

                    buffer_data[i + 2][k] += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                    buffer_data[i + 2][k + 1] += a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
                    buffer_data[i + 2][k + 2] += a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
                    buffer_data[i + 2][k + 3] += a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;

                    buffer_data[i + 3][k] += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
                    buffer_data[i + 3][k + 1] += a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
                    buffer_data[i + 3][k + 2] += a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
                    buffer_data[i + 3][k + 3] += a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;
                }
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
        float a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
        float b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;

        float **buffer_data = ((float **)Buffer.data); // dereference first
        float **A_data = (float **)A.data;
        float **B_data = (float **)B.data;

        gettimeofday(&start_time, 0);
        for (int i = 0; i < A.rows; i += 4)
        {
            for (int j = 0; j < A.columns; j += 4)
            {
                a00 = A_data[i][j], a01 = A_data[i][j + 1], a02 = A_data[i][j + 2], a03 = A_data[i][j + 3];
                a10 = A_data[i + 1][j], a11 = A_data[i + 1][j + 1], a12 = A_data[i + 1][j + 2], a13 = A_data[i + 1][j + 3];
                a20 = A_data[i + 2][j], a21 = A_data[i + 2][j + 1], a22 = A_data[i + 2][j + 2], a23 = A_data[i + 2][j + 3];
                a30 = A_data[i + 3][j], a31 = A_data[i + 3][j + 1], a32 = A_data[i + 3][j + 2], a33 = A_data[i + 3][j + 3];
                for (int k = 0; k < B.columns; k += 4)
                {
                    b00 = B_data[j][k], b01 = B_data[j][k + 1], b02 = B_data[j][k + 2], b03 = B_data[j][k + 3];
                    b10 = B_data[j + 1][k], b11 = B_data[j + 1][k + 1], b12 = B_data[j + 1][k + 2], b13 = B_data[j + 1][k + 3];
                    b20 = B_data[j + 2][k], b21 = B_data[j + 2][k + 1], b22 = B_data[j + 2][k + 2], b23 = B_data[j + 2][k + 3];
                    b30 = B_data[j + 3][k], b31 = B_data[j + 3][k + 1], b32 = B_data[j + 3][k + 2], b33 = B_data[j + 3][k + 3];

                    buffer_data[i][k] += a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
                    buffer_data[i][k + 1] += a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
                    buffer_data[i][k + 2] += a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
                    buffer_data[i][k + 3] += a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;

                    buffer_data[i + 1][k] += a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
                    buffer_data[i + 1][k + 1] += a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
                    buffer_data[i + 1][k + 2] += a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
                    buffer_data[i + 1][k + 3] += a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;

                    buffer_data[i + 2][k] += a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
                    buffer_data[i + 2][k + 1] += a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
                    buffer_data[i + 2][k + 2] += a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
                    buffer_data[i + 2][k + 3] += a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;

                    buffer_data[i + 3][k] += a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
                    buffer_data[i + 3][k + 1] += a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
                    buffer_data[i + 3][k + 2] += a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
                    buffer_data[i + 3][k + 3] += a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;
                }
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
    void *A0, *B0, *C0, *D0;
    void **A, **B, **C, **D;
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
        for (int i = 0; i < m; i++)
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
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                ((float **)D)[i][j] = 0;
            }
        }
    }
    // end of initializing arrays

    matrixMulSequencial(MatA, MatB, MatC, MatD);
    matrixMulUnrolledSequencial(MatA, MatB, MatC, MatD);
    return 0;
}