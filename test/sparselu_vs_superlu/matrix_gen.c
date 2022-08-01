#include "matrix_methods.h"

/**
 * This routine generates a sparse symmetric mxn matrix and writes it to a file
 * stored in Rutherford-Boeing format. The matrix will always be an RSA matrix 
 * (Real, Symmetric, Compressed column form)
 * 
 * Command line: ./matrix_gen -mn [dimension of matrix] -d [density factor] -of [output file location]
 * 
 * Command line args:
 * @arg mn -- Dimension: the generated matrix will be mxn where m = n, or the matrix is square. 
 *                      This argument is the dimension for rows and columns
 * @arg d -- Density factor: the matrix will have a density of 1/d nonzeros. Example: a 2000x2000
 *                      matrix with a density factor of 50 will have, on average,
 *                      (1/50)*2000*2000 = 80000 nonzero entries
 * @arg -of -- Output file: The location and filename for the matrix to be written to. Example:
 *                      ./matrices/Synthetic/test.rb
 **/
int main(int argc, char **argv) {
    int m, n;
    double d;
    char* file;
    FILE *fptr;
    double *A;
    struct MATRIX *matrix;
    int i, j;
    int nonz;

    // Argument checking, set variables
    if ( argc != 7 ||
    strcmp("-mn", argv[1]) != 0 ||
    (m = atoi(argv[2])) == 0 ||
    strcmp("-d", argv[3]) != 0 ||
    (d = atof(argv[4])) == 0 ||
    strcmp("-of", argv[5]) != 0) {
        fprintf(stderr, "Usage is: ./matrix_gen -mn [dimension of matrix] -d [density factor] -of [output file location]\n");
        exit(0);
    }

    if (m <= 0 || d <= 0 || d > 1) {
        fprintf(stderr, "Error: invalid dimensions or density factor given. -mn must be > 0, -d must be between 0 and 1\n");
        exit(0);
    }

    n = m;

    file = argv[6];
    if ((fptr = fopen(file, "w")) == NULL) {
        fprintf(stderr, "Error with file \"%s\": %s\n", file, strerror(errno));
        exit(0);
    }

    // Allocate memory
    A = (double *) malloc(m * n * sizeof(double));
    matrix = (struct MATRIX *) malloc(sizeof(struct MATRIX));

    // Generate the lower triangular of matrix A
    generate_matrix(m, n, d, A);

    // // Print for debugging
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < n; j++) {

    //         if (j > i) {
    //             if (A[j * n + i] != 0) {
    //                 fprintf(stdout, " x");
    //             }
    //             else {
    //                 fprintf(stdout, " .");
    //             }
    //         }
    //         else {
    //             if (A[i * n + j] != 0) {
    //                 fprintf(stdout, " x");
    //             }
    //             else {
    //                 fprintf(stdout, " .");
    //             }
    //         }
    //     }
    //     fprintf(stdout, "\n");
    // }

    // Make it diagonally dominant
    // make_diag_dominant(m, n, A);

    // Count the number of nonzeros
    nonz = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (A[i * n + j] > 0) {
                nonz++;
            }
        }
    }

    // // Print the matrix to file
    // FILE* ex = fopen("example-matrix.txt", "w");
    // for (i = 0; i < m; i++) {
    //     for (j = 0; j < n; j++) {
    //         if (A[i * n + j] != 0) {
    //             fprintf(ex, "%f,", A[i * n + j]);
    //         }
    //         else if (A[j * n + i] != 0) {
    //             fprintf(ex, "%.4f,", A[j * n + i]);
    //         }
    //         else {
    //             fprintf(ex, "0,");
    //         }
    //     }
    //     fprintf(ex, "\n");
    // }
    // fclose(ex);

    // Convert from flat array to compressed sparse column format
    matrix->nonz = nonz;
    matrix->nrow = m;
    matrix->ncol = n;
    matrix->type = RSA;
    matrix->colptr = (int *)malloc(2 * matrix->ncol * sizeof(int));
    matrix->rowind = (int *)malloc(2 * matrix->nonz * sizeof(int));
    matrix->nzval = (double *)malloc(2 * matrix->nonz * sizeof(double));
    flat_to_csc(A, matrix);

    // Write to file
    csc_to_file(matrix, fptr);

    // Close file
    fclose(fptr);

    // Free memory
    free(A);
    free(matrix->nzval);
    free(matrix->rowind);
    free(matrix->colptr);

    return 0;
}
