#include "matrix_methods.h"
#include "math.h"

void generate_matrix(int m, int n, double d, double *A)
{
    srand(time(NULL));
    int nonz = ((m * n * d) - m) / 2;
    // fprintf(stderr, "nons: %d\n", nonz);
    int upper, lower, row, col;
	double rand_val = 1.0;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j) // Diagonal
            {
                do {
                    rand_val = (double)rand()/RAND_MAX; //*2.0-1.0;
                } while (rand_val == 0);
                A[i * n + j] = rand_val;
                // A[i * n + j] = 2;
            }
            else
            {
                A[i * n + j] = 0;
            }
        }
    }

    for (int i = 0; i < nonz; i++)
    {
        do 
        {
            // what row?
            lower = 1;
            upper = m - 1;
            row = (rand() % (upper - lower + 1)) + lower;

            // what col?
            lower = 0;
            upper = row - 1;
            col = (rand() % (upper - lower + 1)) + lower;
        } while (A[row * n + col] != 0);

		do {
            rand_val = (double)rand()/RAND_MAX; //*2.0-1.0;
        } while (rand_val == 0);
        A[row * n + col] = rand_val;
        // A[row * n + col] = 2;
    }
}

void make_diag_dominant(int m, int n, double *A)
{
    double sum;
    int k = MIN(m, n);

    #pragma omp parallel for private(sum) shared(A, m, n, k)
    for (int i = 0; i < m; i++)
    {
        sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum += A[i * n + j];
        }
        if (i < k)
        {
            A[i * n + i] += sum;
        }
    }
}

extern void read_rb_file(struct MATRIX *matrix)
{
    FILE *fp = stdin;
    char *line = NULL;
    char *token = NULL;
    char type[4];
    size_t len = 0;
    ssize_t read;
    int i = 0;
    int lines_pointers = 0, lines_row_indices = 0, lines_values = 0;
    int counter = 0;

    // Read in the first line
    if ((read = getline(&line, &len, fp)) == 0)
    {
        fprintf(stderr, "Error reading file\n");
        exit(0);
    }

    // Get the next line, which contains total lines, lines for pointers, lines for row indices, and lines for numerical values
    if ((read = getline(&line, &len, fp)) != -1)
    {

        // Split the line into tokens
        token = strtok(line, " ");
        // total_lines = atoi(token);
        // fprintf(stderr, "Total Lines (excluding header): %d\n", total_lines);

        token = strtok(NULL, " ");
        lines_pointers = atoi(token);
        // fprintf(stderr, "Lines for pointers: %d\n", lines_pointers);

        token = strtok(NULL, " ");
        lines_row_indices = atoi(token);
        // fprintf(stderr, "Lines for row indeces: %d\n", lines_row_indices);

        token = strtok(NULL, " ");
        lines_values = atoi(token);
        // fprintf(stderr, "Lines for values: %d\n", lines_values);
    }

    // Third line of header
    if ((read = getline(&line, &len, fp)) != -1)
    {

        // Matrix type
        token = strtok(line, " ");
        strcpy(type, token);
        // fprintf(stderr, "Type: %s\n", type);

        // Set matrix type
        if (strcmp("rsa", type) == 0)
        {
            matrix->type = RSA;
        }
        else if (strcmp("isa", type) == 0)
        {
            matrix->type = ISA;
        }
        else if (strcmp("psa", type) == 0)
        {
            matrix->type = PSA;
        }
        else
        {
            fprintf(stderr, "Error: unsupported matrix format. Must be RSA or ISA or PSA\n");
            exit(0);
        }

        // Number of rows
        token = strtok(NULL, " ");
        matrix->nrow = atoi(token);
        // fprintf(stderr, "Number of rows: %d\n", matrix->nrow);

        // Number of columns
        token = strtok(NULL, " ");
        matrix->ncol = atoi(token);
        // fprintf(stderr, "Number of columns: %d\n", matrix->ncol);

        // Number of nonzeros
        token = strtok(NULL, " ");
        matrix->nonz = atoi(token);
        // fprintf(stderr, "Number of nonzero entries: %d\n", matrix->nonz);

        // Unused??
        token = strtok(NULL, " ");
        // unused = atoi(token);
        // fprintf(stderr, "Unused: %d\n", unused);
    }

    // Fourth line (ignore, for FORTRAN use)
    if ((read = getline(&line, &len, fp)) == -1)
    {
        fprintf(stderr, "Error reading file\n");
        exit(0);
    }

    // Allocate memory
    matrix->colptr = (int *)malloc(2 * matrix->ncol * sizeof(int));
    matrix->rowind = (int *)malloc(2 * matrix->nonz * sizeof(int));
    matrix->nzval = (double *)malloc(2 * matrix->nonz * sizeof(double));

    counter = 0;
    // Get column pointers
    for (i = 0; i < lines_pointers; i++)
    {
        // Read the line
        if ((read = getline(&line, &len, fp)) == -1)
        {
            fprintf(stderr, "Error reading file\n");
        }

        // Parse all the values from the line
        // First token
        token = strtok(line, " ");
        matrix->colptr[counter++] = atoi(token) - 1; // from 1-indexed to 0-indexed

        // Rest of the tokens in the line
        while ((token = strtok(NULL, " ")) != NULL)
        {
            matrix->colptr[counter++] = atoi(token) - 1; // from 1-indexed to 0-indexed
        }
    }

    counter = 0;

    // Get row indices
    for (i = 0; i < lines_row_indices; i++)
    {
        // Read the line
        if ((read = getline(&line, &len, fp)) == -1)
        {
            fprintf(stderr, "Error reading file\n");
        }

        // Parse all the values from the line
        // First token
        token = strtok(line, " ");
        matrix->rowind[counter++] = atoi(token) - 1; // from 1-indexed to 0-indexed

        // Rest of the tokens in the line
        while ((token = strtok(NULL, " ")) != NULL)
        {
            matrix->rowind[counter++] = atoi(token) - 1; // from 1-indexed to 0-indexed
        }
    }

    counter = 0; // Reset counter

    // Get matrix values
    for (i = 0; i < lines_values; i++)
    {
        // Read the line
        if ((read = getline(&line, &len, fp)) == -1)
        {
            fprintf(stderr, "Error reading file\n");
        }

        // Parse all the values from the line
        // First token
        token = strtok(line, " ");
        sscanf(token, "%lf", &(matrix->nzval[counter++])); // Read the value as a double

        // Rest of the tokens in the line
        while ((token = strtok(NULL, " ")) != NULL)
        {
            sscanf(token, "%lf", &(matrix->nzval[counter++]));
        }
    }
}

void flat_to_csc(double *A, struct MATRIX *matrix)
{
    int i = 0, j = 0;
    int m = matrix->nrow, n = matrix->ncol;
    int val_ind = 0;
    bool first_in_column = false;

    // Traverse A column by column, setting appropriate values
    for (j = 0; j < n; j++)
    {
        // fprintf(stderr, "j: %d, i: %d\n", j, i);
        first_in_column = true;
        for (i = 0; i < m; i++)
        {
            if (A[i * n + j] != 0 && j <= i) // Since all matrices will be symmetric, we only keep lower daig
            {
                // Set its value in nzval
                matrix->nzval[val_ind] = A[i * n + j];
                // Set its row index
                matrix->rowind[val_ind] = i + 1; // From 0 indexed to 1 indexed
                // If it's the first in the column, set the column ptr
                if (first_in_column)
                {
                    matrix->colptr[j] = val_ind + 1;
                }
                first_in_column = false;
                val_ind++;
            }
        }
    }

    // Do the last col ptr
    matrix->colptr[matrix->ncol] = matrix->nonz + 1;
}

void csc_to_flat(struct MATRIX *matrix, double *A)
{
    int n = matrix->ncol;
    int i;
    int col_strt_val_ind, col_end_val_ind, val_ind, row_pos, col_pos;

    // Convert matrix from CSC format to flat
#pragma omp parallel for default(none) \
    shared(n, matrix, A) private(i, col_strt_val_ind, col_end_val_ind, val_ind, row_pos, col_pos)
    for (i = 0; i < n; i++)
    {
        col_strt_val_ind = matrix->colptr[i];
        col_end_val_ind = matrix->colptr[i + 1];

        for (val_ind = col_strt_val_ind; val_ind < col_end_val_ind; val_ind++)
        {
            row_pos = matrix->rowind[val_ind];
            col_pos = i;
            A[row_pos * n + col_pos] = matrix->nzval[val_ind];

            // Make matrix symmetric (symmetric storage formats only store the bottom triangular matrix)
            if (matrix->type == RSA || matrix->type == CSA || matrix->type == ISA)
            {
                A[col_pos * n + row_pos] = matrix->nzval[val_ind];
            }
        }
    }
}

void csc_to_file(struct MATRIX *matrix, FILE *fptr)
{
    int nlines_total, nlines_pointers, nlines_rowind, nlines_vals;
    int max_items_per_line_colptr = 10;
    int max_items_per_line_rowind = 13;
    int max_items_per_line_val = 4;
    int i, counter, tmp_counter = 0;

    // Number of lines
    if ((matrix->ncol + 1) % max_items_per_line_colptr == 0)
    {
        nlines_pointers = (matrix->ncol + 1) / max_items_per_line_colptr;
    }
    else
    {
        nlines_pointers = ((matrix->ncol + 1) / max_items_per_line_colptr) + 1;
    }

    if (matrix->nonz % max_items_per_line_rowind == 0)
    {
        nlines_rowind = matrix->nonz / max_items_per_line_rowind;
    }
    else
    {
        nlines_rowind = matrix->nonz / max_items_per_line_rowind + 1;
    }

    if (matrix->nonz % max_items_per_line_val == 0)
    {
        nlines_vals = matrix->nonz / max_items_per_line_val;
    }
    else
    {
        nlines_vals = matrix->nonz / max_items_per_line_val + 1;
    }

    nlines_total = nlines_pointers + nlines_rowind + nlines_vals;

    // Title line
    fprintf(fptr, "Pajek/Stranke94; 1994; V. Batagelj; ed: V. Batagelj                    |1524    \n");

    // Line 2: lines excluding header, lines pointers, lines row indices, lines numerical values
    fprintf(fptr, "%14d", nlines_total);
    fprintf(fptr, "%14d", nlines_pointers);
    fprintf(fptr, "%14d", nlines_rowind);
    fprintf(fptr, "%14d", nlines_vals);
    fprintf(fptr, "\n");

    // Line 3: Matrix type, num rows, num columns, number of nonzeros, unused
    fprintf(fptr, "%s", "rsa");
    fprintf(fptr, "%25d", matrix->nrow);
    fprintf(fptr, "%14d", matrix->ncol);
    fprintf(fptr, "%14d", matrix->nonz);
    fprintf(fptr, "%14d", 0);
    fprintf(fptr, "\n");

    // (26I3)          (26I3)          (16I5)

    // Line 4: FORTRAN formatting
    fprintf(fptr, "%-16s", "(10I8)");
    fprintf(fptr, "%-16s", "(13I6)");
    fprintf(fptr, "%-20s", "(4E24.12)");
    fprintf(fptr, "\n");

    // Column pointers
    counter = 0;
    for (i = 0; i < matrix->ncol + 1; i++)
    {

        fprintf(fptr, "%8d", matrix->colptr[i]);

        counter++;
        tmp_counter = counter;
        if (counter == max_items_per_line_colptr)
        {
            fprintf(fptr, "\n");
            counter = 0;
        }
    }
    if (tmp_counter != max_items_per_line_colptr)
    {
        fprintf(fptr, "\n");
    }

    // Row indices
    counter = 0;
    for (i = 0; i < matrix->nonz; i++)
    {

        fprintf(fptr, "%6d", matrix->rowind[i]);

        counter++;
        tmp_counter = counter;
        if (counter == max_items_per_line_rowind)
        {
            fprintf(fptr, "\n");
            counter = 0;
        }
    }
    if (tmp_counter != max_items_per_line_rowind)
    {
        fprintf(fptr, "\n");
    }

    // Numerical values
    counter = 0;
    for (i = 0; i < matrix->nonz; i++)
    {

        fprintf(fptr, "%24.12e", matrix->nzval[i]);

        counter++;
        tmp_counter = counter;
        if (counter == max_items_per_line_val)
        {
            fprintf(fptr, "\n");
            counter = 0;
        }
    }
    if (tmp_counter != max_items_per_line_val)
    {
        fprintf(fptr, "\n");
    }
}