#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "errno.h"

int main (int argc, char** argv) {

    char* matrix_file = argv[1];
    FILE *fptr;
    char* line;
    size_t len = 0;
    ssize_t read;
    char* token;
    char* type;
    int m, n, nonz;
    if ((fptr = fopen(matrix_file, "r")) == NULL){
        fprintf(stderr, "%s\n", strerror(errno));
        exit(0);
    }

    // Firstline
    read = getline(&line, &len, fptr);

    // Second line
    read = getline(&line, &len, fptr);

    // Third line (what we're interested in)
    read = getline(&line, &len, fptr);

    // Type
    type = strtok(line, " ");
    token = strtok(NULL, " ");

    // Number of rows
    m = atoi(token);
    token = strtok(NULL, " ");

    // Number of columns
    n = atoi(token);
    if (m != n) exit(0);
    token = strtok(NULL, " ");
    nonz = atoi(token);

    fclose(fptr);

    // Check for correct matrix properties
    double density = (double) nonz / ( (double) m * (double) n );
    // double Q1 = 0.00015659675; // Q1 for densities in SuiteSparse sparse matrix database
    // double Q3 = 0.0070117045; // Q3 for densities in SuiteSparse sparse matrix database
    // double IQR = Q3 - Q1;

    // if ( m != n ||
    // m > 50000 ||
    // m < 1024 ||
    // density < ( Q1 - 1.5 * IQR ) ||
    // density > ( Q3 + 1.5 * IQR ) ||
    // !(strcmp(type, "rsa") == 0 || strcmp(type, "isa") == 0) ) {
    //     fprintf(stderr, "invalid matrix. Type was %s\n. Density was %f\n", type, density);
    //     // Remove the file
    //     if (remove(filename) == 0) {
    //         fprintf(stdout, "Deleted %s\n", filename);
    //     } else {
    //         fprintf(stderr, "Unable to delete %d\n", filename);
    //     }
    // }

    char *out_file = "/home/80g/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu/out/compare-suitesparse.txt";
    if ((fptr = fopen(out_file, "a")) == NULL){
        fprintf(stderr, "%s\n", strerror(errno));
        exit(0);
    }
    // Filename Size Density
    fprintf(fptr, "%s\t%d\t%f", matrix_file, m, density);
    fclose(fptr);

    return 0;
}