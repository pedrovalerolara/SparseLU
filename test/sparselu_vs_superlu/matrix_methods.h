#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "errno.h"
#include "stdbool.h"
#include "tests.h"

#define SEED 0 // For random generation of matrices

enum TYPE { RSA, RUA, CSA, CUA, ISA, IUA, PSA };

struct MATRIX {
    int nrow;
    int ncol;
    int nonz;
    int *colptr;
    int *rowind;
    double *nzval;
    enum TYPE type;
};

extern void generate_matrix(int m, int n, double d, double *A);

extern void make_diag_dominant(int m, int n, double *A);

extern void read_rb_file(struct MATRIX *);

extern void flat_to_csc(double *A, struct MATRIX *matrix);

extern void csc_to_flat(struct MATRIX *matrix, double *A);

extern void csc_to_file(struct MATRIX *matrix, FILE *fptr);
