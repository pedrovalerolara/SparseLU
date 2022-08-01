#include "tests.h"
#include "matrix_methods.h"

int main(int argc, char *argv[])
{
   // Timing variables
   struct timeval start, end;    // Tracks start and end times of LU factorization
   double flops = 0, flops_dsss; // flops_pdgstrf, flops_dsss; // Flop calculations
   double time;                  // Total time for calculations

   // Routine variables
   int nprocs = 1;
   char *filename = "";
   
   char* output_filename = "/home/80g/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu/out/compare-lass-only.txt";

   if (!(argc == 4 && strcmp("-p", argv[1]) == 0 && (nprocs = atoi(argv[2])) != 0))
   {
      fprintf(stderr, "Usage is: ./main -p [num procs] [-superlu | -lass] < [input file]\n");
      exit(0);
   }

   if (strcmp("-superlu", argv[3]) == 0)

   {
      /* SUPERLU ROUTINE */

      // Run SuperLU pdgstrf routine
      test_superlu(nprocs, &start, &end, &flops);

      // SuperLU pdgstrf execution time
      time = (double)(((end.tv_sec * 1e6 + end.tv_usec) - (start.tv_sec * 1e6 + start.tv_usec)) / 1e6);

      // Print results
       fprintf(stderr, "PDGSTRF TIME(s) = %f\n", time);

      // Print time to file
      FILE *fptr;

      if ((fptr = fopen(output_filename, "a")) == NULL)
      {
         fprintf(stderr, "%s\n", strerror(errno));
         exit(0);
      }

      fprintf(fptr, "%.6f", time);
      fclose(fptr);
   }
   else if (strcmp("-lass", argv[3]) == 0)
   {
      /* SLASS ROUTINE */

      // Run sLASs dsss_dnpgetrf routine
      test_slass(nprocs, filename, &start, &end, &flops);

      // sLASs dsss_dnpgetrf execution time
      time = (double)(((end.tv_sec * 1e6 + end.tv_usec) - (start.tv_sec * 1e6 + start.tv_usec)) / 1e6);

      // FLOPS achieved by the dsss_dnpgetrf routine
      flops_dsss = flops / time;

      // Print results
       fprintf(stderr, "DSSS_DNPGETRF TIME(s) = %f\n", time);
       //fprintf(stderr, "DSSS_DNPGETRF GFLOPS = %1.2f\n", flops_dsss / 1e9);

      // Print time to file
      /*
      FILE *fptr;

      if ((fptr = fopen(output_filename, "a")) == NULL)
      {
         fprintf(stderr, "%s\n", strerror(errno));
         exit(0);
      }

      fprintf(fptr, "%.6f", time);
      fclose(fptr);
      */
      
   }
   else
   {
      fprintf(stderr, "Usage is: ./main -p [num procs] [-superlu | -lass] < [input file]\n");
      exit(0);
   }
}

void test_slass(int nprocs, char *filename,
                struct timeval *start, struct timeval *end,
                double *flops)
{
   // Variables for reading matrix from file
   struct MATRIX *matrix;

   // Variables for sLASs dnpgetrf routine
   double *A;              // A pointer to the start of the flat matrix. Has size m * n * sizeof(double)
   double **TILE_A;        // Tile representation of matrix A
   double **TILE_A_warmup; // Tile representation of matrix A_warmup
   int i, j;               // Iterators
   int_t m = 0;            // Number of rows in the matrix
   int_t n = 0;            // Number of columns in the matrix
   int mt, nt;             // Number of tiles in column and row directions, respectively
   int ret;
   int ldA;

   // Read the input matrix stored in Rutherford-Boeing format.
   matrix = (struct MATRIX *)malloc(sizeof(struct MATRIX));
   read_rb_file(matrix);
   m = matrix->nrow;
   n = matrix->ncol;

   // Checking dimensions
   if (m <= 0)
   {
      fprintf(stderr, "Illegal value of M, M must be >= 0\n");
      exit(0);
   }
   if (n <= 0)
   {
      fprintf(stderr, "Illegal value of N, N must be >= 0\n");
      exit(0);
   }

   ldA = n;

   // Number of tiles
   if (m % TILE_SIZE == 0)
   {
      mt = m / TILE_SIZE;
   }

   else
   {
      mt = (m / TILE_SIZE) + 1;
   }

   if (n % TILE_SIZE == 0)
   {
      nt = n / TILE_SIZE;
   }

   else
   {
      nt = (n / TILE_SIZE) + 1;
   }

   // Allocate memory for matrices
   A = (double *)malloc(m * n * sizeof(double));

   // Allocate memory for tiles
   TILE_A = malloc(nt * mt * sizeof(double *));
   TILE_A_warmup = malloc(nt * mt * sizeof(double *));

   // Initialize matrix to all 0s
#pragma omp parallel for private(i, j) shared(m, n)
   for (i = 0; i < m; i++)
   {
      for (j = 0; j < n; j++)
      {
         A[i * m + j] = 0;
      }
   }

   // Convert matrix from CSC format to flat
   csc_to_flat(matrix, A);
   
   // Free matrix memory
   free(matrix->colptr);
   free(matrix->rowind);
   free(matrix->nzval);
   free(matrix);

   // Make the matrix diagonally dominant
   make_diag_dominant(m, n, A);

   // Print the matrix to file
    FILE* ex = fopen("example-matrix.txt", "w");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (A[i * n + j] != 0) {
                fprintf(ex, "%f,", A[i * n + j]);
            }
            else if (A[j * n + i] != 0) {
                fprintf(ex, "%.4f,", A[j * n + i]);
            }
            else {
                fprintf(ex, "0,");
            }
        }
        fprintf(ex, "\n");
    }
    fclose(ex);

   // Convert matrix to tiles
   dsss_dflat2tiled(m, n, A, ldA, mt, nt, &TILE_A);
   dsss_dflat2tiled(m, n, A, ldA, mt, nt, &TILE_A_warmup);

   *flops = FLOPS_DGETRF(m, n);

   // Warmup
   ret = dsss_dnpgetrf(m, n, mt, nt, TILE_A_warmup);

   // Test
   gettimeofday(start, NULL);
   ret = dsss_dnpgetrf(m, n, mt, nt, TILE_A);
   gettimeofday(end, NULL);

   if (ret == NoSuccess)
   {
      fprintf(stderr, "Test unsuccessful\n");
      exit(0);
   }

   // Free memory
   free(A);
   free(TILE_A);
   free(TILE_A_warmup);
}

extern void test_superlu(int nprocs,
                         struct timeval *start, struct timeval *end,
                         double *flops)
{
   // Variables for making matrix diagonally dominant
   struct MATRIX *matrix;
   double *A_read;
   FILE *fptr;

   SuperMatrix A, AC, L, U, B;
   SCPformat *Lstore;
   NCPformat *Ustore;
   superlumt_options_t superlumt_options;
   fact_t fact;
   trans_t trans;
   yes_no_t refact, usepr;
   double u, drop_tol;
   double *a = NULL;
   int_t *asub = NULL, *xa = NULL;
   int_t *perm_c; /* column permutation vector */
   int_t *perm_r; /* row permutations from partial pivoting */
   void *work;
   int_t info, lwork, nrhs, ldx;
   int_t m = 0, n = 0, nnz = 0, permc_spec, panel_size, relax;
   double *rhsb, *xact;
   Gstat_t Gstat;

   /* Default parameters to control factorization. */
   fact = EQUILIBRATE;
   trans = NOTRANS;
   panel_size = sp_ienv(1);
   relax = sp_ienv(2);
   refact = NO;
   u = 1.0;
   usepr = NO;
   drop_tol = 0.0;
   work = NULL;
   lwork = 0;
   nrhs = 1;

   /* Preliminary read/write of matrix to make it diagonally dominant */

   // Make the matrix diagonally dominant and write it to a temporary file
   matrix = (struct MATRIX *)malloc(sizeof(struct MATRIX));
   read_rb_file(matrix);
   fflush(stdin);
   
   m = matrix->nrow;
   n = matrix->ncol;
   A_read = (double *)malloc(m * n * sizeof(double));
   csc_to_flat(matrix, A_read);
   
   make_diag_dominant(m, n, A_read);
   flat_to_csc(A_read, matrix);
   
   // Free A_read
   free(A_read);

   char *temp_filename = "/home/5pv/SummerIntern21/Cameron-SparseLU/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu/temp.rb";
   if ((fptr = fopen(temp_filename, "w")) == NULL)
   {
      fprintf(stderr, "Error: unable to create/open temporary .rb file for writing\n");
      exit(0);
   }
   csc_to_file(matrix, fptr);

   fclose(fptr);
   
   // Free matrix
   free(matrix->colptr);
   free(matrix->rowind);
   free(matrix->nzval);
   free(matrix);

   // Reopen the file and direct it to stdin
   if ((fptr = fopen(temp_filename, "r")) == NULL)
   {
      fprintf(stderr, "Error: unable to read from temporary .rb file\n");
      exit(0);
   }
   dup2(fileno(fptr), fileno(stdin));
   close(fileno(fptr));

   /* Read the input matrix stored in Rutherford-Boeing format. */
   dreadrb(&m, &n, &nnz, &a, &asub, &xa);

   /* Set up the sparse matrix data structure for A. */
   dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

   if (!(rhsb = doubleMalloc(m * nrhs)))
      SUPERLU_ABORT("Malloc fails for rhsb[].");
   dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
   xact = doubleMalloc(n * nrhs);
   ldx = n;
   dGenXtrue(n, nrhs, xact, ldx);
   dFillRHS(trans, nrhs, xact, ldx, &A, &B);

   if (!(perm_r = intMalloc(m)))
      SUPERLU_ABORT("Malloc fails for perm_r[].");
   if (!(perm_c = intMalloc(n)))
      SUPERLU_ABORT("Malloc fails for perm_c[].");
   if (!(superlumt_options.etree = intMalloc(n)))
      SUPERLU_ABORT("Malloc fails for etree[].");
   if (!(superlumt_options.colcnt_h = intMalloc(n)))
      SUPERLU_ABORT("Malloc fails for colcnt_h[].");
   if (!(superlumt_options.part_super_h = intMalloc(n)))
      SUPERLU_ABORT("Malloc fails for colcnt_h[].");

   /********************************
     * THE FIRST TIME FACTORIZATION *
     ********************************/

   /* ------------------------------------------------------------
       Allocate storage and initialize statistics variables. 
       ------------------------------------------------------------*/
   StatAlloc(n, nprocs, panel_size, relax, &Gstat);
   StatInit(n, nprocs, &Gstat);

   /* ------------------------------------------------------------
       Get column permutation vector perm_c[], according to permc_spec:
       permc_spec = 0: natural ordering 
       permc_spec = 1: minimum degree ordering on structure of A'*A
       permc_spec = 2: minimum degree ordering on structure of A'+A
       permc_spec = 3: approximate minimum degree for unsymmetric matrices
       ------------------------------------------------------------*/
   permc_spec = 1;
   get_perm_c(permc_spec, &A, perm_c);

   /* ------------------------------------------------------------
       Initialize the option structure superlumt_options using the
       user-input parameters;
       Apply perm_c to the columns of original A to form AC.
       ------------------------------------------------------------*/
   // todo: time pdgstrf_init
   refact = NO;
   *flops = FLOPS_DGETRF(m, n);
   pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
                u, usepr, drop_tol, perm_c, perm_r,
                work, lwork, &A, &AC, &superlumt_options, &Gstat);

   /* ------------------------------------------------------------
       Compute the LU factorization of A.
       The following routine will create nprocs threads.
       ------------------------------------------------------------*/
   gettimeofday(start, NULL);
   pdgstrf(&superlumt_options, &AC, perm_r, &L, &U, &Gstat, &info);
   gettimeofday(end, NULL);

   /* ------------------------------------------------------------
       Deallocate storage after factorization.
       ------------------------------------------------------------*/
   pxgstrf_finalize(&superlumt_options, &AC);

   printf("\n** Result of sparse LU **\n");
   dinf_norm_error(nrhs, &B, xact); /* Check inf. norm of the error */

   Lstore = (SCPformat *)L.Store;
   Ustore = (NCPformat *)U.Store;
   printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
   printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
   printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
   fflush(stdout);

   SUPERLU_FREE(rhsb);
   SUPERLU_FREE(xact);
   SUPERLU_FREE(perm_r);
   SUPERLU_FREE(perm_c);

   // Destroy_CompCol_Matrix(&A);
   // Destroy_SuperMatrix_Store(&B);
   // if (lwork == 0)
   // {
   //    Destroy_SuperNode_SCP(&L);
   //    Destroy_CompCol_NCP(&U);
   // }
   // else if (lwork > 0)
   // {
   //    SUPERLU_FREE(work);
   // }
   // StatFree(&Gstat);
}
