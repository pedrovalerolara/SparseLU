#include "../include/lass.h"

int synthetic_matrix(int seed, int m, int n, double d, double *A) {
	
	srand(seed);
	int total_entries = m * n;
	int nonzeros = total_entries * d;
	int lower, upper;
	int m_rand, n_rand;

	// Diagonal
	for (int i = 0; i < m; i++) {
		A[i * n + i] = 2;
		nonzeros--;
	}

	// Lower-Upper
	lower = 0; upper = m - 1;
	while (nonzeros > 0) {
		do {
			m_rand = (rand() % (upper - lower + 1)) + lower;
			n_rand = (rand() % (upper - lower + 1)) + lower;
		} while (A[m_rand * n + n_rand] != 0);

		A[m_rand * n + n_rand] = 2.0;
		A[n_rand * n + m_rand] = 2.0;
		nonzeros -= 2;
	}

	return 0;

    // for (long i = 0; i < m; i++)
    // {
    //     for (long j = 0; j < n; j++)
    //     {
    //         if (i == j) // Diagonal
    //         {
	// 			rand_val = (double)rand()/RAND_MAX;// *2.0-1.0;
    //             A[i * n + j] = rand_val;
    //         }
    //         else
    //         {
    //             A[i * n + j] = 0;
    //         }
    //     }
    // }

    // for (long i = 0; i < nonz_lower_half; i++)
    // {
    //     do 
    //     {
    //         // what row?
    //         lower = 1;
    //         upper = m - 1;
    //         row = (rand() % (upper - lower + 1)) + lower;

    //         // what col?
    //         lower = 0;
    //         upper = row - 1;
    //         col = (rand() % (upper - lower + 1)) + lower;
    //     } while (A[row * n + col] != 0);

	// 	rand_val = (double)rand()/RAND_MAX; //*2.0-1.0;
    //     A[row * n + col] = rand_val;
	// 	A[col * n + row] = rand_val; // To make symmetric
    // }

	// return 0;
}


/**
 *
 * 	@file ddss_test_dnpgetrf.c
 *
 * 	@brief LASs-DDSs ddss_test_dnpgetrf routine.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@author Boro Sofranac boro.sofranac@bsc.es
 * 	@date 2018-04-06
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSs
 *   
 *  main:
 *	Performs the testing of the ddss_dnpgetrf routine. 
 *
**/

/**
 *
 *	@sa ddss_dnpgetrf
 *	@sa kdnpgetrf
 *	@sa dnpgetrf
 *
 **/

int main( int argc, char const *argv[] )
{

	if ( argc < 3 || argc > 3 )
	{ 		
		fprintf(stderr, "Usage: ./test_sparse_dnpgetrf M N\n" );
		return NoSuccess;
	}

	// Local variables
	int i, j;
	int m = 0, n = 0;
	int mt, nt;
	int k;
	int ldA;
	double sum;
	double error;
	double max_error;
	double count_error;	
	double *A;
	double *A_ref;
	double *A_test;
	double *A_warmup;
	double *L;
	double *U;
	struct timeval start, end;
	double time;
	double flops;
	double flops_ref;
	double flops_ddss;
	enum LASS_RETURN ret;
	lapack_int retval;
	int seed[] = {0, 0, 0 , 1};

	m = atoi( argv[1] );
	n = atoi( argv[2] );

	// Row mayor order
	ldA = n; 
	
	max_error = 0.0;
	count_error = 0.0;

	// Checking inputs
	if ( m < 0 )
	{
		fprintf(stderr, "Illegal value of M, M must be >= 0\n");
		return NoSuccess;
	}
	if ( n < 0 )
	{
		fprintf(stderr, "Illegal value of N, N must be >= 0\n");
		return NoSuccess;
	}
	
	k = MIN( m, n );

	// Number of tiles
	if ( m % TILE_SIZE == 0) {
		mt = m / TILE_SIZE;
	}
	
	else {
		mt = ( m / TILE_SIZE ) + 1;
	}

	if ( n % TILE_SIZE == 0) {
		nt = n / TILE_SIZE;
	}

	else {
		nt = ( n / TILE_SIZE ) + 1;
	}

	// Matrices allocation
	A = ( double * ) malloc( sizeof( double ) * m * n );
	A_ref = ( double * ) malloc( sizeof( double ) * m * n );
	A_test = ( double * ) malloc( sizeof( double ) * m * n );
	A_warmup = ( double * ) malloc( sizeof( double ) * m * n );
	L = ( double * ) calloc(  m * k, sizeof( double ) );
	U = ( double * ) calloc(  k * n, sizeof( double ) );

	// Matrix A initialization
	// retval = LAPACKE_dlarnv(1, seed, m * n, A);
	retval = synthetic_matrix(2, m, n, 0.1, A);
	assert( retval == 0 );

	/* DEBUG -- PRINT MATRIX */
        /*
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			fprintf(stdout, "%2.0f", A[i * n + j]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
        */
   
	// Make it diagonally dominant
	// #pragma omp parallel for private(i, m, j, n) shared(A)
	for ( i = 0; i < m; i++ )
	{
		sum = 0.0;
		for ( j = 0; j < n; j++ )
		{
			sum += A[i*n+j];
		}
		if ( i < k )
			A[i*n + i] += sum;
	}

	// Copy matrices for testing purpose
	memcpy( A_ref, A, m * n * sizeof( double ) );
	memcpy( A_test, A, m * n * sizeof( double ) );
	memcpy( A_warmup, A, m * n * sizeof( double ) );
	
	// Make matrix array tiles
	double **TILE_A =
                malloc(nt * mt * sizeof(double*)); // A double pointer to an arrary of double pointers

	double **TILE_A_ref =
                malloc(nt * mt * sizeof(double*)); // A double pointer to an arrary of double pointers

	double **TILE_A_test =
                malloc(nt * mt * sizeof(double*)); // A double pointer to an arrary of double pointers

	double **TILE_A_warmup =
                malloc(nt * mt * sizeof(double*)); // A double pointer to an arrary of double pointers

	// CONVERT MATRICES TO TILES
	dsss_dflat2tiled(m, n, A, ldA, mt, nt, &TILE_A);
	dsss_dflat2tiled(m, n, A, ldA, mt, nt, &TILE_A_ref);
	dsss_dflat2tiled(m, n, A, ldA, mt, nt, &TILE_A_test);
	dsss_dflat2tiled(m, n, A, ldA, mt, nt, &TILE_A_warmup);
	
	// DGETRF FLOPS
	flops = FLOPS_DGETRF( m, n ); 
	
	gettimeofday( &start, NULL );
	
	dnpgetrf( m, n, A_ref, ldA );
	
	gettimeofday( &end, NULL );
	

	flops_ref = ( double ) flops / ( ( ( end.tv_sec * 1e6 + end.tv_usec )
			    - ( start.tv_sec * 1e6 + start.tv_usec ) ) / 1e6 );
	ret = dsss_dnpgetrf( m, n, mt, nt, TILE_A_warmup );

	gettimeofday( &start, NULL );
	
	ret = dsss_dnpgetrf( m, n, mt, nt, TILE_A );	

	gettimeofday( &end, NULL );
	if ( ret == NoSuccess )
		return ret;
	
	// Execution time
	time = ( double ) ( ( ( end.tv_sec * 1e6 + end.tv_usec ) 
			- ( start.tv_sec * 1e6 + start.tv_usec ) ) / 1e6 );

	// FLOPS achieved by the ddss_dnpgetrf routine
	flops_ddss = flops / time;

	// Make the flat matrix again
	dsss_dtiled2flat(m, n, A, ldA, mt, nt, &TILE_A);
	
	// Copy results from A to L and U
	#pragma omp parallel for
	for ( i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++ )
		{
				if ( i > j )
			{
				L[i*k + j] = A[i*n + j];
			}
			else if ( i < j )
			{
				U[i*n + j] = A[i*n+j];
			}
			else if( i == j )
			{
				L[i*k + j] = 1.0;
				U[i*n + j] = A[i*n + j];
			}
		}
	}

	// Multiply A = L * U
	ddss_dgemm( CblasNoTrans, CblasNoTrans,
				m, n, k,
				1.0, L, k,
				     U, n,
				0.0, A, n );

	// Error computation
	for ( int i = 0; i < m; i++ ) 
	{
		for ( int j = 0; j < n; j++ ) 
		{
			error = fabs( A[i * n + j] - A_test[i * n + j] );
			if ( error > max_error )
			{
				#pragma omp atomic write
				max_error = error;
			}
			count_error += error;
		}
	}

	fprintf( stdout, "Max. error = %f\n", max_error );
	fprintf( stdout, "Av. error = %f\n", count_error / ( m * n ) );
	fprintf( stderr, "DDSs_DNPGETRF TIME(s) = %1.2f\n", time );
	fprintf( stderr, "DDSs_DNPGETRF GFLOPS = %1.2f\n", flops_ddss / 1e9 );
	fprintf( stdout, "REF. DNPGETRF GFLOPS = %1.2f\n", flops_ref / 1e9 );
	
	free( A );
	free( A_ref );
	free( A_test );
	free( A_warmup );
	free( U );
	free( L );

	return Success;

}
