#include "../include/lass.h"

/**
 *
 * 	@file ddss_test_dtile.c
 *
 * 	@brief LASs-DDSs ddss_test_dtile routines.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-05-10
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSS
 *   
 *  main:
 *	Performs the testing of the ddss_flat2tiled and ddss_tiled2flat routines. 
 *
**/

/**
 *
 *	@sa ddss_dtile_alloc
 *	@sa ddss_dflat2tiled
 *	@sa ddss_dtiled2flat
 *	@sa ddss_dtile_free
 *
 **/ 

int main()
{

	// Local variables
	int i, j;
	int mt, nt;
	double *A_in;
	double *A_out;

	srand(time(NULL));

	// Initialization of matrix parameters
	int m = rand() % 2047 + 1; // [1, ..., 2048]
	int n = rand() % 2047 + 1;
	int lda = n; // Row-major order

	fprintf( stdout, "MATRIX PARAMETERS:\n" );
	printf( "M=%d\n", m );
	printf( "N=%d\n", n );
	printf( "LDA=%d\n", lda );
	printf( "TILE_SIZE=%d\n", TILE_SIZE );

	// Matrices allocation
	A_in = ( double * ) malloc( sizeof( double ) * m * n );
	A_out = ( double * ) malloc( sizeof( double ) * m * n );

	// Matrix A_in and A_out initialization
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			A_in[ i * n + j ] = ( double ) rand() / (double) rand() 
							  + DBL_MIN;
			A_out[ i * n + j ] = 0.0;
		}
	} 

	// Number of tiles
	if ( m % TILE_SIZE == 0 )
	{
		mt = m / TILE_SIZE;
	}
	else
	{
		mt = ( m / TILE_SIZE ) + 1;
	}

	if ( n % TILE_SIZE == 0 )
	{
		nt = n / TILE_SIZE;
	}
	else
	{
		nt = ( n / TILE_SIZE ) + 1;
	}

	printf( "NUMBER OF TILES IN M ( MT )=%d\n", mt );
	printf( "NUMBER OF TILES IN N ( NT )=%d\n", nt );
	
	// Tile matrix TILE_A 
	// TILE_A allocation
	double (*TILE_A)[nt][TILE_SIZE * TILE_SIZE] = malloc ( mt * nt * 
		TILE_SIZE * TILE_SIZE * sizeof(double) );	

	if ( TILE_A == NULL)
	{
		fprintf( stdout, "Failure in ddss_dtile_alloc for matrix A\n" );
		return NoSuccess;
	}

	// From flat matrix A_in to tiled matrix TILE_A
	ddss_dflat2tiled( m, n, A_in, lda, mt, nt, TILE_A );

	// From tiled matrix TILE_A to flat matrix A_out 
	ddss_dtiled2flat( m, n, A_out, lda, mt, nt, TILE_A );

	#ifdef VERBOSE	
	
	int it, jt;
	int tile_size_m, tile_size_n;
	
	fprintf( stdout, "Matrix A_in\n" );
	
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			printf( "%1.2f ", A_in[ i * n + j ] );
		}
		printf( "\n" );
	} 
		
	fprintf( stdout, "-------\n" );
	fprintf( stdout, "Matrix TILE_A\n" );

	for ( it = 0; it < mt; it++ )
	{
		for ( jt = 0; jt < nt; jt++ )
		{
			printf("\nTILE %d, %d\n", it, jt);
			tile_size_m = ddss_tile_size( m, it );
			tile_size_n = ddss_tile_size( n, jt );
			for ( i = 0; i < tile_size_m; i++ )
			{
				for ( j = 0; j < tile_size_n; j++ )
				{
					printf( "%1.2f ", TILE_A[it][jt][i * tile_size_n + j] );
				}
				printf("\n");
			}
		}
	}
	printf("\n");

	fprintf( stdout, "Matrix A_out\n" );
	
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			printf( "%1.2f ", A_out[ i * n + j ] );
		}
		printf( "\n" );
	}

	#endif 

	// Testing
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			if ( A_in[ i * n + j ] != A_out[ i * n + j ])
			{
				fprintf( stderr, "TEST NOT PASSED\n" );
				return NoSuccess;
			}
		}
	}
    	
	free( A_in );
	free( A_out );
	free( TILE_A );

	fprintf( stdout, "TEST PASSED\n" );
	return Success;

}
