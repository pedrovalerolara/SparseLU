#include "../include/lass.h"

/**
 *
 * 	@file ddss_test_dgemm.c
 *
 * 	@brief LASs-DDSs ddss_test_dgemm routine.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-02-10
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSs
 *   
 *  main:
 *	Performs the testing of the ddss_dgemm routine. 
 *
**/

/**
 *
 *	@sa ddss_dgemm
 *	@sa kdgemm
 *
 **/ 

int main( int argc, char const *argv[] )
{

	if ( argc < 6 || argc > 6 )
	{ 		
		fprintf( stderr, "Usage: ./test_dgemm M N K TRANS_A TRANS_B\n" );
		return NoSuccess;
	}

	// Local variables
	int i;
	int m, n, k;
	int ldA, ldB;
	const char *transA_input = NULL;
	const char *transB_input = NULL;
	enum DDSS_TRANS transA = Trans;
	enum DDSS_TRANS transB = Trans;
	double alpha; 
	double beta;
	double error;
	double max_error;
	double count_error;	
	double *A;
	double *B;
	double *C;
	double *A_test;
	double *B_test;
	double *C_test;
	double *A_ref;
	double *B_ref;
	double *C_ref;
	struct timeval start, end;
	double time;
	double flops;
	double flops_ddss; 
	double flops_ref; 
	enum LASS_RETURN ret;
	lapack_int retval;
	int seed[] = {0, 0, 0 , 1};

	m = atoi( argv[1] );
	n = atoi( argv[2] );
	k = atoi( argv[3] );
	
	if ( strlen( argv[4] ) != 1 ) 
	{
		fprintf( stderr,"Illegal value of TRANS_A, TRANS_A can be T or N\n");
		return NoSuccess;
	}
	transA_input = argv[4];	
	
	if ( strlen( argv[5] ) != 1 ) 
	{
		fprintf( stderr,"Illegal value of TRANS_B, TRANS_B can be T or N\n");
		return NoSuccess;
	}
	transB_input = argv[5];	

	max_error = 0.0;
	count_error = 0.0;

	// Checking inputs
	if ( m < 0 )
	{
		fprintf( stderr, "Illegal value of M, M must be >= 0\n");
		return NoSuccess;
	}
	if ( n < 0 )
	{
		fprintf( stderr, "Illegal value of N, N must be >= 0\n");
		return NoSuccess;
	}
	if ( k < 0 )
	{
		fprintf( stderr, "Illegal value of K, K must be >= 0\n");
		return NoSuccess;
	}

	if ( transA_input[0] == 'T' )
	{
		transA = Trans;
		ldA = m; 
	}
	else if ( transA_input[0] == 'N' )
	{
		transA = NoTrans;
		ldA = k; 
	}
	else
	{
		fprintf( stderr, "Illegal value of TRANS_A, TRANS_A can be T or N\n");
		return NoSuccess;
	}
	
	if ( transB_input[0] == 'T' )
	{
		transB = Trans;
		ldB = k; 
	}
	else if ( transB_input[0] == 'N' )
	{
		transB = NoTrans;
		ldB = n; 
	}
	else
	{
		fprintf( stderr, "Illegal value of TRANS_B, TRANS_B can be T or N\n");
		return NoSuccess;
	}

	// Matrices allocation
	A = ( double * ) malloc( sizeof( double ) * m * k );
	B = ( double * ) malloc( sizeof( double ) * k * n );
	C = ( double * ) malloc( sizeof( double ) * m * n );
	A_test = ( double * ) malloc( sizeof( double ) * m * k );
	B_test = ( double * ) malloc( sizeof( double ) * k * n );
	C_test = ( double * ) malloc( sizeof( double ) * m * n );
	A_ref = ( double * ) malloc( sizeof( double ) * m * k );
	B_ref = ( double * ) malloc( sizeof( double ) * k * n );
	C_ref = ( double * ) malloc( sizeof( double ) * m * n );

	// Alpha and beta initialization
	alpha = ( double ) rand() / (double) rand() + DBL_MIN;
	beta  = ( double ) rand() / (double) rand() + DBL_MIN;
 
	// Matrix A initialization
    retval = LAPACKE_dlarnv( 1, seed, m * k, A );
    assert( retval == 0 );
	
	// Matrix B initialization
    retval = LAPACKE_dlarnv( 1, seed, k * n, B );
    assert( retval == 0 );

	// Matrix C initialization
    retval = LAPACKE_dlarnv( 1, seed, m * n, C );
    assert( retval == 0 );

	// Matrix A_test & A_ref initialization
	memcpy( A_test, A, m * k * sizeof( double ) );
	memcpy( A_ref, A, m * k * sizeof( double ) );
	
	// Matrix B_test & B_ref initialization
	memcpy( B_test, B, k * n * sizeof( double ) );
	memcpy( B_ref, B, k * n * sizeof( double ) );
	
	// Matrix C_test & C_ref initialization
	memcpy( C_test, C, m * n * sizeof( double ) );
	memcpy( C_ref, C, m * n * sizeof( double ) );

	// DGEMM FLOPS
	flops = FLOPS_DGEMM( m, n, k ); 
	
	gettimeofday( &start, NULL );

	cblas_dgemm( CblasRowMajor, 
				 ( CBLAS_TRANSPOSE ) transA,
				 ( CBLAS_TRANSPOSE ) transB,
									 m, n, k,
							 		 alpha,      A, ldA,
							 			         B, ldB,
							 		  beta,  C_ref,   n );

	gettimeofday( &end, NULL );

	// FLOPS achieved by the (reference) cblas_dgemm routine
	flops_ref = (double) flops / (((end.tv_sec * 1e6 + end.tv_usec)
				 - (start.tv_sec * 1e6 + start.tv_usec)) / 1e6);
	
	ret = ddss_dgemm( transA, transB, 
				m, n, k,
                alpha, A_test, ldA,
                       B_test, ldB,
                 beta, C_test,   n );
	
	gettimeofday( &start, NULL );
	
	ret = ddss_dgemm( transA, transB, 
				m, n, k,
                alpha, A, ldA,
                       B, ldB,
                 beta, C,   n );
	
	gettimeofday( &end, NULL );
	if ( ret == NoSuccess )
		return ret;

	// TIME consumed by ddss_dgemm	
	time = ( double ) (((end.tv_sec * 1e6 + end.tv_usec)
				 - (start.tv_sec * 1e6 + start.tv_usec)) / 1e6);
	
	// FLOPS achieved by the ddss_dgemm routine
	flops_ddss = flops / time;
		
	// Error computation
	for ( i = 0; i < m * n; i++ )
	{
		error = fabs( C[ i ] - C_ref[ i ] );
		if ( max_error < error )
			max_error = error;
		count_error += error;
	}

	fprintf( stdout, "Max. error = %1.2f\n", max_error );
	fprintf( stdout, "Av. error = %1.2f\n", count_error / ( m * n ) );
	
	fprintf( stdout, "DDSs_DGEMM TIME(s) = %1.2f\n", time );
	fprintf( stdout, "DDSs_DGEMM GFLOPS = %1.2f\n", flops_ddss / 1e9 );
	fprintf( stdout, "REF. DGEMM GFLOPS = %1.2f\n", flops_ref / 1e9 );
	
	free( A );
	free( B );
	free( C );
	free( A_test );
	free( B_test );
	free( C_test );
	free( A_ref );
	free( B_ref );
	free( C_ref );

	return Success;

}
