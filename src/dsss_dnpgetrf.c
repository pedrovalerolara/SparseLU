#include "../include/lass.h"

/**
 *
 * 	@file dsss_dnpgetrf.c
 *
 * 	@brief sLASs-DDSs dsss_dnpgetrf routine.
 *
 * 	sLASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2019-06-20
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DSSS
 *   
 *	Performs the LU factorization without pivoting of a spare M-by-N matrix A:
 *
 *		A = L * U
 *
 *	where L is a lower triangular ( lower trapezoidal if M > N ) matrix with
 *	unit diagonal elements and U is an upper triangular ( upper trapezoidal 
 *	if M < N ) matrix.
 *
**/

/**
 *	
 *	@param[in]
 *  A  		double *.
 *	    	A is a pointer to a regular matrix of dimension M-by-N.
 *          On exit, if return value is Success, the matrix A is overwriten by 
 *          the factors L and U. The unit diagonal elements of L are not stored.
 *   
 **/

/**
 *
 *	@retval true tile is null
 *	@retval false tile is not null
 *
 **/

/**
 *
 *	@sa dnpgetrf
 *
 **/ 


int is_null( double *A )
{

	if ( A == NULL )
	{
		return 1;
	}
	else
	{
		return 0;
	}	

}


/**
 *	
 *	@param[in]
 *	M       int.
 *          M specifies the number of rows of the matrix A. M >= 0.
 *	
 *	@param[in]
 *	N       int.
 *          N specifies the number of columns of the matrix A. N >= 0.
 *
 *	@param[in,out]
 *  A  		double *.
 *	    	A is a pointer to a regular matrix of dimension M-by-N.
 *          On exit, if return value is Success, the matrix A is overwriten by 
 *          the factors L and U. The unit diagonal elements of L are not stored.
 *   
 **/

/**
 *
 *	@retval Success successful exit
 *	@retval NoSuccess unsuccessful exit
 *
 **/

/**
 *
 *	@sa dnpgetrf
 *
 **/ 

int dsss_dnpgetrf( int M, int N, int MT, int NT, double *TILE_A[ MT * NT ] )
{

	// Local variables
	int mt, nt;
	int mi, mmi, ki, ni, nni, li;
	int tile_size_m;
	int tile_size_mm;
	int tile_size_n;
	int tile_size_nn;
	int tile_size_k;
	int tile_size_km;
	int tile_size_kn;
	int tile_diagonal_null, tile_bellow_null, tile_right_null, tile_null;

	// Argument checking
	if ( M < 0 )
	{
		fprintf( stderr, "Illegal value of M, in dsss_dnpgetrf code\n" );
		return NoSuccess;
	}

	if ( N < 0 )
	{
		fprintf( stderr, "Illegal value of N, in dsss_dnpgetrf code\n" );
		return NoSuccess;
	}

	// Quick return
	if ( MIN( M, N ) == 0 )
	{
		return Success;
	}

	// Number of tiles
	if ( M % TILE_SIZE == 0 )
	{
		mt = M / TILE_SIZE;
	}
	else
	{
		mt = ( M / TILE_SIZE ) + 1;
	}

	if ( N % TILE_SIZE == 0 )
	{
		nt = N / TILE_SIZE;
	}
	else
	{
		nt = ( N / TILE_SIZE ) + 1;
	}

	/****************
	--DNPGETRF tile--
	****************/	
	
	#pragma omp parallel
	#pragma omp master
	{
	for ( ki = 0; ki < MIN( mt, nt ); ki++ )
	{
		tile_size_km = ddss_tile_size( M, ki );
		tile_size_kn = ddss_tile_size( N, ki );

		tile_diagonal_null = is_null( TILE_A[ki * NT + ki] );

		if ( !tile_diagonal_null )
		{
		
			#if defined(LASs_WITH_MKL)
		
			tile_size_k = MIN( tile_size_km, tile_size_kn );
			
			#pragma omp task depend( inout: TILE_A[ki * NT + ki] ) \
                			shared( TILE_A ) \
                			firstprivate( ki, tile_size_kn ) \
					priority( nt )
        		LAPACKE_mkl_dgetrfnpi( CblasRowMajor,
                               		tile_size_km,   tile_size_kn,
                               		tile_size_k,
                               		TILE_A[ki * NT + ki], tile_size_kn );

			#else
			
			#pragma omp task depend( inout: TILE_A[ki * NT + ki] ) \
					shared( TILE_A ) \
					firstprivate( ki, tile_size_kn ) \
					priority( nt )
			dnpgetrf( tile_size_km, tile_size_kn,
				  		TILE_A[ki * NT + ki],
				  		tile_size_kn );

			#endif

		} //end if tile_diagonal_null 
			
		for ( mi = ki + 1; mi < mt; mi++ )
		{
			tile_size_m = ddss_tile_size( M, mi );
		
			tile_bellow_null = is_null( TILE_A[mi * NT + ki] );
			
			if ( !tile_diagonal_null && !tile_bellow_null )
			{
				#pragma omp task depend( in: TILE_A[ki * NT + ki] ) \
						depend( inout: TILE_A[mi * NT + ki] ) \
						shared( TILE_A ) \
						firstprivate( mi, ki, tile_size_kn ) \
						priority( nt )
				cblas_dtrsm( CblasRowMajor, 
				     	 	( CBLAS_SIDE ) Right, ( CBLAS_UPLO ) Upper, 
				     	 	( CBLAS_TRANSPOSE ) NoTrans, ( CBLAS_DIAG ) NonUnit,
					 	 	tile_size_m, tile_size_kn, 
					 	 	1.0, TILE_A[ki * NT + ki], tile_size_kn, 
						  	  	 TILE_A[mi * NT + ki], tile_size_kn );
			}
		}
		for ( ni = ki + 1; ni < nt; ni++ )
		{
			tile_size_n = ddss_tile_size( N, ni );
			
			tile_right_null = is_null( TILE_A[ki * NT + ni] );
			
			if ( !tile_diagonal_null && !tile_right_null )
			{
				#pragma omp task depend( in: TILE_A[ki * NT + ki] ) \
						depend( inout: TILE_A[ki * NT + ni] ) \
						shared( TILE_A ) \
						firstprivate( ni, ki, tile_size_kn ) \
						priority( nt )
				cblas_dtrsm( CblasRowMajor, 
				     	 	( CBLAS_SIDE ) Left, ( CBLAS_UPLO ) Lower, 
				     	 	( CBLAS_TRANSPOSE ) NoTrans, ( CBLAS_DIAG ) Unit, 
					 	 	tile_size_km, tile_size_n, 
					 	 	1.0, TILE_A[ki * NT + ki], tile_size_kn, 
						  	  	 TILE_A[ki * NT + ni], tile_size_n );

			}
		}
		
		for ( nni = ki + 1; nni < nt; nni++ )
		{
			tile_size_nn = ddss_tile_size( N, nni );
			
			tile_right_null = is_null( TILE_A[ki * NT + nni] );

			for ( mmi = ki + 1; mmi < mt; mmi++ )
			{
				tile_size_mm = ddss_tile_size( M, mmi );
			
				tile_bellow_null = is_null( TILE_A[mmi * NT + ki] );
				
				tile_null = is_null( TILE_A[mmi * NT + nni] );

				if ( !tile_bellow_null && !tile_right_null && tile_null )
				{
					TILE_A[mmi * NT + nni] = (double *) malloc ( TILE_SIZE * 
															TILE_SIZE * 
															sizeof( double ) ); 
				}
								
				if ( !tile_bellow_null && !tile_right_null )
				{
					#pragma omp task depend( in: TILE_A[mmi * NT + ki] ) \
							depend( in: TILE_A[ki * NT + nni] ) \
							depend( inout: TILE_A[mmi * NT + nni] ) \
							shared( TILE_A ) \
							firstprivate( nni, mmi, ki, tile_size_kn ) \
							priority( nt - nni )
					{
					cblas_dgemm( CblasRowMajor, 
							 	( CBLAS_TRANSPOSE ) NoTrans,
							 	( CBLAS_TRANSPOSE ) NoTrans,
							 	tile_size_mm, 
							 	tile_size_nn, 
							 	TILE_SIZE,
							 	-1.0, TILE_A[mmi * NT + ki],  tile_size_kn,
							 	   	  TILE_A[ki * NT + nni],  tile_size_nn,
							  	 1.0, TILE_A[mmi * NT + nni], tile_size_nn );
					}
				}
			}
		}
	}
	} // End parallel

	// --Tile A matrix free--
	//free( &TILE_A );

	return Success;

}
