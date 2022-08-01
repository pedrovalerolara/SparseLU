#include "../include/lass.h"

/**
 *
 * 	@file kdnpgetrf.c
 *
 * 	@brief sLASs-DDSs kdnpgetrf routine.
 *
 * 	sLASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2018-07-24
 * 	@reviewer 
 * 	@modified 
 *
 **/

/**
 *  
 *	@ingroup DDSS
 *   
 *	Performs the LU factorization without pivoting of a general M-by-N matrix A:
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
 *	@param[in]
 *  LDA     int.
 * 	    	LDA specifies the number of columns of A ( row-major order ).
 * 	    	LDA must be at least max( 1, N ).
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
 *	@sa tune_dnpgetrf
 *	@sa smart_dnpgetrf
 *	@sa ddss_tile
 *	@sa ddss_flat2tiled
 *	@sa ddss_tiled2flat
 *
 **/ 

enum LASS_RETURN kdnpgetrf( int M, int N, double *A, int LDA )
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

	/***************************
	--Tile A matrix allocation--
	***************************/
		
	double (*TILE_A)[nt][TILE_SIZE * TILE_SIZE] = malloc ( mt * nt * 
		TILE_SIZE * TILE_SIZE * sizeof( double ) );
	
	if ( TILE_A == NULL )
	{
		fprintf( stderr, "Failure in kdnpgetrf for matrix TILE_A\n" );
		return NoSuccess;
	}

	/*********************************************
	--From flat data layout to tiled data layout--	
	*********************************************/

	ddss_dflat2tiled( M, N, A, LDA, mt, nt, TILE_A );

	/****************
	--DNPGETRF tile--
	****************/	
	
	for ( ki = 0; ki < MIN( mt, nt ); ki++ )
	{
		tile_size_km = ddss_tile_size( M, ki );
		tile_size_kn = ddss_tile_size( N, ki );
		
		#if defined(LASs_WITH_MKL)
		
		tile_size_k = MIN( tile_size_km, tile_size_kn );
		
		#pragma oss task inout( TILE_A[ki][ki] ) \
                shared( TILE_A ) \
                firstprivate( ki ) \
				priority( nt ) \
                label( mkl_dgetrfnpi )
        LAPACKE_mkl_dgetrfnpi( CblasRowMajor,
                               tile_size_km, tile_size_kn,
                               tile_size_k,
                               TILE_A[ki][ki], tile_size_kn );

		#else

		#pragma oss task inout( TILE_A[ki][ki] ) \
				shared( TILE_A ) \
				firstprivate( ki ) \
				priority( nt ) \
				label( dnpgetrf )
		dnpgetrf( tile_size_km, tile_size_kn,
				  TILE_A[ki][ki],
				  tile_size_kn );

		#endif
			
		for ( mi = ki + 1; mi < mt; mi++ )
		{
			tile_size_m = ddss_tile_size( M, mi );
			
			#pragma oss task in( TILE_A[ki][ki] ) \
					inout(TILE_A[mi][ki]) \
					shared( TILE_A ) \
					firstprivate( mi, ki ) \
					priority( nt ) \
					label( dtrsm_below )
			cblas_dtrsm( CblasRowMajor, 
				     	 ( CBLAS_SIDE ) Right, ( CBLAS_UPLO ) Upper, 
				     	 ( CBLAS_TRANSPOSE ) NoTrans, ( CBLAS_DIAG ) NonUnit, 
					 	 tile_size_m, tile_size_kn, 
					 	 1.0, TILE_A[ki][ki], tile_size_kn, 
						  	  TILE_A[mi][ki], tile_size_kn );
		}
		for ( ni = ki + 1; ni < nt; ni++ )
		{
			tile_size_n = ddss_tile_size( N, ni );
			
			#pragma oss task in( TILE_A[ki][ki] ) \
					inout(TILE_A[ki][ni]) \
					shared( TILE_A ) \
					firstprivate( ni, ki ) \
					priority( nt ) \
					label( dtrsm_right )
			cblas_dtrsm( CblasRowMajor, 
				     	 ( CBLAS_SIDE ) Left, ( CBLAS_UPLO ) Lower, 
				     	 ( CBLAS_TRANSPOSE ) NoTrans, ( CBLAS_DIAG ) Unit, 
					 	 tile_size_km, tile_size_n, 
					 	 1.0, TILE_A[ki][ki], tile_size_kn, 
						  	  TILE_A[ki][ni], tile_size_n );
		}
		
		for ( nni = ki + 1; nni < nt; nni++ )
		{
			tile_size_nn = ddss_tile_size( N, nni );

			for ( mmi = ki + 1; mmi < mt; mmi++ )
			{
				tile_size_mm = ddss_tile_size( M, mmi );
								
				#pragma oss task in( TILE_A[mmi][ki] ) \
					in( TILE_A[ki][nni] ) \
					inout( TILE_A[mmi][nni] ) \
					shared( TILE_A ) \
					firstprivate( nni, mmi, ki ) \
					priority( nt - nni ) \
					label( dgemm )
				cblas_dgemm( CblasRowMajor, 
							 ( CBLAS_TRANSPOSE ) NoTrans,
							 ( CBLAS_TRANSPOSE ) NoTrans,
							 tile_size_mm, 
							 tile_size_nn, 
							 TILE_SIZE,
							 -1.0, TILE_A[mmi][ki],  tile_size_kn,
							 	   TILE_A[ki][nni],  tile_size_nn,
							  1.0, TILE_A[mmi][nni], tile_size_nn );
			}
		}
	}

	// --From tile data layout to flat data layout--	
	ddss_dtiled2flat( M, N, A, LDA, mt, nt, TILE_A );

	// --Tile A matrix free--
	free( TILE_A );

	return Success;

}
