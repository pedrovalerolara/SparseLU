#include "../include/lass.h"

/**
 *
 * 	@file dnpgetrf.c
 *
 * 	@brief LASs-DDSs dnpgetrf routine.
 *
 * 	LASs-DDSs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@author Boro Sofranac boro.sofranac@bsc.es
 * 	@date 2017-21-11
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
 *	@sa kdnpgetrf
 *
 **/

enum LASS_RETURN dnpgetrf( int M, int N, double *A, int LDA )
{

	// Local variable
    double alpha;
    double sfmin;
    int i, j, k;
    int info = 0;

	// Argument checking
	if ( M < 0 )
	{
		fprintf( stderr, "Illegal value of M, in dnpgetrf code\n" );
		return NoSuccess;
	}

	if ( N < 0 )
	{
		fprintf( stderr, "Illegal value of N, in dnpgetrf code\n" );
		return NoSuccess;
	}

	if ( LDA < MAX( 1, N ) )
	{
		fprintf( stderr, "Illegal value of LDA, in dnpgetrf code\n" );
		return NoSuccess;
	}

	// Qick return
	if ( MIN( M, N ) == 0 )
	{
		return Success;
	}
	
	/****************
	--DNPGETRF tile--
	****************/	

	// Minimum value in double precision	
	sfmin = LAPACKE_dlamch_work( 'S' );
	k = MIN( M, N );
	
	for ( i = 0; i < k; i++ ) 
	{
        alpha = A[i * LDA + i];
		if ( alpha != ( double )0.0 ) 
		{
            // Compute elements from J+1 to M of the J-th column
            if ( i < M ) 
			{	
            	if ( fabs( alpha ) > fabs( sfmin ) )
				{
                    alpha = 1.0 / alpha;  
                    cblas_dscal( M - i - 1, alpha, 
						&( A[( i + 1 ) * LDA + i] ), LDA );
                } 
				else 
				{
                    for( j= i + 1; j < M; j++ )
					{	
                        A[j * LDA + i] = A[j * LDA + i] / alpha;
                	}
            	}
        	}
        } 
		else if ( info == 0 ) 
		{
            info = i;   	
        }
		if ( i < k ) 
		{
			cblas_dger( CblasRowMajor,
                        M - i - 1, N - i - 1, 
						-1.0,
                        &A[( i + 1 ) * LDA + i], LDA,
                        &A[ i * LDA + ( i + 1 )], 1,
                        &A[( i + 1 ) * LDA + ( i + 1 ) ], LDA );		
        }
    }
	
	return Success;

} 
