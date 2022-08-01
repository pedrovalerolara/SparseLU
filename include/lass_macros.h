/**
 *
 * 	@file lass_macros.h
 *
 * 	@brief Macros definition.
 *
 * 	LASs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-01-02
 * 	@reviewer 
 * 	@modified 
 *
 **/

#ifndef DDSS_MACROS_H
#define DDSS_MACROS_H

// VERBOSE (activate prints in testings)
//#define VERBOSE

// Tile size
#define TILE_SIZE 256

// Number of cores
#define NUM_CORES 32

// Return the max/min of two numbers
#define MAX( a, b ) ( ( ( a ) > ( b ) ) ? ( a ) : ( b ) )
#define MIN( a, b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )

// Number of operations 
// GEMM
// Number of multiplications in GEMM
#define FMULS_GEMM(m_, n_, k_) ( (m_) * (n_) * (k_) )
// Number of additions in GEMM
#define FADDS_GEMM(m_, n_, k_) ( (m_) * (n_) * (k_) )
// Flops in DGEMM 
#define FLOPS_DGEMM(m_, n_, k_) ( FMULS_GEMM((double)(m_), (double)(n_), \
	(double)(k_)) + FADDS_GEMM((double)(m_), (double)(n_), (double)(k_)) )

// SYMM
// Number of multiplications in SYMM
#define FMULS_SYMM(side_, m_, n_) ( ( (side_) == Left ) ? FMULS_GEMM((m_), (m_), (n_)) : FMULS_GEMM((m_), (n_), (n_)) )
// Number of additions in SYMM
#define FADDS_SYMM(side_, m_, n_) ( ( (side_) == Left ) ? FADDS_GEMM((m_), (m_), (n_)) : FADDS_GEMM((m_), (n_), (n_)) )
// Flops in DSYMM
#define FLOPS_DSYMM(side_, m_, n_) ( FMULS_SYMM(side_, (double)(m_), (double)(n_)) + FADDS_SYMM(side_, (double)(m_), (double)(n_)) ) 

// TRSM
// Number of multiplications in TRSM
#define FMULS_TRSM_2(m_, n_) ( 0.5 * (n_) * (m_) * ( (m_) + 1 ) )
//Number of additions in TRSM
#define FADDS_TRSM_2(m_, n_) ( 0.5 * (n_) * (m_) * ( (m_) - 1 ) )
// Number of multiplies in TRSM
#define FMULS_TRSM(side_, m_, n_) ( ( (side_) == Left ) ? FMULS_TRSM_2((m_), (n_)) : FMULS_TRSM_2((n_), (m_)) )
// Number of additions in TRSM
#define FADDS_TRSM(side_, m_, n_) ( ( (side_) == Left ) ? FADDS_TRSM_2((m_), (n_)) : FADDS_TRSM_2((n_), (m_)) )
// Flops in DTRSM
#define FLOPS_DTRSM(side_, m_, n_) ( FMULS_TRSM(side_, (double)(m_), (double)(n_)) + FADDS_TRSM(side_, (double)(m_), (double)(n_)) )

// TRMM
// Subfunction: Number of multiplications in TRMM
#define FMULS_TRMM_2(m_, n_) ( 0.5 * (n_) * (m_) * ( (m_) + 1 ) )
// Subfunction: Number of additions in TRMM
#define FADDS_TRMM_2(m_, n_) ( 0.5 * (n_) * (m_) * ( (m_) - 1 ) )
// Number of multiplies in TRMM
#define FMULS_TRMM(side_, m_, n_) ( ( (side_) == Left ) ? FMULS_TRMM_2((m_), (n_)) : FMULS_TRMM_2((n_), (m_)) )
// Number of additions in TRMM
#define FADDS_TRMM(side_, m_, n_) ( ( (side_) == Left ) ? FADDS_TRMM_2((m_), (n_)) : FADDS_TRMM_2((n_), (m_)) )
// Flops in DTRMM
#define FLOPS_DTRMM(side_, m_, n_) (FMULS_TRMM(side_, (double)(m_), (double)(n_)) + FADDS_TRMM(side_, (double)(m_), (double)(n_)) )

// SYRK
// Number of multiplications in SYRK
#define FMULS_SYRK(k_, n_) ( 0.5 * (k_) * (n_) * ( (n_) + 1 ) )
// Number of additions in SYRK
#define FADDS_SYRK(k_, n_) ( 0.5 * (k_) * (n_) * ( (n_)+1 ) )
// Flops in DSYRK
#define FLOPS_DSYRK(k_, n_) ( FMULS_SYRK( (double)(k_), (double)(n_) ) + FADDS_SYRK( (double)(k_), (double)(n_) ) )

//SYR2K
// Number of multiplications in SYR2K
#define FMULS_SYR2K(k_, n_) ((k_) * (n_) * (n_))
// Number of additions in SYR2K
#define FADDS_SYR2K(k_, n_) ((k_) * (n_) * (n_) + (n_))
// Flops in DSYR2K
#define FLOPS_DSYR2K(k_, n_) (FMULS_SYR2K((double)(k_), (double)(n_)) + FADDS_SYR2K((double)(k_), (double)(n_)))

// POTRF
// Number of multiplications in POTRF
#define FMULS_POTRF(n_) ( (1./6.) * (n_) * (n_) * (n_) + (0.5) * (n_) * (n_) \
	+ (1./3.) * n_ )
// Number of additions in POTRF
#define FADDS_POTRF(n_) ( (1./6.) * (n_) * (n_) * (n_) - (1./6.) * (n_))
// Flops in DPOTRF
#define FLOPS_DPOTRF(n_) ( FMULS_POTRF((double)(n_)) + FADDS_POTRF((double)(n_)))

// GETRF
// Number of multiplications in GETRF
#define FMULS_GETRF( m_, n_) ( ( (m_) >= (n_) ) ?  0.5 * (m_) * (n_) * (n_) - 1./6. * (n_) * (n_) * (n_) + 0.5 * (m_) * (n_) - 0.5 * (n_) * (n_) + 2./3. * (n_) : 0.5 * (n_) * (m_) * (m_) - 1./6. * (m_) * (m_) * (m_) + 0.5 * (n_) * (m_) - 0.5 * (m_) * (m_) + 2./3. * (m_) ) 
// Number of additions in GETRF
#define FADDS_GETRF( m_, n_) ( ( (m_) >= (n_) ) ? 0.5 * (m_) * (n_) * (n_) - 1./6. * (n_) * (n_) * (n_) - 0.5 * (m_) * (n_) + 1./6. * (n_) : 0.5 * (n_) * (m_) * (m_) - 1./6. * (m_) * (m_) * (m_) - 0.5 * (n_) * (m_) + 1./6. * (m_) )
// Flops in DGETRF
#define FLOPS_DGETRF(m_, n_) ( FMULS_GETRF((double)(m_), (double)(n_)) + FADDS_GETRF((double)(m_), (double)(n_)) )

#endif

#ifndef DSSS_MACROS_H
#define DSSS_MACROS_H

// Number of operations
// GTSV
// Flops in DGTSV
#define FLOPS_DGTSV( n_ ) ( 8 * (n_) )

// SPMV
// Flops in DSPMV
#define FLOPS_DSPMV( nnz_ ) ( 2 * (nnz_) )

#endif

