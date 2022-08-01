/**
 *
 * 	@file lass.h
 *
 * 	@brief LASs header definition.
 *
 * 	LASs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-09-07
 * 	@reviewer 
 * 	@modified 
 *
 **/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <assert.h> 
#include <float.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#if defined(LASs_WITH_MKL)
#include <mkl.h>
#include <mkl_lapacke.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif
#include "lass_macros.h" 
#include "lass_types.h"

#ifndef DDSS_H
#define DDSS_H

// LASs-DDSs routines

// BLAS-3 routines

int ddss_dgemm( enum DDSS_TRANS TRANS_A, enum DDSS_TRANS TRANS_B,
	    int M, int N, int K,
        double ALPHA, double *A, int LDA,
                      double *B, int LDB,
        double BETA,  double *C, int LDC );

int ddss_dsymm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		int M, int N,
		double ALPHA, double *A, int LDA,
					  double *B, int LDB,
		double BETA,  double *C, int LDC ); 

int ddss_dtrsm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		enum DDSS_TRANS TRANS_A, enum DDSS_DIAG DIAG, 
		int M, int N,
        const double ALPHA, double *A, int LDA,
                            double *B, int LDB );

int ddss_dtrmm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO,
		enum DDSS_TRANS TRANS_A, enum DDSS_DIAG DIAG, 
		int M, int N,
		const double ALPHA, double *A, int LDA, 
							double *B, int LDB ); 

int ddss_dsyrk( enum DDSS_UPLO UPLO, enum DDSS_TRANS TRANS_A,
        int N, int K,
        const double ALPHA, double *A, int LDA,
        const double BETA,  double *C, int LDC );

int ddss_dsyr2k( enum DDSS_UPLO UPLO, enum DDSS_TRANS TRANS, 
		int N, int K,
        const double ALPHA, double *A, int LDA,
						    double *B, int LDB,
    	const double BETA,  double *C, int LDC );

// LAPACK routines

int ddss_dnpgetrf( int M, 
		int N, 
		double *A, int LDA );

int ddss_dnpgesv( int N, int NRHS, 
		double *A, int LDA, 
		double *B, int LDB );

int ddss_dtpgetrf( int M, 
		int N, 
		double *A, int LDA, 
		int *ipiv );

int ddss_dtpgesv( int N, int NRHS, 
		double *A, int LDA, 
		int *IPIV, 
		double *B, int LDB );

int ddss_dpotrf( enum DDSS_UPLO UPLO, 
		int N, 
		double *A, int LDA );

int ddss_dposv( enum DDSS_UPLO UPLO,
        int N, int NRHS,
        double *A, int LDA,
        double *B, int LDB );

// BLAS-3 kernels

enum LASS_RETURN kdgemm( enum DDSS_TRANS TRANS_A, enum DDSS_TRANS TRANS_B,
		int M, int N, int K,	
		const double ALPHA, double *A, int LDA,
							double *B, int LDB,
		const double BETA,  double *C, int LDC );

enum LASS_RETURN kdsymm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO, 
		int M, int N,
		const double ALPHA, double *A, int LDA,
					  		double *B, int LDB,
		const double BETA,  double *C, int LDC );

enum LASS_RETURN kdtrsm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO,
		enum DDSS_TRANS TRANS_A, enum DDSS_DIAG DIAG,
		int M, int N,
		const double ALPHA, double *A, int LDA,
						    double *B, int LDB );

enum LASS_RETURN kdtrmm( enum DDSS_SIDE SIDE, enum DDSS_UPLO UPLO,
		enum DDSS_TRANS TRANS_A, enum DDSS_DIAG DIAG,
		int M, int N,
		const double ALPHA, double *A, int LDA,
						    double *B, int LDB );

enum LASS_RETURN kdsyrk( enum DDSS_UPLO UPLO, enum DDSS_TRANS TRANS_A,
        int N, int K,
        const double ALPHA, double* A, int LDA,
        const double BETA, double *C, int LDC );

enum LASS_RETURN kdsyr2k( enum DDSS_UPLO UPLO, enum DDSS_TRANS TRANS,
		int N, int K,
		const double ALPHA, double *A, int LDA,
							double *B, int LDB,
		const double BETA,  double *C, int LDC ); 

// LAPACK kernels

enum LASS_RETURN kdnpgetrf( int M, 
		int N, 
		double *A, int LDA );

enum LASS_RETURN kdnpgesv( int N, int NRHS, 
		double *A, int LDA, 
		double *B, int LDB );

enum LASS_RETURN kdtpgetrf( int M, 
		int N, 
		double *A, int LDA, 
		int *ipiv );

enum LASS_RETURN kdtpgesv( int N, int NRHS, 
		double *A, int LDA, 
		int *IPIV, 
		double *B, int LDB );

enum LASS_RETURN kdpotrf( enum DDSS_UPLO UPLO, 
		int N, 
		double *A, int LDA );

enum LASS_RETURN kdposv( enum DDSS_UPLO UPLO,
        int N, int NRHS,
        double *A, int LDA,
        double *B, int LDB );

// LAPACK codes

enum LASS_RETURN dnpgetrf( int M, int N, double *A, int LDA );

// ddss_tile.c routines

int ddss_tile_size( int M, int MT );

// ddss_flat2tiled.c routines

void ddss_dflat2tiled( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE] );

void ddss_dsymtiled2flat( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE], enum DDSS_UPLO UPLO );

void ddss_dgather_tile( int M, int N, double *A, int LDA, 
		double  *TILE_A, int MID, int NID );

// ddss_tiled2flat.c routines

void ddss_dtiled2flat( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE] );

void ddss_dtiled2flat_nb( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE] );

void ddss_dsymflat2tiled( int M, int N, double *A, int LDA, int MT, int NT, 
		double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE], enum DDSS_UPLO UPLO );

void ddss_dsymtiled2flat_nb( int M, int N, double *A, int LDA, int MT, int NT,
        double (*TILE_A)[NT][TILE_SIZE * TILE_SIZE], enum DDSS_UPLO UPLO );

void ddss_dscatter_tile( int M, int N, double *A, int LDA, 
		double *TILE_A, int MID, int NID );

// dsss_flat2tiled.c routines

void dsss_dflat2tiled(int M, int N, double *A, int LDA, int MT, int NT,
		double **TILE_A[ MT * NT ]);

void dsss_dgather_tile(int M, int N, double *A, int LDA, int MT, int NT,
		double **TILE_A, int MID, int NID);

// dsss_tiled2flat.c routines

void dsss_dtiled2flat( int M, int N, double *A, int LDA, int MT, int NT,
	double **TILE_A[ MT * NT]);

void dsss_dscatter_tile(int M, int N, double *A, int LDA, int MT, int NT,
	double **TILE_A, int MID, int NID);

#endif

#ifndef DSSS_H
#define DSSS_H

// LASs-DSSs routines

// DGTSV 

int dsss_dgtsv( int N, 
        double *DL, double *D, double *DU, 
        double *RHS ); 

// DSPMV 

int dsss_dspmv( int M, int N, 
		double ALPHA, 
        const double *VAL_A, const int *ROW_PTR_A, const int *COL_IND_A, 
		const double *X, 
        	double BETA, 
			  double *Y );

// DNPGETRF

int dsss_dnpgetrf( int M, int N, 
	int MT, int NT, 
	double *TILE_A[ MT * NT ] ); 

// DGTSV kernel

enum LASS_RETURN kdgtsv( int N, 
        double *DL, double *D, double *DU, 
        double *RHS );

// DSPMV kernel

enum LASS_RETURN kdspmv ( int M, int N, 
		double ALPHA, 
		const double *VAL_A, const int *ROW_PTR_A, const int *COL_IND_A, 
			const double *X, 
				double BETA, 
				  double *Y ); 

// DGTSV code

void dgtsvseq( int N, 
    	double *DL, double *D, double *DU, 
       double *RHS );

// DSPMV code

void dspmvseq( int M, int N, 
		double ALPHA, 
	 	const double *VAL_A, const int *ROW_PTR_A, const int *COL_IND_A, 
		const double *X, 
		    double BETA, 
              double *Y ); 

#endif
