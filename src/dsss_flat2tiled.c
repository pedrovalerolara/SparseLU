#include "../include/lass.h"

/**
*	@file dsss_flat2tiled.c
*
*	@brief LASs-DSSs dsss_flat2tiled routines.
*	
*	LASs-DSSs is a software package provided by:
*	(?)Oak Ridge National Laboratory
*	
*	@author Cameron Greenwalt camerongreenwalt@gmail.com
*	@date dd/mm/yyyy 14/06/2021
*	@reviewer 
*	@modified 
**/

/**
*	@ingroup DDSs
*	 
*	dsss_dflat2tiled:
*	Performs the change of the data layout from flat layout to tiled layout 
*	according to row-major order.
**/

/**
*	@param[in]
*	M		int.	
*	        M specifies the number of rows of the flat matrix.
*	
*	@param[in]
*	N		int.	
*	        N specifies the number of columns of the flat matrix.
*	
*	@param[in]
*	A		double *.	
*	        A is a pointer to the flat matrix.
*	
*	@param[in]
*	LDA		int.
*			LDA specifies the number of columns ( row-major order ) of matrix A.
*	
*	@param[in]
*	MT		int.	
*	        MT specifies the number of rows of the matrix TILE_A.
*			
*	@param[in]
*	NT		int.
*			NT specifies the number of columns of the matrix TILE_A.
*	
*	@param[in, out]
*	TILE_A	double *.
*			TILE_A is a pointer to the tile matrix.
*
**/
	
void dsss_dflat2tiled(int M, int N, double *A, int LDA, int MT, int NT,
		double **TILE_A[ MT * NT ]){

	// Local variables
	int m, n; // indexes to tile  row/col position in matrix


	for (m = 0; m < MT; m++) {
		for (n = 0; n < NT; n++) {
			dsss_dgather_tile( M, N, &A[m * TILE_SIZE * N + n * TILE_SIZE],
				 LDA, MT, NT, *TILE_A, m, n );
		}
	}
}

/**
*	@ingroup DSSs
*	 
*	dsss_dgather_tile:
*	Performs the copy of a tile from the flat matrix A to the tile matrix TILE_A
*	for the MT, NT tile.
*	
**/
	
	
/**	
*	@param[in]
*	M		int.	
*	        M specifies the number of rows of the flat matrix.
*	
*	@param[in]
*	N		int.	
*	        N specifies the number of columns of the flat matrix.
*	
*	@param[in]
*	A		double *.	
*	        A is a pointer to the flat matrix.
*			
*	@param[in]
*	LDA		int.
*			LDA specifies the number of columns ( row-major order ) of matrix A.
*	
*	@param[in, out]
*	TILE_A	double *.
*			TILE_A is a pointer to the tile matrix.
*	
*	@param[in]
*	MID		int.
*			MID specifies the row id of the tile.
*	
*	@param[in]
*	NID		int.
*			NID specifies the column id of the tile. 
**/	
	
	
/**	
*	@sa dsss_tile_size
*	@sa dsss_dflat2tiled
**/

void dsss_dgather_tile(int M, int N, double *A, int LDA, int MT, int NT,
		double **TILE_A, int MID, int NID){

	//Local variables
	int i, j, tile_ind;
	int tile_size_m, tile_size_n;
	int is_empty;

	tile_size_m = ddss_tile_size(M, MID);
	tile_size_n = ddss_tile_size(N, NID);

	tile_ind = MID * NT + NID; // Index of the tile pointer in array

	is_empty = 1;
	for (i = 0; i < tile_size_m; i++) {
		for (j = 0; j < tile_size_n; j++) {
			// The tile isn't empty
			if (A[i * LDA + j] != 0.0) {
				is_empty = 0;
				break;
			}
		}	 
	}

	if (!is_empty) {
		// Allocate memory for the tile
		TILE_A[tile_ind] = malloc(tile_size_m * tile_size_n * sizeof(double));
		
		// Assign values to tile
		for (i = 0; i < tile_size_m; i++) {
			for (j = 0; j < tile_size_n; j++) {
				TILE_A[tile_ind][i * tile_size_n + j] = A[i * LDA + j];
			}	 
		}
	} else {
		// Set tile at that index to null
		TILE_A[tile_ind] = NULL;
	}
	
	/*
	//Local Variables
	int i, j, tile_ind;
	int tile_size_m, tile_size_n;
	int is_empty;

	tile_size_m = ddss_tile_size(M, MID);
	tile_size_n = ddss_tile_size(N, MID);

	tile_ind = MID * NT + NID; // Index of the tile pointer in array
	
	// Allocate memory for the tile
	TILE_A[tile_ind] = malloc(tile_size_m * tile_size_n * sizeof(double));

	// Assume the tile is empty from the beginning
	is_empty = 1;

	for (i = 0; i < tile_size_m; i++) {
			for (j = 0; j < tile_size_n; j++) {
				// The tile isn't empty
			if (A[i * LDA + j] != 0.0) is_empty = 0;
			
			TILE_A[tile_ind][i * tile_size_n + j] = A[i * LDA + j];
		}	 
	}
	
	// The tile is empty. Deallocate memory and set pointer to null.
	if (is_empty) {
			free(TILE_A[tile_ind]);
		TILE_A[tile_ind] = NULL;
	}
	*/
}

	

