#include "../include/lass.h"

void dsss_dtiled2flat( int M, int N, double *A, int LDA, int MT, int NT,
	double **TILE_A[ MT * NT]){
	
	//Local variables
	int m, n;

	for (m = 0; m < MT; m++) {
		for (n = 0; n < NT; n++) {
			dsss_dscatter_tile(M, N,
				&A[m * TILE_SIZE * N + n * TILE_SIZE], LDA, MT, NT,
				*TILE_A, m, n);
		}
	}
}

void dsss_dscatter_tile(int M, int N, double *A, int LDA, int MT, int NT,
	double **TILE_A, int MID, int NID){
	
	//Local variables
	int i, j;
	int tile_size_m, tile_size_n, tile_ind;

	tile_size_m = ddss_tile_size(M, MID);
	tile_size_n = ddss_tile_size(N, NID);

	tile_ind = MID * NT + NID;

	for (i = 0; i < tile_size_m; i++) {
		for (j = 0; j < tile_size_n; j++) {
			if (TILE_A[tile_ind] == NULL) {
				A[i * LDA + j] = 0;
			} else {
				A[i * LDA + j] = TILE_A[tile_ind][i * tile_size_n + j];
			}
		}
	}
}
