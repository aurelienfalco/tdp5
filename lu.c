#include "cblas.h"
#include "mpi.h"
#include "stdlib.h"

// TODO: assigner de vrai tags aux messages

void cblas_lu(const enum CBLAS_ORDER order, const int M, const int N,
			  const double alpha, const double *X, const int incX,
			  const double *Y, const int incY, double *A, const int lda, const int block_size)
{
	MPI_Datatype block_matrix, type_matrix, type_column;
	int blocks_by_line = (M/block_size);
	int nb_blocks = blocks_by_line * blocks_by_line;
	int disps[nb_blocks];
	int counts[nb_blocks];
	int side_size = M, myrank, nb_procs;
	int col_size = side_size * block_size;
	int i, j, next_val = 0, inc = 1, nb_cols = 0;
	double **columns;
	double *cols_ref;
	
	MPI_Type_vector(block_size * side_size, block_size * side_size, block_size * side_size, MPI_DOUBLE, &type_column);
	MPI_Type_commit(&type_column);
	MPI_Type_vector(block_size, block_size, side_size, MPI_DOUBLE, &block_matrix);
	MPI_Type_create_resized(block_matrix, 0, sizeof(double), &type_matrix);
	MPI_Type_commit(&type_matrix);
	
	for (i=0; i<side_size; i++) {
		for (j=0; j<side_size; j++) {
			disps[i*side_size+j] = i*blocks_by_line*block_size+j*block_size;
			counts [i*side_size+j] = 1;
		}
	}
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
	
	// Count the number of columns owned
	next_val = 0;
	for (i = 0; i < blocks_by_line; i++) {
		if (next_val == myrank)
			nb_cols++;
		if (inc) {
			next_val >= nb_procs-1? inc = 0: next_val++;
		}
		else {
			next_val==0? inc = 1: next_val--;
		}
	}
	cols_ref = malloc(nb_cols * sizeof(double));
	columns = malloc(nb_cols * sizeof(double*));
	for (i = 0; i < nb_cols; i++) {
		columns[i] = malloc(col_size * sizeof(double));
	}
	
	// ** Séparer la matrice avec le serpentin et le type_column
	next_val = 0;
	j = 0;
	for (i = 0; i < blocks_by_line; i++) {
		if (next_val == myrank) {
			columns[j] = NULL; // To be completed
			j++; //add column to memory
		}
		else
			MPI_Send(NULL, 1, type_column, next_val, 99, MPI_COMM_WORLD);// send to the right processus
		if (inc) {
			next_val >= nb_procs-1? inc = 0: next_val++;
		}
		else {
			next_val==0? inc = 1: next_val--;
		}
	}
	
	// -- Envoyer la première colonne pour entamer les autres calculs
	//MPI_Bcast(buf, 1, type_column, 0, MPI_COMM_WORLD);
	
	// ++ Chaque processus fait un trsm sur le premier bloc de sa colonne (il faut attendre le bloc de la première colonne)
	// ++ puis fait un gemm sur les autres (il faut aussi attendre le bloc de la première colonne)
	
	
	
	for (i = 0; i < nb_cols; i++) {
		free(columns[i]);
	}
}

#include <assert.h>
/*
void cblas_lu(const enum CBLAS_ORDER order, const int M, const int N,
	const double alpha, const double *X, const int incX,
	const double *Y, const int incY, double *A, const int lda)
{
	assert(order == CblasColMajor);

	for (int i = 0; i < M; i++) {
		// dscal(i);
		// dger(i);
	}
}
*/
