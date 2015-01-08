#include "cblas.h"
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// TODO: assigner de vrai tags aux messages
//--map-by ppr:1:node -by-node

void serpentinSuivant(int *value, int *inc, const int max) {
	if (*inc) {
		if ((*value) >= max-1)
			*inc = 0;
		else
			(*value)++;
	}
	else {
		if (*value == 0)
			*inc = 1;
		else {
			(*value)--;
		}
	}
}

void mpi_cblas_lu(const enum CBLAS_ORDER order, const int m, const int n, double* a, int lda, int block_size)
{
	MPI_Datatype type_column;
	int blocks_per_line = (m/block_size);
	//int disps[nb_blocks];
	//int counts[nb_blocks];
	int side_size = m, myrank, nb_procs;
	int col_size = side_size * block_size;
	int i, j, next_val = 0, inc = 1, nb_cols = 0, *cols_ref;
	double **columns;
	double *left_col;
	MPI_Status status;
	
	MPI_Init(NULL, NULL);
	
	MPI_Type_vector(1, col_size, col_size, MPI_DOUBLE, &type_column);
	MPI_Type_commit(&type_column);
	/*MPI_Type_vector(block_size, block_size, side_size, MPI_DOUBLE, &block_matrix);
	MPI_Type_create_resized(block_matrix, 0, sizeof(double), &type_matrix);
	MPI_Type_commit(&type_matrix);
	
	for (int i = 0; i < side_size; i++) {
		for (int j = 0; j < side_size; j++) {
			printf("i=%d\n", i);
			disps[i * side_size + j] = i * blocks_by_line * block_size + j * block_size;
			counts [i * side_size + j] = 1;
		}
	}*/
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
	
	// Count the number of columns owned
	next_val = 0; inc = 1;
	for (i = 0; i < blocks_per_line; i++) {
		if (next_val == myrank)
			nb_cols++;
		serpentinSuivant(&next_val, &inc, nb_procs);
	}
	cols_ref = malloc(nb_cols * sizeof(double));
	left_col = malloc(col_size * sizeof(double));
	columns = malloc(nb_cols * sizeof(double*));
	for (i = 0; i < nb_cols; i++) {
		columns[i] = malloc(col_size * sizeof(double));
	}
	
	// ** Séparer la matrice avec le serpentin et le type_column par le proc 0, les autres reçoivent les données envoyées
	if (myrank == 0) {
		next_val = 0; inc = 1;
		j = 0;
		for (i = 0; i < blocks_per_line; i++) {
			if (next_val == myrank) {
				memcpy(columns[j], &a[i * col_size], col_size * sizeof(double));
				cols_ref[j++] = i;
			}
			else {
				MPI_Send(&i, 1, MPI_INT, next_val, 99, MPI_COMM_WORLD);
				MPI_Send(&a[i * col_size], 1, type_column, next_val, 99, MPI_COMM_WORLD);// send to the right processus
			}
			serpentinSuivant(&next_val, &inc, nb_procs);
		}
	}
	else {
		for (i = 0; i < nb_cols; i++) {
			MPI_Recv(&cols_ref[i], 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
			MPI_Recv(columns[i], 1, type_column, 0, 99, MPI_COMM_WORLD, &status);
		}
	}
	
	/*next_val = 0; inc = 1;
	for (i = 0; i < blocks_by_line; i++) {
		if (myrank == next_val) {
			// Facto lu ici
			MPI_Send();
		}
		else {
			MPI_Recv();
		}
		serpentinSuivant(&next_val, &inc, nb_procs);
	}*/
	// -- Envoyer la première colonne pour entamer les autres calculs
	//MPI_Bcast(buf, 1, type_column, 0, MPI_COMM_WORLD);
	
	// ++ Chaque processus fait un trsm sur le premier bloc de sa colonne (il faut attendre le bloc de la première colonne)
	// ++ puis fait un gemm sur les autres (il faut aussi attendre le bloc de la première colonne)
	
	if (myrank == 0) {
		next_val = 0; inc = 1; j = 0;
		for (i = 0; i < blocks_per_line; i++) {
			if (next_val == myrank) {
				memcpy(&a[i * col_size], columns[j++], col_size * sizeof(double));
			}
			else {
				MPI_Recv(&a[i * col_size], 1, type_column, next_val, 99, MPI_COMM_WORLD, &status);
			}
			serpentinSuivant(&next_val, &inc, nb_procs);
		}
	}
	else{
		for (i = 0; i < nb_cols; i++) {
			MPI_Send(columns[i], 1, type_column, 0, 99, MPI_COMM_WORLD);
		}
	}
	
	free(cols_ref);
	free(left_col);
	for (i = 0; i < nb_cols; i++) {
		free(columns[i]);
	}
	
	MPI_Finalize();
}
#include "util.h"
#include <assert.h>

void cblas_lu(const enum CBLAS_ORDER order, const int m, const int n, double* a, int lda, int nb)
{
	assert(order == CblasColMajor);
	for (int j = 0; j < MIN(m,n); j+=nb) {
		int jb = MIN(MIN(m,n)-j,nb);
		LAPACKE_dgetf2(order, m-j, jb, &a[j * lda + j], lda, NULL);
		for (int i = j+jb; i < MIN(m,n); i+=jb)
		{
			int ib = MIN(MIN(m,n)-i,jb);
			cblas_dtrsm(order, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit, ib, ib, 1, &a[j + lda * j], lda, &a[j + lda * i], lda );
			cblas_dgemm(order, CblasNoTrans, CblasNoTrans, m-j-jb, ib, ib, -1, &a[j + jb + j * lda], lda,&a[j +  i * lda], lda, 1, &a[j + jb + lda * i], lda);
		}
	}
}
