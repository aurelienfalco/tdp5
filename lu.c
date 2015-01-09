#include "cblas.h"
#include "mpi.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

//--map-by ppr:1:node -by-node
#define TAG_COL_NUM 1
#define TAG_COL_DATA 2

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
	MPI_Datatype type_column, tmp_type_block, type_block;
	int blocks_per_line = (m/block_size);
	//int disps[nb_blocks];
	//int counts[nb_blocks];
	int side_size = m, myrank, nb_procs;
	int col_size = side_size * block_size;
	int i, j, next_val = 0, inc = 1, nb_cols = 0, *cols_ref;
	double **columns;
	double *left_col;
	double initialTime = 0.0;
	MPI_Status status;
	
	MPI_Type_vector(1, col_size, col_size, MPI_DOUBLE, &type_column);
	MPI_Type_commit(&type_column);
	MPI_Type_vector(block_size, block_size, side_size, MPI_DOUBLE, &tmp_type_block);
	MPI_Type_create_resized(tmp_type_block, 0, block_size * sizeof(double), &type_block);
	MPI_Type_commit(&type_block);

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
				MPI_Send(&i, 1, MPI_INT, next_val, TAG_COL_NUM, MPI_COMM_WORLD);
				MPI_Send(&a[i * col_size], 1, type_column, next_val, TAG_COL_DATA, MPI_COMM_WORLD);// send to the right processus
			}
			serpentinSuivant(&next_val, &inc, nb_procs);
		}
	}
	else {
		for (i = 0; i < nb_cols; i++) {
			MPI_Recv(&cols_ref[i], 1, MPI_INT, 0, TAG_COL_NUM, MPI_COMM_WORLD, &status);
			MPI_Recv(columns[i], 1, type_column, 0, TAG_COL_DATA, MPI_COMM_WORLD, &status);
		}
	}
	
	
	if (myrank == 0)
		initialTime = MPI_Wtime();
	// LU decomposition
	
	next_val = 0; inc = 1; j = 0;
	for (int i = 0; i < MIN(m,n); i+=block_size) {
		int ib = MIN(MIN(m,n)-i,block_size);
		if (myrank == next_val) {
			LAPACKE_dgetf2(order, m-i, ib, columns[j]+i, m, NULL);
			memcpy(left_col, columns[j]+i, col_size * sizeof (double));
			MPI_Bcast(left_col, blocks_per_line - (i / block_size), type_block, next_val, MPI_COMM_WORLD);
			j++;
		}
		else {
			MPI_Bcast(left_col, blocks_per_line - (i / block_size), type_block, next_val, MPI_COMM_WORLD);
		}
		
		for(int k = j; k < nb_cols; k++) {
			cblas_dtrsm(order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, ib, ib, 1, left_col, lda, columns[k]+i, lda );
			cblas_dgemm(order, CblasNoTrans, CblasNoTrans, m-i-ib, ib, ib, -1, left_col + ib, lda, columns[k]+i, lda, 1,  columns[k] + i + ib  , lda);
		}

		serpentinSuivant(&next_val, &inc, nb_procs);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == 0)
		fprintf(stdout, "%lf\n", 1000000*(MPI_Wtime() - initialTime));


	// Gather data on proc 0
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
	else {
		for (i = 0; i < nb_cols; i++) {
			MPI_Send(columns[i], 1, type_column, 0, 99, MPI_COMM_WORLD);
		}
	}
	
	free(cols_ref);
	free(left_col);
	for (i = 0; i < nb_cols; i++) {
		free(columns[i]);
	}	
}
