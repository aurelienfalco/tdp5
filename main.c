#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "util.h"
#include "mpi.h"

int main(int argc, char const *argv[])
{
	int myrank, seq = 0;
	int m = 6, n = 6, block_size = 3;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (argc > 1)
		m = n = atoi(argv[1]);
	if (argc > 2)
		seq = atoi(argv[2]);
	
	// Generate L and U
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* A = init_matrix(m,n);
	// Multiply L by U and assigning the result in A
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,A,n);
	if (myrank == 0){
		printf("L\n");
		print_matrix(L,m,n,m);
		printf("U\n");
		print_matrix(U,m,n,m);
	}

	// LU Decomposition of A
	if (seq){
		if (myrank == 0){
			fprintf(stderr, "Sequential LU decomposition of size %d\n", m);
			LAPACKE_dgetrf(CblasColMajor, m, n, A, m, block_size);
		}
	} else {
		if (myrank == 0)
			fprintf(stderr, "Parallel LU decomposition of size %d\n", m);
		mpi_cblas_lu(CblasColMajor, m, n, A, m, block_size);
	}
	
	if (myrank == 0)
		print_matrix(A, m, n, m);
	
	free(A);
	MPI_Finalize();
	return 0;
}
