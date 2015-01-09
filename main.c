#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "util.h"
#include "mpi.h"

int main(int argc, char const *argv[])
{
	int myrank, seq = 0, print = 0;
	int m = 128, n = 128, block_size = 64;
	double initialTime = 0.0;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (argc > 1)
		m = n = atoi(argv[1]);
	if (argc > 2)
		seq = atoi(argv[2]);
	if (argc > 3)
		print = atoi(argv[3]);
	
	// Generate L and U
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* A = init_matrix(m,n);
	// Multiply L by U and assigning the result in A
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,A,n);
	if (myrank == 0 && print){
		printf("L\n");
		print_matrix(L,m,n,m);
		printf("U\n");
		print_matrix(U,m,n,m);
	}
	
	if (myrank == 0)
		initialTime = MPI_Wtime();

	// LU Decomposition of A
	if (seq){
		if (myrank == 0){
			if (print)
				fprintf(stderr, "Sequential LU decomposition of size %d\n", m);
			LAPACKE_dgetrf(CblasColMajor, m, n, A, m, block_size);
		}
	} else {
		if (myrank == 0 && print)
			fprintf(stderr, "Parallel LU decomposition of size %d\n", m);
		mpi_cblas_lu(CblasColMajor, m, n, A, m, block_size);
	}

	if (myrank == 0)
		fprintf(stdout, "%lf\n", 1000000*(MPI_Wtime() - initialTime));
	
	if (myrank == 0){
		if (print) {
			printf("A\n");
			print_matrix(A,m,n,m);
		}
		assert(equal_matrix(EqUpper,A,U,m,n,m) && equal_matrix(EqLower,A,L,m,n,m));
	}
	free(A);
	MPI_Finalize();
	return 0;
}
