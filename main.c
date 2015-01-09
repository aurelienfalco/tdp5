#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "util.h"
#include "mpi.h"

int main(int argc, char const *argv[])
{
	int myrank, scal = 0;
	int m = 16, n = 16, block_size = 4;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (argc > 1)
		m = n = atoi(argv[1]);
	if (argc > 2)
		scal = atoi(argv[2]);
	
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* a = init_matrix(m,n);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,a,n);
	if (myrank == 0){
		printf("L\n");
		print_matrix(L,m,n,m);
		printf("U\n");
		print_matrix(U,m,n,m);
	}
	// double *a = rand_matrix(m, n);
		
	if (scal){
		fprintf(stderr, "Version scalaire %d\n", m);
		if (myrank == 0)
			LAPACKE_dgetrf(CblasColMajor, m, n, a, m, block_size);
	} else {
		fprintf(stderr, "Version parallèle %d\n", m);
		mpi_cblas_lu(CblasColMajor, m, n, a, m, block_size);
	}
	
	if (myrank == 0)
		print_matrix(a, m, n, m);
	
	free(a);
	MPI_Finalize();
	return 0;
}
