#include <stdio.h>
#include <stdlib.h>

#include "cblas.h"
#include "util.h"

int main(int argc, char const *argv[])
{
	(void)argc;
	(void)argv;
	
	int m = 6, n = 6, block_size = 2;
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* a = init_matrix(m,n);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,a,n);
	printf("L\n");
	print_matrix(L,m,n,m);
	printf("U\n");
	print_matrix(U,m,n,m);
	// double *a = rand_matrix(m, n);
	
	// a[0] = 1.1;
	// a[1] = 1.1;
	// a[2] = 1.1;
	// a[3] = 1.1;
	// a[4] = 1.1;
	// a[5] = 1.1;
	// a[6] = 2.2;
	// a[7] = 2.2;
	// a[8] = 2.2;
	// a[9] = 2.2;
	// a[10] = 2.2;
	// a[11] = 2.2;
	// a[12] = 3.3;
	// a[13] = 3.3;
	// a[14] = 3.3;
	// a[15] = 3.3;
	// a[16] = 3.3;
	// a[17] = 3.3;
	// a[18] = 4.4;
	// a[19] = 4.4;
	// a[20] = 4.4;
	// a[21] = 4.4;
	// a[22] = 4.4;
	// a[23] = 4.4;
	// a[24] = 5.5;
	// a[25] = 5.5;
	// a[26] = 5.5;
	// a[27] = 5.5;
	// a[28] = 5.5;
	// a[29] = 5.5;
	// a[30] = 6.6;
	// a[31] = 6.6;
	// a[32] = 6.6;
	// a[33] = 6.6;
	// a[34] = 6.6;
	// a[35] = 6.6;

	mpi_cblas_lu(CblasColMajor, m, n, a, m, block_size);
	
	print_matrix(a, m, n, m);
	
	free(a);
	
    return 0;
}
