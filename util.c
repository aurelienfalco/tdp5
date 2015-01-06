#include "util.h"

double* init_matrix(int m, int n)
{
	double* L = calloc(m*n,sizeof(double));
	return L;		
}

double* rand_matrix(int m, int n, int lda)
{
	double* A = init_matrix(m,n);
	int i,j;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			A[i + j * lda] = rand()%10;
		}
		A[i + i * lda] = 10;
	}
	return A;		
}


void print_matrix(double* A, int m, int n, int lda)
{
	int i,j;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			printf("%g ", A[i + j * lda]);
		}
		printf("\n");
	}

}

double* rand_tri_inf(int size)
{
	double* L = init_matrix(size,size);
	int i,j;
	for (i = 0; i < size; i++){
		for (j = 0; j < i; j++){
			L[i + j * size] = rand()%10;
		}
		L[i + i * size] = 1;
	}
	return L;		
}

double* rand_tri_sup(int size)
{
	double* U = init_matrix(size,size);
	int i,j;
	for (i = 0; i < size; i++){
		for (j = i+1; j < size; j++){
			U[i + j * size] = rand()%10;
		}
		U[i + i * size] = 10;
	}
	return U;		
}
