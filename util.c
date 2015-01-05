#include "util.h"

double* init_matrix(int m, int n)
{
	double* L = calloc(m*n,sizeof(double));
	return L;		
}

double* rand_matrix(int size)
{
	double* A = init_matrix(size,size);
	int i,j;
	for (i = 0; i < size; i++){
		for (j = 0; j < size; j++){
			A[i + j * size] = rand()%10;
		}
		A[i + i * size] = 10;
	}
	return A;		
}


void print_matrix(double* A, int size)
{
	int i,j;
	for (i = 0; i < size; i++){
		for (j = 0; j < size; j++){
			printf("%g ", A[i + j * size]);
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
