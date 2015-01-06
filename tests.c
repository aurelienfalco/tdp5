#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "util.h"

#define N 3

int test_triangular()
{
	int n = N;
	double* L = rand_tri_inf(n);
	printf("matrix INF\n");
	print_matrix(L,n,n,n);
	double* U = rand_tri_sup(n);
	printf("matrix SUP\n");
	print_matrix(U,n,n,n);
	free(L);
	free(U);
	return 0;
}

int test_dgetf2()
{
	int n = N;
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	printf("matrix L\n");
	print_matrix(L,n,n,n);
	printf("matrix U\n");
	print_matrix(U,n,n,n);
	double* A = init_matrix(n,n);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,A,n);
	printf("matrix A\n");
	print_matrix(A,n,n,n);
	LAPACKE_dgetf2(CblasColMajor,n,n,A,n,NULL);
	printf("matrix A = LU\n");
	print_matrix(A,n,n,n);
	free(A);  
	free(L);  
	free(U);
	return 0;
}

int test_dtrsm()
{
	int m = N;
	int n = N + 2;
	double* U = rand_tri_sup(m);
	double* L = rand_tri_inf(m);
	double* B = rand_matrix(m,n,m);
	printf("matrix L\n");
	print_matrix(L,m,m,m);
	printf("matrix U\n");
	print_matrix(U,m,m,m);
	printf("matrix B\n");
	print_matrix(B,m,n,m);
	cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, m, n, 1, L, m, B, m);
	printf("dtrsm L Y = B\n");
	print_matrix(B,m,n,m);
	cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1, U, m, B, m);
	printf("dtrsm U X = Y\n");
	print_matrix(B,m,n,m);
	// vérifier avec dgemm
	free(U);
	free(L);
	free(B);
	return 0;
}

int test_dgetrf()
{
	int n = N;
	double* A = rand_matrix(n,n,n);
	printf("matrix A\n");
	print_matrix(A,n,n,n);
	LAPACKE_dgetrf(CblasColMajor, n, n, A, n, NULL);
	print_matrix(A,n,n,n);
	free(A);
	return 0;
}

int test_dgesv()
{
	int n = N;
	double* A = rand_matrix(n,n,n);
	double* B = rand_matrix(n,n,n);
	printf("matrix A\n");
	print_matrix(A,n,n,n);
	printf("matrix B\n");
	print_matrix(B,n,n,n);
	LAPACKE_dgesv(CblasColMajor, n, 0, A, n, NULL, B, n);
	printf("dgesv A X = B\n");
	print_matrix(B,n,n,n);
	free(A);
	free(B);
	return 0;
}


int main(int argc, char** argv){
	(void)argc;	(void)argv;
	srand(time(NULL));
	/*
	printf("Test random triangular matrix....\n");
	assert(!test_triangular());
	printf("[OK]\n");

	printf("Test dgetf2....\n");
	assert(!test_dgetf2());
	printf("[OK]\n");
	*/
	printf("Test dtrsm....\n");
	assert(!test_dtrsm());
	printf("[OK]\n");

	printf("Test dgetrf....\n");
	assert(!test_dgetrf());
	printf("[OK]\n");


	printf("Test dgesv....\n");
	assert(!test_dgesv());
	printf("[OK]\n");

	// appeler tests

	return 0;
}
