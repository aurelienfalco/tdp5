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

int test_copy()
{
	int n = N;
	int m = N + rand()%5;
	double* A = rand_matrix(m,n);
	double* B = copy_matrix(A, m, n, m);
	int res = equal_matrix(0,A,B,m,n,m);
	free(A);
	free(B);
	return !res;
}

int test_dgetf2()
{
	int n = N;
	int m = n;
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* A = init_matrix(n,n);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,A,n);
	LAPACKE_dgetf2(CblasColMajor,n,n,A,n,NULL);
	int res = !equal_matrix(1,A,U,m,n,m) || !equal_matrix(2,A,L,m,n,m);
	free(A);  
	free(L);  
	free(U);
	return res;
}

int test_dtrsm()
{
	int m = N;
	int n = N + 2;
	double* U = rand_tri_sup(m);
	double* L = rand_tri_inf(m);
	double* A = init_matrix(m,n);
	double* B = rand_matrix(m,n);
	double* C = copy_matrix(B,m,n,m);
	double* E = init_matrix(m,n);
	cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, m, n, 1, L, m, B, m);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m,n,m,1,L,m,B,m,1,A,m);
	int res = !equal_matrix(0,A,C,m,n,m);
	
	double* D = copy_matrix(B,m,n,m);
	cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1, U, m, B, m);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m,n,m,1,U,m,B,m,1,E,m);
	res |= !equal_matrix(0,E,D,m,n,m);
	free(U);
	free(L);
	free(B);
	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	return res;
}

int test_dgetrf()
{
	int n = N;
	int m = N;
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* A = init_matrix(n,n);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,A,n);
	LAPACKE_dgetrf(CblasColMajor, n, n, A, n, NULL);
	int res = !equal_matrix(1,A,U,m,n,m) || !equal_matrix(2,A,L,m,n,m);
	free(L);
	free(U);
	free(A);
	return res;
}

int test_dgesv()
{
	int m = N;
	int n = N;
	double* A = rand_matrix(m,n);
	double* B = rand_matrix(m,n);
	double* C = copy_matrix(A,m,n,m);
	double* D = copy_matrix(B,m,n,m);
	double* E = init_matrix(m,n);
	LAPACKE_dgesv(CblasColMajor, m, 0, A, m, NULL, B, m);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m,n,m,1,C,m,B,m,1,E,m);
	int res = !equal_matrix(0,D,E,m,n,m);
	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	return res;
}

int test_lu_block()
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
	cblas_lu(CblasColMajor, n, n, A, n, 1);
	printf("lu block A = LU\n");
	print_matrix(A,n,n,n);
	// free(L);
	// free(U);
	// free(A);
	return 0;
}


int main(int argc, char** argv){
	(void)argc;	(void)argv;
	srand(time(NULL));
	/*
	printf("Test random triangular matrix....\n");
	assert(!test_triangular());
	printf("[OK]\n");
	*/
	printf("Test dgetf2....\n");
	assert(!test_dgetf2());
	printf("[OK]\n");

	printf("Test copy....\n");
	assert(!test_copy());
	printf("[OK]\n");
	
	printf("Test dtrsm....\n");
	assert(!test_dtrsm());
	printf("[OK]\n");

	printf("Test dgetrf....\n");
	assert(!test_dgetrf());
	printf("[OK]\n");

	printf("Test dgesv....\n");
	assert(!test_dgesv());
	printf("[OK]\n");

	printf("Test lu block....\n");
	assert(!test_lu_block());
	printf("[OK]\n");

	// appeler tests

	return 0;
}
