#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "util.h"

unsigned int N = 3;
extern float epsilon;
int print = 0;

int test_triangular()
{
	printf("Test random triangular matrix....\n");
	int n = N;
	double* L = rand_tri_inf(n);
	printf("matrix INF\n");
	print_matrix(L,n,n,n);
	double* U = rand_tri_sup(n);
	printf("matrix SUP\n");
	print_matrix(U,n,n,n);
	free(L);
	free(U);
	assert(1);
	printf("[OK]\n");
	return 0;
}

int test_copy()
{
	printf("Test copy....\n");
	int m = N + rand()%5;
	int n = N;
	double* A = rand_matrix(m,n);
	double* B = copy_matrix(A, m, n, m);
	int copy = equal_matrix(0,A,B,m,n,m);
	if (print){
		printf("A\n");
		print_matrix(A,m,n,m);
		printf("B\n");
		print_matrix(B,m,n,m);
	}
	free(A);
	free(B);
	assert(copy);
	printf("[OK]\n");
	return 0;
}

int test_dgemm()
{
	printf("Test dgemm....\n");
	int m = N;
	int n = N+1;
	int l = N+2;
	double* A = rand_matrix(m,n);
	double* B = rand_matrix(n,l);
	double* C = init_matrix(m,l);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m,l,n,1,A,m,B,n,1,C,m);
	printf("A\n");
	print_matrix(A,m,n,m);
	printf("B\n");
	print_matrix(B,n,l,n);
	printf("C\n");
	print_matrix(C,m,l,m);
	free(A);  
	free(B);  
	free(C);
	printf("[OK]\n");
	return 0;
}

int test_dgetf2()
{
	printf("Test dgetf2....\n");
	int n = N;
	int m = n;
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* A = init_matrix(n,n);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,A,n);
	LAPACKE_dgetf2(CblasColMajor,n,n,A,n,NULL);
	int equal = equal_matrix(1,A,U,m,n,m) && equal_matrix(2,A,L,m,n,m);
	if (print){
		printf("L\n");
		print_matrix(L,m,n,m);
		printf("U\n");
		print_matrix(U,m,n,m);
		printf("A\n");
		print_matrix(A,m,n,m);
	}
	free(A);
	free(L);  
	free(U);
	assert(equal);
	printf("[OK]\n");
	return 0;
}

int test_dtrsm()
{
	printf("Test dtrsm....\n");
	int m = N;
	int n = N + 2;
	double* U = rand_tri_sup(m);
	if (print){
		printf("U\n");
		print_matrix(U,m,m,m);
	}
	double* L = rand_tri_inf(m);
	if (print){
		printf("L\n");
		print_matrix(L,m,m,m);
	}
	double* A = init_matrix(m,n);
	double* B = rand_matrix(m,n);
	double* C = copy_matrix(B,m,n,m);
	double* E = init_matrix(m,n);
	cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, m, n, 1, L, m, B, m);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m,n,m,1,L,m,B,m,1,A,m);
	int equal = equal_matrix(0,A,C,m,n,m);
	if (print){
		printf("B\n");
		print_matrix(C,m,n,m);
		printf("X in dtrsm L X = B \n");
		print_matrix(B,m,n,m);
		printf("B' (== B?)\n");
		print_matrix(A,m,n,m);
	}
	
	double* D = copy_matrix(B,m,n,m);
	cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1, U, m, B, m);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m,n,m,1,U,m,B,m,1,E,m);
	equal &= equal_matrix(0,E,D,m,n,m);
	if (print){
		printf("Y in dtrsm U Y = X \n");
		print_matrix(B,m,n,m);
		printf("B' (== X ?)\n");
		print_matrix(E,m,n,m);
	}
	free(U);
	free(L);
	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	assert(equal);
	printf("[OK]\n");
	return 0;
}

int test_dgetrf()
{
	printf("Test dgetrf....\n");
	int n = N;
	int m = N;
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* A = init_matrix(n,n);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,A,n);
	LAPACKE_dgetrf(CblasColMajor, n, n, A, n, NULL);
	int equal = equal_matrix(1,A,U,m,n,m) && equal_matrix(2,A,L,m,n,m);
	if (print){
		printf("L\n");
		print_matrix(L,m,n,m);
		printf("U\n");
		print_matrix(U,m,n,m);
		printf("A\n");
		print_matrix(A,m,n,m);
	}
	free(L);
	free(U);
	free(A);
	assert(equal);
	printf("[OK]\n");
	return 0;
}

int test_dgesv()
{
	printf("Test dgesv....\n");
	int m = N;
	int n = N;
	double* A = rand_matrix(m,n);
	double* B = rand_matrix(m,n);
	double* C = copy_matrix(A,m,n,m);
	double* D = copy_matrix(B,m,n,m);
	double* E = init_matrix(m,n);
	LAPACKE_dgesv(CblasColMajor, m, 0, A, m, NULL, B, m);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m,n,m,1,C,m,B,m,1,E,m);
	int equal = equal_matrix(0,D,E,m,n,m);
	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	assert(equal);
	printf("[OK]\n");
	return 0;
}

int test_lu_block()
{
	printf("Test lu block....\n");
	int m = N;
	int n = N;
	double* L = rand_tri_inf(n);
	double* U = rand_tri_sup(n);
	double* A = init_matrix(m,n);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,L,n,U,n,1,A,n);
	if (print){
		printf("A\n");
		print_matrix(A,m,n,m);
	}	

	cblas_lu(CblasColMajor, n, n, A, n, 3);
	int equal = equal_matrix(1,A,U,m,n,m) && equal_matrix(2,A,L,m,n,m);
	if (print){
		printf("L\n");
		print_matrix(L,m,n,m);
		printf("U\n");
		print_matrix(U,m,n,m);
		printf("A\n");
		print_matrix(A,m,n,m);
	}	
	free(L);
	free(U);
	free(A);
	assert(equal);
	printf("[OK]\n");
	return 0;
}


int main(int argc, char** argv){
	(void)argc;	(void)argv;
	srand(time(NULL));
	if (argc > 1){
		N = atoi(argv[1]);
	}
	if (argc > 2){
		epsilon = atof(argv[2]);
	}
	if (argc > 3){
		print = atof(argv[3]);
	}

	// test_triangular();
	// test_dgemm();
	test_dgetf2();
	test_copy();
	test_dtrsm();
	test_dgetrf();
	// test_dgesv();
	test_lu_block();

	return 0;
}
