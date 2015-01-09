#include "util.h"

float epsilon = 0.0001;

double* init_matrix(int m, int n)
{
	double* L = calloc(m * n,sizeof(double));
	return L;		
}

double* rand_matrix(int m, int n)
{
	double* A = init_matrix(m,n);
	int i,j, lda = m;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			A[i + j * lda] = rand()%10;
		}
		if (m <= n)
			A[i + i * lda] = rand()%10 + 10;
	}
	return A;		
}


void print_matrix(double* A, int m, int n, int lda)
{
	int i,j;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			printf("%g  ", A[i + j * lda]);
		}
		printf("\n");
	}
	printf("\n");
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
		U[i + i * size] = rand()%10 + 10;
	}
	return U;		
}

double* copy_matrix(double* B, int m, int n, int lda)
{
	double* A = init_matrix(m,n);
	int i,j;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			A[i + j * lda] = B[i + j * lda];
		}
	}
	return A;
}


int equal_matrix(int UPLO, double* A, double* B, int m, int n, int lda)
{
	switch (UPLO){
		case EqAll: // ordinary matrix
		for (int i = 0; i < m; ++i){
			for (int j = 0; j < n; ++j){
				if (fabs(A[i + j * lda] - B[i + j * lda]) > epsilon){
					fprintf(stderr,"[%d %d]\t %g != %g \t (précision %g)\n",i,j, A[i + j * lda], B[i + j * lda], epsilon);
					return 0;
				}
			}
		}
		break;
		case EqUpper: // UPPER matrix
		for (int i = 0; i < m; ++i){
			for (int j = i; j < n; ++j){
				if (fabs(A[i + j * lda] - B[i + j * lda]) > epsilon){
					fprintf(stderr,"[%d %d]\t %g != %g \t (précision %g)\n",i,j, A[i + j * lda], B[i + j * lda], epsilon);
					return 0;
				}
			}
		}
		break;
		case EqLower: // LOWER matrix
		for (int i = 0; i < m; ++i){
			for (int j = 0; j < i; ++j){
				if (fabs(A[i + j * lda] - B[i + j * lda]) > epsilon){
					fprintf(stderr,"[%d %d]\t %g != %g \t (précision %g)\n",i,j, A[i + j * lda], B[i + j * lda], epsilon);
					return 0;
				}
			}
		}
		break;
	} 
	return 1;
}
