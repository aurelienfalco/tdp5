#include "cblas.h"
#include <assert.h>
#include <stdio.h>
#include "util.h"

void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	const enum CBLAS_DIAG Diag, const int M, const int N,
	const double alpha, const double *A, const int lda,
	double *B, const int ldb)
{
	int i, j, k;
	assert(Order == CblasColMajor);
	assert(TransA == CblasNoTrans);
	assert(Side == CblasLeft);
	assert(alpha);
	if (Uplo == CblasUpper) {
		for (j = 0; j < N; j++){
			for (k = M-1; k >= 0; k--) {
				if (Diag == CblasNonUnit)
					B[k + j * ldb] = B[k + j * ldb] / A[k + k * lda];
				for (i = 0; i < k ; i++) {
					B[i + j * ldb] = B[i + j * ldb] - B[k + j * ldb] * A[i + k * lda];
				}
			}
		}
	}
	else {
		for (j = 0; j < N; j++){
			for (k = 0; k < M-1; k++) {
				if (Diag == CblasNonUnit)
					B[k + j * ldb] = B[k + j * ldb] / A[k + k * lda];
				for (i = k+1; i < M ; i++) {
					B[i + j * ldb] = B[i + j * ldb] - B[k + j * ldb] * A[i + k * lda];
				}
			}
		}
	}
}

