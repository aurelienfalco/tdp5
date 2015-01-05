#include "dgemm_scalaire.h"
#include "assert.h"
#include "util.h"


void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
	const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
	const int K, const double alpha, const double *A,
	const int lda, const double *B, const int ldb,
	const double beta, double *C, const int ldc)
{
	assert(Order == CblasColMajor);
    assert(TransA == CblasNoTrans);
    assert(TransB == CblasNoTrans);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < K; k++) {
				/* C[i, j] += A[i, k] * B[k, j] */
				C[i + j * ldc] += alpha * A[i + k * lda] * B[k + j * ldb];
			}
			C[i + j * ldc] *= beta;
		}
	}
}

void cblas_dgemm_scalaire_ijk(const int M, const double *A, const int lda,
	const double *B, const int ldb, double *C,
	const int ldc)
{
	/* C <- tA * B */

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			for (int k = 0; k < M; k++) {

				/* C[i, j] += A[k, i] * B[k, j] */
				C[i + j * ldc] +=
				A[k + i * lda] * B[k + j * ldb];
			}
		}
	}
}
