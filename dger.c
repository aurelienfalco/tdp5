#include "cblas.h"
#include "util.h"

#include <assert.h>

void cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda)
{
	assert(order == CblasColMajor);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j ++) {
			A[i + j * lda] += alpha * X[i * incX] * Y[j * incY];
		}
	}
}
