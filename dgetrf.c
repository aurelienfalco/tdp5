#include "cblas.h"
#include "util.h"

#include <assert.h>

int LAPACKE_dgetrf(const enum CBLAS_ORDER  matrix_order, int m, int n, double* a, int lda, int* ipiv)
{
	assert(matrix_order == CblasColMajor);
	int nb = 4;
	for (int j = 0; j < MIN(m,n); j+=nb) {
		int jb = MIN(MIN(m,n)-j,nb);
		LAPACKE_dgetf2(matrix_order, m-j, jb, &a[j * lda + j], lda, ipiv);
		cblas_dtrsm(matrix_order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, jb, n-j-jb+1, 1, &a[j + lda * j], lda, &a[j + (j+jb) * lda], lda );
		cblas_dgemm(matrix_order, CblasNoTrans, CblasNoTrans, m-j-jb, n-j-jb, jb, -1, &a[j+jb + j * lda], lda,&a[j +  (j+jb) * lda], lda, 1, &a[j+jb +  (j+jb) * lda], lda);
	}
	return 0;
}
