#include "cblas.h"
#include "util.h"

#include <assert.h>

int LAPACKE_dgesv(const enum CBLAS_ORDER matrix_order, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb)
{
	assert(matrix_order == CblasColMajor);
	(void)nrhs;
	(void)ipiv;

 	LAPACKE_dgetrf(matrix_order, n, n, a, lda, 4);
 	cblas_dtrsm(matrix_order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, n, 1, a, lda, b, ldb);
 	cblas_dtrsm(matrix_order, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1, a, lda, b, ldb);
 	return 0;
}