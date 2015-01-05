#include "cblas.h"
#include "util.h"

#include <assert.h>

void LAPACKE_dgetf2(const enum CBLAS_ORDER matrix_order, int m, int n, double* a, int lda, int* ipiv)
{
	assert(matrix_order == CblasColMajor);
	(void)ipiv;

	for (int i = 0; i < n; i++) {
		cblas_dscal(m-i-1,1/a[i * lda+ i],&a[i * lda + i+1],1);
		cblas_dger(matrix_order,m-i-1,m-i-1,-1,&a[i*lda + i+1],1,&a[(i+1)*lda + i],lda,&a[(i+1)*lda + i+1],lda);
	}
}