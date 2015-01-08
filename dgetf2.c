#include "cblas.h"
#include "util.h"

#include <assert.h>

void LAPACKE_dgetf2(const enum CBLAS_ORDER matrix_order, int m, int n, double* a, int lda, int* ipiv)
{
	assert(matrix_order == CblasColMajor);
	(void)ipiv;
	double epsilon = 0.0001;

	for (int i = 0; i < MIN(m,n); i++) {
		if (a[i * lda+ i] > epsilon)
			cblas_dscal(m-i-1,1/a[i * lda+ i],&a[i * lda + i+1],1);
		else if (a[i * lda+ i] != 0)
			for (int j = 0; j < m-i-1; j++)
				a[i * lda + i+j] /= a[i * lda+ i];
		cblas_dger(matrix_order,m-i-1,n-i-1,-1,&a[i*lda + i+1],1,&a[(i+1)*lda + i],lda,&a[(i+1)*lda + i+1],lda);
	}
}