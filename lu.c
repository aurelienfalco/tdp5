#include "cblas.h"
#include "util.h"

#include <assert.h>

void cblas_lu(const enum CBLAS_ORDER order, const int m, const int n, double* a, int lda, int nb)
{
	assert(order == CblasColMajor);
	for (int j = 0; j < MIN(m,n); j+=nb) {
		int jb = MIN(MIN(m,n)-j,nb);
		LAPACKE_dgetf2(order, m-j, jb, &a[j * lda + j], lda, NULL);
		for (int i = j+jb; i < MIN(m,n); i+=jb)
		{
			int ib = MIN(MIN(m,n)-i,jb);
			cblas_dtrsm(order, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit, ib, ib, 1, &a[j + lda * j], lda, &a[j + lda * i], lda );
			cblas_dgemm(order, CblasNoTrans, CblasNoTrans, m-j-jb, ib, ib, -1, &a[j + jb + j * lda], lda,&a[j +  i * lda], lda, 1, &a[j + jb + lda * i], lda);
		}
	}
}
