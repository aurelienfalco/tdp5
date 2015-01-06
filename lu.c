#include "cblas.h"
#include "util.h"

#include <assert.h>

void cblas_lu(const enum CBLAS_ORDER order, const int M, const int N,
	const double alpha, const double *X, const int incX,
	const double *Y, const int incY, double *A, const int lda, const int bs)
{
	assert(order == CblasColMajor);

	for (int i = 0; i < MIN(M,N); i+=bs) {
		int index = i*(lda + 1); 
		LAPACKE_dgetf2(CblasColMajor,bs,bs,&A[ i*(lda + 1) ],lda,NULL);
		// trsm colonnes
		for (int j = index + bs; j < M; j+=bs){
			cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, bs, bs, 1, &A[ index ], lda, &A[j], lda);
		}
		for (int j = index + bs*lda; j < N; j+=bs*lda){
			cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, bs, bs, 1, &A[ i*(lda + 1) ], lda, &A[ i*(lda + 1) + j], lda);
			for (int k = bs; k < M; k+=bs){
				cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,bs,bs,bs,1,&A[index + k],lda,&A[j],lda,1,&A[j + k],lda);
			}
		}
	}
}
