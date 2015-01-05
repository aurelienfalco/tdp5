#ifndef __DGEMM_SCALAIRE_H
#define __DGEMM_SCALAIRE_H

void cblas_dgemm_scalaire_ijk(const int M, const double *A, const int lda,
			      const double *B, const int ldb, double *C,
			      const int ldc);
#endif /* __DGEMM_SCALAIRE_H */
