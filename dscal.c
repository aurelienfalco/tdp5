#include "cblas.h"

void cblas_dscal(int N, double da, double *dx, int incx)
{
	int i;
	for (i = 0; i < N; i++)
	{
		*dx = *dx * da;
		dx += incx;
	}
}
