#ifndef __UTIL_H__
#define __UTIL_H__

#include "dgemm_scalaire.h"
#include "cblas.h"
#include "time.h"
#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(a,b) ((a)<(b)?(a):(b))

double* init_matrix(int m, int n);
double* rand_matrix(int m, int n);
void print_matrix(double* A, int m, int n, int lda);

double* rand_tri_sup(int size);
double* rand_tri_inf(int size);
double* copy_matrix(double* B, int m, int n, int lda);
int equal_matrix(int UPLO, double* A, double* B, int m, int n, int lda);

#endif