#ifndef __SVD_H__
#define __SVD_H__

int svd(int m,int n,int withu,int withv,double eps,double tol,
		  double **a,double *q,double **u,double **v);

int harper_pinv(int m,int n, double eps,double tol, double **A, double **OUT);

int matVmult(double **OUT, double **U, double *S, int n);
int matTpose(double **INOUT, int m);
int matMatmult(double **OUT, double **USinv, double **Ut, int m, int n);
int invVec(double *S, int m);

#endif
