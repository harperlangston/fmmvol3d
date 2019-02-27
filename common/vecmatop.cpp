/* Kernel Independent Fast Multipole Method

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */
#include "blas.h"
#include "lapack.h"
#include "svdrep.hpp"
//#include "svd.h"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"

#include "vecmatop.hpp"

static double gq(double xx,double yy,double zz) { return sqrt((xx + 2*yy + 3*zz)/6.0); }
static double gx(double xx,double yy,double zz) { return (xx + 2*yy + 3*zz)/18.0; }
static double gy(double xx,double yy,double zz) { return (xx + 2*yy + 3*zz)/9.0; }
static double gz(double xx,double yy,double zz) { return (xx + 2*yy + 3*zz)/6.0; }

using std::cerr;

// ---------------------------------------------------------------------- 
int dscal(double alpha, DblNumVec& X)
{
  dscal( X.m(), alpha, X.data() );
  return 0;
}
// ---------------------------------------------------------------------- 
int dscal(int n, double alpha, double* X)
{
  int incx = 1;
  DSCAL(&n, &alpha, X, &incx);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int daxpy(double a, const DblNumVec& X, DblNumVec& Y)
{
  assert( X.m() == Y.m() );
  daxpy(X.m(), a, X.data(), Y.data());
  return 0;
}
// ---------------------------------------------------------------------- 
int daxpy(double a, const DblNumMat& X, DblNumMat& Y)
{
  assert( X.m() == Y.m() );  assert( X.n() == Y.n() );
  iC( daxpy(X.m()*X.n(), a, X.data(), Y.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int daxpy(int n, double a, double* X, double* Y)
{
  int incx = 1;  int incy = 1;
  DAXPY(&n, &a, X, &incx, Y, &incy);
  return 0;
}

int dcopy(const DblNumVec& X, DblNumVec& Y)
{
  assert( X.m() == Y.m() );
  dcopy(X.m(), X.data(), Y.data());
  return 0;
}
// ---------------------------------------------------------------------- 
int dcopy(const DblNumMat& X, DblNumMat& Y)
{
  assert( X.m() == Y.m() );  assert( X.n() == Y.n() );
  iC( dcopy(X.m()*X.n(), X.data(), Y.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dcopy(int n, double* X, double* Y)
{
  int incx = 1;  int incy = 1;
  DCOPY(&n, X, &incx, Y, &incy);
  return 0;
}

// ---------------------------------------------------------------------- 
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  iC( dgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C)
{
  char transa = 'N';
  char transb = 'N';
  assert(m!=0 && n!=0 && k!=0);
  DGEMM(&transa, &transb, &m, &n, &k,
		  &alpha, A, &m, B, &k, &beta, C, &m);
  return 0;
}
// ---------------------------------------------------------------------- 
int dger(double alpha, const DblNumVec& X, const DblNumVec& Y, DblNumMat& A)
{
  assert(X.m() == A.m());
  assert(Y.m() == A.n());
  iC( dger(A.m(), A.n(), alpha, X.data(), Y.data(), A.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dger(int m, int n, double alpha, double* X, double* Y, double* A)
{
  assert(m!=0 && n!=0);
  int incx = 1;  int incy = 1;
  DGER(&m, &n, &alpha, X, &incx, Y, &incy, A, &m);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  iC( dgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y)
{
  char trans = 'N';
  assert(m!=0 && n!=0);
  int incx = 1;
  int incy = 1;
  DGEMV(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}
//Assume alhpha == 1, beta == 1
int dgemv(const DblNumMat& A, const DblNumVec& X, DblNumVec& Y)
{
  double one = 1.0;
  iC( dgemv(A.m(), A.n(), one, A.data(), X.data(), one, Y.data()));
  return 0;
}
  

// ---------------------------------------------------------------------- 
int tran(const DblNumMat& M, DblNumMat& R)
{
  assert(R.m()==M.n() && R.n()==M.m());  //R.resize(M.n(), M.m());
  for(int i=0; i<M.m(); i++)
	 for(int j=0; j<M.n(); j++)
		R(j,i) = M(i,j);
  return 0;
}

int tReg(const DblNumMat& M, double epsilon, DblNumMat& R){
  DblNumMat MT(M.n(), M.m());
  iC(tran(M, MT));
  DblNumMat I(M.n(), M.n());
  for (int i=0; i<I.m(); i++){ I(i,i) = epsilon; }
  double a = 1.0; 
  dgemm(a, MT, M, a, I);
  DblNumMat Rtmp(I.m(), I.n());
  pinv(I, epsilon, Rtmp, 0);
  dgemm(a, Rtmp, MT, epsilon, R);
  return 0;
}

int pinv(const DblNumMat& M, double epsilon, DblNumMat& R){
  //cerr << "epsilon =  " << epsilon << endl;
  return pinv(M, epsilon, R, 0);
}

// ----------------------------------------------------------------------
int pinv(const DblNumMat& M, double epsilon, DblNumMat& R, int type)
{
  //if (epsilon < pow(0.1, 9)) epsilon = pow(0.1, 9);
  iA(M.m() == R.n());  assert(M.n() == R.m());
  SVDRep svd;
  iC( svd.construct(epsilon, M, type) );
  double cutoff = 0.0;
  if (type == 0){
	 cutoff = svd.S()(0)*epsilon;
  }
  else {
	 iA(type == 1);
	 pinv(M,epsilon,R,0);
	 double mx = (double)(M.m()) > (double)(M.n()) ? (double)(M.m()) : (double)(M.n());
	 cutoff = (svd.S())(0)*mx * 2e-15;
	 cutoff = mx * 2e-15;
	 //if (cutoff > 1e-15) { cutoff = 1e-15; }
	 //cutoff = 1e-10;
  }
  for(int i=0; i<svd.S().m(); i++) {
    if( svd.S()(i) >= cutoff) {
      svd.S()(i) = 1.0/(svd.S()(i));
	 } else {
		//ebiAssert(0);
		svd.S()(i) = 0.0;
	 }
  }

  DblNumMat UT(svd.U().n(),  svd.U().m());
  DblNumMat V( svd.VT().n(), svd.VT().m());
  iC( tran(svd.U(), UT) );
  iC( tran(svd.VT(), V) );
  for(int i=0; i<V.m(); i++)
    for(int j=0; j<V.n(); j++) {
      V(i,j) = V(i,j) * svd.S()(j);
	 }
  char transa = 'N';
  char transb = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int m = V.m();
  int n = UT.n();
  int k = V.n();
  DGEMM(&transa, &transb, &m, &n, &k, &alpha,
		  V.data(), &m, UT.data(), &k, 
		  &beta, R.data(), &m);  

  /*
 //invert Svd
  double cutoff = epsilon * svd.S()(0);  //double cutoff = epsilon;
  //	cerr<<epsilon <<" " <<cutoff << " " <<svd.S()(0)<<" "<<svd.S()(svd.S().m()-1)<<" "<<svd.S()(0)/svd.S()(svd.S().m()-1)<<endl;

  if ((epsilon < pow(0.1, 7))) {
	 for(int i=0; i<svd.S().m(); i++) {
		double t =  svd.S()(i);
		svd.S()(i) = (svd.S()(i)*svd.S()(i)) + 1e-15;
		if( t >= cutoff) {
		  svd.S()(i) = t/(svd.S()(i));
		} else {
		  assert(0);	
		  svd.S()(i) = 0.0;
		}
	 }
  }
  else {
	 for(int i=0; i<svd.S().m(); i++) {
		double t =  svd.S()(i);
		if( t >= cutoff) {
		  svd.S()(i) = 1.0/t;
		} else {
		  assert(0);
		  svd.S()(i) = 0.0;
		}
	 }
  }

  
  DblNumMat UT(svd.U().n(),  svd.U().m());
  DblNumMat V( svd.VT().n(), svd.VT().m());
  iC( tran(svd.U(), UT) );
  iC( tran(svd.VT(), V) );
  for(int i=0; i<V.m(); i++)
    for(int j=0; j<V.n(); j++) {
      V(i,j) = V(i,j) * svd.S()(j);
	 }
  
  char transa = 'N';
  char transb = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int m = V.m();
  int n = UT.n();
  int k = V.n();
  
  DGEMM(&transa, &transb, &m, &n, &k, &alpha,
    V.data(), &m, UT.data(), &k, 
    &beta, R.data(), &m);
  */
	 
  return 0;
}

// ---------------------------------------------------------------------- 
int inv(const DblNumMat& M, DblNumMat& R) //Gaussian Elimination
{
  //OR pinv(M, 0.0, R);
  assert(M.m()==M.n() && R.m()==R.n() && M.m()==R.m());
  memcpy(R.data(), M.data(), M.m()*M.n()*sizeof(double));
  int info;
  int m = M.m();
  int* ipiv = new int[m];
  DGETRF(&m, &m, R.data(), &m, ipiv, &info); assert(info==0);
  int lwork = m;
  double* work = new double[lwork];
  DGETRI(&m, R.data(), &m, ipiv, work, &lwork, &info);  assert(info==0);
  delete [] ipiv;
  delete [] work;
  return 0;
}
// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------- 
int spev1d(int evflag, int dmflag, int dof, double* data, int m, double e, int i, double u, double* res)
{
  int is[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++) is[k]=(i+k-1);
  }
  DblNumMat M(dof,m,false,data); //assert(M.n()==n && M.m()==res.m()); //double dof = M.m();  //int cnt = 0;
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 double scl = 1.0;
	 double us[4]; iC( spcoef(EVFLAG_VL, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++) 
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 double scl = double(m) / e;
	 double us[4]; iC( spcoef(EVFLAG_FD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 double scl = double(m*m)/(e*e);
	 double us[4]; iC( spcoef(EVFLAG_SD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int spev2d(int evflag, int dmflag, int dof, double* data, int* mn, double* ef, int* ij, double* uv, double* res)
{
  int m = mn[0];  int n = mn[1];
  double e = ef[0];  double f = ef[1];
  int i = ij[0];  int j = ij[1];
  double u = uv[0];  double v = uv[1];
  
  int is[4]; int js[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
	 for(int k=0; k<4; k++) js[k]=(j+k-1 + n) % n;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++)	is[k]=(i+k-1);
	 assert(j>=1 && j<=n-3);
	 for(int k=0; k<4; k++) js[k]=(j+k-1);
  }
  DblNumMat M(dof,m*n,false,data);
  double scl;
  double us[4], vs[4];
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 scl = double(m)/e;
	 iC( spcoef(EVFLAG_FD, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b];
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n)/f;
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 scl = double(m*m)/(e*e);
	 iC( spcoef(EVFLAG_SD, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(m*n)/(e*f);
	 iC( spcoef(EVFLAG_FD, u, us) );
	 iC( spcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n*n)/(f*f);
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_SD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int spcoef(int evflag, double u, double* us)
{
  double u1 = u;
  double u2 = u*u;
  double u3 = u*u*u;
  if(       evflag==EVFLAG_VL) {
	 us[0] = (  1 - 3*u1 + 3*u2 -   u3)/6.0;
	 us[1] = (  4        - 6*u2 + 3*u3)/6.0;
	 us[2] = (  1 + 3*u1 + 3*u2 - 3*u3)/6.0;
	 us[3] = (                  +   u3)/6.0;
  } else if(evflag==EVFLAG_FD) {
	 us[0] = (- 3 + 6*u1 - 3*u2)/6.0;
	 us[1] = (    -12*u1 + 9*u2)/6.0;
	 us[2] = (  3 + 6*u1 - 9*u2)/6.0;
	 us[3] = (             3*u2)/6.0;
  } else if(evflag==EVFLAG_SD) {
	 us[0] = (  6 - 6*u1 ) / 6.0;
	 us[1] = (-12 +18*u1 ) / 6.0;
	 us[2] = (  6 -18*u1 ) / 6.0;
	 us[3] = (      6*u1 ) / 6.0;	 //assert(0); //TODO;
  }  
  return 0;
}
// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------- 
int lagev1d(int evflag, int dmflag, int dof, double* data, int m, double e, int i, double u, double* res)
{
  int is[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++) is[k]=(i+k-1);
  }
  DblNumMat M(dof,m,false,data); //assert(M.n()==n && M.m()==res.m()); //double dof = M.m();  //int cnt = 0;
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 double scl = 1.0;
	 double us[4]; iC( lagcoef(EVFLAG_VL, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++) 
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 double scl = double(m) / e;
	 double us[4]; iC( lagcoef(EVFLAG_FD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 double scl = double(m*m)/(e*e);
	 double us[4]; iC( lagcoef(EVFLAG_SD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int lagev2d(int evflag, int dmflag, int dof, double* data, int* mn, double* ef, int* ij, double* uv, double* res)
{
  int m = mn[0];  int n = mn[1];
  double e = ef[0];  double f = ef[1];
  int i = ij[0];  int j = ij[1];
  double u = uv[0];  double v = uv[1];
  
  int is[4]; int js[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
	 for(int k=0; k<4; k++) js[k]=(j+k-1 + n) % n;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++)	is[k]=(i+k-1);
	 assert(j>=1 && j<=n-3);
	 for(int k=0; k<4; k++) js[k]=(j+k-1);
  }
  DblNumMat M(dof,m*n,false,data);
  double scl; 
  double us[4], vs[4];
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( lagcoef(EVFLAG_VL, u, us) );
	 iC( lagcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 scl = double(m)/e;
	 iC( lagcoef(EVFLAG_FD, u, us) );
	 iC( lagcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b];
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n)/f;
	 iC( lagcoef(EVFLAG_VL, u, us) );
	 iC( lagcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 scl = double(m*m)/(e*e);
	 iC( lagcoef(EVFLAG_SD, u, us) );
	 iC( lagcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(m*n)/(e*f);
	 iC( lagcoef(EVFLAG_FD, u, us) );
	 iC( lagcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n*n)/(f*f);
	 iC( lagcoef(EVFLAG_VL, u, us) );
	 iC( lagcoef(EVFLAG_SD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int lagev3d(int evflag, int dmflag, int dof, double* data, int* mno, double* efg, int* ijk, double* uvw, double* res)
{
  assert(evflag==EVFLAG_VL);
  
  int m = mno[0];  int n = mno[1];  int o = mno[2];
  double e = efg[0];  double f = efg[1];  double g = efg[2];
  int i = ijk[0];  int j = ijk[1];  int k = ijk[2];
  double u = uvw[0];  double v = uvw[1];  double w = uvw[2];
  
  int is[4]; int js[4];  int ks[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int h=0; h<4; h++) is[h]=(i+h-1 + m) % m;
	 for(int h=0; h<4; h++) js[h]=(j+h-1 + n) % n;
	 for(int h=0; h<4; h++) ks[h]=(k+h-1 + o) % o;
  } else {
	 assert(i>=1 && i<=m-3);
	 for(int h=0; h<4; h++)	is[h]=(i+h-1);
	 assert(j>=1 && j<=n-3);
	 for(int h=0; h<4; h++) js[h]=(j+h-1);
	 assert(k>=1 && k<=o-3);
	 for(int h=0; h<4; h++) ks[h]=(k+h-1);
  }
  DblNumMat M(dof,m*n*o,false,data);
  double scl; 
  double us[4], vs[4], ws[4];
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( lagcoef(EVFLAG_VL, u, us) );
	 iC( lagcoef(EVFLAG_VL, v, vs) );
	 iC( lagcoef(EVFLAG_VL, w, ws) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++)
		  for(int c=0; c<4; c++) {
			 double coef = us[a]*vs[b]*ws[c]; 
			 for(int d=0; d<dof; d++)
				res[d] += coef * M(d, is[a]+js[b]*m+ks[c]*m*n);
		  }
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int lagcoef(int evflag, double u, double* us)
{
  //double u1 = u;
  double u2 = u*u;
  double u3 = u*u*u;
  if(       evflag==EVFLAG_VL) {
	 us[0] = (u3-3*u2+2*u)/(-6.0);
	 us[1] = (u3-2*u2-u+2)/(2.0);
	 us[2] = (u3-u2-2*u)  /(-2.0);
	 us[3] = (u3-u)       /(6.0);
  } else if(evflag==EVFLAG_FD) {
	 us[0] = (3*u2-6*u+2)/(-6.0);
	 us[1] = (3*u2-4*u-1)/(2.0);
	 us[2] = (3*u2-2*u-2)/(-2.0);
	 us[3] = (3*u2-1)    /(6.0);
  } else if(evflag==EVFLAG_SD) {
	 assert(0); //TODO
  }
  return 0;
}

// ----------------------------------------------------------------------
int shep3d(const int srcdof, const DblNumMat& SPOS, const DblNumVec& SDEN, const DblNumMat& TPOS, DblNumVec& TDEN){
 long int n = SPOS.n();
  iA( srcdof*n == SDEN.m() );
  long int t = TPOS.n();
  iA( srcdof*t == TDEN.m() );

  //iA( n >= 10 );
  //std::cout << SPOS << endl;
  /* Set these parameters to standard values - see qshep3d.f */
  long int nq = 17; long int nw = 32; long int nr = (long int)(pow(((double)n)/3.0,(1.0/3.0))); long int ier;
  if (nr < 3) nr = 3;
  if (n < 33) { nw = n - 1; }
  if (n < 18) { nq = n - 1; }
  long int lout = 6; long int lcell[nr*nr*nr]; long int lnext[n];
  float xyzmin[3]; float xyzdel[3];
  float rsq[n]; float rmax;
  
  float x[n]; float y[n]; float z[n]; float f[n*srcdof];
  float a[9*n];

  //cout <<
  for (int k = 0; k < n; k++) {
	 x[k] = (float)SPOS(0,k); y[k] = (float)SPOS(1,k); z[k] = (float)SPOS(2,k);
  }

  for (int d = 0; d < srcdof; d++){
	 for (int k = 0; k < n; k++) {
		f[k] = (float)SDEN(k*srcdof + d);
	 }
	 
	 QSHEP3(&n, x, y, z, f, &nq, &nw, &nr, lcell, lnext, xyzmin, xyzdel, &rmax, rsq, a, &ier);
	 
	 if ((int)ier != 0){
		//std::cout << SPOS << std::endl;
		std::cout << "ier = " << (int)ier << " " << ier << " " << (short)ier << std::endl;
		std::cout << "nw = " << nw << " nr = " << nr << " nq =  " << nq << " n = " << n << std::endl;
		for (int k = 0; k < n; k++) {
		  std::cout << x[k] << " " << y[k] << " " << z[k] << " " << f[k] << std::endl;
		}
		std::cout << "SPOS = " << SPOS << std::endl;
		std::cout << "TPOS = " << TPOS << std::endl;
		std::cout << "SDEN = " << SDEN << std::endl;
		iA( (int)ier == 0);
	 }
  
	 for (int i = 0; i < t; i++){
		float px = TPOS(0,i); float py = TPOS(1,i); float pz = TPOS(2,i);
		//std::cout << px << " " << py << " " << pz << std::endl;
		TDEN(i*srcdof + d) =  (double)QS3VAL(&px, &py, &pz, &n, x, y, z, f, &nr, lcell, lnext, xyzmin, xyzdel, &rmax, rsq, a);
	 }
	 for (int i = 0; i < nr*nr*nr; i++) lcell[i] = 0;
	 for (int i = 0; i < n; i++) f[i] = 0;
	 
  }
  return 0;
}
