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
/*! \file */
#include "blas.h"
#include "lapack.h"
#include "svdrep.hpp"

using std::min;
using std::max;
using std::abs;

//int    SVDRep::_wssize = 4194304;
//int    SVDRep::_wssize = 4194304;
//double SVDRep::_wsbuf[4194304];

/*! svdrep.cpp provides implementations of int SVDRep::construct(double epsilon, const DblNumMat& K),
 *
 * int SVDRep::dgemv(double alpha, const DblNumVec& X, double beta, DblNumVec& Y, double tol) and
 *
 * int SVDRep::dgemv(double alpha, double* X, double beta, double* Y, double tol)
 *
 * as described in svdrep.hpp */
/* ********************************************************************** */
/*
int SVDRep::construct(double epsilon, const DblNumMat& K)
{
  int m = K.m();
  int n = K.n();
  int k = min(m, n);
  
  DblNumMat tU(m, k);
  DblNumVec tS(k);
  DblNumMat tVT(k, n);
  
  //SVD
  int INFO;
  char JOBU  = 'S';
  char JOBVT = 'S';
  
  int wssize = max(3*min(m,n)+max(m,n), 5*min(m,n));
  wssize *= wssize;
  //wssize *= 10;
  double* wsbuf = new double[wssize];
  DGESVD(&JOBU, &JOBVT, &m, &n, K.data(), &m, tS.data(), tU.data(), &m, tVT.data(), &k, wsbuf, &wssize, &INFO);  iA(INFO==0);
  delete [] wsbuf;
  
  //cutoff
  double cutoff = epsilon*tS(0);
  int cnt=0;
  while(cnt< k)
    if(abs(tS(cnt)) >= cutoff)
      cnt++;
    else
      break;
  _matU.resize(m, cnt);
  _matS.resize(cnt);	
  _matVT.resize(cnt,n);
  
  for(int i=0; i<m; i++)
    for(int j=0; j<cnt; j++)
      _matU(i,j) = tU(i,j);
  for(int i=0; i<cnt; i++)
    _matS(i) = tS(i);
  for(int i=0; i<cnt; i++)
    for(int j=0; j<n; j++)
      _matVT(i,j) = tVT(i,j);
  return 0;
  }
*/

/* ********************************************************************** */
#undef __FUNCT__
#define __FUNCT__ "SVDRep::construct"
int SVDRep::construct(double epsilon, const DblNumMat& K, const int type)
{
  int m = K.m();
  int n = K.n();
  int k = min(m, n);
  
  DblNumMat tU(m, k);
  DblNumVec tS(k);
  DblNumMat tVT(k, n);
  int wssize = 4194304;
  //SVD
  if(1) {
	 int INFO;
	 char JOBU  = 'S';
	 char JOBVT = 'S';
	 iA( wssize >= max(3*min(m,n)+max(m,n), 5*min(m,n)));
	 //double* wsbuf = _wsbuf;
	 double* wsbuf = new double[wssize];
	 DGESVD(&JOBU, &JOBVT, &m, &n, K.data(), &m, tS.data(), tU.data(), &m, tVT.data(), &k, wsbuf, &wssize, &INFO);
	 iA(INFO==0);
	 delete[] wsbuf;
	 
  }
  if(0) {
	 int INFO;
	 char JOBZ = 'S';

#if defined(EBI_linux_mpi) || defined(EBI_linux)
	 int wssize = 4*min(m,n)*min(m,n) + max(m,n) + 9*min(m,n);
#endif
#if defined(EBI_solaris)
	 int wssize = 3*min(m,n)*min(m,n) + max(max(m,n),4*min(m,n)*min(m,n)+4*min(m,n));
#endif
#if defined(EBI_alpha)
	 int wssize = 4*min(m,n)*min(m,n) + max(m,n) + 9*min(m,n);
	 ebiAssert(0);
#endif
	 int wisize = 8*min(m,n);
	 	 
	 double* wsbuf = new double[wssize];
	 int*    wibuf = new int[   wisize];
	 DGESDD(&JOBZ, &m, &n, K.data(), &m, tS.data(), tU.data(), &m, tVT.data(), &k, wsbuf, &wssize, wibuf, &INFO);
	 delete[] wsbuf;
	 delete[] wibuf;
	 
	 iA(INFO==0);
  }
  
  //cutoff
  double cutoff = 0.0;
  if (type == 0){ 
	 double cutoff = epsilon*tS(0);
  }
  else {
	 cutoff = 0.0;
  }
  int cnt=0;
  while(cnt< k)
    if(abs(tS(cnt)) >= cutoff)
      cnt++;
    else
      break;
  
  _matU.resize(m, cnt);
  _matS.resize(cnt);	
  _matVT.resize(cnt,n);
  
  for(int i=0; i<m; i++)
    for(int j=0; j<cnt; j++)
      _matU(i,j) = tU(i,j);
  for(int i=0; i<cnt; i++)
    _matS(i) = tS(i);
  for(int i=0; i<cnt; i++)
    for(int j=0; j<n; j++)
      _matVT(i,j) = tVT(i,j);

  return 0;
}

/* ********************************************************************** */
int SVDRep::dgemv(double alpha, const DblNumVec& X, double beta, DblNumVec& Y, double tol)
{
  iA(Y.m() == _matU.m());
  iA(X.m() == _matVT.n());
  iC( dgemv(alpha, X.data(), beta, Y.data(), tol) );
  iA(0);
  return 0;
}

/* ********************************************************************** */
int SVDRep::dgemv(double alpha, double* X, double beta, double* Y, double tol)
{	
  iA(0);
  int K = 1; //prevent matrix of zero size
  while(K<_matS.m() && _matS(K)>=tol) 
	 K++;
  //buf = VT(1:K,:) * X
  //double* buf = _wsbuf;
  double* buf = new double[K];
  {
	 char TRANS = 'N';
	 int M = _matVT.m();
	 int N = _matVT.n();
	 double ALPHA = 1.0;
	 double BETA = 0.0;
	 int INC = 1;
	 DGEMV(&TRANS, &K, &N, &ALPHA, _matVT.data(), &M, X, &INC, &BETA, buf, &INC); //first K rows
  }
  // buf = S(1:K) .* buf;
  for(int i=0; i<K; i++)
	 buf[i] = buf[i] * _matS(i);
  // y = U(:,1:K) * buf
  {
	 char TRANS = 'N';
	 int M = _matU.m(); //int N = _matU.n();
	 double ALPHA = alpha;
	 double BETA = beta;
	 int INC = 1;
	 DGEMV(&TRANS, &M, &K, &ALPHA, _matU.data(), &M, buf, &INC, &BETA, Y, &INC);	
  }
  delete [] buf;
  
  return 0;
}
