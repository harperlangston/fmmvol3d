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
#include "common/vecmatop.hpp"
#include "fmm3d.hpp"

using std::istringstream;

template <class N>
FMM3d<N>::FMM3d(const string& p):
  KnlMat3d(p), _center(0,0,0), _rootLevel(0),
  _np(6), _let(NULL), _matmgnt(NULL)
{
}

template <class N>
FMM3d<N>::~FMM3d()
{
  if(_let!=NULL){ delete _let; _let = NULL; }
}

// ----------------------------------------------------------------------
template <class N>
int FMM3d<N>::SrcEqu2TrgChk_dgemv(const DblNumMat& srcPos, const DblNumMat& srcNor, const DblNumMat& trgPos, const DblNumVec& srcDen, DblNumVec& trgVal)
{
  int TMAX = 1024;
  if(trgPos.n()<=TMAX) {
	 int MM = trgPos.n() * _knl.trgDOF();
	 int NN = srcPos.n() * _knl.srcDOF();
	 DblNumMat tmp(MM,NN);
	 iC( _knl.kernel(srcPos, srcNor, trgPos, tmp) );
	 iC( dgemv(tmp, srcDen, trgVal) );
  } else {
	 int RUNS = (trgPos.n()-1) / TMAX + 1;
	 for(int r=0; r<RUNS; r++) {
		int stt = r*TMAX;
		int end = min((r+1)*TMAX, trgPos.n());
		int num = end-stt;
		int MM = num * _knl.trgDOF();
		int NN = srcPos.n() * _knl.srcDOF();
		DblNumMat tps(dim(), num, false, trgPos.data() + stt*dim() );
		DblNumVec tvl(num*_knl.trgDOF(), false, trgVal.data() + stt*_knl.trgDOF());
		DblNumMat tmp(MM,NN);
		iC( _knl.kernel(srcPos, srcNor, tps, tmp) );
		iC( dgemv(tmp, srcDen, tvl) );
	 }
  }
  return (0);
}
// ----------------------------------------------------------------------
/*! Used exclusively in the source to multipole computation, so trgRad and srcRad are the same */
template <class N>
int FMM3d<N>::SrcEqu2UpwChk_dgemv(const DblNumMat& srcPos, const DblNumMat& srcNor, Point3 trgCtr, double trgRad, const DblNumVec& srcDen, DblNumVec& trgVal)
{
  DblNumMat trgPos; iC( _matmgnt->localPos(UC, trgCtr, trgRad, trgPos) );
  int MM = trgPos.n() * _knl.trgDOF();
  int NN = srcPos.n() * _knl.srcDOF();
  DblNumMat tmp(MM,NN);
  iC( _knl.kernel(srcPos, srcNor, trgPos, tmp) );
  iC( dgemv(tmp, srcDen, trgVal) );
  return (0);
}

// ---------------------------------------------------------------------- 
template <class N>
int FMM3d<N>::SrcEqu2DwnChk_dgemv(const DblNumMat& srcPos, const DblNumMat& srcNor, Point3 trgCtr, double trgRad, const DblNumVec& srcDen, DblNumVec& trgVal)
{
  DblNumMat trgPos; iC( _matmgnt->localPos(DC, trgCtr, trgRad, trgPos) );
  int MM = trgPos.n() * _knl.trgDOF();
  int NN = srcPos.n() * _knl.srcDOF();
  DblNumMat tmp(MM,NN);
  iC( _knl.kernel(srcPos, srcNor, trgPos, tmp) );
  iC( dgemv(tmp, srcDen, trgVal) );
  return 0;
}
// ---------------------------------------------------------------------- 
template <class N>
int FMM3d<N>::DwnEqu2TrgChk_dgemv(Point3 srcCtr, double srcRad, const DblNumMat& trgPos, const DblNumVec& srcDen, DblNumVec& trgVal)
{
  int TMAX = 1024;
  if(trgPos.n()<=TMAX) {
	 DblNumMat srcPos; iC( _matmgnt->localPos(DE, srcCtr, srcRad, srcPos) );
	 int MM = trgPos.n() * _knl_mm.trgDOF();
	 int NN = srcPos.n() * _knl_mm.srcDOF();
	 DblNumMat tmp(MM,NN);
	 iC( _knl_mm.kernel(srcPos, srcPos, trgPos, tmp) );
	 iC( dgemv(tmp, srcDen, trgVal) );
  } else {
	 DblNumMat srcPos; iC( _matmgnt->localPos(DE, srcCtr, srcRad, srcPos) );
	 int RUNS = (trgPos.n()-1) / TMAX + 1;
	 for(int r=0; r<RUNS; r++) {
		int stt = r*TMAX;
		int end = min((r+1)*TMAX, trgPos.n());
		int num = end-stt;
		int MM = num * _knl_mm.trgDOF();
		int NN = srcPos.n() * _knl_mm.srcDOF();
		DblNumMat tps(dim(), num, false, trgPos.data() + stt*dim());
		DblNumVec tvl(num*_knl_mm.trgDOF(), false, trgVal.data() + stt*_knl_mm.trgDOF());
		DblNumMat tmp(MM, NN);
		iC( _knl_mm.kernel(srcPos, srcPos, tps, tmp) );
		iC( dgemv(tmp, srcDen, tvl) );
	 }
  }
  return 0;
}

// ----------------------------------------------------------------------
template <class N>
int FMM3d<N>::UpwEqu2DwnChk_dgemv(Point3 srcCtr, double srcRad, Point3 trgCtr, double trgRad, const DblNumVec& srcDen, DblNumVec& trgVal)
{
  DblNumMat srcPos; iC( _matmgnt->localPos(UE, srcCtr, srcRad, srcPos) );
  DblNumMat trgPos; iC( _matmgnt->localPos(DC, trgCtr, trgRad, trgPos) );
  int MM = trgPos.n() * _knl_mm.trgDOF();
  int NN = srcPos.n() * _knl_mm.srcDOF();
  DblNumMat tmp(MM,NN);
  iC( _knl_mm.kernel(srcPos, srcPos, trgPos, tmp) );
  iC( dgemv(tmp, srcDen, trgVal) );
  return 0;
}

// ---------------------------------------------------------------------- 
template <class N>
int FMM3d<N>::UpwEqu2TrgChk_dgemv(Point3 srcCtr, double srcRad, const DblNumMat& trgPos, const DblNumVec& srcDen, DblNumVec& trgVal)
{
  int TMAX = 1024;
  if(trgPos.n()<=TMAX) {
	 DblNumMat srcPos; iC( _matmgnt->localPos(UE, srcCtr, srcRad, srcPos) );
	 int MM = trgPos.n() * _knl_mm.trgDOF();
	 int NN = srcPos.n() * _knl_mm.srcDOF();
	 DblNumMat tmp(MM,NN);
	 iC( _knl_mm.kernel(srcPos, srcPos, trgPos, tmp) );
	 iC( dgemv(tmp, srcDen, trgVal) );
  } else {
	 DblNumMat srcPos; iC( _matmgnt->localPos(UE, srcCtr, srcRad, srcPos) );
	 int RUNS = (trgPos.n()-1) / TMAX + 1;
	 for(int r=0; r<RUNS; r++) {
		int stt = r*TMAX;
		int end = min((r+1)*TMAX, trgPos.n());
		int num = end-stt;
		int MM = num * _knl_mm.trgDOF();
		int NN = srcPos.n() * _knl_mm.srcDOF();
		DblNumMat tps(dim(), num, false, trgPos.data() + stt*dim());
		DblNumVec tvl(num*_knl_mm.trgDOF(), false, trgVal.data() + stt*_knl_mm.trgDOF());
		DblNumMat tmp(MM,NN);
		iC( _knl_mm.kernel(srcPos, srcPos, tps, tmp) );
		iC( dgemv(tmp, srcDen, tvl) );
	 }
  }
  return 0;
}

// ---------------------------------------------------------------------- 
template <class N>
DblNumMat FMM3d<N>::srcExaPos(int gNodeIdx) 
{
  iA(_let->terminal(gNodeIdx));
  N& node=_let->node(gNodeIdx);
  int beg = node.srcExaBeg();
  int num = node.srcExaNum();
  return DblNumMat(dim(), num, false, _srcExaPos.data()+beg*dim());
}
template <class N>
DblNumMat FMM3d<N>::srcExaNor(int gNodeIdx)
{
  N& node=_let->node(gNodeIdx);
  int beg = node.srcExaBeg();
  int num = node.srcExaNum();
  return DblNumMat(dim(), num, false, _srcExaNor.data()+beg*dim());
}
template <class N>
DblNumVec FMM3d<N>::srcExaDen(int gNodeIdx)
{
  N& node=_let->node(gNodeIdx);
  int beg = node.srcExaBeg();
  int num = node.srcExaNum();
  return DblNumVec(srcDOF()*num, false, _srcExaDen.data()+beg*srcDOF());
}
template <class N>
DblNumVec FMM3d<N>::srcUpwEquDen(int gNodeIdx)
{
  N& node=_let->node(gNodeIdx);
  int idx = node.srcNodeIdx();
  return DblNumVec(datSze(UE), false, _srcUpwEquDen.data()+idx*datSze(UE));
}
/*
DblNumVec FMM3d::srcUpwChkVal(int gNodeIdx)
{
  letNode& node=_let->node(gNodeIdx);
  int idx = node.srcNodeIdx();
  return DblNumVec(datSze(UC), false, _srcUpwChkVal.data()+idx*datSze(UC));
}
*/
template <class N>
DblNumMat FMM3d<N>::trgExaPos(int gNodeIdx)
{
  N& node=_let->node(gNodeIdx);
  int beg = node.trgExaBeg();
  int num = node.trgExaNum();
  return DblNumMat(dim(), num, false, _trgExaPos.data()+beg*dim());
}
template <class N>
DblNumVec FMM3d<N>::trgExaVal(int gNodeIdx)
{
  N& node=_let->node(gNodeIdx);
  int beg = node.trgExaBeg();
  int num = node.trgExaNum();
  return DblNumVec(trgDOF()*num, false, _trgExaVal.data()+beg*trgDOF());
}
/*
  template <class N>
  DblNumVec FMM3d<N>::trgDwnEquDen(int gNodeIdx)
  {
  N& node=_let->node(gNodeIdx);
  int idx = node.trgNodeIdx();
  return DblNumVec(datSze(DE), false, _trgDwnEquDen.data()+idx*datSze(DE));
}
*/
template <class N>
DblNumVec FMM3d<N>::trgDwnChkVal(int gNodeIdx)
{
  N& node=_let->node(gNodeIdx);
  int idx = node.trgNodeIdx();
  return DblNumVec(datSze(DC), false, _trgDwnChkVal.data()+idx*datSze(DC));
}


