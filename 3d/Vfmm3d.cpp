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
#include "Vfmm3d.hpp"


using std::istringstream;
using namespace std;

template <class VF>
VFMM3d<VF>::VFMM3d(const string& p)
  : FMM3d<VF>(p)
{
  
}

template <class VF>
VFMM3d<VF>::~VFMM3d()
{
  delete _coeffs; _coeffs = NULL;
  delete _tbls; _tbls = NULL;
}

template <class VF>
int VFMM3d<VF>::GrdEqu2UpwChkTbls_dgemv(const int level, const DblNumVec& srcCoeffs, DblNumVec& tmpVal) {
  DblNumMat& _S2M = ((this->matmgnt())->hom()) ? _s2mNbrF[0] : _s2mNbrF[level];
  if (_S2M.m() == 0){
    int NK = vlet()->srcNk();
    _tbls->eval_s2m(level, _S2M);
    iA(_S2M.m() != 0);
  }
  iC( dgemv(_S2M, srcCoeffs, tmpVal) );
  return(0);
}

template <class VF>
int VFMM3d<VF>::GrdEqu2UpwEquTbls_dgemv(const int level, const DblNumVec& srcCoeffs, DblNumVec& tmpVal, DblNumVec& den)
{
  DblNumMat& _S2M = ((this->matmgnt())->hom()) ? _s2mNbrF[0] : _s2mNbrF[level];
  if (_S2M.m() == 0){
    int NK = vlet()->srcNk();
    _tbls->eval_s2m(level, _S2M);
    iA(_S2M.m() != 0);
  }
  iC( dgemv(_S2M, srcCoeffs, tmpVal) );
  iC( (this->matmgnt())->UpwChk2UpwEqu_dgemv(level, tmpVal, den, 2.0) );

  return(0);
}

template <class VF>
int VFMM3d<VF>::GrdCoeffs2TrgVal_dgemv(const int level, const int type, const int nbr, const DblNumVec& srcCoeffs, DblNumVec& trgVal){ 
  map<int, map<int, DblNumMat> >& _SC2TV_levMap =  ((this->matmgnt())->hom()) ? _srcCof2TrgVal[0] : _srcCof2TrgVal[level];
  map<int, DblNumMat>& _SC2TV_type = _SC2TV_levMap[type];

  DblNumMat& _SC2TV =  _SC2TV_type[nbr];
  int NK = vlet()->srcNk();
  //#pragma omp critical 
  if (_SC2TV.m() == 0){
	 _SC2TV.resize(vlet()->trgGrdSze()*(this->srcDOF()),NK*(this->trgDOF()));
	 int bref, symNum, posneg, pnt; Index3 refIdx; Index3 cofIdx;
	 iC( symNbrGetRefs(nbr, type, bref, symNum, posneg, refIdx, cofIdx));

	 int KT = vlet()->kTrgVal();
	 for (int k = 0; k < KT; k++){ for (int j = 0; j < KT; j++){ for (int i = 0; i < KT; i++){
			 pnt = (k*KT + j)*KT + i;	
			 for (int c = 0; c < NK; c++){ for (int ds = 0; ds < (this->srcDOF()); ds++){ for (int dt = 0; dt < (this->trgDOF()); dt++){
					 _SC2TV(pnt*(this->srcDOF()) + ds, c*(this->trgDOF()) + dt) = compSymNbrF(type, level, bref, symNum, posneg, refIdx, i, j, k, c, ds, dt, KT, refIdx, cofIdx);
				  } } } } } }
  }
  
  iA(_SC2TV.m() != 0);
  iC( dgemv(_SC2TV, srcCoeffs, trgVal) );
  return(0);
}

template <class VF>
int VFMM3d<VF>::GrdCoeffs2UnbalTrgVal_dgemv(const int type, const int gNodeIdx, const int nbr, const DblNumVec& coeffs, const Point3 offset, DblNumVec& trgval) {
  iA(type == FINE || type == CRSE);
  int levgni = vlet()->depth(gNodeIdx) + this->rootLevel();
  int levnbr = vlet()->depth(nbr) + this->rootLevel();
  double radgni = vlet()->radius(gNodeIdx);
  double radnbr = vlet()->radius(nbr);
  Point3 ctrgni(vlet()->center(gNodeIdx));
  Point3 ctrnbr(vlet()->center(nbr));
  
  int look_up = _tbls->lookup(type, levnbr, radnbr, ctrnbr, levgni, radgni, ctrgni, offset);
 
  map<int, map<int, DblNumMat> >& _SC2TV_type = (type == FINE) ? _FubSrcCof2TrgVal[vlet()->depth(gNodeIdx)] : _CubSrcCof2TrgVal[vlet()->depth(gNodeIdx)];
  map<int, DblNumMat>& _SC2TV_lev = _SC2TV_type[(vlet()->depth(gNodeIdx) - vlet()->depth(nbr))];
  DblNumMat& _SC2TV =  _SC2TV_lev[look_up];
  if (_SC2TV.m() == 0){
	 iC( _tbls->eval_unbal(type, look_up, levgni, radgni, ctrgni, levnbr, radnbr, ctrnbr, offset, _SC2TV) );
  }
  
  iC( dgemv(_SC2TV, coeffs, trgval) );

  return(0);
}

template <class VF>
int VFMM3d<VF>::DwnEqu2GrdChk_dgemv(const int depth, const DblNumVec& srcDen, DblNumVec& trgVal, bool eval){
  /* Not fixed for scale-varance, but will have to calc so few levels anyway */
  DblNumMat& de2gc  = _dwnEqu2GrdChk[depth];
  
  if (de2gc.m() == 0){
    double rad = (1.0/double(pow2(depth)));
    DblNumMat srcPos; iC( (this->matmgnt())->localPos(DE, Point3(0.0), rad, srcPos) );
    DblNumMat trgPos((this->dim()), vlet()->trgGrdSze());
    daxpy(rad, vlet()->grdTrgSamPos(), trgPos);
    int M = trgPos.n() * (this->_knl_mm).trgDOF();
    int N = srcPos.n() * (this->_knl_mm).srcDOF();
    de2gc.resize(M,N);
    iC( this->_knl_mm.kernel(srcPos, srcPos, trgPos, de2gc) );
  }
  if (eval) { iC( dgemv(de2gc, srcDen, trgVal) ); }

  return(0);
}

/*
  template <class VF>
int VFMM3d<VF>::DwnEqu2GrdChk_dgemv(const int dep, const DblNumVec& srcDen, DblNumVec& trgVal){
  // Not fixed for scale-varance, but will have to calc so few levels anyway 
  DblNumMat& de2gc  = (this->matmgnt())->hom() ? _dwnEqu2GrdChk[0] : _dwnEqu2GrdChk[dep];
  double R = (this->matmgnt())->hom() ? vlet()->radius()/pow2(dep) : 1.0;

  if (de2gc.m() == 0){
  double rad = ((this->matmgnt())->hom() ? vlet()->radius() : vlet()->radius()/double(pow2(dep)));
	 DblNumMat srcPos; iC( (this->matmgnt())->localPos(DE, Point3(0.0), rad, srcPos) );
	 DblNumMat trgPos((this->dim()), vlet()->trgGrdSze());
	 daxpy(rad, grdTrgSamPos(), trgPos);
	 int M = trgPos.n() * _knl_mm.trgDOF();
	 int N = srcPos.n() * _knl_mm.srcDOF();
	 de2gc.resize(M,N);
	 iC( this->_knl_mm.kernel(srcPos, srcPos, trgPos, de2gc) );
  }
  iC( dgemv(1.0/R, de2gc, srcDen, 1.0, trgVal) );
  return(0);
}
*/
/*
  template <class VF>
int VFMM3d<VF>::DwnEqu2GrdChk_dgemv(const int dep, const DblNumVec& srcDen, DblNumVec& trgVal){
  // Not fixed for scale-varance, but will have to calc so few levels anyway 
  
  DblNumMat& de2gc  = (((this->matmgnt())->hom() == true) ? _dwnEqu2GrdChk[0] : _dwnEqu2GrdChk[dep+this->rootLevel()] );
  double R = 1.0/pow(2.0, dep + this->rootLevel());

  if (de2gc.m() == 0){
	 double rad = (((this->matmgnt())->hom() == true) ? vlet()->radius() : vlet()->radius()/double(pow2(dep)));
	 DblNumMat srcPos; iC( (this->matmgnt())->localPos(DE, Point3(0.0), vlet()->radius(), srcPos) );

	 double R = (((this->matmgnt())->hom() == true) ? 1.0 : 1.0/pow(2.0, dep + this->rootLevel()));
	 DblNumMat trgPos((this->dim()), vlet()->trgGrdSze());
	 daxpy(R, grdTrgSamPos(), trgPos);
	 //DblNumMat trgPos(grdTrgSamPos());
	 int MM = trgPos.n() * _knl_mm.trgDOF();
	 int NN = srcPos.n() * _knl_mm.srcDOF();
	 de2gc.resize(MM,NN);
	 iC( this->_knl_mm.kernel(srcPos, srcPos, trgPos, de2gc) );
  }
  iC( dgemv(1.0/R, de2gc, srcDen, 1.0, trgVal) );
  return(0);
}
*/

template <class VF>
int VFMM3d<VF>::WXUpwEqu2DwnChk_dgemv(const int type, const int gNodeIdx, const int gNodeIdx_nbr, const int nbrType, const DblNumVec& sden, DblNumVec& tval, const Point3 offset, bool perdir, int fliptype, bool fullcmpt)
{
  int dep = vlet()->depth(gNodeIdx);
  map<int, map<int, DblNumMat> >& _UE2DC_type = ((type == WNODE) ? _WUpwEqu2DwnChk[dep] : _XUpwEqu2DwnChk[dep]);
  map<int, DblNumMat>& _UE2DC_lev = _UE2DC_type[abs(vlet()->depth(gNodeIdx) - vlet()->depth(gNodeIdx_nbr))];
  DblNumMat& _UE2DC =  _UE2DC_lev[nbrType];

  //if (perdir) { _UE2DC.resize(0,0); } // for now while debugging

  //DblNumMat _UE2DC(0,0);
  
#pragma omp critical	
  if (_UE2DC.m() == 0){
    /* src here is the nbr in the original list, so apply the offset to it */
    Point3 srcCtr(vlet()->center(gNodeIdx));
	 if (vlet()->periodic() && perdir == true){
		if (vlet()->dirichlet()) { vlet()->flpCtrDirNbr(srcCtr, fliptype); }
		for (int d = 0; d < (this->dim()); d++) { srcCtr(d) += offset(d); }
	 }
    Point3 trgCtr(vlet()->center(gNodeIdx_nbr));
	 
    double srcRad = (vlet()->radius(gNodeIdx));
    double trgRad = (vlet()->radius(gNodeIdx_nbr));
    
    DblNumMat srcPos; iC( (this->matmgnt())->localPos(UE, srcCtr, srcRad, srcPos) );
    DblNumMat trgPos; //iC( (this->matmgnt())->localPos(DC, trgCtr, trgRad, trgPos) );
	 
    if (type == XNODE){
		if ((this->matmgnt())->samPos(DC).n() < vlet()->grdTrgSamPos().n()){
		  iC( (this->matmgnt())->localPos(DC, trgCtr, trgRad, trgPos) );
		}
		else{ 
		  //As accurate and may be faster to go to go points from W->X
		  trgPos = vlet()->grdTrgExaPos(gNodeIdx_nbr);
		}
    }
    else {
		//Eventually need to go from coeffs->dwnChk, but for X->W, this seems
		//to be just as accuarte
      iC( (this->matmgnt())->localPos(DC, trgCtr, trgRad, trgPos) );
    }
    
    int MM = trgPos.n() * (this->_knl_mm).trgDOF();
	 int NN = srcPos.n() * (this->_knl_mm).srcDOF();
    
    _UE2DC.resize(MM,NN);
    iC( this->_knl_mm.kernel(srcPos, srcPos, trgPos, _UE2DC) );
  }	

  
  iA( fullcmpt == 0 || fullcmpt == 1);
  if (fullcmpt){
    DblNumVec srcDen(sden.m() == 1 ? this->srcUpwEquDen(gNodeIdx) : sden);

    bool sampos = (_UE2DC.m() == vlet()->grdTrgSamPos().n()*(this->_knl_mm).trgDOF() ? true : false);
    DblNumVec trgVal(tval.m() == 1 ? (((type == XNODE && sampos) ? grdExaVal(gNodeIdx_nbr) : this->trgDwnChkVal(gNodeIdx_nbr))) : tval);
	 
    iA(trgVal.m() == _UE2DC.m() && srcDen.m() == _UE2DC.n());
	 
    if (perdir){
      DblNumVec srcUpwEquDenFlp((this->matmgnt())->plnDatSze(UE));
      DblNumMat& _DEM = ((this->matmgnt())->perdirmaps())->dwnEquGrdMap()[fliptype];
      iC( dgemv(_DEM, srcDen, srcUpwEquDenFlp) );
      iC( dgemv(_UE2DC, srcUpwEquDenFlp, trgVal) );
    }
    else {
		iC( dgemv(_UE2DC, srcDen, trgVal) );
    }
	 
  }
  return 0;
}

template <class VF>
int VFMM3d<VF>::WXcoeffs2TrgVal(const int level, const int type, const int nbr, const DblNumVec& srcCoeffs, DblNumVec& trgVal, bool eval){
  iA(type == WNODE || type == XNODE);
  map<int, map<int, DblNumMat> >& _WXC2TV_levMap =  ((this->matmgnt())->hom()) ? _WXcoeffs2TrgVal[0] : _WXcoeffs2TrgVal[level];
  map<int, DblNumMat>& _WXC2TV_type = _WXC2TV_levMap[type];
  DblNumMat& _WXC2TV =  _WXC2TV_type[nbr];
 
  if (_WXC2TV.m() == 0){
	 int NK = vlet()->srcNk();
	 _WXC2TV.resize(vlet()->trgGrdSze()*(this->srcDOF()),NK*(this->trgDOF()));
	 int bref, symNum, posneg, pnt;
	 Index3 refIdx;
	 Index3 cofIdx;
	 int KT = vlet()->kTrgVal();
	 iC( WXsymNodeGetRefs(nbr, bref, symNum, posneg, refIdx, cofIdx)); 
	 for (int k = 0; k < KT; k++){ for (int j = 0; j < KT; j++){ for (int i = 0; i < KT; i++){
			 pnt = (k*KT + j)*KT + i;
			 for (int c = 0; c < NK; c++){ for (int ds = 0; ds < (this->srcDOF()); ds++){ for (int dt = 0; dt < (this->trgDOF()); dt++){
					 _WXC2TV(pnt*(this->srcDOF()) + ds,c*(this->trgDOF()) + dt) = compSymNbrF(type, level, bref, symNum, posneg, refIdx, i, j, k, c, ds, dt, KT, refIdx, cofIdx);
				  } } }
		  } } }
  }
  iA(_WXC2TV.m() != 0);
  if (eval) { iC( dgemv(_WXC2TV, srcCoeffs, trgVal) ); }
  return(0);
}

// ----------------------------------------------------------------------

/* For use with tables */
template <class VF>
DblNumVec VFMM3d<VF>::srcCoeffs(int gNodeIdx)
{
  VF& node=vlet()->node(gNodeIdx);
  int NK = vlet()->srcNk();
  int num = NK;
  int beg = node.termidx()*num;
  return DblNumVec((this->srcDOF()) * NK, false, _coeffs->data()+beg*(this->srcDOF()));
}

template <class VF>
DblNumVec VFMM3d<VF>::grdExaVal(int gNodeIdx)
{
  VF& node=vlet()->node(gNodeIdx);
  int beg = node.termidx()*vlet()->trgGrdSze();
  int num = vlet()->trgGrdSze();
  return DblNumVec((this->trgDOF())*num, false, _grdExaVal.data()+beg*(this->trgDOF()));
}

/* Effectively deletes this data */
template <class VF>
int VFMM3d<VF>::cleanSrcEqu2UpwChkTbls(){
  int total = ((this->matmgnt())->hom()) ? 1 : vlet()->maxLevel();
  double num;
  for (int t = 0; t < total; t++){
	 map<int, DblNumMat>& _SC2TV_map = _nrNbrPreCompF[t];
	 int mapSze = 10;
	 for (int i = 0; i < mapSze; i++){
		DblNumMat& _SC2TV = _SC2TV_map[i];
		_SC2TV.resize(0,0);
	 }
	 _SC2TV_map.clear();
	 map<int, map<int, DblNumMat> >& _SC2TV_levMap =  ((this->matmgnt())->hom()) ? _srcCof2TrgVal[0] : _srcCof2TrgVal[t];
	 for (int type = 0; type < 3; type++){
		map<int, DblNumMat>& _SC2TV_type = _SC2TV_levMap[type];
		int nbrt = 0;
		if (type == NORM) nbrt = 27;
		else nbrt = 56;
		for (int j = 0; j < nbrt; j++){
		  DblNumMat& _SC2TV =  _SC2TV_type[j];
		  _SC2TV.resize(0,0);
		}
		_SC2TV_type.clear();
	 }
	 _SC2TV_levMap.clear();
  }
  _nrNbrPreCompF.clear();
  
  for (int t = 0; t < total; t++){
	 map<int, DblNumMat>& _WXC2TV_map = _wxNbrPreCompF[t];
	 int mapSze = 12;
	 for (int i = 0; i < mapSze; i++){
		DblNumMat& _WXC2TV = _WXC2TV_map[i];
		_WXC2TV.resize(0,0);
	 }
	 _WXC2TV_map.clear();
	 map<int, map<int, DblNumMat> >& _WXC2TV_levMap =  ((this->matmgnt())->hom()) ? _WXcoeffs2TrgVal[0] : _WXcoeffs2TrgVal[t];
	 for (int type = 3; type < 5; type++){
		map<int, DblNumMat>& _WXC2TV_type = _WXC2TV_levMap[type];
		for (int nbr = 0; nbr < 152; nbr++){
		  DblNumMat& _WXC2TV =  _WXC2TV_type[nbr];
		  _WXC2TV.resize(0,0);
		}
		_WXC2TV_type.clear();
	 }
	 _WXC2TV_levMap.clear();
  }
  _wxNbrPreCompF.clear();
  
  return(0);
}

/*******************************************
 * The following are for symmetries - formerly in VFMM3d_NBRF.cpp file
 *
 *******************************************/
template <class VF>
double VFMM3d<VF>::compSymNbrF(const int type, const int level, const int bref, const int symNum, const int posneg, Index3 &idx, const int i, const int j, const int k, const int c, int ds, int dt, const int base, const Index3 refIdx, const Index3 cofIdx){
  int KVAL = (vlet()->kTrgVal());
  int NK = (vlet()->srcNk());
  iA((0 <= i && KVAL >= i) && (0 <= j && KVAL >= j) && (0 <= k && KVAL >= k));
  iA(0 <= c && NK >= c);
  Index3 pnt = pntRef(idx, i, j, k, base);
  double val = 0.0;
  /* For debugging lookup */
  //int lkp = NK*ds*(this->srcDOF()) + (c*(this->srcDOF()) + dt);
  //cout << "LOKKUP_ = " << lkp << endl;
  int bsr = basSym(bref,c,this->knl().kernelType());
  int dsnew = 0, dtnew = 0;
  if ((this->srcDOF()) > 1){ dsnew = srcTrgSym(bref, ds, this->knl().kernelType()); }
  if ((this->srcDOF()) > 1){ dtnew = srcTrgSym(bref, dt, this->knl().kernelType()); }
  int pnt_ref = (pnt(2)*KVAL + pnt(1))*KVAL + pnt(0);
  int symref = symRef(type, symNum);

  if (type == FINE || type == CRSE || type == NORM){
	 val = NbrF(type, level, symref, pnt_ref, bsr, dsnew, dtnew);
  }
  else if (type == WNODE || type == XNODE){
  	 val = NbrF(type, level, symref, pnt_ref, bsr, dsnew, dtnew);
  }
  else if (type == -1){ /* This option is for S2L for X nodes wherer val is retrieved elsewhere */
	 val = 1.0;
  }
  else { iA(0); }
  val *= posNeg(posneg,c,this->knl().kernelType(), ds, dt);
  return val;
}

template <class VF>
double VFMM3d<VF>::NbrF(const int type, const int level, const int nbr, const int pnt, const int bas, int ds, int dt){
  int NK = vlet()->srcNk();
  iA(type == NORM || type == FINE || type == CRSE || type == WNODE || type == XNODE);
  if (type == NORM) { iA(0 <= nbr < 27); }
  else if (type == FINE || type == CRSE) { iA(0 <= nbr < 56); }
  else { iA( 0 <= nbr < 152); }
  iA( 0<= pnt < vlet()->trgGrdSze()); iA( 0<= bas < vlet()->srcNk());
  iA(level >= 0); /* This is necessary for now - should not be a problem except for large domains */
  double val = -1;

  /* Only Symmetry for now */
  if (type == NORM || type == FINE || type == CRSE) {
	 map<int, DblNumMat>& _SC2TV_map  = ((this->matmgnt())->hom()) ?  _nrNbrPreCompF[0] :  _nrNbrPreCompF[level];
	 DblNumMat& _SC2TV = _SC2TV_map[nbr];	
	 if (_SC2TV.m() == 0){
		iC( _tbls->eval_nbrs(level,_SC2TV_map));
	 }
	 val = _SC2TV(pnt*(this->srcDOF()) + ds, bas*(this->trgDOF()) + dt);
  }
  else if (type == WNODE || type == XNODE) {
	 int nbr_adj = nbr;
	 if (type == XNODE){ nbr_adj += 6; } /* Need to correct this to be less hacky */
	 map<int, DblNumMat>& _WXC2TV_map  = ((this->matmgnt())->hom()) ?  _wxNbrPreCompF[0] :  _wxNbrPreCompF[level];
	 DblNumMat& _WXC2TV = _WXC2TV_map[nbr_adj];
	 
	 if (_WXC2TV.m() == 0){
		iC( _tbls->eval_wx(level,_WXC2TV_map));
	 }
	 
	 val = _WXC2TV(pnt*(this->srcDOF()) + ds,bas*(this->trgDOF()) + dt);
  }
  return val;
}

template <class VF>
int VFMM3d<VF>::nbrType(int me, int you)
{
  vector<VF>&  nodeVec = vlet()->nodeVec();

  //iA( nodeVec[me].tag() & LET_SRCNODE);
  //iA( nodeVec[you].tag() & LET_TRGNODE); 

  int dep = nodeVec[me].depth() - nodeVec[you].depth();
  if (dep == 0) return NORM;
  else if (dep <= -1) return FINE;
  else if (dep >= 1) return CRSE;
  else {
	 std::cerr <<  "Node #1 location = " << vlet()->center(me) << " and Node #2 location = " <<  vlet()->center(you) << endl;
	 std::cerr <<  "Node #1 depth = " << nodeVec[me].depth() << " and Node #2 depth = " <<  nodeVec[you].depth() << endl;
	 std::cerr << "Tree restriction problem with Node # " << me << " with Neighbor # " << you << endl;
	 iA(0);
  }
}

template <class VF>
int VFMM3d<VF>::tabNbrType(int nbrType, int cur, int nbr){
  iA( nbrType == NORM || nbrType == FINE || nbrType == CRSE);
  if (nbrType == NORM) return nrmNbrType(cur, nbr);
  else if (nbrType == FINE) return fnNbrType(cur, nbr);
  else return crsNbrType(cur, nbr);
}

template <class VF>
int VFMM3d<VF>::genNbrTypLsts()
{	
  OffTns<int>& nrmNbrs = _nrmNbrs;
  if (nrmNbrs.m() == 0){
	 nrmNbrs.resize(4,4,4,-2,-2,-2);
	 for (int i = -1; i <= 1; i++){ for (int j = -1; j <= 1; j++){ for (int k = -1; k <= 1; k++){
			 int& nrmNbrsii = nrmNbrs(i,j,k);
			 nrmNbrsii = nrmNbrByDif(Index3(i,j,k));
		  } } }
  }

  OffTns< NumTns<int> >& fnNbrs = _fnNbrs;
  if (fnNbrs.m() == 0){
	 fnNbrs.resize(4,4,4,-2,-2,-2);
	 for (int i = -1; i <= 1; i++){ for (int j = -1; j <= 1; j++){ for (int k = -1; k <= 1; k++){
			 NumTns<int>& fnNbrsii = fnNbrs(i, j, k);
			 fnNbrsii.resize(2,2,2);
			 for (int l = 0; l <= 1; l++){ for (int m = 0; m <= 1; m++){ for (int n = 0; n <= 1; n++){	
					 int& fnNbrsVal = fnNbrsii(l, m, n);
					 fnNbrsVal = fnNbrByDif(Index3(i,j,k), Index3(l,m,n));
				  } } }
		  } } }
  }

  for (int nn = NRM0; nn <= NRM26; nn++){
	 NumTns<int>& Wnodes = _Wnodes[nn];
	 if (Wnodes.m() == 0){
		Wnodes.resize(2,2,2);
		for (int i = 0; i <= 1; i++){ for (int j = 0; j <= 1; j++){ for (int k = 0; k <= 1; k++){
				int& Wnodeii = Wnodes(i,j,k);
				Wnodeii = WNodeByDif(nn,Index3(i,j,k));
			 } } }
	 }
  }
  return (0);
}

  

// ----------------------------------------------------------------------
template <class VF>
int VFMM3d<VF>::nrmNbrType(int cur, int nbr)
{
  vector<VF>&  nodeVec = vlet()->nodeVec();
  iA( nbrType(cur, nbr) == NORM);
  Index3 dif(nodeVec[cur].path2Node() - nodeVec[nbr].path2Node());
  for (int i = 0; i < 3; i++) iA ( dif(i) == 0 || dif(i) == 1 || dif(i) == -1 );

  double D = 2.0 * vlet()->radius(cur);
  Point3 ctrCur = vlet()->center(cur);
  Point3 ctrNbr = vlet()->center(nbr);
  Index3 dif2;
  for(int d=0; d<(this->dim()); d++){
	 dif2(d) = int(round( (ctrCur(d) - ctrNbr(d))/D));
  }

  iA( (dif - dif2) == Index3(0,0,0));
  
  OffTns<int>& nrmNbrs = _nrmNbrs;
  int& nrmNbrsii = nrmNbrs(dif(0), dif(1), dif(2));
  int val = nrmNbrsii;
  iA(val != -1);
  return val;
}

// ----------------------------------------------------------------------
template <class VF>
int VFMM3d<VF>::fnNbrType(int cur, int nbr)
{
  vector<VF>&  nodeVec = vlet()->nodeVec();
  iA( nbrType(cur, nbr) == FINE);
  Index3 dif(nodeVec[cur].path2Node() - nodeVec[nodeVec[nbr].parent()].path2Node());
  Index3 pDif(nodeVec[nbr].path2Node() - 2*nodeVec[nodeVec[nbr].parent()].path2Node());
  for (int i = 0; i < 3; i++) iA ( dif(i) == 0 || dif(i) == 1 || dif(i) == -1 );
  
  OffTns< NumTns<int> >& fnNbrs = _fnNbrs;
  NumTns<int>& fnNbrsii = fnNbrs(dif(0), dif(1), dif(2));
  int& fnNbrsVal = fnNbrsii(pDif(0), pDif(1), pDif(2));
  int val = fnNbrsVal;
  iA( val != -1);

  return val;
}

// ----------------------------------------------------------------------
/* fnNbrType and crsNbrType have an inverse relationship,
 * so call fnNbrType with nodes reversed an take the
 * inverse of the result
 */
template <class VF>
int VFMM3d<VF>::crsNbrType(int cur, int nbr)
{
  int t = fnNbrType(nbr, cur);
  return 55 - t;
}

/* Wnodes and Xnodes have inverse relationships */
template <class VF>
int VFMM3d<VF>::XNodeType(int cur, int xn){
  int t = WNodeType(xn, cur);
  return 151 - t;
  //return t;
}

template <class VF>
int VFMM3d<VF>::WNodeType(int cur, int wn)
{
  vector<VF>&  nodeVec = vlet()->nodeVec();
  int pwn = vlet()->parent(wn);
  //std::cout << wn << " " << pwn << " " << cur << " " << vlet()->radius(cur) << " " << depth(pwn) << " " << depth(cur) << " " << vlet()->center(cur) << " " << vlet()->center(pwn) << 	endl;
  //std::cout << nbrType(cur, pwn) << endl;
  iA( nbrType(cur, pwn) == NORM);
  Index3 pDif(nodeVec[wn].path2Node() - 2*nodeVec[pwn].path2Node());
  int pNrmNum = nrmNbrType(cur, pwn);
  NumTns<int>& Wnodes = _Wnodes[pNrmNum];
  iA (Wnodes.m() != 0);
  int& Wnodeii = Wnodes(pDif(0), pDif(1), pDif(2));
  int val = Wnodeii;
  iA(val != -1);
  return val;
}


// PERIODIC CODE Type Specifications
template <class VF>
int VFMM3d<VF>::perNrmNbrType(const int cur, const int nbr, const Point3 curOffset, const int type){
  //Only for nrm nbrs of same level!
  iA(vlet()->radius(cur) == vlet()->radius(nbr));
  double D = 2.0 * vlet()->radius(cur);
  Point3 ctrCur = vlet()->center(cur);
  if (vlet()->dirichlet() && type != FLN) vlet()->flpCtrDirNbr(ctrCur, type);
  for (int d = 0; d < (this->dim()); d++){ ctrCur(d) += curOffset(d); }
  Point3 ctrNbr = vlet()->center(nbr);
  Index3 dif;
  /*CHECK ME!!*/
  for(int d=0; d<(this->dim()); d++){
	 dif(d) = int(round( (ctrCur(d) - ctrNbr(d))/D));
  }
  //cout << dif << endl;
  for (int i = 0; i < 3; i++) iA ( dif(i) == 0 || dif(i) == 1 || dif(i) == -1 );

  OffTns<int>& nrmNbrs = _nrmNbrs;
  if (nrmNbrs.m() == 0){
	 nrmNbrs.resize(4,4,4,-2,-2,-2);
	 for (int i = -1; i <= 1; i++){
		for (int j = -1; j <= 1; j++){
		  for (int k = -1; k <= 1; k++){
			 int& nrmNbrsii = nrmNbrs(i,j,k);
			 nrmNbrsii = -1;
		  }
		}
	 }
  }

  int& nrmNbrsii = nrmNbrs(dif(0), dif(1), dif(2));
  int val = nrmNbrsii;

  if (val == -1){
	 val = nrmNbrByDif(dif);
	 nrmNbrsii = val;
  }
  iA(val != -1);

  //cout << ctrCur << " " << ctrNbr << " " << dif << " " << val << endl;
  //exit(0);
  
  return val;

}		

/* Most of this code is repeated in let3d_tbls_nbrs.cpp, so I need to combine them */
template <class VF>
int VFMM3d<VF>::perFnNbrType(const int cur, const int nbr, const Point3 offset, const int type){
  //cout << depth(cur) << " " << depth(nbr) << endl;
  //cout << vlet()->radius(cur) << " " << vlet()->radius(nbr) << endl;
  iA( nbrType(cur, nbr) == FINE);
  iA(vlet()->depth(cur) - vlet()->depth(nbr) == -1); 
  vector<VF>& nodeVec = vlet()->nodeVec();

  double D = 2.0 * vlet()->radius(cur);
  int par = nodeVec[nbr].parent();
  Point3 ctrCur(vlet()->center(cur));
  Point3 ctrNbrPar = vlet()->center(par);
  if (vlet()->dirichlet() && type != FLN) vlet()->flpCtrDirNbr(ctrCur, type);
  //if (!(vlet()->dirichlet()) && periodic()) {
  for (int d = 0; d < (this->dim()); d++){ ctrCur(d) += offset(d); }
  Index3 dif;
  for(int d=0; d<(this->dim()); d++){
	 dif(d) = int(round( (-ctrNbrPar(d) + ctrCur(d))/D));
  }
  
  for (int i = 0; i < 3; i++) iA ( dif(i) == 0 || dif(i) == 1 || dif(i) == -1 );

  Index3 pDif(nodeVec[nbr].path2Node() - 2*nodeVec[par].path2Node());
  
  OffTns< NumTns<int> >& fnNbrs = _fnNbrs;
  if (fnNbrs.m() == 0){
	 fnNbrs.resize(4,4,4,-2,-2,-2);
	 for (int i = -1; i <= 1; i++){
		for (int j = -1; j <= 1; j++){
		  for (int k = -1; k <= 1; k++){
			 NumTns<int>& fnNbrsii = fnNbrs(i, j, k);
			 fnNbrsii.resize(2,2,2);
			 for (int l = 0; l <= 1; l++){
				for (int m = 0; m <= 1; m++){
				  for (int n = 0; n <= 1; n++){	
					 int& fnNbrsVal = fnNbrsii(l, m, n);
					 fnNbrsVal = -1;
				  }
				}
			 }
		  }
		}
	 }
  }

  NumTns<int>& fnNbrsii = fnNbrs(dif(0), dif(1), dif(2));
  int& fnNbrsVal = fnNbrsii(pDif(0), pDif(1), pDif(2));
  int val = fnNbrsVal;

  if (val == -1){
	 val = fnNbrByDif(dif, pDif);
	 fnNbrsVal = val;
  }
  iA(val != -1);
  return val;
}

template <class VF>
int VFMM3d<VF>::perCrsNbrType(const int cur, const int nbr, const Point3 offset, const int type)
{

  iA( nbrType(cur, nbr) == CRSE);
  iA(vlet()->depth(cur) - vlet()->depth(nbr) == 1); 
  vector<VF>& nodeVec = vlet()->nodeVec();
  double D = 2.0 * vlet()->radius(nbr);
  int par = nodeVec[cur].parent();
  Point3 ctrCurPar = vlet()->center(par);
  //cout << type << endl;
  if (vlet()->dirichlet() && type != FLN) vlet()->flpCtrDirNbr(ctrCurPar, type);
  for (int d = 0; d < (this->dim()); d++){ ctrCurPar(d) += offset(d); }
  Point3 ctrNbr = vlet()->center(nbr);
  Index3 dif;
  for(int d=0; d<(this->dim()); d++){
	 dif(d) = int(round( (ctrNbr(d) - ctrCurPar(d))/D));
  }
  /* Need to negate dif here - not anympre*/
  //dif = -dif;
  for (int i = 0; i < 3; i++) iA ( dif(i) == 0 || dif(i) == 1 || dif(i) == -1 );
  Index3 pDif(nodeVec[cur].path2Node() - 2*nodeVec[par].path2Node());
  if (vlet()->dirichlet()) vlet()->flpPDif(pDif, type);
  
  /* Adjust for COARSE: */
  int val = 55 - fnNbrByDif(dif, pDif);
  iA(val != -1);
  return val;
}

template <class VF>
int VFMM3d<VF>::perWnodeType(const int cur, const int wn, const Point3 offset, const int type)
{
  vector<VF>&  nodeVec = vlet()->nodeVec();
  int pwn = vlet()->parent(wn);
  iA( nbrType(cur, pwn) == NORM);
  Index3 pDif(nodeVec[wn].path2Node() - 2*nodeVec[pwn].path2Node());
  Point3 ctrPwn = vlet()->center(pwn);
  Point3 ctrCur = vlet()->center(cur);
  if (vlet()->dirichlet() && type != FLN) vlet()->flpCtrDirNbr(ctrCur, type);
  for (int d = 0; d < (this->dim()); d++){ ctrCur(d) += offset(d); }
  double D = 2.0 * vlet()->radius(cur);
  Index3 difint;
  for(int d=0; d<(this->dim()); d++){
	 difint(d) = int(round( (ctrCur(d) -  ctrPwn(d))/D));
  }
  
  OffTns<int>& nrmNbrs = _nrmNbrs;
  int& nrmNbrsii = nrmNbrs(difint(0), difint(1), difint(2));
  int pNrmNum = nrmNbrsii;

  NumTns<int>& Wnodes = _Wnodes[pNrmNum];
  if (Wnodes.m() == 0){
	 Wnodes.resize(2,2,2);
	 for (int i = 0; i <= 1; i++){
		for (int j = 0; j <= 1; j++){
		  for (int k = 0; k <= 1; k++){
			 int& Wnodeii = Wnodes(i,j,k);
			 Wnodeii = -1;
		  }
		}
	 }
  }
  int& Wnodeii = Wnodes(pDif(0), pDif(1), pDif(2));
  int val = Wnodeii;
    
  if(val == -1){
	 val =  WNodeByDif(pNrmNum, pDif);
	 Wnodeii = val;
  }
  iA(val != -1);
  return val;
}

template <class VF>
int VFMM3d<VF>::perXnodeType(const int cur, const int xn, const Point3 offset, const int type)
{
  iA( nbrType(vlet()->parent(cur), xn) == NORM);
  vector<VF>&  nodeVec = vlet()->nodeVec();
  double D = 2.0 * vlet()->radius(xn);
  int par = nodeVec[cur].parent();
  Point3 ctrCurPar = vlet()->center(par);
  if (vlet()->dirichlet() && type != FLN) vlet()->flpCtrDirNbr(ctrCurPar, type);
  for (int d = 0; d < (this->dim()); d++){ ctrCurPar(d) += offset(d); }
  Point3 ctrX = vlet()->center(xn);
  Index3 dif;
  /* This is reversed to conform with the way it is being done in general
	* The problem was with the periodicity */
  for(int d=0; d<(this->dim()); d++){
	 dif(d) = int(round( (ctrX(d) - ctrCurPar(d))/D));
  }
  
  Index3 pDif(nodeVec[cur].path2Node() - 2*nodeVec[par].path2Node());
  if (vlet()->dirichlet()) vlet()->flpPDif(pDif, type); 

  
  OffTns<int>& nrmNbrs = _nrmNbrs;
  iA (nrmNbrs.m() != 0);
  int& nrmNbrsii = nrmNbrs(dif(0), dif(1), dif(2));
  int pNrmNum = nrmNbrsii;
  
  if (pNrmNum == -1){
	 pNrmNum = nrmNbrByDif(dif);
	 nrmNbrsii = pNrmNum;
  }
  iA( pNrmNum != -1);
  
  int val =  151 - WNodeByDif(pNrmNum, pDif);
  return val;
}

/* Type refers to the specific type of Dirichlet boundary reflection.  These follow the same
 * rules as the posNeg array in Syms.hpp.  Specifically:
 * type 0 = ( X, Y, Z), type 1 = ( X, Y,-Z), type 2 = ( X,-Y, Z)
 * type 3 = ( X,-Y,-Z), type 4 = (-X, Y, Z), type 5 = (-X, Y,-Z)
 * type 6 = (-X,-Y, Z), type 7 = (-X,-Y,-Z)
 * We use posneg to change the sign of the coefficients only; however,
 * due to the fact that we are reflecting across boundaries as well as the sign
 * of the forces in some cases, we must negate the changes (hence, the "flip")
 * Specifically, across bdry types 1, 2, 4, and 7, we negate the posneg results
 * For periodic conditions, type should only be zero always
 */
template <class VF>
int VFMM3d<VF>::reflectCoeffs(const int type, const DblNumVec& incoeffs, DblNumVec& outcoeffs){
  iA( type >= 0 && type <= 7 );
  int kt = 111; /* kt fixed at 111 as we cycle through the coefficients */
  double flip = (type == 1 || type == 2 || type ==4 || type == 7) ? -1.0 : 1.0;
  for (int c = 0; c < vlet()->srcNk(); c++){
	 for (int ds = 0; ds < (this->srcDOF()); ds++){
		int cnt = c*(this->srcDOF()) + ds;
		outcoeffs(cnt) = flip*(incoeffs(cnt)*posNeg(type,c,kt,0,0));
	 }
  }
  return(0);
}

template <class VF>
int VFMM3d<VF>::PerCleanup(){
  
  int l = 0; /* Only one level for now - otherwise iterate */
  NumTns<DblNumMat>& _UE2UC = (this->knl().homogeneous()==true ) ? ((this->matmgnt())->perdirmaps())->PerUpwEqu2UpwChk()[0] : ((this->matmgnt())->perdirmaps())->PerUpwEqu2UpwChk()[l];
  if(_UE2UC.m()==0)	 _UE2UC.resize(3,3,3);
  for (int i = 0; i <= 2; i++){
	 for (int j = 0; j <= 2; j++){
		for (int k = 0; k <= 2; k++){
		  DblNumMat& _UE2UCii = _UE2UC(i,j,k);
		  _UE2UCii.resize(0,0);
		}
	 }
  }
  _UE2UC.resize(0,0,0);

  DblNumMat& _UC2UE = (this->knl().homogeneous()==true ) ? ((this->matmgnt())->perdirmaps())->PerUpwChk2UpwEqu()[0] : ((this->matmgnt())->perdirmaps())->PerUpwChk2UpwEqu()[l];
  _UC2UE.resize(0,0);

  
  NumTns<DblNumMat>& _DE2DC = (this->knl().homogeneous()==true ) ? ((this->matmgnt())->perdirmaps())->PerDwnEqu2DwnChk()[0] : ((this->matmgnt())->perdirmaps())->PerDwnEqu2DwnChk()[l];
  DblNumMat& _DE2DCii = _DE2DC(0, 0, 0);
  _DE2DCii.resize(0,0);
  _DE2DC.resize(0,0,0);


  DblNumMat& _DC2DE = (this->knl().homogeneous()==true ) ? ((this->matmgnt())->perdirmaps())->PerDwnChk2DwnEqu()[0]: ((this->matmgnt())->perdirmaps())->PerDwnChk2DwnEqu()[l];
  _DC2DE.resize(0,0);
  return(0);
}
