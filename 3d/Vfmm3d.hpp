#ifndef _VFMM3D_HPP_
#define _VFMM3D_HPP_

#include "common/nummat.hpp"
#include "common/numtns.hpp"
#include "common/offtns.hpp"
#include "common/vec3t.hpp"
#include "common/CmptTbls.hpp"
#include "fmm3d.hpp"
#include "Vlet3d.hpp"


template <class VF>
class VFMM3d: public FMM3d<VF>
{
private:

protected:
  CmptTbls* _tbls;
  
  /* Whether to do adaptive Gaussian quadrature at all */
  
  int _nonzero;
  int _nonzeropts;

  //Lookups
  OffTns<int> _nrmNbrs;
  OffTns< NumTns<int> > _fnNbrs;
  map<int, NumTns<int> > _Wnodes;

  double _lambda;
  int _eqnType;
  
  /*! Grid Exact Values */
  DblNumVec _grdExaVal;

  /*! Coefficients - referenced by terminal indices, which should change soon */
  DblNumVec* _coeffs;
  
  map<int, map<int, map<int, DblNumMat> > >_pinv;
  map<int, DblNumMat> _coeff2GrdVal;
  /* First is for level for scale-invariant kernels and
	* then for value of k and then nk for par-chicoeff building */
  map<int, DblNumMat> _dwnEqu2GrdChk;

  

  /*! For use with tables */
  map<int, map<int, DblNumMat> > _nrNbrPreCompF;
  map<int, map<int, DblNumMat> > _wxNbrPreCompF;
  map<int, DblNumMat> _s2mNbrF;
  //map<int, map<int, DblNumMat> > _dwnChkNbrF;

  /*! Downward Grid Check To Downward Equivalent - different than matmgnt version */
  //map<int, DblNumMat> _dwnGrdChk2DwnEqu;
  
  map<int, map<int, map<int, DblNumMat> > > _srcCof2TrgVal;
  map<int, map<int, map<int, DblNumMat> > > _FubSrcCof2TrgVal;
  map<int, map<int, map<int, DblNumMat> > > _CubSrcCof2TrgVal;

  map<int, map<int, map<int, DblNumMat> > >_WXcoeffs2TrgVal;
  //map<int, map<int, map<int, DblNumMat> > >_WXUpwEqu2DwnChk;
  map<int, map<int, map<int, DblNumMat> > >_WUpwEqu2DwnChk;
  map<int, map<int, map<int, DblNumMat> > >_XUpwEqu2DwnChk;

  //map<int, map<int, DblNumMat> > _srcCof2DwnChkVal;

public:
  VFMM3d(const string& p);
  ~VFMM3d();

  VLet3d<VF>* vlet() { return ((VLet3d<VF>*)(this->_let)); }
  CmptTbls* tbls() { return _tbls; }
  
  double compSymNbrF(const int type, const int level, const int bref, const int symNum, const int posneg, Index3 &idx, const int i, const int j, const int k, const int c, int ds, int dt, const int base, const Index3 refIdx, const Index3 cofIdx);
  double NbrF(const int type, const int level, const int nbr, const int pnt, const int bas, int ds, int dt);
  int BldNbrFonFly(const int type, const int level);

  int genNbrTypLsts();
  int nbrType(int me, int you);
  int nrmNbrType(int cur, int nbr);
  int fnNbrType(int cur, int nbr);

  int crsNbrType(int cur, int nbr);
  int tabNbrType(int nbrTyp, int cut, int nbr);
  double nbrF(int type, int nbr, int i, int j);

  int WNodeType(int cur, int wNode);
  int XNodeType(int cur, int xNode);
  double WXNodeF(int type, int nde, int i, int j);
  

  /* symmetry test */
  int symstest();
  
  DblNumVec*& coeffs() { return _coeffs; }
  /* Set coef at termIdx from termVec and 0 <= j < Nk() to c */
  int bldCoeffs(int type);
  virtual int bldSrcCoeffsGrd(const bool over, const int gNodeIdx);
  int prntCoeffs(int gNodeIdx);

  int bldBasPolyMat(const bool scale, const bool chebyshev, const int type, const int gni, const int NK, const DblNumMat pos, DblNumMat& bp, bool dep=false);
  int bldBasPolyMat(const bool scale, const bool chebyshev, const int type, const int gni, const int KVAL, const int NK, const DblNumMat pos, DblNumMat& bp, bool dep=false);
  int bldDOFMat(const DblNumMat Ain, DblNumMat& Aout, const int s, const int t);
  int pinType(const int NK);

  /* Returns the gni for a target pos based on original position ordering (before rearrangement by let) */
  //int trgPosGni(int tposidx) { iA(tposidx >= 0 && tposidx < (*trgPos()).n()); return vlet()->trgPosTrmGni()[tposidx]; }

  int potCoeffs2TrgVal();
  int potCoeffs2TrgVal(int gNodeIdx);

  int grdVals2TrgVal(DblNumVec& trgval);

  double& lambda() { return _lambda; }
  int& eqnType() { return _eqnType; }
  
  bool hom() { return (this->_matmgnt)->hom(); }

  virtual int setFromOptions(map<string,string>& opts);
  virtual int setup();
  virtual int evaluate(const DblNumVec& srcDen, DblNumVec& trgVal);
  
    
  int vcheck(double& rerr, double& inferr);
  
protected:

  int GrdEqu2UpwChkTbls_dgemv(const int level, const DblNumVec& srcoeffs, DblNumVec& tmpVal);

  int GrdEqu2UpwEquTbls_dgemv(const int level, const DblNumVec& srcoeffs, DblNumVec& tmpVal, DblNumVec& den);
  int GrdEqu2UpwEquTbls_dgemv(int level, int gNodeIdx, const DblNumVec& srcoeffs, DblNumVec& trgDen);
  
  int GrdCoeffs2TrgVal_dgemv(const int level, const int type, const int nbr, const DblNumVec& srcCoeffs, DblNumVec& trgVal);
  
  int GrdCoeffs2UnbalTrgVal_dgemv(const int type, const int gNodeIdx, const int nbr, const DblNumVec& coeffs, const Point3 offset, DblNumVec& trgval);
  int GrdCoeffs2UnbalTrgVal_dgemv(const int type, const int gNodeIdx, const int nbr, const DblNumVec& coeffs, DblNumVec& trgval) { return GrdCoeffs2UnbalTrgVal_dgemv(type, gNodeIdx, nbr, coeffs, Point3(0.0,0.0,0.0), trgval); }

  int WXUpwEqu2DwnChk_dgemv(const int type, const int gNodeIdx, const int gNodeIdx_nbr, const int nbr_type, const DblNumVec& sden, DblNumVec& tval, const Point3 offset=Point3(0.0,0.0,0.0), bool perdir=0, int fliptype=0, bool fullcmpt=true);
  int WXUpwEqu2DwnChk_dgemv(const int type, const int gNodeIdx, const int gNodeIdx_nbr, const int nbr_type, const Point3 offset=Point3(0.0,0.0,0.0), bool perdir=0, int fliptype=0, bool fullcmpt=true) { DblNumVec s(1); DblNumVec t(1); return WXUpwEqu2DwnChk_dgemv(type, gNodeIdx, gNodeIdx_nbr, nbr_type, s, t, offset, perdir, fliptype, fullcmpt); }
  int WXUpwEqu2DwnChk_dgemv(const int type, const int gNodeIdx, const int gNodeIdx_nbr, const int nbr_type, bool fullcmpt) { DblNumVec s(1); DblNumVec t(1); return WXUpwEqu2DwnChk_dgemv(type, gNodeIdx, gNodeIdx_nbr, nbr_type, s, t, Point3(0.0,0.0,0.0), false, 0, fullcmpt); }
  //  int WXUpwEqu2DwnChk_dgemv(const int type, const int gNodeIdx, const int gNodeIdx_nbr, const int nbr_type, const DblNumVec& sden, DblNumVec& tval, bool fullcmpt=true) { return  WXUpwEqu2DwnChk_dgemv(type, gNodeIdx, gNodeIdx_nbr, nbr_type, sden, tval, Point3(0.0, 0.0, 0.0), 0, 0, fullcmpt); }
  
  int DwnEqu2GrdChk_dgemv(const int dep, const DblNumVec& srcDen, DblNumVec& trgVal, bool eval=true);

  int WXcoeffs2TrgVal(const int level, const int type, const int nbr, const DblNumVec& srcCoeffs, DblNumVec& trgVal, bool eval=true);  
  int cleanSrcEqu2UpwChkTbls();

  DblNumVec srcCoeffs(int gNodeIdx);
  
  /*! compute trg data - different slightly for Vfmm id periodicity */
  int grdData();
  int grdTolData(bool descend_all);
  
  DblNumVec grdExaVal(int gNodeIdx);
   
  int S2Mclean();
  
  int bldDirMaps();

  int perNrmNbrType(const int cur, const int nbr, const Point3 offset, const int type);
  int perFnNbrType(const int cur, const int nbr, const Point3 offset, const int type);
  int perCrsNbrType(const int cur, const int nbr, const Point3 offset, const int type);
  int perWnodeType(const int cur, const int wn, const Point3 offset, const int type);
  int perXnodeType(const int cur, const int xn, const Point3 offset, const int type);

  int evaluate_far_BCs(const int gni, DblNumVec& FFtrgDwnEquDen);
  int evaluateUnodes_BCs(const int gNodeIdx);
  int evaluateVnodes_BCs(const int gNodeIdx, DblNumVec& effVal);
  int evaluateWnodes_BCs(const int gNodeIdx);
  int evaluateXnodes_BCs(const int gNodeIdx);
  int evalPerDirVNodesRoot(const int gNodeIdx);
  
  int PerUpwChk2UpwEqu_dgemv (int level,             const DblNumVec&, DblNumVec&);
  int PerUpwEqu2UpwChk_dgemv (int level, Index3 ii,  const DblNumVec&, DblNumVec&);
  int PerDwnChk2DwnEqu_dgemv (int level,             const DblNumVec&, DblNumVec&);
  int PerDwnEqu2DwnChk_dgemv (int level, Index3 ii,  const DblNumVec&, DblNumVec&);
  int PerUpwEqu2DwnChk_dgemv (int level, Index3 ii,  const DblNumVec& effDen, DblNumVec& effVal);
  int reflectCoeffs(const int type, const DblNumVec& incoeffs, DblNumVec& outcoeffs);
  int PerCleanup();
  


  static double _wsbuf[];
  
};

#endif
