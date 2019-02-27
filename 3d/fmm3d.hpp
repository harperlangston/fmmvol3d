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
#ifndef _FMM3D_HPP_
#define _FMM3D_HPP_

#include "knlmat3d.hpp"
#include "let3d.hpp"
#include "matmgnt3d.hpp"


//-------------------------------------------
//! FMM3d sequential class.
/*! FMM3d implements KnlMat3d
 */

template <class N>
class FMM3d: public KnlMat3d
{
public:
  typedef pair<int,int> intpair;
  /*! UE = Upper Equivalent
	*  UC = Upper Check
	*  DE = Downward Equivalent
	*  DC = Downward Check
	*/
  enum {	 UE=0,	 UC=1,	 DE=2,	 DC=3,  };
  //------------------------------------
    
protected:
  //PARAMS (REQ)
  /*! Center of the toplevel box */
  Point3 _center;
  /*! The level of the root box, radius of the box is 2^(-_rootlvl)*/
  int    _rootLevel;
  //PARAMS (OPT)
  int _np;
  //COMPONENTS, local member and data  vector<int> _stphds, _evlhds;
  Let3d<N>* _let;
  MatMgnt3d* _matmgnt;

  /*! Source Exact Positions */
  DblNumMat _srcExaPos;
  DblNumMat _srcExaNor;
  /*! Source Exa Densities */
  DblNumVec _srcExaDen;
  /*! Source Upward Equivalent Densities */
  DblNumVec _srcUpwEquDen;
  /*! Source Upward Check Values - no longer used */
  //DblNumVec _srcUpwChkVal;
  /*! Target Exact Positions */
  DblNumMat _trgExaPos;
  /*! Target Exact Values */
  DblNumVec _trgExaVal;
  /*! Target Downward Equivalent Densities - no longer used */
  //DblNumVec _trgDwnEquDen;
  /*! Target Downward Check Values */ 
  DblNumVec _trgDwnChkVal;
  
  Kernel3d _knl_mm; //elq used in matmgnt
  int _mul_mm; //mul used in matmgnt  //FUNCTIONS: eq, lt, qt, sdof, tdof
public:
  FMM3d(const string& p);
  ~FMM3d();

  /*! Return center of toplevel box */
  Point3& center() { return _center; }
  /*! Return the rootlevel.   2^(-rootlvl) is the radius of the toplevel box */
  int& rootLevel() { return _rootLevel; }
  int& np() { return _np; }
  
  /*! Setup function */
  virtual int setFromOptions(map<string,string>& opts);
  virtual int setup();
  /*! Evaluate - in fmm3d_eval.cpp */
  virtual int evaluate(const DblNumVec& srcDen, DblNumVec& trgVal);
  /*! Check relative error - in fmm_check.cpp */
  virtual int check(const DblNumVec& srcDen, DblNumVec& trgVal, int numChk, double& relativeErr);  
  
  Let3d<N>* let() { return _let; }
  MatMgnt3d* matmgnt() { return _matmgnt; }
  //vector<N>& nodeVec() { let()->nodeVec(); }
  N& node(int gNodeIdx) { return _let->nodeVec()[gNodeIdx]; }

public:
  /*! Data Size */
  int datSze(int tp) { return _matmgnt->plnDatSze(tp); }
  
  /*! Source Equivaluent To Target Check Multiplication */
  int SrcEqu2TrgChk_dgemv(const DblNumMat& srcPos, const DblNumMat& srcNor, const DblNumMat& trgPos, const DblNumVec& srcDen, DblNumVec& trgVal);
  /*! Source Equivalent To Upward Check Multiplication */
  int SrcEqu2UpwChk_dgemv(const DblNumMat& srcPos, const DblNumMat& srcNor, Point3 trgCtr, double trgRad, const DblNumVec& srcDen, DblNumVec& trgVal);
  /*! Source Equivalent To Downward Check Multiplcation */
  int SrcEqu2DwnChk_dgemv(const DblNumMat& srcPos, const DblNumMat& srcNor, Point3 trgCtr, double trgRad, const DblNumVec& srcDen, DblNumVec& trgVal);
  /*! Downward Equivalent To Target Check */
  int DwnEqu2TrgChk_dgemv(Point3 srcCtr, double srcRad, const DblNumMat& trgPos, const DblNumVec& srcDen, DblNumVec& trgVal);
  /* Upward Equivalent To Target Check */
  int UpwEqu2TrgChk_dgemv(Point3 srcCtr, double srcRad, const DblNumMat& trgPos, const DblNumVec& srcDen, DblNumVec& trgVal);

  int UpwEqu2DwnChk_dgemv(Point3 srcCtr, double srcRad, Point3 trgCtr, double trgRad, const DblNumVec& srcDen, DblNumVec& trgVal);

  DblNumVec& srcUpwEquDen()   { return _srcUpwEquDen;}
  DblNumVec& trgDwnChkVal()   { return _trgDwnChkVal; }
  DblNumVec& trgExaVal()   { return _trgExaVal; }
  
  //contributor data
  DblNumMat srcExaPos(int gNodeIdx);
  DblNumMat srcExaNor(int gNodeIdx);
  DblNumVec srcExaDen(int gNodeIdx);
  DblNumVec srcUpwEquDen(int gNodeIdx);
  //DblNumVec srcUpwChkVal(int gNodeIdx);
  //evaluator data
  DblNumMat trgExaPos(int gNodeIdx);
  DblNumVec trgExaVal(int gNodeIdx);
  //DblNumVec trgDwnEquDen(int gNodeIdx);
  DblNumVec trgDwnChkVal(int gNodeIdx);
  
  /*! compute source data */
  int srcData();
  /*! compute target data */
  int trgData();

};

#endif
