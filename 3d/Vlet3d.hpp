#ifndef _VLET3D_HPP_
#define _VLET3D_HPP_

#include "let3d.hpp"
#include "common/numtns.hpp"
#include "common/vecmatop.hpp"
#include <iostream>
#include "common/syms.hpp"
#include "exsol3d.hpp"

#define FRC 0
#define POT 1

/* "Flip" conditions */
#define FLN 0 /* No Flip */
#define FLZ 1 /* Flip in z-direction */
#define FLY 2 /* etc... see notes */
#define FLYZ 3
#define FLX 4
#define FLXZ 5
#define FLXY 6
#define FLXYZ 7

class VolNode: public Node {
public:
  VolNode(int p, int c, Index3 t, int d) : Node(p,c,t,d), _prmVltr(0), _errcheck(0), _scdVltr(0), _dscPrmVltr(0), _termidx(-1) {;}  
  bool& errcheck() { return _errcheck; }
  
  bool& prmVltr()          { return _prmVltr; }
  bool& scdVltr()          { return _scdVltr; }
  bool& dscPrmVltr()       { return _dscPrmVltr; } 
  
  int& termidx()        { return _termidx; }
  
protected:
  bool _errcheck;
  
  bool _grd;
  int _termidx;
  
  bool _prmVltr; /* primary violator for tree balancing - also used for tolerance testing */
  bool _scdVltr; /* secondary violator for tree balancing */
  bool _dscPrmVltr; /* Descendant of Primary Violator */
  
private:
  
};

//---------------------------------------
class BdryNode{
protected:
  int _gni;
  Point3 _bdryOffset;
  int _type; /* Only used if Dirichlet */
public:
  //BdryNode(){;}
  
  BdryNode(int gni, Point3 off):
	 _gni(gni), _bdryOffset(off) { _type = FLN; }
  
  BdryNode(int gni, Point3 off, int type):
	 _gni(gni), _bdryOffset(off), _type(type) { ; }
  
  //~BdryNode(){ ; }
  
  int gni() { return _gni; }
  Point3& bdryOffset() {  return _bdryOffset; }
  int& type() { return _type; }
};
class PerNode {
public:
  PerNode() :
	 _bdry(false){ }
  
  
  bool& bdry()          { return _bdry; }
  Point3 & bdryLoc()    { return _bdryLoc; }
  
  vector<BdryNode>& bdryNbrs() { return _bdryNbrs; }
  vector<BdryNode>& bdryUnodes() { return _bdryUnodes; }
  vector<BdryNode>& bdryVnodes() { return _bdryVnodes; }
  vector<BdryNode>& bdryWnodes() { return _bdryWnodes; }
  vector<BdryNode>& bdryXnodes() { return _bdryXnodes; }
  
  vector<int>& dscList() { return _dscList;} 
protected:
  bool _bdry;
  Point3 _bdryLoc;
  
  // These are used only if bdry is set to true 
  vector<BdryNode> _bdryNbrs;	
  vector<BdryNode> _bdryUnodes;
  vector<BdryNode> _bdryVnodes;
  vector<BdryNode> _bdryWnodes;
  vector<BdryNode> _bdryXnodes;
  vector<int> _dscList;  
};	 

template <class V>
class VLet3d: public Let3d<V>
{
public:
  //---------------------------------------
  

  //----------------------------------------------
protected:
  int _kSrcVal; //default 4
  int _kTrgVal; //default 4
  int _srcNk;
  int _trgNk;
  int _rhs;

  bool _adaptive;
  bool _balance;

  bool _periodic;
  bool _dirichlet;
 
  int _numViolators;
  int _sdof;
  
  //COMPONENTS
  int _trmnodecnt;
  int _grdNodeCnt;

  /* Used when periodic or dirichlet solver present */
  vector<PerNode> _perNodeVec;

  Exsol3d* _exsol3d;
  Kernel3d _knl; //barely used by let - just for sdof, tdof, etc.

  DblNumMat _grdSrcSamPos;
  DblNumMat _grdTrgSamPos;
  DblNumMat _grdOverSamPos;
  DblNumMat _grdDblSrcSamPos;
  
public:
  VLet3d(const string& p);
  ~VLet3d();

  //virtual int upwOrderCollect(vector<int>&);
  
  Exsol3d*& exsol3d() { return _exsol3d; }
  Kernel3d& knl() { return _knl; }

  const DblNumMat& grdOverSamPos() { return _grdOverSamPos; }
  const DblNumMat& grdSrcSamPos() { return _grdSrcSamPos; }
  const DblNumMat& grdTrgSamPos() { return _grdTrgSamPos; }
  const DblNumMat& grdDblSrcSamPos() { return _grdDblSrcSamPos; }

  DblNumMat grdOverExaPos(int gNodeIdx, bool depth=false);
  DblNumMat grdDblSrcExaPos(int gNodeIdx, bool depth=false);
  DblNumMat grdSrcExaPos(int gNodeIdx, bool depth=false);
  DblNumMat grdTrgExaPos(int gNodeIdx, bool depth=false);

  int setFromOptions(map<string,string>& optionsMap);
  int setup();

  int build(); /* separated for tolerance test purposes */
  int bal();
  int subVltrs();
  int initTree();
  int srcBalData();
  int srcTrgBalData();
  int nbrsBalBld();
  int nbrsBalBld(int g);

  bool cheb(int kval);
  
  /* Extra infor for periodic or dirichelt solver */
  vector<PerNode>& perNodeVec() { return _perNodeVec; }
  PerNode& pernode(int gNodeIdx) { return _perNodeVec[gNodeIdx]; }

  int grdNodeCnt()   { return _grdNodeCnt; }
  int grdExaCnt()  { return trmNodeCnt()*trgGrdSze(); }
  int trmNodeCnt() { return _trmnodecnt; }
    
  int tblsSrcTrgData(int req_levels);
  int tblsSrcTrgData();

  double symNbrF(int type, int bref, int symNum, int posneg, Index3 &idx, int i, int j, int k, int c, int base);
  double dpScl(int gNodeIdx, int j);

  //int srcNodeCnt()   { return Let3d::srcNodeCnt(); }  
  int& kSrcVal()      { return _kSrcVal; }
  int& kTrgVal()      { return _kTrgVal; }
  int srcGrdSze()     { return _kSrcVal*_kSrcVal*_kSrcVal; }
  int trgGrdSze()     { return _kTrgVal*_kTrgVal*_kTrgVal; }
  int& srcNk()        { return _srcNk; }
  int& trgNk()        { return _trgNk; }
  int& rhs() { return _rhs; }
  double eps_rhs() { return pow(0.1, (double)(_rhs)); }

  int termidx(int gNodeIdx) { return this->node(gNodeIdx).termidx(); } //no child

  bool& balance() { return _balance; }
  bool& adaptive() { return _adaptive; }

  bool& periodic() { return _periodic; }	
  bool& dirichlet() { return _dirichlet; }

  int& numVltrs() { return _numViolators; }
  int& sdof() { return _sdof; }

  int bldPerNbrs();
  int bldPerUWXnodes();
  int bldPerUnodes(); /* Just U nodes */
  int bldPerVnodes();
  int setBdryLoc(int gni, int chiNum);

  int bldDirNbrs();
  int bldDirVnodes();
  int flpCtrDirNbr(Point3 &ctrDirNbr, const int flpTyp);
  int flpPDif(Index3& pDif, const int flpTyp);

  int bldDscList(int gni);
  
};

#endif
