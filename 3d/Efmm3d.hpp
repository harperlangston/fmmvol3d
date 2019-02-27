#ifndef __EFMM3D_HPP__
#define __EFMM3D_HPP__

#include "Vfmm3d.hpp"

class EbiNode: public VolNode {
public:
  EbiNode(int p, int c, Index3 t, int d) : VolNode(p,c,t,d), _color(-1), _grdMrkd(false), _dist_from_bdry(-1.0), _numIntPts(0), _numIntDblPts(0), _srcnode(false) { }
  bool& bndCrv() { return _bndCrv; }
  vector<int>& bndCrvPts() { return _bndCrvPts; }
  
  char& color() { return _color; }//Indicates if in or out
  double& dist() { return _dist_from_bdry; } //closest distance to closest bdry point - set in bisdov.cpp::markgrid
  
  vector<int>& isInt() { return _isInt; }
  vector<int>& isIntDbl() { return _isIntDbl; }
  
  int& isInt(int i) { return _isInt[i]; }
  int& isIntDbl(int i) { return _isIntDbl[i]; }
  
  int& numIntPts() { return _numIntPts; }
  int& numIntDblPts() { return _numIntDblPts; }
  
  bool& grdMrkd() { return _grdMrkd; }
  bool& srcnode() { return _srcnode; }
  
protected:
  bool _bndCrv; // Whether bdry curve intersects this node - for leaves only
  // for use with the Ebi code for curves 
  
  vector<int> _bndCrvPts ; // A listing of points which are on the
  //  curve (from Ebi) which are inside of this
  // node - for leaves only - initiated outside
  // in Ebi code 	 
  char _color; /* Indicates if in or out - for Ebi code */
  
  vector<int> _isInt;
  vector<int> _isIntDbl;
  
  int _numIntPts;
  int _numIntDblPts;
  
  bool _srcnode;
  bool _grdMrkd;
  double _dist_from_bdry; //closest distance ANY corner of the box is from the closest bdry point;
  
private:
  
};


//---------------------------------------
class EFMM3d: public VFMM3d<EbiNode>
{
private:
  int _ptsMin; // Used if Vfmm3d uses interpolation 
  int _interp;
  bool _adaptive;
  bool _ebs;
  bool _tolbsd;

public:
  //---------------------------------------
  

  //----------------------------------------------
protected:
  
  
public:
  EFMM3d(const string& p);
  ~EFMM3d();

  //int setFromOptions(map<string,string>& optionsMap);
  //virtual int setup();

  bool& ebs() { return _ebs; }
  bool& tolbsd() { return _tolbsd; }
  int& ptsMin()    { return _ptsMin; }
  int& interp()    { return _interp; }

  bool& adaptive() { return _adaptive; }
  
};

#endif
