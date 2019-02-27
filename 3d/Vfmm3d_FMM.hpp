#ifndef _VFMM3D_FMM_HPP_
#define _VFMM3D_FMM_HPP_

#include "Vfmm3d.hpp"

class VFMM3d_FMM: public VFMM3d<VolNode>
{
private:
  
protected:
  FMM3d<Node>* _fmm;
  DblNumVec _grdSrcDen_FMM;
  
public:
  VFMM3d_FMM(const string& p);
  ~VFMM3d_FMM();

  FMM3d<Node>*& fmm() { return _fmm; }

  int setup();

  int evaluate(const DblNumVec& srcDen, DblNumVec& trgVal);

  int bldSrcCoeffsGrd_notused(const bool over, const int gNodeIdx);

  DblNumVec grdExaDen(int gNodeIdx);
  
protected:

  
};

#endif
