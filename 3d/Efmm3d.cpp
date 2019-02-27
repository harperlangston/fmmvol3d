#include "Efmm3d.hpp"

using std::istringstream;
using namespace std;


EFMM3d::EFMM3d(const string& p)
  : VFMM3d<EbiNode>(p)
{
  
}

EFMM3d::~EFMM3d()
{

}

/*
// ---------------------------------------------------------------------- 
int EFmm3d::setFromOptions(map<string,string>& optionsMap)
{
  //------------
  map<string,string>::iterator mapindex;
  mapindex = optionsMap.find("-" + prefix() + "maxLevel"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_maxLevel; }

  mapindex = optionsMap.find("-" + prefix() + "periodic"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_periodic; }
  mapindex = optionsMap.find("-" + prefix() + "dirichlet"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_dirichlet; }

  mapindex = optionsMap.find("-" + prefix() + "ksrcval"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_kSrcVal; }
  mapindex = optionsMap.find("-" + prefix() + "ktrgval"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_kTrgVal; }

  mapindex = optionsMap.find("-" + prefix() + "rhs"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_rhs; }

  mapindex = optionsMap.find("-" + prefix() + "balance"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_balance; }

  mapindex = optionsMap.find("-" + prefix() + "ptsMin"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_ptsMin; }

  mapindex = optionsMap.find("-" + prefix() + "interp"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_interp; }

  mapindex = optionsMap.find("-" + prefix() + "adaptive"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_adaptive; }

  mapindex = optionsMap.find("-" + prefix() + "ebs"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_ebs; }

  mapindex = optionsMap.find("-" + prefix() + "tolbsd"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_tolbsd; }
  return(0);
}



*/
