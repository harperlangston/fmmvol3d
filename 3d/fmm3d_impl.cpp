#include "Vfmm3d_check.cpp"
#include "Vfmm3d_setup.cpp"
#include "Vfmm3d_eval.cpp"
#include "Vfmm3d_Coeffs.cpp"
#include "Vfmm3d.cpp"

//#include "Efmm3d.cpp"
#include "Vlet3d.cpp"

#include "fmm3d.cpp"
#include "fmm3d_setup.cpp"
#include "fmm3d_eval.cpp"
#include "fmm3d_check.cpp"
#include "let3d.cpp"

//template class FMM3d<Node>;
template class FMM3d<VolNode>;
template class VFMM3d<VolNode>;
//template class VFMM3d<EbiNode>;

//template class Let3d<Node>;
template class Let3d<VolNode>;
//template class Let3d<EbiNode>;
template class VLet3d<VolNode>;
//template class VLet3d<EbiNode>;
