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
#include <omp.h>

using std::istringstream;
using namespace std;

template <class VF>
int VFMM3d<VF>::bldCoeffs(int type){
  int sum = 0;
  int srcDOF = this->srcDOF(); int trgDOF = this->trgDOF();
  iA(type == FRC);
  map<int, vector<int> > lvlOrdVec; iC( vlet()->revLvlOrderCollect(lvlOrdVec) ); //BOTTOM UP
  vector<bool> lvsAtLev;
  vector<int> srcTrmVec;
  for (int j = 0; j < lvlOrdVec.size(); j++){
	 lvsAtLev.push_back(false);
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
	 for (int i = 0; i < thisLevBoxes.size(); i++){
		int gNodeIdx = thisLevBoxes[i];
		if (vlet()->terminal(gNodeIdx) == true) { lvsAtLev[j] = true; }
		if (vlet()->terminal(gNodeIdx) && vlet()->tag(gNodeIdx) & LET_SRCNODE) { //evaluator
		  srcTrmVec.push_back(gNodeIdx);
		}
	 }
  }

  bool oversamp = true;
  int levnum = ((this->matmgnt())->hom() ? 1 : lvlOrdVec.size());
  for (int j = 0; j < levnum; j++){
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
	 if (lvsAtLev[j] == true || (this->matmgnt())->hom()){
		DblNumMat srcPos(oversamp ? vlet()->grdOverExaPos(j,true) : vlet()->grdSrcExaPos(j,true)); /* Use oversampling to get good coeff approx */

		map<int, map<int, DblNumMat> >& pin_LevMap = _pinv[j];
		map<int, DblNumMat>& pin_map =  pin_LevMap[1];
		int NK = vlet()->srcNk();
		DblNumMat& pinDOF = pin_map[pinType(NK)];
		if (pinDOF.m() == 0){
		  DblNumMat polys(srcPos.n(), NK);
		  DblNumMat pin(NK, srcPos.n());
		  /* Build polys from depth = j */
		  iC( bldBasPolyMat((this->matmgnt())->hom(), vlet()->cheb(vlet()->kSrcVal()), FRC, j, NK, srcPos, polys, true));
		  iC( pinv(polys, 1e-12, pin));
		  pinDOF.resize(NK*srcDOF, srcPos.n() * trgDOF);
		  iC( bldDOFMat(pin, pinDOF, srcDOF, trgDOF));
		}
	 }
  }
  
#pragma omp parallel for
  for(int i=0; i<srcTrmVec.size(); i++) {
	 int gNodeIdx = srcTrmVec[i];
	 bldSrcCoeffsGrd(oversamp, gNodeIdx);
  }
	
  for (int i = 0; i < vlet()->maxLevel(); i++){
	 map<int, map<int, DblNumMat> >& pin_LevMap = ((this->matmgnt())->hom()) ? _pinv[0] : _pinv[i];
	 map<int, DblNumMat>& pin_map =  pin_LevMap[1];
	 DblNumMat& pinDOF = pin_map[pinType(vlet()->srcNk())];
	 pinDOF.resize(0,0);
  }
  //_grdOverSamPos.resize(0,0); 
  return(0);
}


/* This is for testing straight from the grid ONLY where the function input is known in functions.h */
template <class VF>
int VFMM3d<VF>::bldSrcCoeffsGrd(const bool over, const int gNodeIdx)
{
  iA(vlet()->terminal(gNodeIdx));
  int srcDOF = this->srcDOF(); int trgDOF = this->trgDOF();
  int NK = vlet()->srcNk();
 
  /* If the sources come from the forces on the GRID, then we are computing
	* a polynomial approximation to the force distribution.  Otherwise, we are computing
	* an approximation to the potential values on the GRID
	*/
  DblNumVec coeffVec(srcCoeffs(gNodeIdx));   setvalue(coeffVec,0.0);
  DblNumMat srcPos(over ? vlet()->grdOverExaPos(gNodeIdx) : vlet()->grdSrcExaPos(gNodeIdx)); /* Use oversampling to get good coeff approx */
  DblNumVec srcDen(srcPos.n() * srcDOF); /* Use oversampling to get good coeff approx */

  (this->exsol3d())->quantity(QNT_RHS, srcPos, srcDen);
  
  if (srcDen.linfty() > 0){
	 double scale = ((this->matmgnt())->hom()) ? pow(pow(0.5, (vlet()->depth(gNodeIdx)+(this->rootLevel()))),2.0) : 1.0;
	 map<int, map<int, DblNumMat> >& pin_LevMap = ((this->matmgnt())->hom()) ? _pinv[0] : _pinv[vlet()->depth(gNodeIdx)];
	 map<int, DblNumMat>& pin_map =  pin_LevMap[1];
	 DblNumMat& pinDOF = pin_map[pinType(NK)];
	 if (pinDOF.m() == 0){
		DblNumMat polys(srcPos.n(), NK);
		DblNumMat pin(NK, srcPos.n());
		/* Build off of the node at this depth */
		iC( bldBasPolyMat((this->matmgnt())->hom(), vlet()->cheb(vlet()->kSrcVal()), FRC, gNodeIdx, NK, srcPos, polys));
		iC( pinv(polys, 1e-12, pin));
		pinDOF.resize(NK*srcDOF, srcPos.n() * trgDOF);
		iC( bldDOFMat(pin, pinDOF, srcDOF, trgDOF));
	 }
	 DblNumVec resid(NK*srcDOF);
	 iC( dgemv(scale, pinDOF, srcDen, 1.0, resid));
	 for (int i = 0; i < NK*srcDOF; i++) { 		coeffVec(i) = resid(i); }

	 bool all_zero = true;
	 for (int i = 0; i < NK*srcDOF; i++){ if (abs(coeffVec(i)) != 0.0) all_zero = false; }
	 iA(!all_zero);
  }

  return(0);
}

/* The next two functions take potential values as computed by the fmm algorithm
 * onto a grid on each leaf and treats that solution as coefficients as computed by
 * bldCoeffs up above.  These coefficients are then translated to the actual target
 * positions
 */
template <class VF>
int VFMM3d<VF>::potCoeffs2TrgVal(){
  
  vector<int> ordVec; iC( vlet()->upwOrderCollect(ordVec) );
  for (int i = 0; i < ordVec.size(); i++){
	 int gNodeIdx = ordVec[i];
	 if (vlet()->terminal(gNodeIdx)) {
		iC( potCoeffs2TrgVal(gNodeIdx));
	 }
  }
  return(0);
}

template <class VF>
int VFMM3d<VF>::potCoeffs2TrgVal(int gNodeIdx){
  iA(0);
  /* build coefficients for potential - come straight from grd vals */
  int srcDOF = this->srcDOF(); int trgDOF = this->trgDOF();
  int NK = vlet()->trgNk();
  DblNumVec srcVec(grdExaVal(gNodeIdx));
  DblNumVec coeffVec(NK);
  
  /* Clear coeffVec to make sure the values are zero */
  setvalue(coeffVec,0.0);
  /* Doesn't matter if scale-invariant as building from potential */
  map<int, map<int, DblNumMat> >& pin_LevMap = _pinv[0]; 
  map<int, DblNumMat>& pin_map = pin_LevMap[vlet()->grdTrgSamPos().n()];
  DblNumMat& pinDOF = pin_map[pinType(NK)];
  double scale = pow(pow(0.5,(vlet()->depth(gNodeIdx)+(this->rootLevel()))),2);
  if (pinDOF.m() == 0){
	 DblNumMat bp(vlet()->trgGrdSze(), NK); DblNumMat pin(NK, vlet()->trgGrdSze());
	 iA(vlet()->grdTrgSamPos().n() == vlet()->trgGrdSze());
	 /* Build off of the root node - potential values ALL scale */
	 iA(vlet()->root(0));
	 iC( bldBasPolyMat(true, vlet()->cheb(vlet()->kTrgVal()), POT, 0, vlet()->kTrgVal(), NK, vlet()->grdTrgSamPos(), bp));
	 pinv(bp, 1e-12, pin);
	 pinDOF.resize(NK * srcDOF, vlet()->trgGrdSze()*trgDOF);
	 iC( bldDOFMat(pin, pinDOF, srcDOF, trgDOF));
  }
  iC( dgemv(scale, pinDOF, srcVec, 1.0, coeffVec));
 
  DblNumMat trgExaPos(this->trgExaPos(gNodeIdx));
  DblNumVec trgExaVal(this->trgExaVal(gNodeIdx));
  if ( trgExaPos.n() == 0) return 0;

  DblNumMat c2v(trgExaPos.n(),NK);
  double scalef = 1.0/pow(pow(0.5,(vlet()->depth(gNodeIdx)+(this->rootLevel()))),2);
  setvalue(trgExaVal,0.0);
  iC( bldBasPolyMat(true, vlet()->cheb(vlet()->kTrgVal()), POT, gNodeIdx, vlet()->kTrgVal(), NK, trgExaPos, c2v));
  DblNumMat c2vDOF(trgExaPos.n() * trgDOF,NK * srcDOF);
  iC( bldDOFMat(c2v, c2vDOF, trgDOF, srcDOF));
  
  iC( dgemv(scalef, c2vDOF, coeffVec, 1.0, trgExaVal));
  return(0);
}

/* Take grd values and interpolate to random grid locations - not fully tested */
template <class VF>
int VFMM3d<VF>::grdVals2TrgVal(DblNumVec& trgval){
  int trgDOF = this->trgDOF();
  
  iC( potCoeffs2TrgVal() );

  vector<int> ordVec; iC( vlet()->upwOrderCollect(ordVec) ); //BOTTOM UP
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( vlet()->tag(gNodeIdx) & LET_TRGNODE ) {
		if( vlet()->terminal(gNodeIdx)==true) {
		  DblNumVec grdExaVal(this->grdExaVal(gNodeIdx));
		  vector<int>& curVecIdxs = vlet()->node(gNodeIdx).trgOwnVecIdxs();
		  for(int k=0; k<curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k];
			 for(int d=0; d<trgDOF; d++) {
				trgval(poff*trgDOF+d) = grdExaVal(k*trgDOF+d);
			 }	
		  }
		}
	 }
  }
  return(0);
}


// ----------------------------------------------------------------------
template <class VF>
int VFMM3d<VF>::prntCoeffs(int gNodeIdx)
{
  int srcDOF = this->srcDOF(); int trgDOF = this->trgDOF();
  int NK = vlet()->srcNk();
  DblNumVec& coeffs = *(_coeffs);  iA( coeffs.m()==NK*(vlet()->trmNodeCnt()));
  DblNumVec coeffVec(srcCoeffs(gNodeIdx));
  for (int i = 0; i < NK * srcDOF; i++){
	 std::cout << coeffVec(i) << " ";
  }
  std::cout << endl;
  return(0);
}

/* Build DOF Mat from poly or pinv matrices */
template <class VF>
int VFMM3d<VF>::bldDOFMat(const DblNumMat Ain, DblNumMat& Aout, const int sdof, const int tdof){	
  for (int i = 0; i < Ain.m(); i++){
	 for (int j = 0; j < Ain.n(); j++){ 	
		for (int ds = 0; ds < sdof; ds++){
		  for (int dt = 0; dt < tdof; dt++){
			 /* CHECKME */
			 if (ds == dt)
				Aout(i*sdof + ds,j*tdof + dt) = Ain(i,j);
		  }	
		}
	 }
  }
  return 0;
}

/* For lookup in pinv array only */
template <class VF>
int VFMM3d<VF>::pinType(const int NK){
  int pin_nk;
  if (NK == 1) pin_nk = 0;
  else if (NK == 4) pin_nk = 1;
  else if (NK == 10) pin_nk = 2;
  else if (NK == 20) pin_nk = 3;
  else if (NK == 35) pin_nk = 4;
  else if (NK == 56) pin_nk = 5;
  else if (NK == 84) pin_nk = 6;
  else if (NK == 120) pin_nk = 7;
  else if (NK == 165) pin_nk = 8;
  else if (NK == 220) pin_nk = 9;
  else if (NK == 286) pin_nk = 10;
  else if (NK == 364) pin_nk = 11;
  //else if (NK == 455) pin_nk = 12;
  //else if (NK == 560) pin_nk = 13;
  //else if (NK == 680) pin_nk = 14;
  //else if (NK == 816) pin_nk = 15;
  
  else { iA(0); }
  return pin_nk;
}

template <class VF>
int VFMM3d<VF>::bldBasPolyMat(const bool scale, const bool cheb, const int type, const int gni, const int NK, const DblNumMat pos, DblNumMat& bp, bool dep){
  return bldBasPolyMat(scale, cheb, type, gni, vlet()->kSrcVal(), NK, pos, bp, dep);
}		
/* Build the basis polynomial matrices */
/* bool dep allows us to build polynomial matrices from a depth value directly */
template <class VF>
int VFMM3d<VF>::bldBasPolyMat(const bool scale, const bool cheb, const int type, const int gni, const int KVAL, const int NK, const DblNumMat pos, DblNumMat& bp, bool dep){

  iA(bp.m() == pos.n()); iA(bp.n() >= 0);
  iA(type == POT || type == FRC);
  //int max_nk = bp.n(); /* In case not enough */
  Point3 ctr(dep ? vlet()->center(gni) : Point3(0,0,0));

  /* If dep is true, gni is the depth of the node; otherwise, get the depth of gni */
  double dd = (dep ? (double)((gni) + (this->rootLevel())) : (double)(vlet()->depth(gni) + (this->rootLevel())));
  
  if (bp.n() == pos.n()){
	 for (int i = 0; i < bp.m(); i++){
		for (int j = 0; j < bp.n(); j++){
		  if (i == j) bp(i,j) = 1.0;
		}
	 }
  }

  int K = KVAL;	 
  iA(NK == K*(K+1)*(K+2)/6);

  for (int idx = 0; idx < pos.n(); idx++){
	 double x = pos(0,idx) - ctr(0);
	 double y = pos(1,idx) - ctr(1);
	 double z = pos(2,idx) - ctr(2);
	 if (scale) { /* Not called when !_hom && type == FRC*/
		x = pow(2.0,dd)*x;
		y = pow(2.0,dd)*y;
		z = pow(2.0,dd)*z;
	 }
		 
	 if (!cheb){ /* 4th or 6th order approximations */
		bp(idx,0) = 1.0;
		int ii=0,jj=0;
		for (ii = 1; ii < K; ii++){
		  for (jj = (1 + ii*(ii+1)*(ii+2)/6); jj <= ((ii+1)*(ii+2)*(ii+3)/6 - (ii+1)); jj++){
			 bp(idx,jj-1) = x*bp(idx,jj-(ii*(ii+1)/2) -1);
		  }
		  for (jj = ((ii+1)*(ii+2)*(ii+3)/6 - ii); jj <= ((ii+1)*(ii+2)*(ii+3)/6 - 1); jj++){
			 bp(idx,jj-1) = y*bp(idx,jj-(ii*(ii+3)/2)-1);
		  }
		  jj = (ii+1)*(ii+2)*(ii+3)/6;
		  bp(idx,jj-1) = z*bp(idx,(ii*(ii+1)*(ii+2)/6)-1);
		}
	 }
	 else {
		double T[3*K];
		T[0] = 1.0; T[K] = 1.0;   T[2*K] = 1.0;
		T[1] = x;   T[K + 1] = y; T[2*K + 1] = z;
		for (int i=2; i < K; i++){
		  for (int j = 0; j < 3; j++){
			 int tmp = j*K + i;
			 T[tmp] = 2.0*T[j*K + 1]*T[tmp -1] - T[tmp - 2];
		  }
		}
		int cnt=0;
		//Degree is J
		for (int J = 0; J < K; J++){
		  for (int i = J; i >= 0; i--){
			 for (int j = (J-i); j >= 0; j--){
				int k = (J-(i+j));
				bp(idx,cnt) = T[i]*T[K + j]*T[2*K +k];
				cnt++;
			 }
		  }
		}
		iA( cnt == NK );
	 }
  }
  return(0);
}
