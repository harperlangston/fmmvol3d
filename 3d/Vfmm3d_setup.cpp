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
#include "Vfmm3d.hpp"
#include "common/vecmatop.hpp"

#include "common/memAlloc.hpp"

using namespace std;

using std::cerr;
using std::endl;

using std::istringstream;

// ----------------------------------------------------------------------
template <class VF>
int VFMM3d<VF>::setFromOptions(map<string,string>& opts)
{  
  //-----------------------------------------------------
  map<string,string>::iterator mapIdx;
  mapIdx = opts.find("-" + this->prefix() + "np"); iA(mapIdx!=opts.end());
  { istringstream ss((*mapIdx).second); ss>>(this->np()); }
  mapIdx = opts.find("-" + this->prefix() + "lambda"); assert(mapIdx!=opts.end());
  { istringstream ss((*mapIdx).second);  ss>>_lambda; }
  cerr <<_lambda << endl;
  cerr << (double)(((this->knl()).coefs())[1]) << endl;
  iA(_lambda == (double)(((this->knl()).coefs())[1]));

  double rval; //mostly for Colella test
  mapIdx = opts.find("-" + this->prefix() + "rval"); assert(mapIdx!=opts.end());
  { istringstream ss((*mapIdx).second);  ss>>rval; }
  cerr <<_lambda << " " << (double)(((this->knl()).coefs())[0]) << endl;
  iA(rval == (double)(((this->knl()).coefs())[0]));
  
  mapIdx = opts.find("-" + this->prefix() + "eqnType"); assert(mapIdx!=opts.end());
  { istringstream ss((*mapIdx).second);  ss>>_eqnType; }
  
  //2. decide _eq_mm and _mul_mm, and get matmgnt based on that
  switch((this->knl()).kernelType()) {
	 //laplace kernels
  case KNL_LAP_S_U: (this->_knl_mm) = Kernel3d(KNL_LAP_S_U, (this->knl()).coefs()); break;
  case KNL_LAP_D_U: (this->_knl_mm) = Kernel3d(KNL_LAP_S_U, (this->knl()).coefs()); break;
	 //Mod Lap kernels
  case KNL_MODHEL_S_U: (this->_knl_mm) = Kernel3d(KNL_MODHEL_S_U, (this->knl()).coefs()); break;
  case KNL_MODHEL_D_U: (this->_knl_mm) = Kernel3d(KNL_MODHEL_S_U, (this->knl()).coefs()); break;
	 //stokes kernels
  case KNL_STK_S_U: (this->_knl_mm) = Kernel3d(KNL_STK_F_U, (this->knl()).coefs()); break;
  case KNL_STK_S_P: (this->_knl_mm) = Kernel3d(KNL_LAP_S_U, vector<double>()); break;
  case KNL_STK_D_U: (this->_knl_mm) = Kernel3d(KNL_STK_F_U, (this->knl()).coefs()); break;
  case KNL_STK_D_P: (this->_knl_mm) = Kernel3d(KNL_LAP_S_U, vector<double>()); break;
	 //Unsteady stokes kernels
  case KNL_UNSTK_S_U: (this->_knl_mm) = Kernel3d(KNL_UNSTK_F_U, (this->knl()).coefs()); break;
  case KNL_UNSTK_D_U: (this->_knl_mm) = Kernel3d(KNL_UNSTK_F_U, (this->knl()).coefs()); break;
	 //navier kernels
  case KNL_NAV_S_U: (this->_knl_mm) = Kernel3d(KNL_NAV_S_U, (this->knl()).coefs()); break;
  case KNL_NAV_D_U: (this->_knl_mm) = Kernel3d(KNL_NAV_S_U, (this->knl()).coefs()); break;
	 //others
  case KNL_SQRTLAP: (this->_knl_mm) = Kernel3d(KNL_SQRTLAP, (this->knl()).coefs()); break;
  case KNL_EXP    : (this->_knl_mm) = Kernel3d(KNL_EXP    , (this->knl()).coefs()); break;
  default: iA(0);
  }

  this->_mul_mm = 1; //for the time being
  this->_matmgnt = new MatMgnt3d((this->_knl_mm), (this->np()));
  //matmgnt()  = MatMgnt3d::getmmptr((this->_knl_mm), (this->np()));

  
   this->_let = new VLet3d<VF>(this->prefix()+"vlet3d_");
  vlet()->srcPos()=this->_srcPos;
  vlet()->trgPos()=this->_trgPos;
  vlet()->center()=this->_center;  vlet()->rootLevel()=this->_rootLevel;
  vlet()->srcDen()=this->_srcDen;
  vlet()->exsol3d() =this->_exsol3d;
  vlet()->knl() = this->knl();

  iC( vlet()->setFromOptions(opts) );

  iC( vlet()->setup() );
  vlet()->numVltrs() = -1;
  iC( vlet()->initTree() );


  _tbls = new CmptTbls(this->prefix()+"tbls_");
  _tbls->setup(opts);
  iA(_tbls->ksrcval() == vlet()->kSrcVal());
  iA(_tbls->ktrgval() == vlet()->kTrgVal());
  iA(_tbls->kt() == (this->knl()).kernelType());
  iA(_tbls->np() == (this->np()));
  iA(_tbls->lambda() == _lambda);
        
  return(0);
}

template <class VF>
int VFMM3d<VF>::setup()
{

  bool descend_all = true;
  vlet()->ptsMax() = 1;
  while(vlet()->numVltrs() != 0 && vlet()->nodeVec().size() <= 2250000){ //Reasonable size max
	 cout << "nodevec size = " << vlet()->nodeVec().size() << endl;
	 iC(grdTolData(descend_all) );
	 descend_all = false;
	 iC( vlet()->subVltrs() );
  }
  
  {
	 int cnt = 0;
	 vector<int> ordVec; iC( vlet()->upwOrderCollect(ordVec) );
	 for(int i=0; i<ordVec.size(); i++) {
		int gNodeIdx = ordVec[i];
		if (vlet()->terminal(gNodeIdx)){ cnt++; }
	 }
	 (this->_srcPos)->resize(0,0);
	 (this->_srcPos)->resize(3,cnt);
	 cnt = 0;
	 for(int i=0; i<ordVec.size(); i++) {
		int gNodeIdx = ordVec[i];
		if (vlet()->terminal(gNodeIdx)){ 
		  /* JUST put center in to save space for large grid tests */
		  Point3 ctr(vlet()->center(gNodeIdx));
		  DblNumMat spos(this->dim(), 1, false, (*(this->_srcPos)).data()+cnt*this->dim());
		  for (int d = 0; d < this->dim(); d++){
			 spos(d,0) = ctr(d);
		  }
		  cnt++;
		}
	 }	 
  }

  iC( vlet()->build() );
  iC( vlet()->tblsSrcTrgData() );
  //3. self setup
  iC( grdData() ); //VFMM3d - adds some info for periodicity


  map<int, map<int, DblNumMat> >& pin_LevMap = _pinv[0];
  map<int, DblNumMat>& pin_map =  pin_LevMap[1];
  DblNumMat& pinDOF = pin_map[pinType(vlet()->srcNk())];
  if (pinDOF.m() != 0) pinDOF.resize(0,0);
  DblNumMat& c2vDOF = _coeff2GrdVal[vlet()->srcNk()];
  c2vDOF.resize(0,0);

  
  iC( genNbrTypLsts());  
  if (vlet()->dirichlet() || vlet()->periodic()){
	 iA (this->srcDOF() == 1 && this->trgDOF() == 1); //Only sdof=tdof=1 for now
	 (this->matmgnt())->pdMapsAlloc();
	 iC( bldDirMaps() );
  }


  
  cout << "Building Coefficients" << endl;
  _coeffs = new DblNumVec((vlet()->srcNk() * (this->srcDOF())) * vlet()->trmNodeCnt());
  
  iC( bldCoeffs(FRC) );

  this->_srcExaDen.resize(0);
  this->_srcPos->resize(0,0);
  this->_srcDen->resize(0);
  
  return (0);
}

template <class VF>
int VFMM3d<VF>::grdData()
{
  //1. create vecs
  int grdExaCnt = vlet()->grdExaCnt();
  (this->_trgPos)->resize(3,grdExaCnt/vlet()->trgGrdSze());
  _grdExaVal.resize(grdExaCnt * this->trgDOF());
  
  cout << "Number of source points = " << (*(this->_srcPos)).n() << endl;
  cout << "Number of target points = " << (*(this->_trgPos)).n() << endl;
  int cnt = 0;
  vector<int> ordVec; iC( vlet()->upwOrderCollect(ordVec) );
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(vlet()->tag(gNodeIdx) & LET_SRCNODE) {
		VF& node=vlet()->node(gNodeIdx);
		node.srcNodeIdx() = cnt;
		cnt++;
		if(vlet()->terminal(gNodeIdx)==true) {
		  int beg = node.termidx();
		  int num = 1;
		  /* JUST put center in to save space for large grid tests */
		  Point3 ctr(vlet()->center(gNodeIdx));
		  DblNumMat tpos(this->dim(), num, false, (*(this->_trgPos)).data()+beg*this->dim());
		  for (int d = 0; d < this->dim(); d++){
			 tpos(d,0) = ctr(d);
		  }
		}
	 }
  }
	
  /* let positions */
  vlet()->trgPos()=(this->_trgPos);
  vlet()->trgVal()=(this->_trgVal);


  /* let trg data */
  iC( vlet()->trgData()); cout << "let trgdata" << endl;

  /*
  if (vlet()->periodic()){
	 for(int i=0; i<ordVec.size(); i++) {
		int gNodeIdx = ordVec[i];
		if (vlet()->tag(gNodeIdx) & LET_TRGNODE) {
		  //Let3d::Node& gg = vlet()->node(gNodeIdx);
		  PerNode& gg = vlet()->pernode(gNodeIdx);
		  _nodeVec[gNodeIdx].VinNum() += gg.bdryVnodes().size();
		  for(int a = 0; a < gg.bdryVnodes().size(); a++) {
			 int v = gg.bdryVnodes()[a].gni();
			 _nodeVec[v].VotNum() ++;
		  }
		}
	 }	 		
  }
  */
  
//{ (this->_trgPos)->resize(0,0); } /* Don't need it anymore for grid tests */

  return(0);
}

template <class VF>
int VFMM3d<VF>::grdTolData(bool dsc)
{

  int sdof = this->srcDOF();
  int tdof = this->trgDOF();
  int NK = vlet()->srcNk();
  vector<int> ordVec;  iC( vlet()->upwOrderCollect(ordVec) );

  vlet()->numVltrs() = 0; /* reset number of vltrsop ! */
  double _maxval = 0.0;
  vector<int> vltrVec;
  for (int i = 0; i < ordVec.size(); i++){
    int gNodeIdx = ordVec[i];
    VF& node=vlet()->node(gNodeIdx);
    node.prmVltr() = false;
    node.scdVltr() = false;
    bool test = ((dsc || node.dscPrmVltr()) ? true : false);
    if((test && vlet()->terminal(gNodeIdx)) && (vlet()->depth(gNodeIdx)) < vlet()->maxLevel()) {
      //if (let()->terminal(gNodeIdx) && let()->tag(gNodeIdx) & LET_TRGNODE) { //evaluator
      vltrVec.push_back(gNodeIdx);
      //}
    }
  }

#pragma omp parallel for
  for(int i=0; i<vltrVec.size(); i++) {
    int gNodeIdx = vltrVec[i];
	 VF& node=vlet()->node(gNodeIdx);
	 node.prmVltr() = false;
	 node.scdVltr() = false;
	 if (vlet()->adaptive()){
		bool test = ((dsc || node.dscPrmVltr()) ? true : false);

		if((test && vlet()->terminal(gNodeIdx)) && (vlet()->depth(gNodeIdx)) < vlet()->maxLevel()) {
		  int num = (vlet()->grdSrcSamPos().n());
		  num = NK;
		
		  DblNumMat srcPos(vlet()->grdSrcExaPos(gNodeIdx));
		  DblNumVec coeffs(num*sdof);

		  DblNumVec srcDen(vlet()->grdSrcSamPos().n() * sdof);
		  (this->exsol3d())->quantity(QNT_RHS, srcPos, srcDen);

		  double scale = pow(pow(0.5,(vlet()->depth(gNodeIdx)+(this->_rootLevel))),2.0);
		  DblNumMat pinDOF(num*sdof, srcPos.n() * this->trgDOF());
		  DblNumMat polys(srcPos.n(), num);
		  DblNumMat pin(num, srcPos.n());

		  
		  /* Build off of the node at this depth */
		  iC( bldBasPolyMat(true, vlet()->cheb(vlet()->kSrcVal()), POT, 0, NK, vlet()->grdSrcSamPos(), polys));

		  //iC( pinv(polys, 1e-14, pin));
		  iC( pinv(polys, 0, pin, 1));
		  iC( bldDOFMat(pin, pinDOF, sdof, tdof));

		  DblNumVec resid(num*sdof);
		  iC( dgemv(1.0, pinDOF, srcDen, 1.0, resid));
		  for (int i = 0; i < num*sdof; i++) coeffs(i) = resid(i);

		  DblNumMat gniPos(vlet()->grdDblSrcExaPos(gNodeIdx));
		  DblNumMat c2vDOF(gniPos.n()* this->trgDOF(), num*sdof);
		  DblNumMat c2v(gniPos.n(),num);
		  iC( bldBasPolyMat(true, vlet()->cheb(vlet()->kSrcVal()), POT, 0, NK, vlet()->grdDblSrcSamPos(), c2v));
		  iC( bldDOFMat(c2v, c2vDOF, tdof, sdof));				
		  DblNumVec val(gniPos.n() * sdof); setvalue(val,0.0);
		  iC( dgemv(1.0, c2vDOF, coeffs, 1.0, val));
		  DblNumVec realVal(val.m());
		  (this->exsol3d())->quantity(QNT_RHS, gniPos, realVal);
	 
		  node.prmVltr() = false;
		  if (1 || !(srcDen.linfty() == 0.0 && realVal.linfty() == 0.0)){
			 double diff[sdof];
			 double ldiff[sdof];
			 double sum[sdof];
			 double mainf[sdof];
			 double maxval[sdof];
			 for (int s = 0; s < sdof; s++){
				diff[s] = 0.0;
				ldiff[s] = 0.0;
				sum[s] = 0.0;
				mainf[s] = 1e-16;
				maxval[s] = 0.0;
			 }
			 double h = vlet()->radius(gNodeIdx)*2.0/vlet()->kSrcVal();
			 for (int i = 0; i < gniPos.n(); i++){
				for (int s = 0; s < sdof; s++){
				  diff[s] += (realVal(i*sdof + s) - val(i*sdof + s))*(realVal(i*sdof + s) - val(i*sdof + s));
				  ldiff[s] += diff[s]*h*h*h;
				  sum[s] += realVal(i*sdof + s)*realVal(i*sdof + s);
				  mainf[s] = max(mainf[s], abs(realVal(i*sdof + s) - val(i*sdof + s)));
				  maxval[s] = max(maxval[s], abs(realVal(i*sdof + s)));
				}
			 }
		  
		
			 double volume = (2.0*vlet()->radius(gNodeIdx))*(2.0*vlet()->radius(gNodeIdx))*(2.0*vlet()->radius(gNodeIdx));

			 DblNumVec comp_val(sdof); setvalue(comp_val,0.0);
			 for (int s = 0; s < sdof; s++){

				comp_val(s) = max(mainf[s],sqrt(diff[s]))*(h*h*h);
			 }
			 _maxval = max(_maxval,maxval[0]);
			 Point3 cdif(vlet()->center(gNodeIdx) - (Point3(1.0,1.0,1.0)));
			 if (comp_val.linfty() > (vlet()->eps_rhs()) || vlet()->depth(gNodeIdx) <= 1) {
				node.prmVltr() = true;
				vlet()->numVltrs()++;
			 }
			 else {
				node.prmVltr() = false;
			 }
		  }
		  else {
			 //cerr << gNodeIdx << " = 0" << endl;
		  }
		}
	 }
	 else {
		if (vlet()->depth(gNodeIdx) + this->rootLevel() < (vlet()->rhs())){
		  node.prmVltr() = true;
		  vlet()->numVltrs()++;
		}
		else {
		  node.prmVltr() = false;
		}
	 }
  }
  return(0);
}

template <class VF>
int VFMM3d<VF>::bldDirMaps(){
  
  string str = "../include/DIRICHLET/DwnEquGrdMap_";
  stringstream converter;
  converter << (this->np());
  str.append(converter.str());
  ifstream indata; // indata is like cin
  
  {
	 indata.open(str.c_str()); // opens the file
	 if(!indata) { // file couldn't be opened
		cerr << "Error: Dwn Equ Grd Map file could not be opened" << endl;
		exit(1);
	 }
  }

  int num;
  /* Need to fix these to just be arrays as opposed to matrices to save expense */
  for (int i = 0; i < 8; i++){
	 DblNumMat& _DEM = ((this->matmgnt())->perdirmaps())->dwnEquGrdMap()[i]; /* For Dwn Equ Mapping for Dirichlet */
	 _DEM.resize(this->datSze(DE), this->datSze(DE));
	 double flip;
	 if ( i == 1 || i == 2 || i == 4 || i == 7) flip = -1.0;
	 else flip = 1.0;
	 for (int j = 0; j < this->datSze(DE); j++){
		indata >> num; 
		_DEM(j, num) = flip*1.0;
	 }
  }

  indata.close();
  return(0);
}

